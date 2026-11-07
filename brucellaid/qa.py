"""
The facts a quality assurance report rests on, in a form the HTML and the PDF renderers share.

A result is only as good as what can be checked about it afterwards. That means, for every sample: which
files went in and what they were (size, date, checksum), every distance MASH measured rather than only the
winning one, what the sequence type contributed to the species call, and how far each allele call was ahead
of the allele it beat. And for the run: who ran what, where, with which versions, and against which
reference data, identified by checksum rather than by name.

Nothing here renders anything. The numbers are gathered once and the two report formats lay the same numbers
out differently.
"""

import getpass
import hashlib
import os
import platform
import socket
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from .errors import BrucellaIDError
from .mlst import NEW_PROFILE, UNKNOWN
from .references import EXTENDED_SKETCH_FILE, data_path, lookup
from .species import CONFLICT, FROM_ST, UNRESOLVED

BLOCK = 1 << 20  # Read a file for its checksum a megabyte at a time: an input file can be gigabytes

OK = 'ok'
WARNING = 'warning'
FAILED = 'failed'

# How many reference genomes a sample block of a report lists before it stops: enough to show the call
# and what it beat, few enough that a sample stays on one page.
TOP_HITS = 5


def file_md5(path):
    """
    MD5 of a file.

    usedforsecurity=False because this is a checksum of a file a reviewer may want to match against the
    one the sequencer wrote, not a security measure; without it, MD5 is unavailable on a FIPS enabled host.
    """
    digest = hashlib.md5(usedforsecurity=False)
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(BLOCK), b''):
            digest.update(block)
    return digest.hexdigest()


def human_size(size):
    """A file size a person can read: 1.2 GB."""
    for unit in ('B', 'kB', 'MB', 'GB', 'TB'):
        if size < 1024 or unit == 'TB':
            return '{:.0f} {}'.format(size, unit) if unit == 'B' else '{:.1f} {}'.format(size, unit)
        size /= 1024.0
    return '{:.1f} TB'.format(size)  # pragma: no cover - the loop returns at TB


@dataclass(frozen=True)
class FileFacts:
    """What one input file was when the report was written."""
    path: str  # Absolute, with symbolic links left alone: that is what the operator typed
    size: int
    modified: str
    md5: str = ''  # Empty with --no-md5
    target: str = ''  # The real file, when path is a symbolic link

    @classmethod
    def describe(cls, path, md5=True):
        stat = os.stat(path)
        absolute = os.path.abspath(path)
        real = os.path.realpath(path)
        modified = datetime.fromtimestamp(stat.st_mtime).astimezone().strftime('%Y-%m-%d %H:%M:%S %Z')
        return cls(absolute, stat.st_size, modified, file_md5(path) if md5 else '',
                   real if real != absolute else '')

    @property
    def name(self):
        return os.path.basename(self.path)

    @property
    def size_text(self):
        return human_size(self.size)


def describe_inputs(results, md5=True):
    """
    {sample: [FileFacts]} for every sample of a run.

    A file that has gone away since the run is skipped rather than raised over: the report is written after
    the analysis, and a missing input should not cost the operator the whole report.
    """
    described = {}
    for result in results:
        facts = []
        for path in result.files:
            try:
                facts.append(FileFacts.describe(path, md5=md5))
            except OSError:
                continue
        described[result.sample] = facts
    return described


def fasta_count(path):
    """How many sequences a fasta file holds, for the line that describes it in a report."""
    count = 0
    with open(path) as f:
        for line in f:
            count += line.startswith('>')
    return count


def user_name():
    """The login name, or the numeric user id when there is no name (some containers have none)."""
    try:
        return getpass.getuser()
    except (KeyError, OSError):
        return str(os.getuid()) if hasattr(os, 'getuid') else 'unknown'


def environment(operator=None):
    """Who ran it, on what, and where."""
    user = user_name()
    return {'operator': operator or user, 'user': user, 'host': socket.gethostname(),
            'system': platform.platform(), 'working_directory': os.getcwd()}


def is_packaged(path):
    """
    Whether a reference data file is the one that ships with the package.

    Any packaged file would do to find the folder; the extended sketch is only the one it is cheapest to
    name. A run that never needs that file -- "--reference-set vsnp", or a fully custom pair -- must not
    fail here because it is missing, so a missing data file means "not packaged" rather than an error.
    """
    try:
        return Path(path).resolve().parent == Path(data_path(EXTENDED_SKETCH_FILE)).resolve().parent
    except (OSError, BrucellaIDError):
        return False


def data_facts(files, md5=True):
    """
    The reference data a run measured against, each file identified by its checksum.

    A report that names "the MLST profile table" proves nothing; one that gives its MD5 lets a reviewer show
    that two runs a year apart used the same table, or that they did not.

    :param files: [(name, path, detail)], the files the run actually used
    :return: [{'name', 'path', 'detail', 'md5', 'packaged'}], skipping the files that cannot be read
    """
    facts = []
    for name, path, detail in files:
        if not path or not os.path.isfile(path):
            continue
        facts.append({'name': name, 'path': os.path.abspath(path), 'detail': detail,
                      'md5': file_md5(path) if md5 else '', 'packaged': is_packaged(path)})
    return facts


def status(result):
    """
    Whether a sample needs a second look: "ok", "warning" or "failed".

    A warning is for the things a reviewer has to resolve, not for every note: a read set that was not
    identified, a sequence type that could not be called or is not in the table, a species the two kinds
    of evidence disagree about, and a species that rests on the sequence type alone because the distance
    had none to offer. A sample whose second reference was nearly as close is normal, and says so in its
    notes without being flagged.
    """
    if result.error:
        return FAILED
    if not result.identified:
        return WARNING
    if result.mlst is not None and result.sequence_type in (UNKNOWN, NEW_PROFILE):
        return WARNING
    if result.call is not None and result.call.verdict in (CONFLICT, UNRESOLVED):
        return WARNING
    # Only one of the two kinds of evidence spoke. A species from the ST while MASH still named a
    # reference is the ordinary arbitration and stays "ok"; a species from the ST with no reference behind
    # it means the sketch and the table do not match, or the distance separated nothing, and half the
    # evidence this report exists to show is missing.
    if result.call is not None and result.call.source == FROM_ST and result.call.reference is None:
        return WARNING
    return OK


def counted(results):
    """{status: how many samples}, for the summary line of a report."""
    counts = {OK: 0, WARNING: 0, FAILED: 0}
    for result in results:
        counts[status(result)] += 1
    return counts


def groups(results):
    """
    The samples of a run, grouped by what they were called: the same species and the same ST together.

    The largest group first, so a report opens on the bulk of the run, and the failed samples last, where a
    reviewer looks for them. Within a group the samples keep the order they were analysed in.

    :return: [((species, ST), [SampleResult])], the failed samples under ('', '')
    """
    buckets = {}
    for result in results:
        key = ('', '') if result.error else (result.species, result.sequence_type)
        buckets.setdefault(key, []).append(result)
    failed = buckets.pop(('', ''), None)
    ordered = sorted(buckets.items(), key=lambda item: (-len(item[1]), item[0]))
    if failed:
        ordered.append((('', ''), failed))
    return ordered


def margin(call):
    """
    How far an allele call was ahead of the allele it had to beat: (difference, ratio).

    The ratio is what says whether a call is comfortable. A locus whose runner-up is one read behind is a
    call, but it is the one to look at when a profile comes back as "new".

    :return: (reads ahead, ratio) or (None, None) when there was no call or no runner-up
    """
    if not call.allele or not call.runner_up:
        return None, None
    if not call.runner_up_reads:
        return call.reads, None
    return call.reads - call.runner_up_reads, call.reads / call.runner_up_reads


def close_hits(result, references):
    """
    The MASH hits that were in the running: the closest, and the ones near enough to be reported with it.

    :return: [(MashHit, Reference or None, is_close)] in the order MASH ranked them, closest first
    """
    if result.identification is None:
        return []
    close = {id(hit) for hit in result.identification.close_matches}
    rows = []
    for hit in result.identification.hits:
        reference = lookup(references, hit.reference)
        rows.append((hit, reference, id(hit) in close))
    return rows


def hits_caption(rows, hidden):
    """
    What to say under the table of closest reference genomes.

    "Closest first", and "the rest are further away", are untrue whenever more references tie at the best
    distance than the table has room for. That happens two ways, both of them from a coarse distance: a
    read set sharing no hashes sits at distance 1 from all 432 of the default set, and one sharing a few
    hashes ties every reference that shares the same few. Either way the rows shown are an arbitrary slice
    of the tie, and the ones left out are not further away. Saying so is the point of the table.

    :param rows: every row from close_hits(), not only the ones being shown
    :param hidden: how many top_hits() left out
    """
    if not rows:
        return 'Closest first.'
    tied = sum(1 for hit, _, _ in rows if hit.distance == rows[0][0].distance)
    if tied > len(rows) - hidden:
        lead = ('Every one of these {} reference genomes is'.format(len(rows)) if tied == len(rows)
                else '{} of these reference genomes are'.format(tied))
        caption = ('{} the same distance from the reads, so the rows below are an arbitrary few of a tie '
                   'rather than the closest.'.format(lead))
        return caption + ' brucellaid_mash.tsv has every distance.' if hidden else caption
    caption = 'Closest first.'
    if any(close for _, _, close in rows):
        caption += ' A shaded row is too close to the first to be separated by distance alone.'
    if hidden:
        caption += (' The other {} reference(s) are further away; brucellaid_mash.tsv has every '
                    'distance.'.format(hidden))
    return caption


def top_hits(rows, top=TOP_HITS):
    """
    The hits a report shows: the closest few, and any other that is near enough to be in the running.

    A report that lists every reference would run to pages of rows saying the sample is not Brucella ovis,
    and the default set has 432 of them. The ones that matter are the closest handful and anything inside
    the close-match
    margin, which is the evidence for the species call; the rest are in brucellaid_mash.tsv, which holds
    every distance.

    :param rows: [(hit, reference, is_close)] from close_hits()
    :return: (rows to show, how many were left out)
    """
    shown = [row for index, row in enumerate(rows) if index < top or row[2]]
    return shown, len(rows) - len(shown)


def warnings(results):
    """
    [(sample, status, text)] of everything a reviewer has to read: the errors first, then the notes.

    One row per note rather than one row per sample, because a sample with four notes is exactly the sample
    whose notes should not be run together into a paragraph.
    """
    rows = []
    for result in results:
        if result.error:
            rows.append((result.sample, FAILED, result.error))
    for result in results:
        for note in result.notes:
            if note:
                rows.append((result.sample, status(result), note))
    return rows
