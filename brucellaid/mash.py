"""Find the closest Brucella reference genome of a read set with MASH."""

import contextlib
import gzip
import logging
import shlex
import shutil
import subprocess
import sys
import zlib
from dataclasses import dataclass
from functools import cache
from pathlib import Path

from .errors import BrucellaIDError

log = logging.getLogger(__name__)

# Minimum number of the reference's 1000 min-hashes the read set must share for the identification to be
# reported. The first version of BrucellaID used the same threshold: below it the read set is not Brucella, or
# holds too little of it to tell.
DEFAULT_MIN_HASHES = 150
# A second reference within this distance of the best one is reported as a close match: the Brucella genomes
# differ by little, and the biovar references in particular are nearly identical.
CLOSE_MATCH_MARGIN = 0.0005
# More references than this within the margin means the margin is not separating anything: see
# Identification.undiscriminating. The densest real case over the 432 packaged genomes is 15.
MAX_CLOSE_MATCHES = 40
# Beyond this distance the closest reference is reported as a loose match. It says nothing about which
# species is right: it says the reads are not close to anything in the sketch, which happens when the
# library holds more than one organism and when the sample belongs to a lineage the sketch has no near
# relative for.
#
# The rule is the 95th percentile of the distance of calls that AGREE with the ENA submitter, taken to a
# round number: past that, the closest reference is further away than a correctly identified library
# normally sits. Over the public sweep, re-measured under the copy threshold above, that percentile is
# 0.00066 and this threshold flags 4.0% of the 4,020 identified runs -- as rare as the 4.9% the 2.0
# calibration accepted, which had a percentile of 0.00464 and a threshold of 0.005. The old threshold was
# not wrong for the sketches it was measured on; it was calibrated on distances the old copy threshold had
# inflated about sevenfold, and left here it would have missed the thin-reference cases it exists to
# catch -- a B. pituitosa run at 0.0029, for one.
#
# Two things say the choice is not on a knife edge. Everything from 0.0005 to 0.002 flags 2.8% to 6.6% of
# runs with 96% to 98% of them still agreeing with the submitter, so the behaviour barely moves; and the
# ordinary libraries are far below it -- the 3,226 runs already within 0.001 before now sit within
# 0.00046, so this clears the clean bulk by better than a factor of two. Rounding the percentile to one
# significant figure would give 0.0007 and behave the same; 0.001 is preferred for the headroom and
# because it is a value a user may reasonably type.
#
# What justifies the note NOT changing the species is a separate finding, and the only non-circular one
# here: agreement with the submitter is flat in distance (99.8%, 100%, 95.3%, 100%, 100% across the bands
# in the data file), so a loose match does not predict a wrong species. Note that "the agreeing calls past
# the threshold all shared few hashes", which earlier versions of this comment offered as corroboration,
# is no evidence at all: mash derives the distance from the shared-hash count, so the two are one
# measurement in different units and the claim is entailed by the threshold.
#
# scripts/loose_distance_calibration.tsv holds the measurements and
# scripts/recalibrate_loose_distance.py prints the tables above from it.
LOOSE_MATCH_DISTANCE = 0.001
# How many copies of a k-mer the read set must hold before MASH sketches it. The point is to drop the
# k-mers that only exist because of a sequencing error, and the right number grows with the coverage: a
# k-mer seen once in a 5x library is almost certainly an error, while at 70x the errors reach several
# copies on their own and crowd the true k-mers out of the bottom-1000 sketch. So the threshold is derived
# from the sample rather than fixed.
#
# Re-derived by measurement. Three real libraries were subsampled by striding -- which separates depth from
# library quality, as comparing two different libraries cannot -- and swept over every threshold from 1 to
# 20 at each depth, against the reference each one matches. Recovery against depth is a single-peaked
# curve, so the target is the PLATEAU: every threshold within 1% of the best recovery at that depth. The
# plateau is narrow when the library is shallow and wide when it is deep, and it moves up with depth:
#
# Each cell below is the UNION of the plateaus of the curves at that depth, which is the rule used for
# both tables here: the widest span any curve there allowed.
#
#     depth     5x   9x  14x  20x  28x  41x  46x  54x  68x
#     plateau  2-2  2-2  2-3  4-4  3-9 6-11 4-13 6-15 7-20
#
# Those three libraries stop at 70x, and the median public Brucella run is 141x, so a second set of 60 ENA
# runs was swept at their own depth from 51x to 1903x, thresholds to 64. 81 curves from 63 libraries in
# all. One more copy per 5x with a ceiling of 48 lands inside the plateau for 80 of the 81, recovering at
# least 99.6% of the best for nine curves in ten and 100% at the median.
#
#     depth      51-54x  101-112x  198-250x  420-535x  700-780x  1218-1903x
#     plateau       6-24      3-48      6-64     32-64     12-64        20-64
#
# What this replaces was one more copy per 20x with a ceiling of 20, which lands outside for 55 of the 81
# and as low as 63.6%. It was about three times too cautious in the middle of the range, so a deep library
# was measured with a sketch half full of error k-mers and its distance came out inflated -- ERR2136546 at
# 47x read 0.00652 against 773/1000 hashes where it is really 0.00012 against 995. Per 4x and per 6x both
# land outside more often than per 5x, and ceilings of 24, 32 and 40 all land outside more often than 48,
# so neither number is a point on a flat surface. Simulated reads cannot settle any of this: uniform error
# and no duplicates make their plateau 3-20 at every depth above 50x, which chooses nothing.
#
# The single curve the rule misses is SRR24420453, and it is not a Brucella library: at a nominal 1482x it
# never shares more than 396 of 1000 hashes with anything, its best distance is 0.027, and its recovery
# peaks at 20 copies and then falls -- the signature of a minority Brucella population being filtered out
# along with the errors. Its distance means nothing at any threshold, and both the shared-hash count and
# the loose-match note say so.
MIN_KMER_COPIES = 2             # Floor: the plateau is exactly 2-2 at every depth up to about 10x
MAX_KMER_COPIES = 48            # Ceiling: it binds past 240x, where the plateau runs to 64
DEPTH_PER_KMER_COPY = 5         # One more required copy per this much coverage
KMER_COPIES_AUTO = 0            # "Derive it from the coverage", the default
BRUCELLA_GENOME_SIZE = 3_300_000  # Used only to turn a number of bases into a coverage
BASE_FRACTION_RECORDS = 1000    # Records read from the head of a file to measure its bases per byte
RATIO_SAMPLE = 4 * 1024 * 1024  # Compressed bytes read to measure a file's compression ratio
READ_CHUNK = 1024 * 1024


class MashError(BrucellaIDError):
    pass


@cache
def executable():
    """
    Path of mash: the one in the PATH, otherwise the one installed next to the running Python.
    The fallback makes "/path/to/env/bin/brucellaid" work without activating the conda environment.

    :return: path to mash, or None if it cannot be found
    """
    found = shutil.which('mash')
    if found:
        return found
    # Do not resolve symlinks: a virtualenv's python is a symlink to the base interpreter
    candidate = Path(sys.executable).parent / 'mash'
    if candidate.is_file():
        return str(candidate)
    return None


def check_mash():
    """Make sure mash is installed and return its path."""
    path = executable()
    if path is None:
        raise MashError('"mash" was not found in your PATH or in "{}". Install MASH in the environment with '
                        '"conda install -c bioconda mash".'.format(Path(sys.executable).parent))
    return path


def versions():
    """{"MASH": version} of the installation, "unknown" when it cannot be read."""
    path = executable()
    if path is None:
        return {'MASH': 'unknown'}
    try:
        out = subprocess.run([path, '--version'], capture_output=True, text=True, timeout=60)
        version = out.stdout.strip().splitlines()
        return {'MASH': version[0].strip() if version else 'unknown'}
    except (OSError, subprocess.SubprocessError):
        return {'MASH': 'unknown'}


def parse_info(output, sketch='the sketch'):
    """
    Read the table "mash info -t" writes, and return the ID of each genome of the sketch.

    #Hashes	Length	ID	Comment
    1000	3286445	NC_006932-NC_006933.fasta	[2 seqs] NC_006932.1 Brucella abortus biovar 1 ...

    :return: list of IDs, in the order of the sketch
    """
    ids = []
    for line in output.splitlines():
        fields = line.rstrip('\n').split('\t')
        if not line.strip() or line.startswith('#'):
            continue
        if len(fields) < 3:
            raise MashError('Unexpected line in the "mash info" output for {}: {}'.format(sketch, line.strip()))
        ids.append(fields[2])
    if not ids:
        raise MashError('The reference sketch {} holds no genome.'.format(sketch))
    return ids


@cache
def sketch_references(sketch):
    """
    The IDs of the genomes a sketch holds, read with "mash info -t".

    The ID of a sketch is the name of the file it was built from, which is what "mash dist" reports and what
    the reference table keys on.

    :return: list of IDs, in the order of the sketch
    """
    cmd = [check_mash(), 'info', '-t', str(sketch)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise MashError('The reference sketch {} cannot be read by MASH: {}'.format(
            sketch, (proc.stderr or proc.stdout).strip() or 'no output'))
    return parse_info(proc.stdout, sketch)


@dataclass
class MashHit:
    """One line of "mash dist": how close the read set is to one reference genome."""
    reference: str  # The sketch ID of the reference, which is the name of the file it was built from
    distance: float  # MASH distance, an estimate of the fraction of positions that differ
    p_value: float
    shared_hashes: int
    total_hashes: int

    @property
    def identity(self):
        """The MASH distance as a percentage of identity, as MASH's documentation reads it: (1 - distance)."""
        return 100.0 * (1.0 - self.distance)


def parse_dist(output):
    """
    Read the table "mash dist" writes, closest reference first.

    Bceti1Cudo.fasta	stdin	0.0294416	0	493/1000
    """
    hits = []
    for line in output.splitlines():
        fields = line.rstrip('\n').split('\t')
        if not line.strip():
            continue
        if len(fields) < 5:
            raise MashError('Unexpected line in the MASH dist output: {}'.format(line.strip()))
        shared, _, total = fields[4].partition('/')
        try:
            hits.append(MashHit(fields[0], float(fields[2]), float(fields[3]), int(shared), int(total)))
        except ValueError:
            raise MashError('Unexpected line in the MASH dist output: {}'.format(line.strip())) from None
    # Closest first: by distance, then by shared hashes, then by name, so that the order never depends on the
    # order mash happened to write
    hits.sort(key=lambda hit: (hit.distance, -hit.shared_hashes, hit.reference))
    return hits


@contextlib.contextmanager
def opened_reads(path):
    """Open a read file for binary reading, decompressing it when it is gzipped."""
    with open(path, 'rb') as f:
        gzipped = f.read(2) == b'\x1f\x8b'
    with (gzip.open(path, 'rb') if gzipped else open(path, 'rb')) as reads:
        yield reads


def write_reads(files, stream):
    """
    Copy the reads of one sample into stream, decompressing the gzipped files.

    The files of a sample are piped into one mash sketch: mash sketches each of its input files on its own, and
    the R1 and R2 files of a sample are one read set.
    """
    for path in files:
        with opened_reads(path) as reads:
            shutil.copyfileobj(reads, stream, READ_CHUNK)


def fastq_bases_per_byte(path, records=BASE_FRACTION_RECORDS):
    """
    Bases per uncompressed byte of a FASTQ file, measured on its first records.

    A FASTQ record spends its bytes on the sequence, an equally long quality string, two description lines
    and four newlines, so the fraction is around 0.44 but depends on the read length and on how long the
    descriptions are. Measuring it on the head of the file is exact enough to pick a k-mer copy threshold
    and costs one short read.

    :return: bases per byte, 0.0 when the file holds no readable record
    """
    bases = total = 0
    try:
        with opened_reads(path) as reads:
            for number, line in enumerate(reads):
                if number >= records * 4:
                    break
                total += len(line)
                if number % 4 == 1:
                    bases += len(line.rstrip(b'\r\n'))
    except (OSError, EOFError, zlib.error):
        return 0.0
    return bases / total if total else 0.0


def gzip_sample(path, sample=RATIO_SAMPLE):
    """
    (uncompressed bytes per compressed byte, gzip members seen) on the head of a gzipped file.

    A .gz file may be several gzip members concatenated -- lane files cat'ed together, or the 64 KB blocks
    bgzip writes -- so the sample is decompressed member by member: a ratio taken on the first member alone
    would say nothing about the rest, and the member count tells the caller whether the file's own trailer
    can be trusted at all.

    :return: (ratio, members); (0.0, 0) when the file cannot be read or holds no gzip stream
    """
    try:
        with open(path, 'rb') as f:
            chunk = f.read(sample)
    except OSError:
        return 0.0, 0
    if not chunk:
        return 0.0, 0
    out, members, data = 0, 0, chunk
    try:
        while data:
            decompressor = zlib.decompressobj(47)   # 47: gzip or zlib, whichever the member is
            out += len(decompressor.decompress(data))
            members += 1
            if not decompressor.eof:   # The sample ends inside this member
                break
            data = decompressor.unused_data
    except zlib.error:
        if not out:   # Nothing decompressed at all: not a gzip stream after all
            return 0.0, 0
    return out / len(chunk), members


def uncompressed_size(path):
    """
    Size the file would have uncompressed, without decompressing it.

    A gzip member ends with the uncompressed size of that member modulo 4 GiB, which is exact for the files
    most samples are made of. But the trailer describes only the LAST member, so for a file of several
    members -- cat'ed lane files, or bgzip's 64 KB blocks -- and for a member past 4 GiB, the file's own
    measured compression ratio decides instead: outright when the head sample saw more than one member, and
    as the arbiter whenever the trailer disagrees with it by more than the head of a file can drift from
    the rest.

    :return: number of bytes, 0 when the file cannot be read
    """
    path = Path(path)
    try:
        size = path.stat().st_size
        if size < 4:
            return size
        with open(path, 'rb') as f:
            gzipped = f.read(2) == b'\x1f\x8b'
            if not gzipped:
                return size
            f.seek(-4, 2)
            trailer = int.from_bytes(f.read(4), 'little')
    except OSError:
        return 0
    ratio, members = gzip_sample(path)
    if ratio <= 0:
        return trailer
    estimate = round(size * ratio)
    if members > 1:
        return estimate   # Concatenated members: the trailer describes only the last one
    wraps = round((size * ratio - trailer) / 2 ** 32)
    exact = trailer + max(0, wraps) * 2 ** 32
    if not 0.7 * estimate <= exact <= 1.4 * estimate:
        return estimate   # A trailer the measured ratio contradicts: members the sample never reached
    return exact


def estimate_bases(files):
    """
    Number of bases in the read files of one sample, without reading them through.

    :return: number of bases, 0 when none of the files can be read
    """
    total = 0.0
    for path in files:
        total += uncompressed_size(path) * fastq_bases_per_byte(path)
    return int(total)


def min_kmer_copies(bases, genome_size=BRUCELLA_GENOME_SIZE):
    """
    How many copies of a k-mer to require of a read set holding this many bases.

    One more copy per DEPTH_PER_KMER_COPY of coverage, never below MIN_KMER_COPIES, never above
    MAX_KMER_COPIES. See the comment on MIN_KMER_COPIES for where the numbers come from.
    """
    if not bases or genome_size <= 0:
        return MIN_KMER_COPIES
    depth = bases / genome_size
    return int(min(MAX_KMER_COPIES, max(MIN_KMER_COPIES, round(depth / DEPTH_PER_KMER_COPY))))


def run_dist(files, sketch, threads=1, min_copies=MIN_KMER_COPIES):
    """
    Compare the reads of one sample to every reference genome of the sketch.

    :param files: one fastq file, or two paired-end fastq files
    :param sketch: path of the MASH sketch of the reference genomes
    :return: list of MashHit, closest reference first
    """
    # No -k or -s: mash refuses them when the reference is a sketch, and the read set inherits them from it
    cmd = [check_mash(), 'dist', '-p', str(threads), '-r', '-m', str(min_copies), str(sketch), '-']
    log.debug('Running: %s < %s', shlex.join(cmd), ' '.join(str(f) for f in files))
    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=False)
    try:
        write_reads(files, proc.stdin)
    except BrokenPipeError:  # mash stopped reading: it failed, and its error is the useful one
        pass
    except (OSError, EOFError, zlib.error) as e:  # An unreadable or truncated file, or a broken gzip stream
        proc.kill()
        proc.wait()
        raise MashError('The reads of the sample cannot be read: {}'.format(e)) from None
    finally:
        with contextlib.suppress(BrokenPipeError):
            proc.stdin.close()
        proc.stdin = None  # communicate() would flush a pipe that is already closed
    stdout, stderr = proc.communicate()
    stdout = stdout.decode('utf-8', 'replace')
    stderr = stderr.decode('utf-8', 'replace').strip()
    if proc.returncode != 0:
        raise MashError('MASH failed (exit code {}): {}\nCommand: {}'.format(
            proc.returncode, stderr.splitlines()[-1].strip() if stderr else 'no output', shlex.join(cmd)))
    log.debug('MASH output:\n%s', stderr)
    hits = parse_dist(stdout)
    if not hits:
        raise MashError('MASH compared the reads to no reference: the sketch {} is empty.'.format(sketch))
    return hits


@dataclass
class Identification:
    """What MASH says about one sample."""
    hits: list  # list of MashHit, closest reference first
    min_hashes: int = DEFAULT_MIN_HASHES
    min_copies: int = MIN_KMER_COPIES  # k-mer copies required of the read set, see min_kmer_copies()
    coverage: float | None = None      # Estimated coverage the threshold was derived from, None when given

    @property
    def best(self):
        return self.hits[0]

    @property
    def identified(self):
        """Whether the closest reference shares enough hashes with the read set to be reported."""
        return self.best.shared_hashes >= self.min_hashes

    @property
    def within_margin(self):
        """Every other reference within CLOSE_MATCH_MARGIN of the closest one, closest first."""
        return [hit for hit in self.hits[1:] if hit.distance - self.best.distance <= CLOSE_MATCH_MARGIN]

    @property
    def undiscriminating(self):
        """
        How many references sit within the margin when that number is too large to be evidence, else 0.

        A close match is meant to say "MASH could not choose between these two". It can only say that
        while the distance is resolved finely enough to separate them, and at a low shared-hash count it
        is not: mash derives the distance from the shared count, so every reference sharing the same few
        hashes gets exactly the same distance and ties with the best one. Sharing 4 of 1000 hashes has
        been seen to tie 127 references across twelve species, which is not ambiguity between two
        candidates but an absence of measurement.

        The cap is set against the densest real case: over the 432 genomes of the default set, comparing
        each with all the others, at most 15 others ever fall within the margin of a genome's closest
        hit, and the 99th percentile is 12. (Counting the closest hit itself gives 16 and 13; this
        property counts the others, as within_margin does.) Anything past 40 is the degenerate kind.
        """
        count = len(self.within_margin)
        return count if count > MAX_CLOSE_MATCHES else 0

    @property
    def close_matches(self):
        """
        The other references within CLOSE_MATCH_MARGIN of the closest one, closest first.

        Empty when the sample was not identified, because a read set that shares no hashes sits at
        distance 1 from every reference and would tie with all 431 others; and empty when the margin is
        not discriminating, which `undiscriminating` explains. Either way the report says so in words
        instead of printing a reference table, and brucellaid_mash.tsv still holds every distance.
        """
        if not self.identified or self.undiscriminating:
            return []
        return self.within_margin


def coverage_text(coverage):
    """The estimated coverage as the log and the reports write it: never a "0x" next to a real base count."""
    if coverage >= 10:
        return '{:.0f}x'.format(coverage)
    if coverage >= 0.1:
        return '{:.1f}x'.format(coverage)
    return 'under 0.1x'


def kmer_copies_note(identification):
    """How the k-mer copy threshold of one identification was arrived at, for the reports."""
    if identification.coverage is None:
        return '{} (given)'.format(identification.min_copies)
    return '{} (from the estimated {} coverage)'.format(identification.min_copies,
                                                        coverage_text(identification.coverage))


def identify(files, sketch, min_hashes=DEFAULT_MIN_HASHES, threads=1, min_copies=KMER_COPIES_AUTO,
             sample=None):
    """
    Compare one sample to the reference genomes.

    :param min_copies: k-mer copies the read set must hold for a k-mer to be sketched. The default,
                       KMER_COPIES_AUTO, derives it from the coverage of the sample: see min_kmer_copies().
    :param sample: name to prefix the coverage log line with, so that parallel samples stay apart
    :return: Identification
    """
    coverage = None
    if not min_copies:
        bases = estimate_bases(files)
        coverage = bases / BRUCELLA_GENOME_SIZE if bases else None
        min_copies = min_kmer_copies(bases)
        log.info('%s %s bases, roughly %s coverage: sketching the k-mers seen at least %d time(s)',
                 '{}: about'.format(sample) if sample else 'About',
                 '{:,}'.format(bases) if bases else 'an unknown number of',
                 coverage_text(coverage) if coverage else 'unknown',
                 min_copies)
    return Identification(run_dist(files, sketch, threads=threads, min_copies=min_copies),
                          min_hashes=min_hashes, min_copies=min_copies, coverage=coverage)
