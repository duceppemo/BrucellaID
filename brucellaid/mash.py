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
# How many copies of a k-mer the read set must hold before MASH sketches it. The point is to drop the
# k-mers that only exist because of a sequencing error, and the right number depends on the coverage: an
# error k-mer appears about depth * error_rate times, so at 30x a single copy is almost certainly an error
# and two copies almost certainly real, while at 400x erroneous k-mers reach two copies routinely and flood
# the sketch. Measured on 225 real runs plus one library strided from 5x to 92x, a fixed 2 starts losing the
# true reference above about 240x (the shared hashes fall below the --min-hashes threshold and nothing is
# reported) and a fixed 10 destroys a 5x sketch, so the value is derived from the coverage instead of being
# a constant. One extra required copy per 20x fits every one of those measurements.
MIN_KMER_COPIES = 2             # Floor, and the right value below about 50x
MAX_KMER_COPIES = 20            # Ceiling: beyond about 400x the sketch is already stable
DEPTH_PER_KMER_COPY = 20        # One more required copy per this much coverage
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
    def close_matches(self):
        """The other references within CLOSE_MATCH_MARGIN of the closest one, closest first."""
        return [hit for hit in self.hits[1:] if hit.distance - self.best.distance <= CLOSE_MATCH_MARGIN]


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
