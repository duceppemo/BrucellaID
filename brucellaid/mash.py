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
MIN_KMER_COPIES = 2  # Ignore the k-mers seen once in the read set: most of them are sequencing errors
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


def write_reads(files, stream):
    """
    Copy the reads of one sample into stream, decompressing the gzipped files.

    The files of a sample are piped into one mash sketch: mash sketches each of its input files on its own, and
    the R1 and R2 files of a sample are one read set.
    """
    for path in files:
        with open(path, 'rb') as f:
            gzipped = f.read(2) == b'\x1f\x8b'
        with (gzip.open(path, 'rb') if gzipped else open(path, 'rb')) as reads:
            shutil.copyfileobj(reads, stream, READ_CHUNK)


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


def identify(files, sketch, min_hashes=DEFAULT_MIN_HASHES, threads=1):
    """
    Compare one sample to the reference genomes.

    :return: Identification
    """
    return Identification(run_dist(files, sketch, threads=threads), min_hashes=min_hashes)
