"""Fixtures: synthetic reads, a MASH sketch of a synthetic genome, and the MLST data of the package."""

import gzip
import random
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from brucellaid import mash, references, seal
from brucellaid.htmlreport import logo
from brucellaid.mlst import alleles_fasta, load_profiles
from brucellaid.species import load_st_species

QUALITY = 'I'  # Phred 40, so that no tool trims the synthetic reads


def tool(name):
    """Path of an external tool, in the PATH or next to the running Python."""
    found = shutil.which(name)
    if found:
        return found
    candidate = Path(sys.executable).parent / name
    return str(candidate) if candidate.exists() else None


requires_mash = pytest.mark.skipif(tool('mash') is None, reason='MASH is not installed')
requires_seal = pytest.mark.skipif(tool('seal.sh') is None, reason='Seal (BBTools) is not installed')


def read_fasta(path):
    """{name: sequence} of a fasta file."""
    records, name = {}, None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                name = line[1:].split()[0]
                records[name] = []
            elif line and name:
                records[name].append(line.upper())
    return {name: ''.join(parts) for name, parts in records.items()}


def write_fastq(path, reads, gzipped=None):
    """Write (name, sequence) pairs as a fastq file, gzipped when the name ends with ".gz"."""
    path = Path(path)
    gzipped = path.name.endswith('.gz') if gzipped is None else gzipped
    opener = gzip.open if gzipped else open
    with opener(path, 'wt') as f:
        for name, sequence in reads:
            f.write('@{}\n{}\n+\n{}\n'.format(name, sequence, QUALITY * len(sequence)))
    return path


def tile(sequence, length=150, step=25, prefix='read'):
    """
    Reads covering a sequence from end to end.

    The last window is always included: two alleles of a locus sometimes differ at one position only, and a
    position no read covers cannot tell them apart.
    """
    starts = list(range(0, max(1, len(sequence) - length + 1), step))
    last = max(0, len(sequence) - length)
    if last not in starts:
        starts.append(last)
    return [('{}_{}'.format(prefix, start), sequence[start:start + length]) for start in starts]


def random_sequence(length, seed):
    generator = random.Random(seed)
    return ''.join(generator.choice('ACGT') for _ in range(length))


def random_reads(sequence, count, length=150, seed=1, prefix='read'):
    """Reads taken at random positions of a sequence, on either strand."""
    generator = random.Random(seed)
    complement = str.maketrans('ACGT', 'TGCA')
    reads = []
    for i in range(count):
        start = generator.randrange(0, len(sequence) - length)
        read = sequence[start:start + length]
        if generator.random() < 0.5:
            read = read.translate(complement)[::-1]
        reads.append(('{}_{}'.format(prefix, i), read))
    return reads


@pytest.fixture
def mlst_alleles():
    """{allele name: sequence} of the MLST alleles of the package."""
    return read_fasta(alleles_fasta())


@pytest.fixture
def mlst_loci():
    """The nine loci of the scheme, in the order of the profile table."""
    return load_profiles()[0]


@pytest.fixture
def mlst_profiles():
    """{(allele, ...): ST} of the profile table of the package."""
    return load_profiles()[1]


@pytest.fixture
def st_reads(tmp_path, mlst_alleles, mlst_loci, mlst_profiles):
    """
    A function writing a fastq file of reads covering the alleles of one sequence type.

    reads = st_reads('1') -> a fastq file whose sample must be typed as ST 1.
    """
    by_st = {st: alleles for alleles, st in mlst_profiles.items()}

    def make(st='1', name=None, step=25, length=150, gzipped=True, folder=None):
        reads = []
        for locus, allele in zip(mlst_loci, by_st[st], strict=True):
            sequence = mlst_alleles['{}_{}'.format(locus, allele)]
            reads += tile(sequence, length=length, step=step, prefix='{}_{}'.format(locus, allele))
        folder = Path(folder or tmp_path)
        folder.mkdir(parents=True, exist_ok=True)
        path = folder / '{}_R1.fastq{}'.format(name or 'ST{}'.format(st), '.gz' if gzipped else '')
        return write_fastq(path, reads, gzipped=gzipped)

    return make


@pytest.fixture
def genome(tmp_path):
    """A synthetic genome, standing in for a reference: a 60 kb sequence written as a fasta file."""
    sequence = random_sequence(60000, seed=7)
    path = tmp_path / 'test_genome.fasta'
    with open(path, 'w') as f:
        f.write('>chromosome test genome\n')
        for i in range(0, len(sequence), 70):
            f.write(sequence[i:i + 70] + '\n')
    return path, sequence


@pytest.fixture
def sketch(tmp_path, genome):
    """A MASH sketch of the synthetic genome, and a reference table naming it."""
    path, _ = genome
    prefix = tmp_path / 'test_references'
    subprocess.run([tool('mash'), 'sketch', '-o', str(prefix), str(path)], check=True, capture_output=True)
    table = tmp_path / 'test_references.tsv'
    table.write_text('\t'.join(references.COLUMNS) + '\n'
                     + '\t'.join([path.name, 'test', 'Brucella testii', '9', 'T1', 'NC_000000.1']) + '\n')
    return prefix.with_suffix('.msh'), table


@pytest.fixture
def genome_reads(tmp_path, genome):
    """Reads of the synthetic genome, deep enough for MASH to recognize it."""
    _, sequence = genome
    return write_fastq(tmp_path / 'genome_R1.fastq.gz', random_reads(sequence, 6000, seed=3))


@pytest.fixture(autouse=True)
def clear_caches():
    """The cached readers and tool lookups must not leak between tests."""
    cached_readers = (mash.executable, mash.sketch_references, seal.executable, references.load_references,
                      load_profiles, load_st_species, logo)
    for cached in cached_readers:
        cached.cache_clear()
    yield
    for cached in cached_readers:
        cached.cache_clear()
