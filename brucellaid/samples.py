"""Find the samples in an input folder: single-end or paired-end fastq files."""

import logging
import os
import re
from dataclasses import dataclass
from pathlib import Path

from .errors import BrucellaIDError

log = logging.getLogger(__name__)

FASTQ_EXTENSIONS = ('.fastq', '.fq')
FASTA_EXTENSIONS = ('.fasta', '.fa', '.fna', '.fas', '.fsa')
# Read suffixes of fastq file names, e.g. "S1_R1", "S1_S12_L001_R2_001", "SRR123_1". Group 1 is the mate number.
READ_SUFFIX = re.compile(r'_R?([12])(?:_001)?$')

ASSEMBLY_MESSAGE = ('BrucellaID reads fastq files, not assemblies: both the MASH sketch and the MLST allele '
                    'calling count reads. Give the reads the assembly was built from.')


@dataclass
class Sample:
    name: str
    files: list  # One file, or R1 and R2

    @property
    def paired(self):
        return len(self.files) == 2


def split_extension(filename):
    """(stem, "fastq" or "fasta") for a sequence file, gzipped or not, or None for any other file."""
    lower = filename.lower()
    if lower.endswith('.gz'):
        lower = lower[:-3]
    for file_type, extensions in (('fastq', FASTQ_EXTENSIONS), ('fasta', FASTA_EXTENSIONS)):
        for extension in extensions:
            if lower.endswith(extension) and len(lower) > len(extension):
                return filename[:len(lower) - len(extension)], file_type
    return None


def sample_name(path):
    """
    Sample name from a fastq file name: the extension and the read suffix are removed.
    "S1_R1.fastq.gz" -> "S1", "SRR123_1.fq" -> "SRR123".
    """
    name = Path(path).name
    parsed = split_extension(name)
    if parsed is not None:
        name = parsed[0]
    name = READ_SUFFIX.sub('', name) or name
    return re.sub(r'[\t\r\n]', '_', name)  # Tabs and line breaks would break the table


def find_samples(folder, exclude=()):
    """
    Look for fastq files in a folder and its subfolders, and group R1 and R2 files per sample.

    Name collisions are errors rather than guesses: two files giving the same sample name would otherwise be
    merged or overwritten silently.

    :param exclude: folders to skip, e.g. the output folder when it is inside the input folder
    :return: list of Sample, sorted by name
    """
    if not Path(folder).is_dir():
        raise BrucellaIDError('Input folder not found: {}'.format(folder))
    exclude = {Path(p).resolve() for p in exclude}
    found = {}  # {sample name: [(mate number or None, path)]}
    assemblies = []
    visited = set()  # Real paths of the folders searched: symbolic links to folders are followed, loops are not
    for root, dirs, filenames in os.walk(folder, followlinks=True):
        real = Path(root).resolve()
        if real in visited:
            dirs[:] = []
            continue
        visited.add(real)
        dirs[:] = sorted(d for d in dirs if Path(root, d).resolve() not in exclude and not d.startswith('.'))
        for filename in sorted(filenames):
            parsed = split_extension(filename)
            if parsed is None or filename.startswith('.'):
                continue
            stem, file_type = parsed
            if file_type == 'fasta':
                assemblies.append(os.path.join(root, filename))
                continue
            mate = None
            match = READ_SUFFIX.search(stem)
            if match and match.start() > 0:
                mate, stem = int(match.group(1)), stem[:match.start()]
            stem = re.sub(r'[\t\r\n]', '_', stem)  # Tabs and line breaks would break the table
            found.setdefault(stem, []).append((mate, os.path.join(root, filename)))

    samples, errors = [], []
    for name, files in sorted(found.items()):
        if len(files) == 1:
            mate, path = files[0]
            if mate == 2:
                log.warning('%s: R2 file without R1, used as single-end reads: %s', name, path)
            samples.append(Sample(name, [path]))
        elif len(files) == 2 and {f[0] for f in files} == {1, 2}:
            samples.append(Sample(name, [f[1] for f in sorted(files, key=lambda f: f[0])]))
        else:
            errors.append('"{}": {}'.format(name, ', '.join(f[1] for f in files)))
    if errors:
        raise BrucellaIDError('Several files give the same sample name. Rename, move or remove them:\n  '
                              + '\n  '.join(errors))
    if assemblies:
        log.warning('%d fasta file(s) in %s were skipped. %s', len(assemblies), folder, ASSEMBLY_MESSAGE)
    if not samples:
        raise BrucellaIDError('No fastq file found in {}.{}'.format(
            folder, ' ' + ASSEMBLY_MESSAGE if assemblies else ''))
    return samples


def single_sample(r1, r2=None, name=None):
    """The Sample for one input given on the command line, with an optional sample name."""
    for path in (r1, r2):
        if path is not None and not Path(path).is_file():
            raise BrucellaIDError('Input file not found: {}'.format(path))
    for path in (r1, r2):
        parsed = split_extension(Path(path).name) if path is not None else None
        if parsed is not None and parsed[1] == 'fasta':
            raise BrucellaIDError('"{}" is a fasta file. {}'.format(path, ASSEMBLY_MESSAGE))
    files = [os.path.abspath(r1)] + ([os.path.abspath(r2)] if r2 else [])
    if len(files) == 2 and files[0] == files[1]:
        raise BrucellaIDError('-r1 and -r2 are the same file: {}'.format(files[0]))
    return Sample(name or sample_name(r1), files)
