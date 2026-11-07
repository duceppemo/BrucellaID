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
    # Tabs and line breaks would break the table. Surrounding spaces are dropped because report.one_cell()
    # strips the TSV cell while the JSON and the MultiQC file keep it, so the four output files would
    # disagree about what the sample is called; a name that is nothing but whitespace is no name at all,
    # and --sample already refuses one.
    name = re.sub(r'[\t\r\n]', '_', name).strip()
    return name or Path(path).name.strip()


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
            path = os.path.join(root, filename)
            if os.path.exists(path) and not os.path.isfile(path):
                # A FIFO blocks the whole run on open() with no writer, and a symbolic link to a
                # device reads without end. single_sample() already refuses anything that is not a
                # regular file; a folder of samples has to do the same, and name what it skipped.
                # exists(), not isfile() alone: a broken symbolic link does not exist, and it is a
                # missing input rather than a strange one, so it falls through to be reported as a
                # sample that failed. Dropping it here lost a named sample from a signed report.
                log.warning('not a regular file, skipped: %s', path)
                continue
            if file_type == 'fasta':
                assemblies.append(path)
                continue
            mate = None
            match = READ_SUFFIX.search(stem)
            if match and match.start() > 0:
                mate, stem = int(match.group(1)), stem[:match.start()]
            # The key is the stem as the file names give it, with only tabs and line breaks taken
            # out: it is what pairs R1 with R2, so two stems that differ must stay two keys. The
            # display name is normalised below instead.
            found.setdefault(re.sub(r'[\t\r\n]', '_', stem), []).append((mate, path))

    # One display name per group of files, and two groups that would answer to the same name are a
    # collision rather than a merge. Surrounding spaces go, because report.one_cell() strips the TSV
    # cell while the JSON and the MultiQC file keep them, and a name of nothing but whitespace falls
    # back to the folder, which is the same for both mates of a pair where the stem itself is not.
    by_name = {}
    for stem, files in sorted(found.items()):
        # The folder name needs the same cleaning as the stem, or it reintroduces the disagreement
        # between the TSV and the JSON that this fallback exists to prevent.
        folder_name = re.sub(r'[\t\r\n]', '_', Path(files[0][1]).parent.name).strip()
        name = stem.strip() or folder_name or 'sample'
        by_name.setdefault(name, []).append(files)

    samples, errors = [], []
    for name, groups in sorted(by_name.items()):
        files = [one for group in groups for one in group]
        if len(groups) > 1:
            errors.append('"{}": {}'.format(name, ', '.join(f[1] for f in files)))
        elif len(files) == 1:
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
