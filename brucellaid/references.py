"""The reference genomes of the MASH sketch, and the species they stand for."""

import logging
from dataclasses import dataclass
from functools import cache
from importlib import resources
from pathlib import PurePosixPath

from .errors import BrucellaIDError

log = logging.getLogger(__name__)

REFERENCES_FILE = 'references.tsv'
SKETCH_FILE = 'references.msh'
EXTENDED_REFERENCES_FILE = 'references_extended.tsv'
EXTENDED_SKETCH_FILE = 'references_extended.msh'
COLUMNS = ('reference', 'code', 'species', 'biovar', 'strain', 'accessions')


def data_path(name):
    """Absolute path of a file of the package's data folder."""
    path = resources.files('brucellaid') / 'data' / name
    if not path.is_file():
        raise BrucellaIDError('The data file "{}" is missing from the installed brucellaid package. Reinstall '
                              'it.'.format(name))
    return str(path)


@dataclass(frozen=True)
class Reference:
    """One genome of the MASH sketch."""
    reference: str  # The sketch ID, which is the name of the file the sketch was built from
    code: str  # The vSNP species code, e.g. "ab1"
    species: str  # e.g. "Brucella abortus"
    biovar: str = ''
    strain: str = ''
    accessions: str = ''  # NCBI nucleotide accessions, separated by ";", or ""

    @property
    def label(self):
        """The species and, when there is one, the biovar: "Brucella abortus bv. 1"."""
        return '{} bv. {}'.format(self.species, self.biovar) if self.biovar else self.species

    @property
    def description(self):
        """The species, biovar and strain: "Brucella abortus bv. 1 str. 9-941"."""
        return '{} str. {}'.format(self.label, self.strain) if self.strain else self.label


def parse_references(lines):
    """
    Read the reference table: a commented header, then one tab separated row per reference.

    Trailing empty columns may be missing from a row: a trailing tab is whitespace, which the repository's
    pre-commit hooks strip.

    :return: {reference: Reference}
    """
    references = {}
    header = None
    for number, line in enumerate(lines, start=1):
        line = line.rstrip('\n')
        if not line.strip() or line.startswith('#'):
            continue
        fields = line.split('\t')
        if header is None:
            header = fields
            if tuple(header) != COLUMNS:
                raise BrucellaIDError('The reference table must start with the columns {}, not {}.'.format(
                    ', '.join(COLUMNS), ', '.join(header)))
            continue
        if len(fields) > len(COLUMNS):
            raise BrucellaIDError('Line {} of the reference table has {} columns, more than the {} of the '
                                  'header.'.format(number, len(fields), len(COLUMNS)))
        fields += [''] * (len(COLUMNS) - len(fields))
        reference = Reference(*(f.strip() for f in fields))
        if not reference.reference or not reference.code or not reference.species:
            raise BrucellaIDError('Line {} of the reference table has no reference, code or species.'.format(number))
        if reference.reference in references:
            raise BrucellaIDError('The reference "{}" appears twice in the reference table.'.format(
                reference.reference))
        references[reference.reference] = reference
    if not references:
        raise BrucellaIDError('The reference table holds no reference.')
    return references


@cache
def load_references(path=None):
    """
    The reference table, read once.

    :param path: a reference table of your own, instead of the one of the package
    :return: {reference: Reference}
    """
    path = path or data_path(REFERENCES_FILE)
    try:
        with open(path) as f:
            return parse_references(f)
    except OSError as e:
        raise BrucellaIDError('The reference table {} cannot be read: {}'.format(path, e)) from None


def lookup(references, sketch_id):
    """
    The Reference a MASH sketch ID stands for, or None.

    "mash dist" reports a genome under the name the sketch was built from, which is a path when the sketch was
    built from one ("/data/refs/NC_006932-NC_006933.fasta"), so the file name is tried as well.

    :param references: {reference: Reference}
    :param sketch_id: the ID "mash dist" reported
    """
    return references.get(sketch_id) or references.get(PurePosixPath(sketch_id).name)


def reference_sketch(path=None):
    """Path of the MASH sketch of the reference genomes."""
    return path or data_path(SKETCH_FILE)


def compare_with_sketch(sketch_ids, references):
    """
    What the MASH sketch and the reference table disagree on.

    Neither is an error: a sketch may hold a genome the table does not name (its sample is then reported
    without a species), and the table may name a genome a custom sketch left out.

    :param sketch_ids: IDs of the genomes of the sketch, as "mash info" reports them
    :return: list of messages, empty when the two match
    """
    missing = [name for name in sketch_ids if lookup(references, name) is None]
    named = {lookup(references, name).reference for name in sketch_ids if lookup(references, name) is not None}
    unused = [name for name in references if name not in named]
    messages = []
    if missing:
        messages.append('{} genome(s) of the MASH sketch are not in the reference table, and a sample closest to '
                        'one of them is reported without a species: {}'.format(len(missing), ', '.join(missing)))
    if unused:
        messages.append('{} genome(s) of the reference table are not in the MASH sketch: {}'.format(
            len(unused), ', '.join(unused)))
    return messages
