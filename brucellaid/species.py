"""
Decide the species of a sample: what the MASH distances and the sequence type together say.

MASH compares the whole read set to the reference genomes and is right about the species nearly always, but
not always: some Brucella species are closer to each other than the variation within either of them. In the
sketch of the package, B. canis ATCC 23365 and B. suis bv. 4 strain 40 are 0.00024 apart -- the closest pair
of the sixteen, and closer than two B. suis references are to each other -- so which of the two a canis or a
suis bv. 4 read set lands on first is noise. Every other pair of species in the sketch is at least 0.0007
apart, so MASH decides those on its own.

The sequence type decides the pair MASH cannot. The nine-locus scheme gives B. canis ST 20, 21 and 106 and
gives B. suis fourteen others, and not one ST is shared between them: an ST is a different kind of evidence,
not a sharper version of the same one.

So the rule is: MASH names the species, and where MASH cannot choose, the ST chooses for it among the
species MASH put within reach. The ST is never allowed to overrule a species MASH picked cleanly -- a hit on
nine loci is thinner evidence than a hit on a thousand min-hashes -- but a disagreement is reported, because
it means one of the two is looking at something unexpected.

The ST speaks only when its isolates agree on a species: the labels of mlst_st_species.tsv are the
submitters' and some of them are wrong, which is the same problem that showed up in the NCBI assemblies. An
ST whose isolates disagree, or that rests on one or two of them, is left out of the decision rather than
settled by majority.
"""

import logging
from dataclasses import dataclass
from functools import cache

from .errors import BrucellaIDError
from .mlst import UNKNOWN
from .references import data_path, lookup

log = logging.getLogger(__name__)

ST_SPECIES_FILE = 'mlst_st_species.tsv'
COLUMNS = ('ST', 'Species', 'Purity', 'Isolates', 'Biovars', 'OtherSpecies')
# Share of the isolates of an ST that must carry the same species for the ST to be used, and how many
# isolates it must rest on. 0.95 keeps the STs whose only dissenters are one or two isolates out of dozens --
# a wrong label rather than a species boundary -- and drops the few that are really split.
DEFAULT_MIN_PURITY = 0.95
DEFAULT_MIN_ISOLATES = 3
# The species field of an isolate whose submitter did not name one. Nineteen STs of the table are this and
# nothing else, six of them with enough isolates to look reliable, so it has to be read as the absence of a
# species rather than as one: an ST of unnamed isolates cannot confirm, choose or contradict anything.
UNNAMED = 'Brucella sp.'

# What the Species column rests on.
FROM_MASH = 'MASH'
FROM_BOTH = 'MASH+ST'
FROM_ST = 'ST'

# What the two kinds of evidence did, reported in the SpeciesSource column's company and in the notes.
AGREE = 'agree'  # Both name the same species, and MASH could have named it alone
ARBITRATED = 'arbitrated'  # MASH left several species within reach and the ST chose between them
CONFLICT = 'conflict'  # The two name species that are not even close to each other
UNRESOLVED = 'unresolved'  # MASH left several species within reach and the ST could not choose


class SpeciesError(BrucellaIDError):
    pass


@dataclass(frozen=True)
class StSpecies:
    """One row of the ST to species table: the species the isolates of an ST carry."""
    st: str
    species: str
    purity: float  # Share of the isolates of this ST that carry this species, 1.0 when they all do
    isolates: int  # Isolates of this ST in the PubMLST isolate database
    biovars: str = ''  # The biovars recorded for this species, "1:86;2:16", or ""
    other: str = ''  # The other species of this ST, "Brucella melitensis:3", or ""

    @property
    def agreeing(self):
        """How many isolates of the ST carry the species."""
        return int(round(self.purity * self.isolates))

    def reliable(self, min_purity=DEFAULT_MIN_PURITY, min_isolates=DEFAULT_MIN_ISOLATES):
        """Whether the ST rests on enough isolates, agreeing enough, to be used in a decision."""
        return self.isolates >= min_isolates and self.purity >= min_purity

    @property
    def evidence(self):
        """The row as a note reads it: "ST 20 is Brucella canis in 70 of 70 isolate(s)"."""
        text = 'ST {} is {} in {} of {} isolate(s)'.format(self.st, self.species, self.agreeing, self.isolates)
        return '{} (the rest: {})'.format(text, self.other.replace(';', ', ')) if self.other else text


def parse_st_species(lines):
    """
    Read the ST to species table: comments, a header, then one tab separated row per ST.

    ST	Species	Purity	Isolates	Biovars	OtherSpecies
    20	Brucella canis	1.000	70

    Trailing empty columns may be missing from a row: a trailing tab is whitespace, which the repository's
    pre-commit hooks strip.

    :return: {ST: StSpecies}
    """
    table = {}
    header = None
    for number, line in enumerate(lines, start=1):
        line = line.rstrip('\n')
        if not line.strip() or line.startswith('#'):
            continue
        fields = line.split('\t')
        if header is None:
            header = fields
            if tuple(header) != COLUMNS:
                raise SpeciesError('The ST to species table must start with the columns {}, not {}.'.format(
                    ', '.join(COLUMNS), ', '.join(header)))
            continue
        if len(fields) > len(COLUMNS):
            raise SpeciesError('Line {} of the ST to species table has {} columns, more than the {} of the '
                               'header.'.format(number, len(fields), len(COLUMNS)))
        fields += [''] * (len(COLUMNS) - len(fields))
        st, species, purity, isolates, biovars, other = (f.strip() for f in fields)
        if not st or not species:
            raise SpeciesError('Line {} of the ST to species table has no ST or no species.'.format(number))
        try:
            purity, isolates = float(purity), int(isolates)
        except ValueError:
            raise SpeciesError('Line {} of the ST to species table has a purity or an isolate count that is not a '
                               'number: "{}", "{}".'.format(number, purity, isolates)) from None
        if not 0.0 < purity <= 1.0 or isolates < 1:
            raise SpeciesError('Line {} of the ST to species table has a purity outside 0-1 or no isolate: {}, '
                               '{}.'.format(number, purity, isolates))
        if st in table:
            raise SpeciesError('ST {} appears twice in the ST to species table.'.format(st))
        table[st] = StSpecies(st, species, purity, isolates, biovars, other)
    if not table:
        raise SpeciesError('The ST to species table holds no ST.')
    return table


@cache
def load_st_species(path=None):
    """
    The ST to species table, read once.

    :param path: a table of your own, instead of the one of the package
    :return: {ST: StSpecies}
    """
    path = path or data_path(ST_SPECIES_FILE)
    try:
        with open(path) as f:
            return parse_st_species(f)
    except OSError as e:
        raise SpeciesError('The ST to species table {} cannot be read: {}'.format(path, e)) from None


def candidates_of(identification, references):
    """
    The species MASH put within reach: the species of the closest reference genome, and of the ones nearly as
    close.

    More than one means MASH cannot choose. The genomes nearly as close are the ones within
    mash.CLOSE_MATCH_MARGIN of the closest, which is wide enough to hold the canis/suis bv. 4 pair and
    nothing else of the sketch of the package.

    :return: tuple of species, the closest first, without repetition
    """
    names = []
    for hit in [identification.best, *identification.close_matches]:
        reference = lookup(references, hit.reference)
        if reference is not None and reference.species not in names:
            names.append(reference.species)
    return tuple(names)


def closest_of(identification, references, species):
    """
    The closest reference genome of one species, which is where its biovar and its vSNP code come from.

    When the ST moves the call to another species, the biovar of the genome MASH landed on belongs to the
    species it named, not to the one being reported, so it cannot be carried over.

    :return: references.Reference, or None when the sketch holds no genome of that species
    """
    for hit in identification.hits:
        reference = lookup(references, hit.reference)
        if reference is not None and reference.species == species:
            return reference
    return None


@dataclass(frozen=True)
class SpeciesCall:
    """The species reported for one sample, and what it rests on."""
    species: str = UNKNOWN
    reference: object = None  # references.Reference the biovar and the vSNP code are taken from, or None
    source: str = ''  # FROM_MASH, FROM_BOTH, FROM_ST, or "" when no species could be called
    verdict: str = ''  # AGREE, ARBITRATED, CONFLICT, UNRESOLVED, or "" when only MASH had anything to say
    candidates: tuple = ()  # The species MASH put within reach, the closest first
    st: str = ''  # The sequence type the ST evidence came from, or ""
    entry: object = None  # The StSpecies row used, or None when the ST said nothing usable
    st_reason: str = ''  # Why the ST was not used, when it was not: a report has to be able to say
    notes: tuple = ()

    @property
    def biovar(self):
        return self.reference.biovar if self.reference else ''

    @property
    def code(self):
        """The vSNP species code, written to the "species" column of the legacy xlsx report."""
        return self.reference.code if self.reference else UNKNOWN

    @property
    def st_species(self):
        """The species the ST stands for, or "" when the ST was not usable."""
        return self.entry.species if self.entry else ''


def st_evidence(mlst, table, min_purity=DEFAULT_MIN_PURITY, min_isolates=DEFAULT_MIN_ISOLATES):
    """
    What the sequence type of a sample says about its species.

    :param mlst: mlst.MlstResult, or None when the MLST was not run
    :param table: {ST: StSpecies}, or None when the table is not to be used
    :return: (StSpecies or None, why it cannot be used); the reason is "" when it can
    """
    if table is None:
        return None, 'the ST to species table is not in use'
    if mlst is None:
        return None, 'the MLST was not run'
    if mlst.st in (UNKNOWN, '') or not mlst.st:
        return None, 'the sample has no ST'
    entry = table.get(mlst.st)
    if entry is None:
        reason = ('this allele profile has no ST yet' if mlst.st == 'new'
                  else 'no isolate of ST {} has a species in the table'.format(mlst.st))
        return None, reason
    if entry.species == UNNAMED:
        return None, 'the isolates of ST {} are recorded as {}, with no species named'.format(entry.st, UNNAMED)
    if entry.isolates < min_isolates:
        return None, 'ST {} rests on {} isolate(s), fewer than the {} required'.format(
            entry.st, entry.isolates, min_isolates)
    if entry.purity < min_purity:
        return None, '{}, below the {:.0%} agreement required'.format(entry.evidence, min_purity)
    return entry, ''


def call_species(identification, references, mlst=None, table=None, min_purity=DEFAULT_MIN_PURITY,
                 min_isolates=DEFAULT_MIN_ISOLATES):
    """
    The species of one sample, from the MASH hits and, where it decides anything, from the sequence type.

    See the module docstring for the rule. In short: MASH names the species; where MASH left more than one
    species within reach, a reliable ST chooses between them; an ST that contradicts a species MASH picked
    cleanly is reported but not obeyed.

    :param identification: mash.Identification, or None when MASH did not run
    :param references: {reference: references.Reference}
    :param mlst: mlst.MlstResult, or None
    :param table: {ST: StSpecies}, or None not to use the ST at all
    :return: SpeciesCall
    """
    entry, why_not = st_evidence(mlst, table, min_purity=min_purity, min_isolates=min_isolates)
    st = mlst.st if mlst else ''
    if identification is None:
        return SpeciesCall(st=st, entry=entry, st_reason=why_not)

    best = lookup(references, identification.best.reference)
    candidates = candidates_of(identification, references) if identification.identified else ()
    mash_species = best.species if identification.identified and best else ''

    # MASH found too little Brucella to name a species. A clean ST is worth repeating -- nine loci were
    # called, so there are Brucella reads -- but not worth calling a species on: the two statements are
    # about different amounts of the genome, and this one is the weaker.
    if not identification.identified:
        notes = ['{}, but the reads are not identified as Brucella'.format(entry.evidence)] if entry else []
        return SpeciesCall(st=st, entry=entry, st_reason=why_not, notes=tuple(notes))

    # A custom sketch whose closest genome the reference table does not name: MASH has no species to offer,
    # so a reliable ST is all there is.
    if not mash_species:
        if entry:
            return SpeciesCall(species=entry.species, source=FROM_ST, verdict=ARBITRATED, st=st, entry=entry,
                               notes=('the closest reference has no species in the reference table, so the species '
                                      'is the ST\'s: {}'.format(entry.evidence),))
        return SpeciesCall(st=st, entry=entry, st_reason=why_not)

    if entry is None:
        # Nothing to arbitrate with. Worth saying so only when MASH actually needed it.
        if len(candidates) > 1:
            return SpeciesCall(species=mash_species, reference=best, source=FROM_MASH, verdict=UNRESOLVED,
                               candidates=candidates, st=st, st_reason=why_not,
                               notes=('MASH cannot separate {}: the closest is {}, and the ST cannot choose '
                                      'because {}'.format(' and '.join(candidates), mash_species, why_not),))
        return SpeciesCall(species=mash_species, reference=best, source=FROM_MASH, candidates=candidates,
                           st=st, st_reason=why_not)

    if entry.species == mash_species:
        # The ST agrees. Where MASH had a choice to make, that agreement is the thing that settles it.
        if len(candidates) > 1:
            return SpeciesCall(species=mash_species, reference=best, source=FROM_BOTH, verdict=ARBITRATED,
                               candidates=candidates, st=st, entry=entry,
                               notes=('MASH cannot separate {}: {}, which agrees with the closest '
                                      'reference'.format(' and '.join(candidates), entry.evidence),))
        return SpeciesCall(species=mash_species, reference=best, source=FROM_BOTH, verdict=AGREE,
                           candidates=candidates, st=st, entry=entry)

    if entry.species in candidates:
        # The two species are inside the distance MASH cannot resolve, and the ST names one of them. This is
        # the canis/suis case the table was built for.
        reference = closest_of(identification, references, entry.species)
        return SpeciesCall(species=entry.species, reference=reference, source=FROM_ST, verdict=ARBITRATED,
                           candidates=candidates, st=st, entry=entry,
                           notes=('MASH cannot separate {}: the closest is {}, and {}, so the species is the '
                                  'ST\'s'.format(' and '.join(candidates), mash_species, entry.evidence),))

    # The ST names a species that is not even among the ones MASH put within reach. MASH keeps the call, on a
    # thousand min-hashes against nine loci, but something is wrong with the sample or with one of the two
    # databases and the report has to say so.
    return SpeciesCall(species=mash_species, reference=best, source=FROM_MASH, verdict=CONFLICT,
                       candidates=candidates, st=st, entry=entry,
                       notes=('the ST disagrees with the MASH identification: the closest reference is {} and {}. '
                              'The species reported is the MASH one; check the sample'.format(
                                  best.description, entry.evidence),))
