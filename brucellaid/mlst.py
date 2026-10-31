"""Call the nine MLST alleles of a read set with Seal, and look up the sequence type."""

import logging
from dataclasses import dataclass, field
from functools import cache

from .errors import BrucellaIDError
from .references import data_path
from .seal import DEFAULT_KMER_SIZE, DEFAULT_MISMATCHES, run_seal

log = logging.getLogger(__name__)

ALLELES_FILE = 'mlst_alleles.fasta'
PROFILES_FILE = 'mlst_profiles.tsv'
# Reads behind an allele call for it to be kept
DEFAULT_MIN_READS = 5
UNKNOWN = 'ND'  # The value the first version of BrucellaID wrote for a sequence type it could not call
NEW_PROFILE = 'new'


class MlstError(BrucellaIDError):
    pass


def parse_allele(name):
    """
    (locus, allele number) of an allele name: "gap_2" -> ("gap", "2"), "int_hyp_12" -> ("int_hyp", "12").

    :return: (locus, allele), or None when the name is not "<locus>_<number>"
    """
    locus, _, allele = name.rpartition('_')
    if not locus or not allele.isdigit():
        return None
    return locus, allele


def parse_profiles(lines):
    """
    Read the profile table: a header of "ST" and the loci, then one row per sequence type.

    ST	gap	aroA	glk	dnaK	gyrB	trpE	cobQ	int_hyp	omp25
    1	2	1	1	2	1	3	1	1	1

    :return: (loci, {(allele, ...): ST}) with the alleles in the order of the loci
    """
    loci, profiles = None, {}
    for number, line in enumerate(lines, start=1):
        fields = line.rstrip('\n').rstrip('\t').split('\t')
        if not line.strip() or line.startswith('#'):
            continue
        if loci is None:
            if fields[0] != 'ST' or len(fields) < 2:
                raise MlstError('The profile table must start with a header of "ST" and the loci, not "{}".'.format(
                    line.strip()))
            loci = tuple(fields[1:])
            continue
        if len(fields) != len(loci) + 1:
            raise MlstError('Line {} of the profile table has {} columns, and the header {}.'.format(
                number, len(fields), len(loci) + 1))
        key = tuple(fields[1:])
        if key in profiles:
            raise MlstError('ST {} and ST {} of the profile table have the same alleles.'.format(
                profiles[key], fields[0]))
        profiles[key] = fields[0]
    if loci is None or not profiles:
        raise MlstError('The profile table holds no sequence type.')
    return loci, profiles


@cache
def load_profiles(path=None):
    """
    The profile table, read once.

    :param path: a profile table of your own, instead of the one of the package
    :return: (loci, {(allele, ...): ST})
    """
    path = path or data_path(PROFILES_FILE)
    try:
        with open(path) as f:
            return parse_profiles(f)
    except OSError as e:
        raise MlstError('The profile table {} cannot be read: {}'.format(path, e)) from None


def alleles_fasta(path=None):
    """Path of the fasta file of the MLST alleles."""
    return path or data_path(ALLELES_FILE)


@dataclass
class AlleleCall:
    """The allele called at one locus, and what the call rests on."""
    locus: str
    allele: str | None = None  # None when no allele could be called
    reads: int = 0  # Reads whose best match among the alleles of the locus is the called allele
    runner_up: str | None = None  # The next allele of the locus, which the call must beat
    runner_up_reads: int = 0
    note: str = ''

    @property
    def name(self):
        """The allele as it is named in the fasta file: "gap_2"."""
        return '{}_{}'.format(self.locus, self.allele) if self.allele else ''


@dataclass
class MlstResult:
    """The sequence type of one sample, and the allele calls it comes from."""
    st: str = UNKNOWN
    calls: list = field(default_factory=list)
    assigned_reads: int = 0  # Reads matching an allele of any locus
    total_reads: int | None = None  # Reads in the input
    note: str = ''

    @property
    def alleles(self):
        """{locus: allele or None}, in the order of the loci."""
        return {call.locus: call.allele for call in self.calls}

    @property
    def profile(self):
        """The allele numbers in the order of the loci: "2-1-1-2-1-3-1-1-1", with "-" for a locus not called."""
        return '-'.join(call.allele or '-' for call in self.calls)

    @property
    def complete(self):
        return bool(self.calls) and all(call.allele for call in self.calls)


def counts_by_locus(counts, loci):
    """
    Group Seal's per-allele read counts by locus.

    :param counts: {allele name: reads}, as Seal reports them
    :return: {locus: {allele: reads}} for the loci of the scheme
    """
    by_locus = {locus: {} for locus in loci}
    unknown = set()
    for name, reads in counts.items():
        parsed = parse_allele(name)
        if parsed is None or parsed[0] not in by_locus:
            unknown.add(name)
            continue
        locus, allele = parsed
        by_locus[locus][allele] = reads
    if unknown:
        raise MlstError('The allele fasta file holds {} sequence(s) that do not belong to a locus of the profile '
                        'table ({}): {}'.format(len(unknown), ', '.join(loci), ', '.join(sorted(unknown)[:5])))
    empty = [locus for locus, alleles in by_locus.items() if not alleles]
    if empty:
        raise MlstError('The allele fasta file holds no allele for the locus/loci {} of the profile '
                        'table.'.format(', '.join(empty)))
    return by_locus


def call_locus(locus, alleles, min_reads=DEFAULT_MIN_READS):
    """
    The allele of one locus: the one the most reads match best.

    Seal is run with "ambiguous=all" and "clearzone=0", so a read is counted for every allele it shares the
    highest number of k-mers with, and for no other. The allele the sample really carries is therefore counted
    for every read of the locus, and another allele only for the reads that fall where the two agree. The
    allele with the highest count wins, and a tie is reported rather than broken: at a low depth, two alleles
    can agree everywhere the reads happen to cover.

    :param alleles: {allele: reads} for this locus
    :return: AlleleCall
    """
    ranked = sorted(alleles.items(), key=lambda item: (-item[1], int(item[0])))
    best, best_reads = ranked[0]
    second, second_reads = ranked[1] if len(ranked) > 1 else (None, 0)
    call = AlleleCall(locus, runner_up=second, runner_up_reads=second_reads)
    if best_reads < min_reads:
        call.note = ('no allele reached {} read(s) (best: {}_{} with {})'.format(min_reads, locus, best, best_reads)
                     if best_reads else 'no read matched this locus')
        return call
    if second_reads == best_reads:
        tied = [allele for allele, reads in ranked if reads == best_reads]
        call.note = 'tied at {} read(s) between {}'.format(
            best_reads, ', '.join('{}_{}'.format(locus, allele) for allele in tied))
        return call
    call.allele, call.reads = best, best_reads
    return call


def call_alleles(counts, loci, min_reads=DEFAULT_MIN_READS):
    """The allele of each locus of the scheme, in the order of the loci."""
    by_locus = counts_by_locus(counts, loci)
    return [call_locus(locus, by_locus[locus], min_reads=min_reads) for locus in loci]


def lookup_st(calls, profiles):
    """
    The sequence type of a set of allele calls.

    :return: (ST, note): the ST of the profile table, "new" for a profile it does not hold, or "ND" when a
        locus could not be called
    """
    missing = [call.locus for call in calls if not call.allele]
    if missing:
        return UNKNOWN, 'no allele at {}'.format(', '.join(missing))
    st = profiles.get(tuple(call.allele for call in calls))
    if st is None:
        return NEW_PROFILE, 'this allele profile is not in the profile table'
    return st, ''


def run_mlst(files, alleles=None, profiles=None, threads=1, memory='1g', kmer_size=DEFAULT_KMER_SIZE,
             mismatch=DEFAULT_MISMATCHES, min_reads=DEFAULT_MIN_READS):
    """
    Call the MLST alleles of one sample and look up its sequence type.

    :param files: one fastq file, or two paired-end fastq files
    :param alleles: fasta file of the alleles, instead of the one of the package
    :param profiles: profile table, instead of the one of the package
    :return: MlstResult
    """
    loci, table = load_profiles(profiles)
    stats = run_seal(files, alleles_fasta(alleles), threads=threads, memory=memory, k=kmer_size, hdist=mismatch,
                     ambiguous='all', clearzone=0)
    calls = call_alleles(stats.counts, loci, min_reads=min_reads)
    st, note = lookup_st(calls, table)
    return MlstResult(st=st, calls=calls, assigned_reads=stats.matched or 0, total_reads=stats.reads, note=note)
