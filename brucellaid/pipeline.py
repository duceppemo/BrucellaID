"""Run the identification and the typing on each sample."""

import logging
from concurrent import futures
from dataclasses import dataclass, field

from .errors import BrucellaIDError
from .mash import DEFAULT_MIN_HASHES, identify
from .mlst import DEFAULT_MIN_READS, UNKNOWN, run_mlst
from .references import load_references, lookup, reference_sketch
from .seal import DEFAULT_KMER_SIZE, DEFAULT_MISMATCHES
from .species import DEFAULT_MIN_ISOLATES, DEFAULT_MIN_PURITY, call_species, load_st_species

log = logging.getLogger(__name__)


@dataclass
class SampleResult:
    """What one sample gave: the closest reference genome, the sequence type, or an error."""
    sample: str
    files: list
    identification: object = None  # mash.Identification
    reference: object = None  # references.Reference of the closest genome, when it is in the reference table
    mlst: object = None  # mlst.MlstResult
    call: object = None  # species.SpeciesCall: the species reported, and what it rests on
    notes: list = field(default_factory=list)
    error: str | None = None

    @property
    def identified(self):
        return self.identification is not None and self.identification.identified

    @property
    def species(self):
        """The species reported, which the closest reference genome and the ST decided between them, or "ND"."""
        return self.call.species if self.call else UNKNOWN

    @property
    def closest_reference(self):
        """The name of the closest reference genome, as the reference table names it when it is in it."""
        if self.identification is None:
            return ''
        return self.reference.reference if self.reference else self.identification.best.reference

    @property
    def code(self):
        """The vSNP species code, written to the "species" column of the legacy xlsx report."""
        return self.call.code if self.call else UNKNOWN

    @property
    def sequence_type(self):
        return self.mlst.st if self.mlst else UNKNOWN

    @property
    def note(self):
        return '; '.join(note for note in self.notes if note)


def close_label(references, hit):
    """How a nearly as close reference is named in a note: its species, or that it is not in the table."""
    reference = lookup(references, hit.reference)
    return reference.label if reference else 'not in the reference table'


def species_notes(identification, references):
    """
    What is worth saying about an identification: a read set that is not Brucella, and the references that are
    nearly as close as the closest one.

    :return: (Reference of the closest genome or None, list of notes)
    """
    best = identification.best
    reference = lookup(references, best.reference)
    notes = []
    if not identification.identified:
        notes.append('not identified: the closest reference ({}) shares only {}/{} hashes with the reads, fewer '
                     'than the {} required. The sample may not be Brucella, or may hold too few Brucella '
                     'reads'.format(best.reference, best.shared_hashes, best.total_hashes,
                                    identification.min_hashes))
        return reference, notes
    if reference is None:
        notes.append('the closest reference ({}) is not in the reference table, so no species is '
                     'reported'.format(best.reference))
    close = identification.close_matches
    if close:
        notes.append('nearly as close: {}'.format(', '.join(
            '{} ({}, {}/{} hashes)'.format(hit.reference, close_label(references, hit),
                                           hit.shared_hashes, hit.total_hashes) for hit in close)))
    return reference, notes


def analyse_sample(sample, sketch=None, references=None, min_hashes=DEFAULT_MIN_HASHES,
                   min_reads=DEFAULT_MIN_READS, kmer_size=DEFAULT_KMER_SIZE, mismatch=DEFAULT_MISMATCHES,
                   threads=1, memory='1g', alleles=None, profiles=None, mlst=True, st_table=None,
                   use_st_species=True, min_st_purity=DEFAULT_MIN_PURITY,
                   min_st_isolates=DEFAULT_MIN_ISOLATES):
    """
    Identify and type one sample.

    An error on one sample is kept in its SampleResult rather than raised, so that a run over a folder reports
    the samples that worked.

    :param st_table: {ST: species.StSpecies}, instead of the table of the package
    :param use_st_species: whether the sequence type may take part in the species call at all
    """
    references = load_references() if references is None else references
    if use_st_species and st_table is None:
        st_table = load_st_species()
    sketch = reference_sketch(sketch)
    result = SampleResult(sample.name, list(sample.files))
    log.info('%s: identifying the species with MASH (%d file(s))', sample.name, len(sample.files))
    try:
        result.identification = identify(sample.files, sketch, min_hashes=min_hashes, threads=threads)
        result.reference, result.notes = species_notes(result.identification, references)
        best = result.identification.best
        if result.identified:
            log.info('%s: %s (%s, MASH distance %.5f, %d/%d hashes)', sample.name,
                     result.reference.description if result.reference else best.reference, best.reference,
                     best.distance, best.shared_hashes, best.total_hashes)
        else:
            log.warning('%s: %s', sample.name, result.notes[0])
        if mlst:
            log.info('%s: calling the MLST alleles with Seal', sample.name)
            result.mlst = run_mlst(sample.files, alleles=alleles, profiles=profiles, threads=threads,
                                   memory=memory, kmer_size=kmer_size, mismatch=mismatch, min_reads=min_reads)
            log.info('%s: ST %s (profile %s)', sample.name, result.mlst.st, result.mlst.profile)
            for call in result.mlst.calls:
                if call.note:
                    result.notes.append('{}: {}'.format(call.locus, call.note))
            if result.mlst.note:
                result.notes.append('ST {}: {}'.format(result.mlst.st, result.mlst.note))
        # The species is decided once both are in: where MASH left more than one species within reach, the
        # ST is what chooses between them.
        result.call = call_species(result.identification, references, mlst=result.mlst,
                                   table=st_table if use_st_species else None, min_purity=min_st_purity,
                                   min_isolates=min_st_isolates)
        result.notes.extend(result.call.notes)
        if result.call.verdict:
            log.info('%s: %s (%s, %s)', sample.name, result.call.species, result.call.source,
                     result.call.verdict)
    except BrucellaIDError as e:
        log.error('%s: %s', sample.name, e)
        result.error = str(e)
        if result.call is None and result.identification is not None:
            # MASH had already answered when the typing failed: report the species it gave.
            result.call = call_species(result.identification, references)
    return result


def analyse_samples(samples, jobs=1, threads=1, **options):
    """
    Identify and type several samples, up to jobs samples at a time.

    The threads are shared between the jobs, as MASH and Seal use them within one sample.

    :return: list of SampleResult, in the order of samples
    """
    jobs = max(1, min(jobs, len(samples)))
    per_job = max(1, threads // jobs)
    if jobs > 1:
        log.info('%d sample(s), %d at a time, %d thread(s) each', len(samples), jobs, per_job)
    with futures.ThreadPoolExecutor(max_workers=jobs) as executor:
        return list(executor.map(lambda s: analyse_sample(s, threads=per_job, **options), samples))
