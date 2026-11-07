"""What the pipeline makes of the MASH and MLST results of a sample."""

import pytest

from brucellaid import pipeline
from brucellaid.errors import BrucellaIDError
from brucellaid.mash import Identification, parse_dist
from brucellaid.mlst import UNKNOWN, AlleleCall, MlstResult
from brucellaid.pipeline import SampleResult, analyse_sample, analyse_samples, close_label, species_notes
from brucellaid.references import COLUMNS, parse_references
from brucellaid.samples import Sample

REFERENCES = parse_references([
    '\t'.join(COLUMNS) + '\n',
    'NC_006932-NC_006933.fasta\tab1\tBrucella abortus\t1\t9-941\tNC_006932.1\n',
    'NZ_CP007682-NZ_CP007683.fasta\tab3\tBrucella abortus\t3\tBER\tNZ_CP007682.1\n',
])


def identification(*rows, min_hashes=150):
    return Identification(parse_dist(''.join(rows)), min_hashes=min_hashes)


def test_species_notes_of_a_clear_match():
    reference, notes = species_notes(identification('NC_006932-NC_006933.fasta\tq\t0.001\t0\t950/1000\n'),
                                     REFERENCES)
    assert reference.code == 'ab1'
    assert notes == []


def test_species_notes_say_when_the_distance_cannot_separate_the_references():
    """
    Past MAX_CLOSE_MATCHES the margin is measuring nothing, so the names are replaced by the count.

    Sharing a few hashes of 1000 makes mash give many references exactly the same distance, and they all
    tie with the best one. Listing them put 127 genomes and 10 kB of text in one Notes cell.
    """
    hits = 'a.fasta\tq\t0.2\t1\t4/1000\n' + ''.join(
        'ref{}.fasta\tq\t0.2\t1\t4/1000\n'.format(i) for i in range(60))
    found = identification(hits, min_hashes=1)
    assert found.identified
    assert len(found.within_margin) == 60
    assert found.undiscriminating == 60
    assert found.close_matches == []       # not offered as evidence
    _, notes = species_notes(found, REFERENCES)
    note = next(n for n in notes if 'cannot separate' in n)
    assert '60 of them sit within 0.0005 of the closest' in note
    assert '4/1000 hashes' in note
    assert 'ref1.fasta' not in note
    assert not any('nearly as close' in n for n in notes)


def test_species_notes_still_name_a_handful_of_close_references():
    """Below the cap the names are the evidence for the call, and are kept."""
    hits = 'a.fasta\tq\t0.001\t0\t950/1000\n' + ''.join(
        'ref{}.fasta\tq\t0.001\t0\t949/1000\n'.format(i) for i in range(3))
    found = identification(hits)
    assert found.undiscriminating == 0
    assert len(found.close_matches) == 3
    _, notes = species_notes(found, REFERENCES)
    assert any('nearly as close' in n for n in notes)
    assert not any('cannot separate' in n for n in notes)


def test_species_notes_call_a_distant_closest_reference_a_loose_match():
    found = identification('NC_006932-NC_006933.fasta\tq\t0.00900\t0\t700/1000\n')
    reference, notes = species_notes(found, REFERENCES)
    assert reference.code == 'ab1'
    assert len(notes) == 1
    assert notes[0].startswith('loose match: the closest reference is 0.00900 away and shares 700/1000 hashes')
    assert 'the species named is still the closest one' in notes[0]


def test_species_notes_say_nothing_about_a_distance_just_inside_the_threshold():
    found = identification('NC_006932-NC_006933.fasta\tq\t0.00500\t0\t870/1000\n')
    assert species_notes(found, REFERENCES)[1] == []


def test_species_notes_take_a_loose_distance_of_their_own():
    found = identification('NC_006932-NC_006933.fasta\tq\t0.00100\t0\t950/1000\n')
    assert species_notes(found, REFERENCES, loose_distance=0.0005)[1][0].startswith('loose match:')
    assert species_notes(found, REFERENCES, loose_distance=0.01)[1] == []


def test_species_notes_of_a_sample_that_is_not_brucella_do_not_mention_a_loose_match():
    found = identification('NC_006932-NC_006933.fasta\tq\t0.30000\t0\t3/1000\n')
    notes = species_notes(found, REFERENCES)[1]
    assert len(notes) == 1
    assert notes[0].startswith('not identified:')


def test_species_notes_count_a_nearly_as_close_genome_of_the_same_species():
    found = identification('NC_006932-NC_006933.fasta\tq\t0.00100\t0\t950/1000\n',
                           'NZ_CP007682-NZ_CP007683.fasta\tq\t0.00110\t0\t940/1000\n')
    reference, notes = species_notes(found, REFERENCES)
    assert reference.code == 'ab1'
    assert notes == ['nearly as close: 1 other Brucella abortus genome(s)']


def test_species_notes_list_another_species_but_count_the_same_one():
    references = parse_references([
        '\t'.join(COLUMNS) + '\n',
        'NC_006932-NC_006933.fasta\tab1\tBrucella abortus\t1\t9-941\tNC_006932.1\n',
        'NZ_CP007682-NZ_CP007683.fasta\tab3\tBrucella abortus\t3\tBER\tNZ_CP007682.1\n',
        'NC_010169-NC_010167.fasta\tsuis2\tBrucella suis\t2\tATCC 23445\tNC_010169.1\n',
    ])
    found = identification('NC_006932-NC_006933.fasta\tq\t0.00100\t0\t950/1000\n',
                           'NC_010169-NC_010167.fasta\tq\t0.00105\t0\t945/1000\n',
                           'NZ_CP007682-NZ_CP007683.fasta\tq\t0.00110\t0\t940/1000\n',
                           'mystery.fasta\tq\t0.00115\t0\t938/1000\n')
    reference, notes = species_notes(found, references)
    assert reference.code == 'ab1'
    assert notes == ['nearly as close: NC_010169-NC_010167.fasta (Brucella suis bv. 2, 945/1000 hashes), '
                     'mystery.fasta (not in the reference table, 938/1000 hashes), '
                     '1 other Brucella abortus genome(s)']


def test_species_notes_below_the_threshold():
    found = identification('NC_006932-NC_006933.fasta\tq\t0.2\t0\t20/1000\n')
    reference, notes = species_notes(found, REFERENCES)
    assert reference.code == 'ab1'  # The closest reference is known, it is simply too far away
    assert 'not identified' in notes[0]
    assert '20/1000' in notes[0]


def test_species_notes_of_a_reference_that_is_not_in_the_table():
    found = identification('unknown.fasta\tq\t0.001\t0\t950/1000\n',
                           'NC_006932-NC_006933.fasta\tq\t0.0010\t0\t949/1000\n')
    reference, notes = species_notes(found, REFERENCES)
    assert reference is None
    assert 'not in the reference table' in notes[0]
    assert close_label(REFERENCES, found.best) == 'not in the reference table'


def test_sample_result_without_any_result():
    result = SampleResult('s', ['s_R1.fastq.gz'])
    assert not result.identified
    assert (result.species, result.code) == (UNKNOWN, UNKNOWN)
    assert result.sequence_type == UNKNOWN
    assert result.closest_reference == ''
    assert result.note == ''


def test_sample_result_of_a_sample_whose_reference_is_not_in_the_table():
    result = SampleResult('s', [], identification=identification('unknown.fasta\tq\t0.001\t0\t950/1000\n'))
    assert result.identified
    assert result.species == UNKNOWN
    assert result.closest_reference == 'unknown.fasta'


def _raise(error):
    raise error


def fake_tools(monkeypatch, hits=None, mlst=None, error=None):
    """Replace MASH and Seal by what they would have returned."""
    def fake_identify(files, sketch, min_hashes=150, threads=1, min_copies=0, sample=None):
        if error:
            raise BrucellaIDError(error)
        return identification(*(hits or ['NC_006932-NC_006933.fasta\tq\t0.001\t0\t950/1000\n']),
                              min_hashes=min_hashes)

    def fake_run_mlst(files, **options):
        return mlst or MlstResult(st='1', calls=[AlleleCall('gap', '2', 40)], total_reads=1000)

    monkeypatch.setattr(pipeline, 'identify', fake_identify)
    monkeypatch.setattr(pipeline, 'run_mlst', fake_run_mlst)


def test_analyse_sample(monkeypatch):
    fake_tools(monkeypatch)
    result = analyse_sample(Sample('s', ['s_R1.fastq.gz']), references=REFERENCES)
    assert result.species == 'Brucella abortus'
    assert result.code == 'ab1'
    assert result.sequence_type == '1'
    assert result.error is None


def test_analyse_sample_collects_the_notes_of_the_loci(monkeypatch):
    mlst = MlstResult(st=UNKNOWN, calls=[AlleleCall('gap', None, note='no read matched this locus')],
                      note='no allele at gap')
    fake_tools(monkeypatch, mlst=mlst)
    result = analyse_sample(Sample('s', ['s_R1.fastq.gz']), references=REFERENCES)
    assert 'gap: no read matched this locus' in result.note
    assert 'ST ND: no allele at gap' in result.note


def test_analyse_sample_without_mlst(monkeypatch):
    fake_tools(monkeypatch)
    result = analyse_sample(Sample('s', ['s_R1.fastq.gz']), references=REFERENCES, mlst=False)
    assert result.mlst is None
    assert result.sequence_type == UNKNOWN


def test_analyse_sample_keeps_the_error_of_a_sample(monkeypatch, caplog):
    fake_tools(monkeypatch, error='MASH failed')
    result = analyse_sample(Sample('s', ['s_R1.fastq.gz']), references=REFERENCES)
    assert result.error == 'MASH failed'
    assert result.identification is None
    assert 's: MASH failed' in caplog.text


def test_analyse_sample_keeps_the_species_when_the_typing_fails(monkeypatch):
    """MASH had already answered: a Seal failure costs the sample its ST, not its species."""
    fake_tools(monkeypatch)
    monkeypatch.setattr(pipeline, 'run_mlst', lambda files, **options: _raise(BrucellaIDError('Seal failed')))
    result = analyse_sample(Sample('s', ['s_R1.fastq.gz']), references=REFERENCES)
    assert result.error == 'Seal failed'
    assert (result.species, result.code) == ('Brucella abortus', 'ab1')
    assert result.mlst is None


def test_analyse_sample_reads_the_reference_table_of_the_package(monkeypatch):
    """With no table passed, the extended one of the package is read, as it is from 2.0."""
    fake_tools(monkeypatch, hits=['GCA_000163115.1_ASM16311v1_genomic.fna.gz\tq\t0.001\t0\t950/1000\n'])
    result = analyse_sample(Sample('s', ['s_R1.fastq.gz']), mlst=False)
    assert result.species == 'Brucella abortus'
    assert result.code == 'abortus5'  # The extended set codes the epithet and the biovar, not vSNP's "ab5"


def test_an_unexpected_error_is_reported_against_its_sample_and_not_raised(monkeypatch, caplog):
    """
    Anything escaping analyse_sample() ends the run before a single report is written.

    analyse_samples() maps it over a thread pool, so a MemoryError from one pathological input file lost
    every other sample's result. It is now that sample's error, with its type named so a bug still reads
    as a bug.
    """
    fake_tools(monkeypatch)

    def boom(*args, **options):
        raise MemoryError('out of memory')

    monkeypatch.setattr(pipeline, 'identify', boom)
    samples = [Sample(name, ['{}_R1.fastq.gz'.format(name)]) for name in ('a', 'b')]
    with caplog.at_level('ERROR'):
        results = pipeline.analyse_samples(samples, jobs=2, threads=1, references=REFERENCES, mlst=False)
    assert [r.sample for r in results] == ['a', 'b']
    assert all(r.error == 'MemoryError: out of memory' for r in results)
    assert 'MemoryError' in caplog.text


def test_analyse_samples_keeps_the_order_of_the_samples(monkeypatch):
    fake_tools(monkeypatch)
    samples = [Sample(name, ['{}_R1.fastq.gz'.format(name)]) for name in ('c', 'a', 'b')]
    results = analyse_samples(samples, jobs=3, threads=4, references=REFERENCES, mlst=False)
    assert [r.sample for r in results] == ['c', 'a', 'b']


def test_analyse_samples_of_one_sample(monkeypatch, caplog):
    fake_tools(monkeypatch)
    results = analyse_samples([Sample('a', ['a_R1.fastq.gz'])], jobs=4, threads=1, references=REFERENCES,
                              mlst=False)
    assert len(results) == 1
    assert 'at a time' not in caplog.text  # One sample: nothing to say about sharing the threads


@pytest.mark.parametrize('min_hashes, identified', [(150, True), (1000, False)])
def test_analyse_sample_follows_the_threshold(monkeypatch, min_hashes, identified):
    fake_tools(monkeypatch, hits=['NC_006932-NC_006933.fasta\tq\t0.001\t0\t950/1000\n'])
    result = analyse_sample(Sample('s', ['s_R1.fastq.gz']), references=REFERENCES, min_hashes=min_hashes,
                            mlst=False)
    assert result.identified is identified
