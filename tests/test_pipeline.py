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


def test_species_notes_of_two_nearly_identical_references():
    found = identification('NC_006932-NC_006933.fasta\tq\t0.00100\t0\t950/1000\n',
                           'NZ_CP007682-NZ_CP007683.fasta\tq\t0.00110\t0\t940/1000\n')
    reference, notes = species_notes(found, REFERENCES)
    assert reference.code == 'ab1'
    assert len(notes) == 1
    assert 'nearly as close' in notes[0]
    assert 'Brucella abortus bv. 3' in notes[0]


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
    def fake_identify(files, sketch, min_hashes=150, threads=1):
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
    fake_tools(monkeypatch)
    result = analyse_sample(Sample('s', ['s_R1.fastq.gz']), mlst=False)
    assert result.species == 'Brucella abortus'  # ab1 is in the table of the package too


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
