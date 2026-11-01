"""The facts a quality assurance report rests on: checksums, status, grouping and margins."""

import os
import socket

import pytest

from brucellaid import qa
from brucellaid.mash import Identification, parse_dist
from brucellaid.mlst import UNKNOWN, AlleleCall, MlstResult
from brucellaid.pipeline import SampleResult, species_notes
from brucellaid.references import COLUMNS, parse_references
from brucellaid.species import call_species, load_st_species

REFERENCES = parse_references([
    '\t'.join(COLUMNS) + '\n',
    'NC_006932-NC_006933.fasta\tab1\tBrucella abortus\t1\t9-941\tNC_006932.1\n',
    'NC_010103-NC_010104.fasta\tcanis\tBrucella canis\t\tATCC 23365\tNC_010103.1\n',
    'B-REF-BS4-40.fasta\tsuis4\tBrucella suis\t4\t40\t\n',
])


def sample(name='s', hits='NC_006932-NC_006933.fasta\tq\t0.0012\t0\t950/1000\n', st='1', error=None,
           min_hashes=150, files=()):
    found = Identification(parse_dist(hits), min_hashes=min_hashes)
    reference, notes = species_notes(found, REFERENCES)
    mlst = MlstResult(st=st, calls=[AlleleCall('gap', '2', 40, '3', 12)], total_reads=1000,
                      assigned_reads=900) if st else None
    call = call_species(found, REFERENCES, mlst=mlst, table=load_st_species())
    return SampleResult(name, list(files), identification=found, reference=reference, mlst=mlst, call=call,
                        notes=notes, error=error)


def test_file_md5_and_size(tmp_path):
    path = tmp_path / 'reads.fastq'
    path.write_bytes(b'hello')
    # The MD5 of "hello", so that a reviewer matching this report against another tool gets the same string
    assert qa.file_md5(path) == '5d41402abc4b2a76b9719d911017c592'
    facts = qa.FileFacts.describe(path)
    assert (facts.name, facts.size, facts.size_text) == ('reads.fastq', 5, '5 B')
    assert facts.md5 == '5d41402abc4b2a76b9719d911017c592'
    assert facts.target == ''  # Not a symbolic link
    assert qa.FileFacts.describe(path, md5=False).md5 == ''


def test_file_facts_of_a_symbolic_link(tmp_path):
    real = tmp_path / 'real.fastq'
    real.write_bytes(b'x' * 2048)
    link = tmp_path / 'link.fastq'
    link.symlink_to(real)
    facts = qa.FileFacts.describe(link)
    assert facts.path == str(link)  # The path the operator gave, not the one it points at
    assert facts.target == str(real)
    assert facts.size_text == '2.0 kB'


@pytest.mark.parametrize('size, text', [
    (0, '0 B'), (512, '512 B'), (1024, '1.0 kB'), (1536, '1.5 kB'), (1048576, '1.0 MB'),
    (3 * 1024 ** 3, '3.0 GB'), (2 * 1024 ** 4, '2.0 TB'), (5 * 1024 ** 5, '5120.0 TB'),
])
def test_human_size(size, text):
    assert qa.human_size(size) == text


def test_describe_inputs_skips_a_file_that_has_gone(tmp_path):
    there = tmp_path / 'a_R1.fastq'
    there.write_bytes(b'reads')
    result = sample(files=[str(there), str(tmp_path / 'gone_R2.fastq')])
    inputs = qa.describe_inputs([result])
    assert [f.name for f in inputs['s']] == ['a_R1.fastq']


def test_fasta_count(tmp_path):
    path = tmp_path / 'alleles.fasta'
    path.write_text('>gap_1\nACGT\nACGT\n>gap_2\nACGT\n')
    assert qa.fasta_count(path) == 2


def test_environment_and_user_name(monkeypatch):
    info = qa.environment()
    assert info['operator'] == info['user'] == qa.user_name()
    assert info['host'] == socket.gethostname()
    assert info['working_directory'] == os.getcwd()
    assert qa.environment(operator='A Reviewer')['operator'] == 'A Reviewer'
    assert qa.environment(operator='A Reviewer')['user'] == info['user']


def test_user_name_without_a_login_name(monkeypatch):
    """Some containers have no passwd entry for the user running the process."""
    monkeypatch.setattr('brucellaid.qa.getpass.getuser', lambda: (_ for _ in ()).throw(KeyError('no name')))
    monkeypatch.setattr('brucellaid.qa.os.getuid', lambda: 1234, raising=False)
    assert qa.user_name() == '1234'


def test_user_name_without_a_login_name_or_a_uid(monkeypatch):
    monkeypatch.setattr('brucellaid.qa.getpass.getuser', lambda: (_ for _ in ()).throw(OSError('no name')))
    monkeypatch.delattr('brucellaid.qa.os.getuid', raising=False)
    assert qa.user_name() == 'unknown'


def test_data_facts(tmp_path):
    custom = tmp_path / 'my_sketch.msh'
    custom.write_bytes(b'sketch')
    facts = qa.data_facts([('MASH sketch', str(custom), '1 genome'),
                           ('Missing', str(tmp_path / 'nothing.tsv'), ''),
                           ('Nothing at all', '', '')])
    assert len(facts) == 1
    assert facts[0]['name'] == 'MASH sketch'
    assert facts[0]['md5'] == qa.file_md5(custom)
    assert facts[0]['packaged'] is False
    assert qa.data_facts([('MASH sketch', str(custom), '')], md5=False)[0]['md5'] == ''


def test_data_facts_knows_the_packaged_files():
    from brucellaid.references import SKETCH_FILE, data_path
    facts = qa.data_facts([('MASH sketch', data_path(SKETCH_FILE), '16 genomes')])
    assert facts[0]['packaged'] is True


def test_is_packaged_of_an_unreadable_path(monkeypatch):
    monkeypatch.setattr('brucellaid.qa.data_path', lambda name: (_ for _ in ()).throw(OSError('gone')))
    assert qa.is_packaged('/anything') is False


@pytest.mark.parametrize('kwargs, expected', [
    ({}, qa.OK),
    ({'error': 'Seal failed'}, qa.FAILED),
    ({'min_hashes': 1000}, qa.WARNING),  # Not identified
    ({'st': UNKNOWN}, qa.WARNING),  # No sequence type
    ({'st': 'new'}, qa.WARNING),  # A profile the table does not hold
    ({'st': '20'}, qa.WARNING),  # ST 20 is B. canis, and the closest reference is B. abortus: a conflict
])
def test_status(kwargs, expected):
    assert qa.status(sample(**kwargs)) == expected


def test_status_of_a_sample_whose_st_was_not_asked_for():
    """--no-mlst leaves the ST at ND, which is not a warning: nobody asked for it."""
    result = sample(st=None)
    assert result.mlst is None
    assert qa.status(result) == qa.OK


def test_counted():
    assert qa.counted([sample('a'), sample('b'), sample('c', error='x')]) == {
        qa.OK: 2, qa.WARNING: 0, qa.FAILED: 1}


def test_groups_puts_the_largest_first_and_the_failed_last():
    results = [sample('a'), sample('b', st='2'), sample('c'), sample('d', error='boom'), sample('e')]
    ordered = qa.groups(results)
    assert [key for key, _ in ordered] == [('Brucella abortus', '1'), ('Brucella abortus', '2'), ('', '')]
    assert [r.sample for r in ordered[0][1]] == ['a', 'c', 'e']  # The order they were analysed in
    assert [r.sample for r in ordered[-1][1]] == ['d']


def test_groups_without_a_failed_sample():
    assert [key for key, _ in qa.groups([sample('a')])] == [('Brucella abortus', '1')]


@pytest.mark.parametrize('call, expected', [
    (AlleleCall('gap', '2', 40, '3', 10), (30, 4.0)),
    (AlleleCall('gap', '2', 40, '3', 0), (40, None)),  # A runner-up with no reads
    (AlleleCall('gap', '2', 40), (None, None)),  # No runner-up at all
    (AlleleCall('gap', None, 0, '3', 2), (None, None)),  # No call
])
def test_margin(call, expected):
    assert qa.margin(call) == expected


def test_close_hits_marks_the_ones_in_the_running():
    result = sample(hits='NC_010103-NC_010104.fasta\tq\t0.000239\t0\t990/1000\n'
                         'B-REF-BS4-40.fasta\tq\t0.000420\t0\t980/1000\n'
                         'NC_006932-NC_006933.fasta\tq\t0.0029\t0\t888/1000\n', st='20')
    rows = qa.close_hits(result, REFERENCES)
    assert [(reference.code, close) for _, reference, close in rows] == [
        ('canis', False), ('suis4', True), ('ab1', False)]


def test_close_hits_of_a_reference_the_table_does_not_name():
    result = sample(hits='/data/refs/unknown.fasta\tq\t0.001\t0\t950/1000\n')
    assert qa.close_hits(result, REFERENCES)[0][1] is None


def test_close_hits_without_an_identification():
    assert qa.close_hits(SampleResult('s', []), REFERENCES) == []


def test_warnings_puts_the_errors_first():
    results = [sample('a', hits='NC_010103-NC_010104.fasta\tq\t0.000239\t0\t990/1000\n'
                                'B-REF-BS4-40.fasta\tq\t0.000420\t0\t980/1000\n', st='20'),
               sample('b', error='Seal failed')]
    rows = qa.warnings(results)
    assert rows[0] == ('b', qa.FAILED, 'Seal failed')
    assert all(row[0] == 'a' for row in rows[1:])
    assert any('nearly as close' in row[2] for row in rows[1:])
