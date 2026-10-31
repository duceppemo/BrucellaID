"""The MASH identification: the command, its output, and what is made of it."""

import gzip

import pytest

from brucellaid import mash
from brucellaid.errors import BrucellaIDError
from brucellaid.mash import (
    CLOSE_MATCH_MARGIN,
    DEFAULT_MIN_HASHES,
    Identification,
    check_mash,
    executable,
    identify,
    parse_dist,
    parse_info,
    run_dist,
    sketch_references,
    versions,
    write_reads,
)
from tests.conftest import requires_mash, write_fastq

DIST_OUTPUT = (
    'NC_006932-NC_006933.fasta\tstdin\t0.00123\t0\t900/1000\n'
    'NZ_CP007682-NZ_CP007683.fasta\tstdin\t0.0011\t0\t910/1000\n'
    '\n'
    'Bceti1Cudo.fasta\tstdin\t0.0294416\t1.5e-20\t493/1000\n'
)


def test_parse_dist_sorts_the_closest_reference_first():
    hits = parse_dist(DIST_OUTPUT)
    assert [hit.reference for hit in hits] == ['NZ_CP007682-NZ_CP007683.fasta', 'NC_006932-NC_006933.fasta',
                                               'Bceti1Cudo.fasta']
    assert hits[0].shared_hashes == 910
    assert hits[0].total_hashes == 1000
    assert hits[2].p_value == 1.5e-20
    assert hits[0].identity == pytest.approx(99.89)


def test_parse_dist_breaks_a_tie_on_the_hashes_then_the_name():
    hits = parse_dist('b.fasta\tq\t0.001\t0\t500/1000\n'
                      'a.fasta\tq\t0.001\t0\t500/1000\n'
                      'c.fasta\tq\t0.001\t0\t600/1000\n')
    assert [hit.reference for hit in hits] == ['c.fasta', 'a.fasta', 'b.fasta']


@pytest.mark.parametrize('output', ['a.fasta\tq\t0.001\t0\n', 'a.fasta\tq\tnot a number\t0\t1/1000\n'])
def test_parse_dist_errors(output):
    with pytest.raises(BrucellaIDError, match='Unexpected line'):
        parse_dist(output)


def test_identification_of_a_good_match():
    hits = parse_dist(DIST_OUTPUT)
    identification = Identification(hits)
    assert identification.min_hashes == DEFAULT_MIN_HASHES
    assert identification.identified
    assert identification.best.reference == 'NZ_CP007682-NZ_CP007683.fasta'
    assert [hit.reference for hit in identification.close_matches] == ['NC_006932-NC_006933.fasta']


def test_identification_below_the_threshold():
    identification = Identification(parse_dist('a.fasta\tq\t0.3\t0\t100/1000\n'))
    assert not identification.identified
    assert identification.close_matches == []


def test_close_match_margin_is_a_distance_not_a_ratio():
    hits = parse_dist('a.fasta\tq\t0.001\t0\t900/1000\n'
                      'b.fasta\tq\t{}\t0\t800/1000\n'.format(0.001 + CLOSE_MATCH_MARGIN * 2))
    assert Identification(hits).close_matches == []


def test_write_reads_decompresses_the_gzipped_files(tmp_path):
    plain = write_fastq(tmp_path / 'a.fastq', [('r1', 'ACGT')], gzipped=False)
    gzipped = write_fastq(tmp_path / 'b.fastq.gz', [('r2', 'TGCA')])
    out = tmp_path / 'joined.fastq'
    with open(out, 'wb') as f:
        write_reads([plain, gzipped], f)
    assert out.read_text() == '@r1\nACGT\n+\nIIII\n@r2\nTGCA\n+\nIIII\n'


def test_write_reads_of_a_broken_gzip_file(tmp_path):
    broken = tmp_path / 'broken.fastq.gz'
    broken.write_bytes(b'\x1f\x8b\x08\x00broken')
    with open(tmp_path / 'out', 'wb') as f, pytest.raises((OSError, EOFError, gzip.BadGzipFile)):
        write_reads([broken], f)


def test_check_mash_when_it_is_not_installed(monkeypatch):
    monkeypatch.setattr(mash, 'executable', lambda: None)
    with pytest.raises(BrucellaIDError, match='was not found in your PATH'):
        check_mash()
    assert versions() == {'MASH': 'unknown'}


def test_executable_falls_back_next_to_python(monkeypatch, tmp_path):
    monkeypatch.setattr(mash.shutil, 'which', lambda name: None)
    monkeypatch.setattr(mash.sys, 'executable', str(tmp_path / 'bin' / 'python'))
    (tmp_path / 'bin').mkdir()
    assert executable() is None
    executable.cache_clear()
    (tmp_path / 'bin' / 'mash').write_text('#!/bin/sh\n')
    assert executable() == str(tmp_path / 'bin' / 'mash')


@requires_mash
def test_versions_of_the_installation():
    assert versions()['MASH'] != 'unknown'


@requires_mash
def test_sketch_references_of_the_package():
    ids = sketch_references(mash.Path(__file__).parent.parent / 'brucellaid' / 'data' / 'references.msh')
    assert 'NC_006932-NC_006933.fasta' in ids
    assert sketch_references.cache_info().currsize == 1  # Read once per sketch


@requires_mash
def test_sketch_references_of_a_file_that_is_not_a_sketch(tmp_path):
    not_a_sketch = tmp_path / 'broken.msh'
    not_a_sketch.write_text('not a sketch')
    with pytest.raises(BrucellaIDError, match='cannot be read by MASH'):
        sketch_references(not_a_sketch)


@requires_mash
def test_identify_a_read_set_of_the_reference_genome(sketch, genome_reads):
    msh, _ = sketch
    identification = identify([genome_reads], msh, threads=2)
    assert identification.identified
    assert identification.best.reference.endswith('test_genome.fasta')
    assert identification.best.shared_hashes > 800


@requires_mash
def test_run_dist_reports_what_mash_wrote_when_it_fails(tmp_path, genome_reads):
    not_a_sketch = tmp_path / 'broken.msh'
    not_a_sketch.write_bytes(b'\x00\x01\x02not a sketch')
    with pytest.raises(BrucellaIDError, match='MASH'):
        run_dist([genome_reads], not_a_sketch)


@requires_mash
def test_run_dist_of_an_unreadable_read_file(tmp_path, sketch):
    msh, _ = sketch
    broken = tmp_path / 'broken.fastq.gz'
    broken.write_bytes(b'\x1f\x8b\x08\x00not a gzip stream')
    with pytest.raises(BrucellaIDError, match='cannot be read'):
        run_dist([broken], msh)


def test_parse_info():
    assert parse_info('#Hashes\tLength\tID\tComment\n'
                      '1000\t3286445\ta.fasta\t[2 seqs] ...\n'
                      '\n'
                      '1000\t3286445\tb.fasta\n') == ['a.fasta', 'b.fasta']


@pytest.mark.parametrize('output, message', [
    ('1000\t3286445\n', 'Unexpected line'),
    ('#Hashes\tLength\tID\n', 'holds no genome'),
])
def test_parse_info_errors(output, message):
    with pytest.raises(BrucellaIDError, match=message):
        parse_info(output, 'refs.msh')


def test_versions_when_mash_cannot_be_run(monkeypatch):
    monkeypatch.setattr(mash, 'executable', lambda: '/bin/mash')

    def fails(*args, **options):
        raise OSError('no such file')

    monkeypatch.setattr(mash.subprocess, 'run', fails)
    assert versions() == {'MASH': 'unknown'}


def test_versions_when_mash_prints_nothing(monkeypatch):
    monkeypatch.setattr(mash, 'executable', lambda: '/bin/mash')
    monkeypatch.setattr(mash.subprocess, 'run',
                        lambda *args, **options: type('R', (), {'stdout': '  \n', 'stderr': ''})())
    assert versions() == {'MASH': 'unknown'}


@requires_mash
def test_run_dist_when_mash_compares_nothing(monkeypatch, sketch, genome_reads):
    msh, _ = sketch
    monkeypatch.setattr(mash, 'parse_dist', lambda output: [])
    with pytest.raises(BrucellaIDError, match='the sketch .* is empty'):
        run_dist([genome_reads], msh)
