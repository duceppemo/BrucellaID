"""The Seal wrapper: the command it builds, the paths it has to work around, and its output."""

import tempfile
from pathlib import Path

import pytest

from brucellaid import seal
from brucellaid.errors import BrucellaIDError
from brucellaid.seal import (
    SealStats,
    auto_memory,
    check_seal,
    content_extension,
    executable,
    extensions,
    failure_reason,
    parse_stats,
    run_seal,
    safe_paths,
    safe_temporary_folder,
    seal_command,
    system_memory_mb,
    versions,
)
from tests.conftest import requires_seal, write_fastq

STATS = ('#File\tsample_R1.fastq.gz\n'
         '#Total\t822714\t123407100\n'
         '#Matched\t799\t0.09712%\n'
         '#Name\tReads\tReadsPct\n'
         'gap_2 the description is dropped\t62\t0.00754%\n'
         '\n'
         'aroA_1\t0\t0.00000%\n')


def test_seal_command_of_paired_reads(monkeypatch):
    monkeypatch.setattr(seal, 'executable', lambda: '/bin/seal.sh')
    cmd = seal_command(['a_R1.fastq.gz', 'a_R2.fastq.gz'], 'alleles.fasta', 'stats.tsv', 4, '2g')
    assert cmd[:2] == ['/bin/seal.sh', '-Xmx2g']
    assert 'in2=a_R2.fastq.gz' in cmd
    assert 'int=f' not in cmd
    assert 'ambiguous=toss' in cmd and 'clearzone=0' in cmd and 'hdist=0' in cmd
    assert 'threads=4' in cmd


def test_seal_command_of_single_reads_and_several_references(monkeypatch):
    monkeypatch.setattr(seal, 'executable', lambda: '/bin/seal.sh')
    cmd = seal_command(['a.fastq'], ['one.fasta', 'two.fasta'], 'stats.tsv', 1, '1g', k=31, hdist=1,
                       ambiguous='toss', clearzone=5)
    assert 'int=f' in cmd
    assert 'ref=one.fasta,two.fasta' in cmd
    assert 'k=31' in cmd and 'hdist=1' in cmd and 'ambiguous=toss' in cmd and 'clearzone=5' in cmd


@pytest.mark.parametrize('name, expected', [
    ('a.fastq.gz', '.fastq.gz'),
    ('a.fastq', '.fastq'),
    ('a.b.fastq', '.b.fastq'),
    ('dataset_1.dat', '.dat'),
    ('weird name.fastq', '.fastq'),
    ('a.xmx.fastq', '.fastq'),  # BBTools would read ".xmx" as a memory setting
])
def test_extensions(name, expected):
    assert extensions(name) == expected


def test_content_extension(tmp_path):
    assert content_extension(write_fastq(tmp_path / 'a.dat', [('r', 'ACGT')], gzipped=False)) == '.fastq'
    assert content_extension(write_fastq(tmp_path / 'b.dat', [('r', 'ACGT')], gzipped=True)) == '.fastq.gz'
    fasta = tmp_path / 'c.dat'
    fasta.write_text('>r\nACGT\n')
    assert content_extension(fasta) == '.fasta'
    empty = tmp_path / 'd.dat'
    empty.write_text('')
    assert content_extension(empty) == '.fasta'
    broken = tmp_path / 'e.dat'
    broken.write_bytes(b'\x1f\x8b\x08\x00broken')
    assert content_extension(broken) == '.fasta.gz'  # A broken stream: Seal reports the real error
    assert content_extension(tmp_path / 'missing.dat') == '.fastq'


def test_safe_paths_links_what_seal_cannot_read(tmp_path):
    folder = tmp_path / 'tmp'
    folder.mkdir()
    plain = write_fastq(tmp_path / 'a_R1.fastq.gz', [('r', 'ACGT')])
    spaces = write_fastq(tmp_path / 'with space.fastq.gz', [('r', 'ACGT')])
    mislabelled = write_fastq(tmp_path / 'c.fasta', [('r', 'ACGT')], gzipped=False)
    safe = safe_paths([plain, spaces, mislabelled], folder, 'in', check_extension=True)
    assert safe[0] == str(plain)  # Already safe
    assert safe[1] == str(folder / 'in1.fastq.gz')
    assert safe[2] == str(folder / 'in2.fastq')  # A fastq file named ".fasta" is linked under the right name
    assert all(Path(path).exists() for path in safe)  # Every link resolves


def test_safe_temporary_folder():
    folder = safe_temporary_folder()
    try:
        assert 'xmx' not in folder.lower()
        assert Path(folder).is_dir()
    finally:
        seal.shutil.rmtree(folder, ignore_errors=True)


def test_safe_temporary_folder_of_a_temporary_path_seal_cannot_read(monkeypatch, tmp_path):
    made = []

    def mkdtemp(prefix):
        folder = tmp_path / '{}xmx{}'.format(prefix, len(made))
        folder.mkdir()
        made.append(folder)
        return str(folder)

    monkeypatch.setattr(seal.tempfile, 'mkdtemp', mkdtemp)
    with pytest.raises(BrucellaIDError, match='Set TMPDIR'):
        safe_temporary_folder(attempts=2)
    assert len(made) == 2
    assert not any(folder.exists() for folder in made)  # The folders it had to reject are removed


def test_parse_stats():
    stats = parse_stats_from_text(STATS)
    assert stats.reads == 822714
    assert stats.bases == 123407100
    assert stats.matched == 799
    assert stats.counts == {'gap_2': 62, 'aroA_1': 0}


def parse_stats_from_text(text, tmp=None):
    folder = Path(tmp or tempfile.mkdtemp())
    path = folder / 'stats.tsv'
    path.write_text(text)
    return parse_stats(path)


def test_parse_stats_of_an_unexpected_line():
    with pytest.raises(BrucellaIDError, match='Unexpected line in Seal stats file'):
        parse_stats_from_text('gap_2\tnot a number\n')


def test_seal_stats_defaults():
    assert SealStats().counts == {}
    assert SealStats().reads is None


@pytest.mark.parametrize('lines, expected', [
    (['java.lang.OutOfMemoryError: Java heap space'], 'java.lang.OutOfMemoryError: Java heap space'),
    (['Exception in thread "main" java.lang.RuntimeException:', '\tat jgi.Seal.main()', 'Bad input'],
     'Exception in thread "main" java.lang.RuntimeException: Bad input'),
    (['Exception in thread "main" x', 'java.io.FileNotFoundException: missing.fastq'],
     'java.io.FileNotFoundException: missing.fastq'),
    (['nothing to report here', 'the last line'], 'the last line'),
    (['java.lang.RuntimeException:'], 'java.lang.RuntimeException:'),  # Nothing follows the message
])
def test_failure_reason(lines, expected):
    assert failure_reason(lines) == expected


def test_check_seal_when_it_is_not_installed(monkeypatch):
    monkeypatch.setattr(seal, 'executable', lambda: None)
    with pytest.raises(BrucellaIDError, match='was not found in your PATH'):
        check_seal()
    assert versions() == {'BBTools': 'unknown', 'Java': 'unknown'}


def test_executable_falls_back_next_to_python(monkeypatch, tmp_path):
    monkeypatch.setattr(seal.shutil, 'which', lambda name: None)
    monkeypatch.setattr(seal.sys, 'executable', str(tmp_path / 'bin' / 'python'))
    (tmp_path / 'bin').mkdir()
    assert executable() is None
    executable.cache_clear()
    script = tmp_path / 'bin' / 'seal.sh'
    script.write_text('#!/bin/sh\n')
    script.chmod(0o755)
    assert executable() == str(script)


def test_system_memory_mb():
    assert system_memory_mb() is None or system_memory_mb() > 0


def test_auto_memory_rounds_up_and_stays_within_the_machine(monkeypatch):
    monkeypatch.setattr(seal, 'system_memory_mb', lambda: 64000)
    assert auto_memory(0) == '1024m'
    assert auto_memory(100_000_000) == '4608m'
    monkeypatch.setattr(seal, 'system_memory_mb', lambda: 2000)
    assert auto_memory(100_000_000, jobs=2) == '850m'  # 85% of the machine, shared between the jobs
    monkeypatch.setattr(seal, 'system_memory_mb', lambda: None)
    assert auto_memory(0) == '1024m'


@requires_seal
def test_versions_of_the_installation():
    assert versions()['BBTools'] != 'unknown'


@requires_seal
def test_run_seal_counts_the_reads_of_each_reference(tmp_path):
    refs = tmp_path / 'refs.fasta'
    refs.write_text('>one\n{}\n>two\n{}\n'.format('ACGT' * 20, 'TTGCA' * 16))
    reads = write_fastq(tmp_path / 'reads.fastq', [('r1', 'ACGT' * 20), ('r2', 'TTGCA' * 16), ('r3', 'A' * 80)],
                        gzipped=False)
    stats = run_seal([str(reads)], refs)
    assert stats.counts == {'one': 1, 'two': 1}
    assert stats.reads == 3
    assert stats.matched == 2


@requires_seal
def test_run_seal_reports_what_seal_wrote_when_it_fails(tmp_path):
    reads = write_fastq(tmp_path / 'reads.fastq', [('r1', 'ACGT' * 20)], gzipped=False)
    with pytest.raises(BrucellaIDError, match='Seal failed'):
        run_seal([str(reads)], tmp_path / 'missing_refs.fasta')


@requires_seal
def test_run_seal_of_an_input_file_that_is_not_there(tmp_path):
    refs = tmp_path / 'refs.fasta'
    refs.write_text('>one\n{}\n'.format('ACGT' * 20))
    with pytest.raises(BrucellaIDError, match='Seal failed'):
        run_seal([str(tmp_path / 'missing.fastq')], refs)


class FakeRun:
    """What subprocess.run returns, for the version lookups."""
    def __init__(self, stdout='', stderr=''):
        self.stdout, self.stderr = stdout, stderr


def test_versions_when_seal_prints_no_version(monkeypatch, tmp_path):
    monkeypatch.setattr(seal, 'executable', lambda: str(tmp_path / 'seal.sh'))
    monkeypatch.setattr(seal.shutil, 'which', lambda name: None)  # No java either
    monkeypatch.setattr(seal.subprocess, 'run', lambda *args, **options: FakeRun('no version here'))
    assert versions() == {'BBTools': 'unknown', 'Java': 'unknown'}


def test_versions_reads_the_java_next_to_seal(monkeypatch, tmp_path):
    (tmp_path / 'seal.sh').write_text('#!/bin/sh\n')
    (tmp_path / 'java').write_text('#!/bin/sh\n')
    monkeypatch.setattr(seal, 'executable', lambda: str(tmp_path / 'seal.sh'))
    runs = [FakeRun('BBMap version 39.01'), FakeRun(stderr='openjdk version "21"\nmore lines')]
    monkeypatch.setattr(seal.subprocess, 'run', lambda *args, **options: runs.pop(0))
    assert versions() == {'BBTools': '39.01', 'Java': 'openjdk version "21"'}


def test_versions_when_java_prints_nothing(monkeypatch, tmp_path):
    (tmp_path / 'seal.sh').write_text('#!/bin/sh\n')
    (tmp_path / 'java').write_text('#!/bin/sh\n')
    monkeypatch.setattr(seal, 'executable', lambda: str(tmp_path / 'seal.sh'))
    runs = [FakeRun('BBTools version 40.05'), FakeRun('')]
    monkeypatch.setattr(seal.subprocess, 'run', lambda *args, **options: runs.pop(0))
    assert versions() == {'BBTools': '40.05', 'Java': 'unknown'}


def test_versions_when_seal_cannot_be_run(monkeypatch, tmp_path):
    monkeypatch.setattr(seal, 'executable', lambda: str(tmp_path / 'seal.sh'))

    def fails(*args, **options):
        raise OSError('no such file')

    monkeypatch.setattr(seal.subprocess, 'run', fails)
    assert versions() == {'BBTools': 'unknown', 'Java': 'unknown'}


def test_system_memory_mb_when_the_machine_does_not_say(monkeypatch):
    def fails(name):
        raise ValueError(name)

    monkeypatch.setattr(seal.os, 'sysconf', fails)
    assert system_memory_mb() is None
