"""The command line: the arguments it accepts, and what it does with a failure."""

import json

import pytest

from brucellaid import __version__, cli
from brucellaid.errors import BrucellaIDError
from brucellaid.mlst import MlstResult
from brucellaid.pipeline import SampleResult
from brucellaid.references import COLUMNS
from brucellaid.samples import Sample
from tests.conftest import write_fastq


def test_version_and_help(capsys):
    with pytest.raises(SystemExit) as exit_info:
        cli.main(['--version'])
    assert exit_info.value.code == 0
    assert __version__ in capsys.readouterr().out
    with pytest.raises(SystemExit):
        cli.main(['--help'])
    assert 'brucellaid' in capsys.readouterr().out


@pytest.mark.parametrize('argv, message', [
    (['-i', 'reads', '-o', 'out', '-r2', 'b.fastq'], 'only used with -r1'),
    (['-i', 'reads', '-o', 'out', '-s', 'name'], 'only used with -r1'),
    (['-r1', 'a.fastq', '-o', 'out', '-p', '2'], '--parallel is only used with -i'),
    (['-r1', 'a.fastq', '-o', 'out', '--min-hashes', '0'], '--min-hashes must be >= 1'),
    (['-r1', 'a.fastq', '-o', 'out', '--loose-distance', '0'], '--loose-distance must be greater than 0'),
    (['-r1', 'a.fastq', '-o', 'out', '--min-kmer-copies', '-1'], '--min-kmer-copies cannot be negative'),
    (['-r1', 'a.fastq', '-o', 'out', '--min-reads', '0'], '--min-reads must be >= 1'),
    (['-r1', 'a.fastq', '-o', 'out', '-m', '-1'], 'cannot be negative'),
    (['-r1', 'a.fastq', '-o', 'out', '-m', '25'], 'must be smaller than --kmer-size'),
    (['-r1', 'a.fastq', '-o', 'out', '-s', 'a/b'], 'must be a name, not a path'),
    (['-r1', 'a.fastq', '-o', 'out', '-s', ' '], 'must be a name, not a path'),
    (['-r1', 'a.fastq', '-o', 'out', '-s', 'a\tb'], 'cannot contain tabs'),
    (['-o', 'out'], 'one of the arguments'),  # Neither -r1 nor -i
    (['-r1', 'a.fastq', '-i', 'reads', '-o', 'out'], 'not allowed with argument'),
    (['-r1', 'a.fastq'], 'required'),  # No -o
    (['-r1', 'a.fastq', '-o', 'out', '-t', '0'], 'must be >= 1'),
    (['-r1', 'a.fastq', '-o', 'out', '-t', 'two'], 'not an integer'),
    (['-r1', 'a.fastq', '-o', 'out', '--memory', '2k'], 'not a Java memory size'),
    (['-r1', 'a.fastq', '-o', 'out', '--min-st-purity', '1.5'], 'must be a fraction between 0 and 1'),
    (['-r1', 'a.fastq', '-o', 'out', '--min-st-purity', '0'], 'must be a fraction between 0 and 1'),
    (['-r1', 'a.fastq', '-o', 'out', '--min-st-isolates', '0'], 'must be >= 1'),
    (['-r1', 'a.fastq', '-o', 'out', '--st-species', 'st.tsv', '--no-st-species'], 'pick one'),
])
def test_rejected_arguments(argv, message, capsys):
    with pytest.raises(SystemExit) as exit_info:
        cli.main(argv)
    assert exit_info.value.code == 2
    assert message in capsys.readouterr().err


@pytest.mark.parametrize('seconds, expected', [
    (0, '0s'), (0.4, '0s'), (59.6, '1m'), (61, '1m1s'), (3661, '1h1m1s'), (90061, '1d1h1m1s'),
])
def test_format_elapsed(seconds, expected):
    assert cli.format_elapsed(seconds) == expected


def test_available_cpus(monkeypatch):
    assert cli.available_cpus() >= 1
    monkeypatch.delattr(cli.os, 'sched_getaffinity', raising=False)
    assert cli.available_cpus() >= 1


def test_java_memory():
    assert cli.java_memory('2g') == '2g'
    assert cli.java_memory('500M') == '500M'


def test_file_bases(tmp_path):
    path = tmp_path / 'refs.fasta'
    path.write_text('>a\nACGT\n')
    assert cli.file_bases(path) == 8
    assert cli.file_bases(tmp_path / 'missing.fasta') == 0


def fake_run(monkeypatch, results=None, error=None, sketch_ids=None):
    """
    Replace the tool checks and the analysis, so that main() can be run without MASH or Seal.

    sketch_ids has to name a genome of whichever reference table the run will load, because main() refuses
    a sketch none of whose genomes the table names. The default is a genome of the extended table, which
    is the default set from 2.0.
    """
    monkeypatch.setattr(cli, 'check_mash', lambda: '/bin/mash')
    monkeypatch.setattr(cli, 'check_seal', lambda: '/bin/seal.sh')
    ids = ['GCA_000742195.1_ASM74219v1_genomic.fna.gz'] if sketch_ids is None else list(sketch_ids)
    monkeypatch.setattr(cli, 'sketch_references', lambda sketch: ids)
    monkeypatch.setattr(cli, 'compare_with_sketch', lambda ids, references: [])

    def analyse(samples, **options):
        if error:
            raise error
        return results if results is not None else [
            SampleResult(sample.name, list(sample.files), mlst=MlstResult(st='1')) for sample in samples]

    monkeypatch.setattr(cli, 'analyse_samples', analyse)


def test_main_writes_the_reports_and_prints_the_table(tmp_path, monkeypatch, capsys):
    fake_run(monkeypatch)
    reads = write_fastq(tmp_path / 'a_R1.fastq.gz', [('r', 'ACGT')])
    cli.main(['-r1', str(reads), '-o', str(tmp_path / 'out'), '--verbose'])
    out = capsys.readouterr()
    assert out.out.splitlines()[0].startswith('Sample\t')
    assert '\ta\t' in out.out.replace('\n', '\t')
    assert (tmp_path / 'out' / 'brucellaid.tsv').exists()
    assert 'Done in' in out.err


def test_main_without_the_st_species_table(tmp_path, monkeypatch, capsys):
    """--no-st-species: the table is never read, and the species is the MASH one alone."""
    fake_run(monkeypatch)
    reads = write_fastq(tmp_path / 'a_R1.fastq.gz', [('r', 'ACGT')])
    cli.main(['-r1', str(reads), '-o', str(tmp_path / 'out'), '--no-st-species'])
    out = capsys.readouterr()
    assert 'sequence type(s) with a species' not in out.err
    report = json.loads((tmp_path / 'out' / 'brucellaid.json').read_text())
    assert report['parameters']['ST to species table'] is None


def test_main_writes_the_quality_assurance_reports(tmp_path, monkeypatch, capsys):
    """--pdf and --operator: the HTML, the MultiQC file and the PDF, with the operator's name on them."""
    pytest.importorskip('reportlab')
    fake_run(monkeypatch)
    reads = write_fastq(tmp_path / 'a_R1.fastq.gz', [('r', 'ACGT')])
    out = tmp_path / 'out'
    cli.main(['-r1', str(reads), '-o', str(out), '--pdf', '--operator', 'A Reviewer', '--no-md5'])
    markup = (out / 'brucellaid_report.html').read_text(encoding='utf-8')
    assert 'A Reviewer' in markup
    assert 'not computed (--no-md5)' in markup
    assert (out / 'brucellaid_report.pdf').stat().st_size > 1000
    assert json.loads((out / 'brucellaid_mqc.json').read_text())['id'] == 'brucellaid'
    report = json.loads((out / 'brucellaid.json').read_text())
    assert report['operator'] == 'A Reviewer'
    assert report['parameters']['Input checksums'] == 'none'
    assert [item['name'] for item in report['reference_data']] == [
        'MASH sketch', 'Reference table', 'MLST alleles', 'MLST profiles', 'ST to species']
    assert all(item['md5'] == '' for item in report['reference_data'])  # --no-md5 covers these too


def test_main_without_the_html_report(tmp_path, monkeypatch, capsys):
    fake_run(monkeypatch)
    reads = write_fastq(tmp_path / 'a_R1.fastq.gz', [('r', 'ACGT')])
    out = tmp_path / 'out'
    cli.main(['-r1', str(reads), '-o', str(out), '--no-html', '--no-mlst'])
    assert not (out / 'brucellaid_report.html').exists()
    assert (out / 'brucellaid_mqc.json').exists()  # The MultiQC file is small and always written
    report = json.loads((out / 'brucellaid.json').read_text())
    # Without the MLST, the allele and profile files are not part of what the run measured against
    assert [item['name'] for item in report['reference_data']] == [
        'MASH sketch', 'Reference table', 'ST to species']


def test_main_reports_a_user_error_without_a_traceback(tmp_path, monkeypatch, capsys):
    fake_run(monkeypatch, error=BrucellaIDError('the sketch is empty'))
    reads = write_fastq(tmp_path / 'a_R1.fastq.gz', [('r', 'ACGT')])
    with pytest.raises(SystemExit) as exit_info:
        cli.main(['-r1', str(reads), '-o', str(tmp_path / 'out')])
    assert exit_info.value.code == 1
    err = capsys.readouterr().err
    assert 'the sketch is empty' in err
    assert 'Traceback' not in err


def test_main_exits_with_130_when_interrupted(tmp_path, monkeypatch, capsys):
    fake_run(monkeypatch, error=KeyboardInterrupt())
    reads = write_fastq(tmp_path / 'a_R1.fastq.gz', [('r', 'ACGT')])
    with pytest.raises(SystemExit) as exit_info:
        cli.main(['-r1', str(reads), '-o', str(tmp_path / 'out')])
    assert exit_info.value.code == 130
    assert 'Interrupted' in capsys.readouterr().err


def test_main_exits_with_1_when_a_sample_failed(tmp_path, monkeypatch, capsys):
    fake_run(monkeypatch, results=[SampleResult('a', ['a_R1.fastq.gz'], error='Seal failed')])
    reads = write_fastq(tmp_path / 'a_R1.fastq.gz', [('r', 'ACGT')])
    with pytest.raises(SystemExit) as exit_info:
        cli.main(['-r1', str(reads), '-o', str(tmp_path / 'out')])
    assert exit_info.value.code == 1
    assert '1 of 1 sample(s) failed: a' in capsys.readouterr().err


def test_main_warns_about_a_sketch_that_does_not_match_the_table(tmp_path, monkeypatch, capsys):
    fake_run(monkeypatch)
    monkeypatch.setattr(cli, 'compare_with_sketch', lambda ids, references: ['a genome is missing'])
    reads = write_fastq(tmp_path / 'a_R1.fastq.gz', [('r', 'ACGT')])
    cli.main(['-r1', str(reads), '-o', str(tmp_path / 'out'), '--no-mlst'])
    assert 'a genome is missing' in capsys.readouterr().err


def test_main_of_a_folder_of_samples(tmp_path, monkeypatch, capsys):
    fake_run(monkeypatch)
    reads = tmp_path / 'reads'
    reads.mkdir()
    write_fastq(reads / 'a_R1.fastq.gz', [('r', 'ACGT')])
    write_fastq(reads / 'b_R1.fastq.gz', [('r', 'ACGT')])
    cli.main(['-i', str(reads), '-o', str(tmp_path / 'out'), '-p', '2'])
    assert '2 sample(s) found' in capsys.readouterr().err


def test_setup_logging_does_not_add_a_handler_twice():
    cli.setup_logging()
    cli.setup_logging(verbose=True)
    assert len(cli.log.handlers) == 1


def test_sample_fixture_names_the_sample_from_the_file_name():
    assert Sample('a', ['a_R1.fastq.gz']).paired is False


def reference_paths(out):
    report = json.loads((out / 'brucellaid.json').read_text())
    paths = {item['name']: item['path'] for item in report['reference_data']}
    detail = next(item['detail'] for item in report['reference_data'] if item['name'] == 'Reference table')
    return paths, detail


@pytest.mark.parametrize('argv', [[], ['--reference-set', 'extended']])
def test_main_identifies_against_the_extended_set_by_default(tmp_path, monkeypatch, argv):
    """From 2.0 the extended set is what a run with no --reference-set uses."""
    fake_run(monkeypatch)
    reads = write_fastq(tmp_path / 'a_R1.fastq.gz', [('r', 'ACGT')])
    cli.main(['-r1', str(reads), '-o', str(tmp_path / 'out'), *argv])
    paths, detail = reference_paths(tmp_path / 'out')
    assert paths['MASH sketch'].endswith('references_extended.msh')
    assert paths['Reference table'].endswith('references_extended.tsv')
    assert '432 genome(s)' in detail


def test_main_with_the_vsnp_reference_set(tmp_path, monkeypatch):
    fake_run(monkeypatch, sketch_ids=['NC_006932-NC_006933.fasta'])
    reads = write_fastq(tmp_path / 'a_R1.fastq.gz', [('r', 'ACGT')])
    cli.main(['-r1', str(reads), '-o', str(tmp_path / 'out'), '--reference-set', 'vsnp'])
    paths, detail = reference_paths(tmp_path / 'out')
    assert paths['MASH sketch'].endswith('references.msh')
    assert not paths['MASH sketch'].endswith('references_extended.msh')
    assert '16 genome(s)' in detail


@pytest.mark.parametrize('which', ['extended', 'vsnp'])
@pytest.mark.parametrize('extra', [['--sketch', 'mine.msh'], ['--references', 'mine.tsv']])
def test_main_rejects_a_packaged_set_next_to_a_custom_file(tmp_path, monkeypatch, capsys, which, extra):
    fake_run(monkeypatch)
    reads = write_fastq(tmp_path / 'a_R1.fastq.gz', [('r', 'ACGT')])
    with pytest.raises(SystemExit) as exit_info:
        cli.main(['-r1', str(reads), '-o', str(tmp_path / 'out'), '--reference-set', which, *extra])
    assert exit_info.value.code == 1
    assert 'cannot be combined' in capsys.readouterr().err


@pytest.mark.parametrize(('argv', 'warned'), [
    ([], True),                                  # the default set: epithet codes, not vSNP's
    (['--reference-set', 'extended'], True),
    (['--reference-set', 'vsnp'], False),         # vSNP's own codes, which is what the report was for
])
def test_main_warns_that_xlsx_codes_changed_with_the_default_set(tmp_path, monkeypatch, capsys, argv, warned):
    """The legacy report exists for pipelines reading ab1/mel1, and the default table has no such codes."""
    ids = ['NC_006932-NC_006933.fasta'] if 'vsnp' in argv else None
    fake_run(monkeypatch, sketch_ids=ids)
    reads = write_fastq(tmp_path / 'a_R1.fastq.gz', [('r', 'ACGT')])
    cli.main(['-r1', str(reads), '-o', str(tmp_path / 'out'), '--xlsx', '--no-mlst', *argv])
    assert ('Add --reference-set vsnp for the 1.x codes' in capsys.readouterr().err) is warned


def test_main_does_not_warn_about_xlsx_with_a_custom_table(tmp_path, monkeypatch, capsys):
    """A custom table's codes are the caller's own business, so there is nothing to warn about."""
    fake_run(monkeypatch, sketch_ids=['mine.fasta'])
    reads = write_fastq(tmp_path / 'a_R1.fastq.gz', [('r', 'ACGT')])
    table = tmp_path / 'mine.tsv'
    table.write_text('\t'.join(COLUMNS) + '\nmine.fasta\tab1\tBrucella abortus\t\t\t\n')
    cli.main(['-r1', str(reads), '-o', str(tmp_path / 'out'), '--xlsx', '--no-mlst',
              '--references', str(table)])
    assert 'Add --reference-set vsnp' not in capsys.readouterr().err


def test_main_rejects_a_sketch_the_table_cannot_name_when_the_st_is_off_too(tmp_path, monkeypatch, capsys):
    """
    A 1.x command line that passed --sketch alone used to pair with the 16-genome table.

    From 2.0 it pairs with the extended one, whose IDs are GCA file names, so nothing in the sketch can
    be named. With the sequence type switched off as well there is no species to be had from anywhere,
    and every sample would come back ND from a run that exited 0.
    """
    fake_run(monkeypatch, sketch_ids=['NC_006932-NC_006933.fasta'])
    reads = write_fastq(tmp_path / 'a_R1.fastq.gz', [('r', 'ACGT')])
    sketch = tmp_path / 'mine.msh'
    sketch.write_bytes(b'')
    with pytest.raises(SystemExit) as exit_info:
        cli.main(['-r1', str(reads), '-o', str(tmp_path / 'out'), '--sketch', str(sketch), '--no-mlst'])
    assert exit_info.value.code == 1
    error = capsys.readouterr().err
    assert 'is in the reference table' in error
    assert '--reference-set vsnp' in error


@pytest.mark.parametrize('off', [['--no-mlst'], ['--no-st-species']])
def test_main_rejects_it_for_either_way_of_switching_the_st_off(tmp_path, monkeypatch, capsys, off):
    fake_run(monkeypatch, sketch_ids=['NC_006932-NC_006933.fasta'])
    reads = write_fastq(tmp_path / 'a_R1.fastq.gz', [('r', 'ACGT')])
    with pytest.raises(SystemExit):
        cli.main(['-r1', str(reads), '-o', str(tmp_path / 'out'), '--sketch', 'mine.msh', *off])
    assert 'no sample could be given a species at all' in capsys.readouterr().err


def test_main_allows_a_sketch_of_your_own_genomes_when_the_st_can_name_them(tmp_path, monkeypatch, capsys):
    """
    A sketch of your own assemblies is a supported way to run: the species comes from the ST.

    species.call_species has a branch for exactly this -- "a custom sketch whose closest genome the
    reference table does not name" -- so the guard must warn rather than refuse while an ST table is
    loaded, or that branch is unreachable from the command line.
    """
    fake_run(monkeypatch, sketch_ids=['my_assembly_1.fna.gz', 'my_assembly_2.fna.gz'])
    reads = write_fastq(tmp_path / 'a_R1.fastq.gz', [('r', 'ACGT')])
    cli.main(['-r1', str(reads), '-o', str(tmp_path / 'out'), '--sketch', 'mine.msh'])
    error = capsys.readouterr().err
    assert 'the distance has no species to offer' in error
    assert 'species of its sequence type' in error


def test_main_takes_a_custom_sketch_without_a_reference_set(tmp_path, monkeypatch):
    """--sketch alone must not trip the guard: no --reference-set was asked for."""
    fake_run(monkeypatch, sketch_ids=['mine.fasta'])
    reads = write_fastq(tmp_path / 'a_R1.fastq.gz', [('r', 'ACGT')])
    sketch = tmp_path / 'mine.msh'
    sketch.write_bytes(b'')
    table = tmp_path / 'mine.tsv'
    table.write_text('\t'.join(COLUMNS) + '\nmine.fasta\tx\tBrucella abortus\t\t\t\n')
    cli.main(['-r1', str(reads), '-o', str(tmp_path / 'out'), '--sketch', str(sketch),
              '--references', str(table)])
    paths, _ = reference_paths(tmp_path / 'out')
    assert paths['MASH sketch'] == str(sketch.resolve())
