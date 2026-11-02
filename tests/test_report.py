"""The reports: the tables, the JSON report and the legacy xlsx report."""

import json
import os
import subprocess
import sys
from datetime import datetime

import pytest

from brucellaid import htmlreport, report
from brucellaid.errors import BrucellaIDError
from brucellaid.mash import Identification, parse_dist
from brucellaid.mlst import UNKNOWN, AlleleCall, MlstResult
from brucellaid.pipeline import SampleResult, species_notes
from brucellaid.references import COLUMNS, parse_references
from brucellaid.report import (
    LEGACY_HEADER,
    MAIN_HEADER,
    MLST_HEADER,
    legacy_filename,
    main_table,
    mash_table,
    mlst_table,
    one_cell,
    print_table,
    run_info,
    write_legacy_xlsx,
    write_reports,
)
from brucellaid.species import call_species, load_st_species

REFERENCES = parse_references([
    '\t'.join(COLUMNS) + '\n',
    'NC_006932-NC_006933.fasta\tab1\tBrucella abortus\t1\t9-941\tNC_006932.1;NC_006933.1\n',
    'NZ_CP007682-NZ_CP007683.fasta\tab3\tBrucella abortus\t3\tBER\n',
])
STARTED = datetime(2026, 10, 31, 9, 15, 0)


def good_sample():
    found = Identification(parse_dist('NC_006932-NC_006933.fasta\tq\t0.0012345\t0\t950/1000\n'
                                      'NZ_CP007682-NZ_CP007683.fasta\tq\t0.0013\t0\t940/1000\n'))
    reference, notes = species_notes(found, REFERENCES)
    mlst = MlstResult(st='1', calls=[AlleleCall('gap', '2', 40, 'gap_3', 12), AlleleCall('aroA', '1', 38, '2', 10)],
                      assigned_reads=900, total_reads=1000)
    call = call_species(found, REFERENCES, mlst=mlst, table=load_st_species())
    return SampleResult('good', ['good_R1.fastq.gz', 'good_R2.fastq.gz'], identification=found,
                        reference=reference, mlst=mlst, notes=notes, call=call)


def failed_sample():
    return SampleResult('failed', ['failed_R1.fastq.gz'], error='Seal failed (exit code 1):\tno\noutput')


def test_one_cell_keeps_a_cell_on_one_line():
    assert one_cell('a\tb\nc\r\nd ') == 'a b c d'
    assert one_cell(42) == '42'


def test_main_table():
    rows = main_table([good_sample(), failed_sample()])
    assert rows[0] == MAIN_HEADER
    good = dict(zip(MAIN_HEADER, rows[1], strict=True))
    assert good['Species'] == 'Brucella abortus'
    assert 'Biovar' not in good
    # ST 1 is B. abortus in the isolate database too, and MASH had no second species to choose from
    assert good['SpeciesSource'] == 'MASH+ST (agree)'
    assert good['STSpecies'] == 'Brucella abortus'
    assert good['ClosestReference'] == 'NC_006932-NC_006933.fasta'
    assert good['MashDistance'] == '0.001234'
    assert good['SharedHashes'] == '950/1000'
    assert good['SequenceType'] == '1'
    assert good['AlleleProfile'] == '2-1'
    assert good['Reads'] == '1000'
    assert 'nearly as close' in good['Notes']
    failed = dict(zip(MAIN_HEADER, rows[2], strict=True))
    assert failed['Species'] == UNKNOWN
    assert failed['SequenceType'] == UNKNOWN
    assert '\t' not in failed['Error'] and '\n' not in failed['Error']


def test_mlst_table_leaves_out_the_samples_without_mlst():
    rows = mlst_table([good_sample(), failed_sample()])
    assert rows[0] == MLST_HEADER
    assert len(rows) == 3
    assert rows[1] == ['good', 'gap', '2', '40', 'gap_3', '12', '']


def test_mlst_table_of_a_locus_without_an_allele():
    result = SampleResult('s', [], mlst=MlstResult(calls=[AlleleCall('gap', None, note='tied')]))
    assert mlst_table([result])[1] == ['s', 'gap', UNKNOWN, '0', '', '0', 'tied']


def test_mash_table_names_the_species_of_every_reference():
    rows = mash_table([good_sample(), failed_sample()], REFERENCES)
    assert rows[1] == ['good', 'NC_006932-NC_006933.fasta', 'Brucella abortus', '1', '0.001234', '0', '950', '1000']
    assert rows[2][1] == 'NZ_CP007682-NZ_CP007683.fasta'
    assert len(rows) == 3  # The failed sample has no MASH result


def test_write_reports(tmp_path):
    results = [good_sample(), failed_sample()]
    run = run_info('brucellaid -i reads/ -o out/', {'Threads': 4}, STARTED)
    rows, paths = write_reports(results, REFERENCES, tmp_path / 'out', STARTED, command='brucellaid -i reads/',
                               run=run)
    assert [path.name for path in paths] == [report.MAIN_FILE, report.MLST_FILE, report.MASH_FILE,
                                            report.JSON_FILE, report.MQC_FILE, htmlreport.REPORT_FILE]
    assert rows[0] == MAIN_HEADER
    assert (tmp_path / 'out' / report.MAIN_FILE).read_text().startswith('#brucellaid -i reads/\n')

    data = json.loads((tmp_path / 'out' / report.JSON_FILE).read_text())
    assert data['version'] and data['tool'] == 'brucellaid'
    assert data['command'] == 'brucellaid -i reads/ -o out/'
    assert data['parameters'] == {'Threads': 4}
    assert len(data['references']) == 2
    assert data['references'][0]['accessions'] == ['NC_006932.1', 'NC_006933.1']
    assert data['references'][1]['accessions'] == []
    good, failed = data['samples']
    assert good['closest_reference']['accessions'] == ['NC_006932.1', 'NC_006933.1']
    assert good['closest_reference']['identity_percent'] == pytest.approx(99.87655)
    assert good['close_matches'][0]['description'] == 'Brucella abortus bv. 3 str. BER'
    assert good['mlst']['alleles'] == {'gap': '2', 'aroA': '1'}
    assert good['mlst']['loci'][0]['runner_up'] == 'gap_3'
    assert failed['mlst'] is None
    assert failed['closest_reference'] is None
    assert failed['error']


def test_write_reports_of_a_run_without_mlst(tmp_path):
    result = good_sample()
    result.mlst = None
    _, paths = write_reports([result], REFERENCES, tmp_path / 'out', STARTED)
    assert [path.name for path in paths] == [report.MAIN_FILE, report.MASH_FILE, report.JSON_FILE,
                                            report.MQC_FILE, htmlreport.REPORT_FILE]


def test_write_table_is_utf8_whatever_the_locale(tmp_path):
    """A sample name with an accent must not lose the reports of a finished run under LC_ALL=C."""
    out = tmp_path / 'r.tsv'
    code = ("from brucellaid.report import write_table; "
            "write_table([['Sample'], ['\\u00e9chantillon']], {!r})".format(str(out)))
    env = {**os.environ, 'LC_ALL': 'C', 'PYTHONCOERCECLOCALE': '0', 'PYTHONUTF8': '0'}
    env.pop('PYTHONIOENCODING', None)
    done = subprocess.run([sys.executable, '-c', code], env=env, capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
    assert '\u00e9chantillon' in out.read_bytes().decode('utf-8')


def test_run_info():
    info = run_info('brucellaid -h', {'Threads': 1}, STARTED, finished=datetime(2026, 10, 31, 9, 20, 0))
    assert info['started'] == '2026-10-31T09:15:00'
    assert info['finished'] == '2026-10-31T09:20:00'
    assert set(info['software']) == {'brucellaid', 'Python', 'MASH', 'BBTools', 'Java'}


def test_print_table(capsys):
    print_table([['a', 'b'], ['1', '2']])
    assert capsys.readouterr().out == 'a\tb\n1\t2\n'


def test_legacy_filename():
    assert legacy_filename(STARTED) == 'stat_alignment_summary_2026-10-31_09-15-00.xlsx'


def test_legacy_xlsx(tmp_path):
    pytest.importorskip('xlsxwriter')
    openpyxl = pytest.importorskip('openpyxl', reason='openpyxl reads the report back')
    path = write_legacy_xlsx([good_sample(), failed_sample()], tmp_path / legacy_filename(STARTED), STARTED)
    sheet = openpyxl.load_workbook(path).active
    values = [[cell.value for cell in row] for row in sheet.iter_rows()]
    assert values[0] == LEGACY_HEADER
    assert values[1] == ['2026-10-31_09-15-00', 'good', 'ab1', '1']
    assert values[2] == ['2026-10-31_09-15-00', 'failed', UNKNOWN, UNKNOWN]


def test_legacy_xlsx_without_xlsxwriter(tmp_path, monkeypatch):
    monkeypatch.setitem(__import__('sys').modules, 'xlsxwriter', None)
    with pytest.raises(BrucellaIDError, match='needs the xlsxwriter package'):
        write_legacy_xlsx([good_sample()], tmp_path / 'report.xlsx', STARTED)
