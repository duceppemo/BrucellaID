"""Write the reports of a run."""

import json
import logging
import platform
import re
import sys
from datetime import datetime
from pathlib import Path

from . import __version__, htmlreport, qa
from .errors import BrucellaIDError
from .mash import versions as mash_versions
from .mlst import UNKNOWN
from .references import lookup
from .seal import versions as seal_versions

log = logging.getLogger(__name__)

MAIN_FILE = 'brucellaid.tsv'
MLST_FILE = 'brucellaid_mlst.tsv'
MASH_FILE = 'brucellaid_mash.tsv'
JSON_FILE = 'brucellaid.json'
MQC_FILE = 'brucellaid_mqc.json'
PDF_FILE = 'brucellaid_report.pdf'

MAIN_HEADER = ['Sample', 'Species', 'SpeciesSource', 'ClosestReference', 'MashDistance',
               'SharedHashes', 'SequenceType', 'AlleleProfile', 'STSpecies', 'Reads', 'Notes', 'Error']
MLST_HEADER = ['Sample', 'Locus', 'Allele', 'Reads', 'RunnerUp', 'RunnerUpReads', 'Note']
MASH_HEADER = ['Sample', 'Reference', 'Species', 'Biovar', 'MashDistance', 'PValue', 'SharedHashes', 'TotalHashes']

# The name and the columns of the report the first version of BrucellaID wrote, kept for the pipelines that read it
LEGACY_PREFIX = 'stat_alignment_summary'
LEGACY_HEADER = ['time_stamp', 'sample_name', 'species', 'mlst_type']


def one_cell(text):
    """
    Text that cannot break the table it is written into: no tab, no line break.

    Seal's error output is tab aligned and several lines long, and a file name may hold anything, so neither can
    be written to a tab separated file as it comes.
    """
    return re.sub(r'[\t\r\n]+', ' ', str(text)).strip()


def species_source(call):
    """The SpeciesSource cell: "MASH", or "MASH+ST (agree)" when the two kinds of evidence had a tie to break."""
    if call is None or not call.source:
        return ''
    return '{} ({})'.format(call.source, call.verdict) if call.verdict else call.source


def main_row(result):
    """
    One row of the main table: what one sample is.

    "SpeciesSource" says which evidence the species rests on -- MASH, the ST, or both agreeing -- and, when
    the two had something to settle, how it was settled: "MASH+ST (agree)", "ST (arbitrated)",
    "MASH (conflict)". "STSpecies" is the species the ST stands for, empty when the ST was not usable; the
    note says why.
    """
    best = result.identification.best if result.identification else None
    mlst = result.mlst
    call = result.call
    return [result.sample,
            result.species,
            species_source(call),
            result.closest_reference,
            '{:.6f}'.format(best.distance) if best else '',
            '{}/{}'.format(best.shared_hashes, best.total_hashes) if best else '',
            mlst.st if mlst else UNKNOWN,
            mlst.profile if mlst else '',
            call.st_species if call else '',
            '' if mlst is None or mlst.total_reads is None else str(mlst.total_reads),
            one_cell(result.note),
            one_cell(result.error or '')]


def main_table(results):
    """The table written to brucellaid.tsv and printed: one row per sample."""
    return [MAIN_HEADER] + [[one_cell(cell) for cell in main_row(result)] for result in results]


def write_table(rows, path, command=None):
    """Write a table, tab separated, with the command line as a first comment line."""
    with open(path, 'w', encoding='utf-8') as f:
        if command:
            f.write('#{}\n'.format(command))
        for row in rows:
            f.write('\t'.join(row) + '\n')
    return rows


def mlst_table(results):
    """One row per sample and locus: the allele, and what the call rests on."""
    rows = [MLST_HEADER]
    for result in results:
        if result.mlst is None:
            continue
        for call in result.mlst.calls:
            rows.append([result.sample, call.locus, call.allele or UNKNOWN, str(call.reads),
                         call.runner_up or '', str(call.runner_up_reads), one_cell(call.note)])
    return rows


def mash_table(results, references):
    """One row per sample and reference genome: every MASH distance, closest reference first."""
    rows = [MASH_HEADER]
    for result in results:
        if result.identification is None:
            continue
        for hit in result.identification.hits:
            reference = lookup(references, hit.reference)
            rows.append([result.sample, reference.reference if reference else hit.reference,
                         reference.species if reference else '',
                         reference.biovar if reference else '',
                         '{:.6f}'.format(hit.distance), '{:.3g}'.format(hit.p_value),
                         str(hit.shared_hashes), str(hit.total_hashes)])
    return rows


def json_report(results, references, run=None, inputs=None):
    """Everything the run produced, for a pipeline: parameters, versions, and the result of each sample."""
    samples = []
    for result in results:
        best = result.identification.best if result.identification else None
        reference = result.reference
        samples.append({
            'sample': result.sample,
            'files': result.files,
            'species': result.species,
            'species_code': result.code,
            'identified': result.identified,
            'species_call': None if result.call is None else {
                'species': result.call.species,
                'source': result.call.source,
                'verdict': result.call.verdict,
                'mash_candidates': list(result.call.candidates),
                'sequence_type': result.call.st,
                'st_species': result.call.st_species,
                'st_purity': result.call.entry.purity if result.call.entry else None,
                'st_isolates': result.call.entry.isolates if result.call.entry else None,
                'st_not_used': result.call.st_reason or None,
                'notes': list(result.call.notes),
            },
            'closest_reference': None if best is None else {
                'reference': result.closest_reference,
                'sketch_id': best.reference,
                'description': reference.description if reference else None,
                'accessions': reference.accessions.split(';') if reference and reference.accessions else [],
                'mash_distance': best.distance,
                'p_value': best.p_value,
                'shared_hashes': best.shared_hashes,
                'total_hashes': best.total_hashes,
                'identity_percent': best.identity,
                'min_kmer_copies': result.identification.min_copies,
                'estimated_coverage': result.identification.coverage,
            },
            'close_matches': [] if result.identification is None else [
                {'reference': hit.reference,
                 'description': lookup(references, hit.reference).description
                 if lookup(references, hit.reference) else None,
                 'mash_distance': hit.distance, 'shared_hashes': hit.shared_hashes}
                for hit in result.identification.close_matches],
            'mlst': None if result.mlst is None else {
                'sequence_type': result.mlst.st,
                'profile': result.mlst.profile,
                'alleles': result.mlst.alleles,
                'complete': result.mlst.complete,
                'assigned_reads': result.mlst.assigned_reads,
                'reads': result.mlst.total_reads,
                'loci': [{'locus': call.locus, 'allele': call.allele, 'reads': call.reads,
                          'runner_up': call.runner_up, 'runner_up_reads': call.runner_up_reads,
                          'note': call.note}
                         for call in result.mlst.calls],
                'note': result.mlst.note,
            },
            'notes': list(result.notes),
            'status': qa.status(result),
            'input_files': [{'path': f.path, 'size': f.size, 'modified': f.modified, 'md5': f.md5 or None,
                             'symlink_to': f.target or None}
                            for f in (inputs or {}).get(result.sample, [])],
            'error': result.error,
        })
    report = {'tool': 'brucellaid', 'version': __version__,
              'references': [{'reference': r.reference, 'code': r.code, 'species': r.species, 'biovar': r.biovar,
                              'strain': r.strain, 'accessions': r.accessions.split(';') if r.accessions else []}
                             for r in references.values()],
              'samples': samples}
    report.update(run or {})
    return report


def write_json(results, references, path, run=None, inputs=None):
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(json_report(results, references, run=run, inputs=inputs), f, indent=2, sort_keys=False)
        f.write('\n')


MQC_COLUMNS = ('Species', 'SpeciesSource', 'STSpecies', 'SequenceType', 'SharedHashes',
               'ClosestReference', 'Status')


def write_multiqc(results, path):
    """
    A MultiQC custom content file: run "multiqc" on the output folder and the run appears as a
    "BrucellaID" table.

    JSON rather than TSV, because a sequence type is a label and not a number: MultiQC would right-align
    "20" and happily average a column of them.
    """
    data = {}
    for result in results:
        values = (result.species, species_source(result.call),
                  result.call.st_species if result.call else '', result.sequence_type,
                  '{}/{}'.format(result.identification.best.shared_hashes,
                                 result.identification.best.total_hashes) if result.identification else '',
                  result.closest_reference, qa.status(result))
        data[result.sample] = {column: value or '-' for column, value in zip(MQC_COLUMNS, values, strict=True)}
    link = '<a href="https://github.com/duceppemo/BrucellaID">BrucellaID</a> {}'.format(__version__)
    content = {
        'id': 'brucellaid',
        'section_name': 'BrucellaID',
        'description': 'Brucella species from the MASH distance to a set of reference genomes, and the nine-locus '
                       'sequence type, by {}. "SpeciesSource" says whether the species rests on the distance, '
                       'on the sequence type, or on both.'.format(link),
        'plot_type': 'table',
        'pconfig': {'id': 'brucellaid_table', 'namespace': 'BrucellaID', 'title': 'BrucellaID'},
        'headers': {column: {'title': column, 'scale': False} for column in MQC_COLUMNS},
        'data': data,
    }
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(content, f, indent=2)
        f.write('\n')
    return path


def legacy_filename(started):
    """The name the first version of BrucellaID gave its report: stat_alignment_summary_2026-10-31_09-15-00.xlsx."""
    return '{}_{:%Y-%m-%d_%H-%M-%S}.xlsx'.format(LEGACY_PREFIX, started)


def write_legacy_xlsx(results, path, started):
    """
    The xlsx report of the first version of BrucellaID: the time stamp, the sample, the vSNP species code and
    the sequence type, in Courier New 10, written with xlsxwriter.
    """
    try:
        import xlsxwriter
    except ImportError:
        raise BrucellaIDError('--xlsx needs the xlsxwriter package: install it with "pip install xlsxwriter" '
                              '(or "pip install brucellaid[xlsx]").') from None
    stamp = '{:%Y-%m-%d_%H-%M-%S}'.format(started)
    workbook = xlsxwriter.Workbook(str(path))
    try:
        worksheet = workbook.add_worksheet()
        bold = workbook.add_format({'bold': True, 'font_name': 'Courier New', 'font_size': 10})
        courier = workbook.add_format({'font_name': 'Courier New', 'font_size': 10})
        courier.set_align('top')
        for column, header in enumerate(LEGACY_HEADER):
            worksheet.write(0, column, header, bold)
        for row, result in enumerate(results, start=1):
            for column, value in enumerate([stamp, result.sample, result.code, result.sequence_type]):
                worksheet.write(row, column, value, courier)
    finally:
        workbook.close()
    return path


def run_info(command, parameters, started, finished=None, operator=None, data=None):
    """
    The run metadata of the JSON report and of the quality assurance report: what was run, by whom, on what,
    with which versions and parameters, and against which reference data.

    :param operator: the name to put on the report, when it is not the user account's
    :param data: reference data facts from qa.data_facts(), identifying each file by its checksum
    """
    software = {'brucellaid': __version__, 'Python': platform.python_version()}
    software.update(mash_versions())
    software.update(seal_versions())
    info = {'command': command, 'started': started.isoformat(timespec='seconds'),
            'finished': (finished or datetime.now().astimezone()).isoformat(timespec='seconds'),
            'software': software, 'parameters': parameters}
    info.update(qa.environment(operator))
    info['reference_data'] = list(data or [])
    return info


def print_table(rows, stream=None):
    """Print the main table, tab separated: "brucellaid ... > results.tsv" works."""
    stream = stream or sys.stdout
    for row in rows:
        print('\t'.join(row), file=stream)
    stream.flush()


def write_reports(results, references, output, started, command=None, run=None, xlsx=False, html=True,
                  pdf=False, md5=True):
    """
    Write the reports in the output folder and return (main table, paths), the main table first.

    :param html: write the self-contained quality assurance report
    :param pdf: write the same report as a PDF, which needs reportlab
    :param md5: checksum the input files, for the reports
    """
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    inputs = qa.describe_inputs(results, md5=md5)
    rows = main_table(results)
    write_table(rows, output / MAIN_FILE, command=command)
    paths = [output / MAIN_FILE]
    mlst_rows = mlst_table(results)
    if len(mlst_rows) > 1:
        write_table(mlst_rows, output / MLST_FILE, command=command)
        paths.append(output / MLST_FILE)
    mash_rows = mash_table(results, references)
    if len(mash_rows) > 1:
        write_table(mash_rows, output / MASH_FILE, command=command)
        paths.append(output / MASH_FILE)
    write_json(results, references, output / JSON_FILE, run=run, inputs=inputs)
    paths.append(output / JSON_FILE)
    write_multiqc(results, output / MQC_FILE)
    paths.append(output / MQC_FILE)
    if html:
        htmlreport.write(results, references, run or {}, output / htmlreport.REPORT_FILE, inputs=inputs,
                         data=(run or {}).get('reference_data'))
        paths.append(output / htmlreport.REPORT_FILE)
    if pdf:
        from . import pdfreport  # Imported here: it needs reportlab, which only --pdf requires
        pdfreport.write(results, references, run or {}, output / PDF_FILE, inputs=inputs,
                        data=(run or {}).get('reference_data'))
        paths.append(output / PDF_FILE)
    if xlsx:
        paths.append(write_legacy_xlsx(results, output / legacy_filename(started), started))
    return rows, paths
