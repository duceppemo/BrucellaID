"""End to end runs of the command line, with MASH and Seal."""

import json

import pytest

from brucellaid import cli, report
from tests.conftest import requires_mash, requires_seal, write_fastq


def read_table(path):
    """[rows] of a report, without the "#command line" comment."""
    rows = [line.split('\t') for line in path.read_text().splitlines() if not line.startswith('#')]
    return rows


def columns(rows):
    """[{column: value}] of a table whose first row is the header."""
    return [dict(zip(rows[0], row, strict=True)) for row in rows[1:]]


@requires_mash
@requires_seal
def test_one_sample_identified_and_typed(tmp_path, st_reads, sketch, genome, genome_reads):
    """A sample holding the reads of the reference genome and of the ST 1 alleles is identified and typed."""
    msh, table = sketch
    _, sequence = genome
    reads = st_reads('1')
    combined = tmp_path / 'sample_R1.fastq.gz'
    with open(reads, 'rb') as f1, open(genome_reads, 'rb') as f2, open(combined, 'wb') as out:
        out.write(f1.read())
        out.write(f2.read())
    output = tmp_path / 'out'
    cli.main(['-r1', str(combined), '-o', str(output), '--sketch', str(msh), '--references', str(table),
              '-t', '2'])

    rows = columns(read_table(output / report.MAIN_FILE))
    assert len(rows) == 1
    row = rows[0]
    assert row['Sample'] == 'sample'
    assert row['Species'] == 'Brucella testii'
    assert 'Biovar' not in row  # BrucellaID does not claim a biovar for a sample
    assert row['ClosestReference'] == 'test_genome.fasta'
    assert row['SequenceType'] == '1'
    assert row['AlleleProfile'] == '2-1-1-2-1-3-1-1-1'
    assert row['Error'] == ''
    assert int(row['SharedHashes'].split('/')[0]) >= 150

    loci = columns(read_table(output / report.MLST_FILE))
    assert [locus['Locus'] for locus in loci] == ['gap', 'aroA', 'glk', 'dnaK', 'gyrB', 'trpE', 'cobQ', 'int_hyp',
                                                  'omp25']
    assert [locus['Allele'] for locus in loci] == ['2', '1', '1', '2', '1', '3', '1', '1', '1']
    assert all(int(locus['Reads']) > int(locus['RunnerUpReads']) for locus in loci)

    hits = columns(read_table(output / report.MASH_FILE))
    assert hits[0]['Reference'] == 'test_genome.fasta'

    data = json.loads((output / report.JSON_FILE).read_text())
    assert data['tool'] == 'brucellaid'
    assert data['samples'][0]['mlst']['sequence_type'] == '1'
    assert data['samples'][0]['closest_reference']['description'] == 'Brucella testii bv. 9 str. T1'
    assert data['samples'][0]['mlst']['complete'] is True
    assert data['software']['MASH'] != 'unknown'
    assert data['parameters']['Threads'] == 2


@requires_mash
@requires_seal
def test_folder_of_samples_paired_and_single(tmp_path, st_reads, sketch, mlst_alleles, mlst_loci, mlst_profiles):
    """A folder gives one row per sample, paired-end and single-end alike."""
    msh, table = sketch
    reads = tmp_path / 'reads'
    st_reads('1', name='paired', folder=reads)
    (reads / 'paired_R1.fastq.gz').rename(reads / 'paired_R1.fastq.gz')
    # The R2 file of the pair holds the same reads: both are the same read set to Seal and MASH
    st_reads('1', name='mate', folder=reads)
    (reads / 'mate_R1.fastq.gz').rename(reads / 'paired_R2.fastq.gz')
    st_reads('6', name='single', folder=reads)
    output = tmp_path / 'out'
    cli.main(['-i', str(reads), '-o', str(output), '--sketch', str(msh), '--references', str(table), '-p', '2'])

    rows = {row['Sample']: row for row in columns(read_table(output / report.MAIN_FILE))}
    assert sorted(rows) == ['paired', 'single']
    assert rows['paired']['SequenceType'] == '1'
    assert rows['single']['SequenceType'] == '6'
    assert rows['paired']['Species'] == 'ND'  # Allele reads only: no genome, so no species


@requires_mash
@requires_seal
def test_no_mlst_and_legacy_xlsx(tmp_path, sketch, genome_reads):
    """--no-mlst identifies the species only, and --xlsx writes the report of the first version."""
    msh, table = sketch
    output = tmp_path / 'out'
    cli.main(['-r1', str(genome_reads), '-o', str(output), '--sketch', str(msh), '--references', str(table),
              '--no-mlst', '--xlsx', '-s', 'legacy'])
    rows = columns(read_table(output / report.MAIN_FILE))
    assert rows[0]['Species'] == 'Brucella testii'
    assert rows[0]['SequenceType'] == 'ND'
    assert not (output / report.MLST_FILE).exists()
    xlsx = list(output.glob('{}_*.xlsx'.format(report.LEGACY_PREFIX)))
    assert len(xlsx) == 1

    openpyxl = pytest.importorskip('openpyxl', reason='openpyxl reads back the xlsx report')
    sheet = openpyxl.load_workbook(xlsx[0]).active
    values = [[cell.value for cell in row] for row in sheet.iter_rows()]
    assert values[0] == report.LEGACY_HEADER
    assert values[1][1:] == ['legacy', 'test', 'ND']


@requires_mash
@requires_seal
def test_reads_that_are_not_brucella(tmp_path, sketch):
    """A read set sharing nothing with the references is reported as not identified, and the run still succeeds."""
    from tests.conftest import random_reads, random_sequence
    msh, table = sketch
    unrelated = random_sequence(20000, seed=99)
    reads = write_fastq(tmp_path / 'other_R1.fastq.gz', random_reads(unrelated, 2000, seed=5))
    output = tmp_path / 'out'
    cli.main(['-r1', str(reads), '-o', str(output), '--sketch', str(msh), '--references', str(table), '--no-mlst'])
    row = columns(read_table(output / report.MAIN_FILE))[0]
    assert row['Species'] == 'ND'
    assert 'not identified' in row['Notes']
    assert row['Error'] == ''


@requires_mash
@requires_seal
def test_failed_sample_exits_with_an_error(tmp_path, sketch):
    """A sample whose reads cannot be read fails on its own, and the run exits with 1."""
    msh, table = sketch
    reads = tmp_path / 'broken_R1.fastq.gz'
    reads.write_bytes(b'\x1f\x8b\x08\x00not a gzip file')
    output = tmp_path / 'out'
    with pytest.raises(SystemExit) as exit_info:
        cli.main(['-r1', str(reads), '-o', str(output), '--sketch', str(msh), '--references', str(table),
                  '--no-mlst'])
    assert exit_info.value.code == 1
    row = columns(read_table(output / report.MAIN_FILE))[0]
    assert row['Error']


@requires_mash
@pytest.mark.parametrize(('which', 'genomes'), [('EXTENDED', 432), ('VSNP', 16)])
def test_the_sketch_and_the_reference_table_of_the_package_match(which, genomes):
    """Every genome of each packaged sketch is named in its reference table, and the other way round."""
    from brucellaid import references as refs
    from brucellaid.mash import sketch_references
    sketch = refs.data_path(getattr(refs, '{}_SKETCH_FILE'.format(which)))
    table = refs.data_path(getattr(refs, '{}_REFERENCES_FILE'.format(which)))
    ids = sketch_references(refs.reference_sketch(sketch))
    assert len(ids) == genomes
    assert refs.compare_with_sketch(ids, refs.load_references(table)) == []
