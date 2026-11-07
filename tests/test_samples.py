"""Finding the samples of a folder, and naming them."""

import pytest

from brucellaid.errors import BrucellaIDError
from brucellaid.samples import Sample, find_samples, sample_name, single_sample, split_extension


def touch(folder, *names):
    folder.mkdir(parents=True, exist_ok=True)
    for name in names:
        (folder / name).write_text('@r\nACGT\n+\nIIII\n')
    return folder


@pytest.mark.parametrize('name, expected', [
    ('S1.fastq', ('S1', 'fastq')),
    ('S1.fastq.gz', ('S1', 'fastq')),
    ('S1.FQ.GZ', ('S1', 'fastq')),
    ('S1.fq', ('S1', 'fastq')),
    ('genome.fasta', ('genome', 'fasta')),
    ('genome.fna.gz', ('genome', 'fasta')),
    ('notes.txt', None),
    ('.fastq', None),  # Nothing left of the name
    ('reads', None),
])
def test_split_extension(name, expected):
    assert split_extension(name) == expected


@pytest.mark.parametrize('name, expected', [
    ('S1_R1.fastq.gz', 'S1'),
    ('S1_S12_L001_R2_001.fastq.gz', 'S1_S12_L001'),
    ('SRR123_1.fq', 'SRR123'),
    ('sample.fastq', 'sample'),
    ('_R1.fastq', '_R1'),  # The whole name is the read suffix: it is kept
    ('weird\tname_R1.fastq', 'weird_name'),
    ('reads_R1', 'reads'),  # No extension to remove
])
def test_sample_name(name, expected):
    assert sample_name(name) == expected


def test_find_samples_pairs_and_singles(tmp_path):
    reads = touch(tmp_path / 'reads', 'A_R1.fastq.gz', 'A_R2.fastq.gz', 'B.fastq', 'C_1.fq.gz', 'C_2.fq.gz')
    touch(reads / 'sub', 'D_R1.fastq')  # Subfolders are searched
    samples = find_samples(reads)
    assert [s.name for s in samples] == ['A', 'B', 'C', 'D']
    assert [s.paired for s in samples] == [True, False, True, False]
    assert [f.endswith('A_R1.fastq.gz') for f in samples[0].files] == [True, False]


def test_find_samples_skips_what_is_not_a_regular_file(tmp_path, caplog):
    """
    A FIFO blocks the run forever on open(), and a link to a device reads without end.

    single_sample() has always refused anything that is not a regular file; a folder of samples did not,
    so one such file in a folder took the whole run down or hung it, losing every other sample.
    """
    import os
    reads = touch(tmp_path / 'reads', 'A_R1.fastq.gz')
    os.mkfifo(reads / 'pipe_R1.fastq')
    (reads / 'dir_R1.fastq').mkdir()
    with caplog.at_level('WARNING'):
        samples = find_samples(reads)
    assert [s.name for s in samples] == ['A']
    assert 'not a regular file, skipped' in caplog.text
    assert 'pipe_R1.fastq' in caplog.text


@pytest.mark.parametrize(('filename', 'expected'), [
    ('   _R1.fastq.gz', '_R1.fastq.gz'),   # Nothing but spaces is not a name
    (' spaced _R1.fastq.gz', 'spaced'),    # Surrounding spaces would differ between the TSV and the JSON
])
def test_find_samples_never_names_a_sample_blank(tmp_path, filename, expected):
    reads = touch(tmp_path / 'reads', filename)
    assert [s.name for s in find_samples(reads)] == [expected]


def test_find_samples_skips_hidden_and_excluded(tmp_path):
    reads = touch(tmp_path / 'reads', 'A_R1.fastq.gz', '.hidden.fastq')
    touch(reads / '.git', 'B_R1.fastq')
    touch(reads / 'out', 'C_R1.fastq')
    samples = find_samples(reads, exclude=[reads / 'out'])
    assert [s.name for s in samples] == ['A']


def test_find_samples_warns_about_an_r2_without_an_r1(tmp_path, caplog):
    reads = touch(tmp_path / 'reads', 'A_R2.fastq.gz')
    samples = find_samples(reads)
    assert [s.name for s in samples] == ['A']
    assert 'R2 file without R1' in caplog.text


def test_find_samples_rejects_a_name_collision(tmp_path):
    reads = touch(tmp_path / 'reads', 'A_R1.fastq.gz', 'A_R1.fq')
    with pytest.raises(BrucellaIDError, match='same sample name'):
        find_samples(reads)


def test_find_samples_rejects_three_files_for_one_sample(tmp_path):
    reads = touch(tmp_path / 'reads', 'A_R1.fastq', 'A_R2.fastq', 'A_R1.fastq.gz')
    with pytest.raises(BrucellaIDError, match='same sample name'):
        find_samples(reads)


def test_find_samples_skips_assemblies_with_a_warning(tmp_path, caplog):
    reads = touch(tmp_path / 'reads', 'A_R1.fastq.gz', 'assembly.fasta')
    samples = find_samples(reads)
    assert [s.name for s in samples] == ['A']
    assert 'fasta file(s)' in caplog.text


def test_find_samples_without_a_fastq_file(tmp_path):
    reads = touch(tmp_path / 'reads', 'assembly.fasta')
    with pytest.raises(BrucellaIDError, match='No fastq file found'):
        find_samples(reads)


def test_find_samples_in_a_missing_folder(tmp_path):
    with pytest.raises(BrucellaIDError, match='Input folder not found'):
        find_samples(tmp_path / 'nowhere')


def test_find_samples_follows_a_link_once(tmp_path):
    reads = touch(tmp_path / 'reads', 'A_R1.fastq.gz')
    (reads / 'loop').symlink_to(reads)  # A folder linking to itself is searched once, not forever
    assert [s.name for s in find_samples(reads)] == ['A']


def test_single_sample(tmp_path):
    reads = touch(tmp_path, 'A_R1.fastq.gz', 'A_R2.fastq.gz')
    sample = single_sample(str(reads / 'A_R1.fastq.gz'), str(reads / 'A_R2.fastq.gz'))
    assert sample == Sample('A', [str(reads / 'A_R1.fastq.gz'), str(reads / 'A_R2.fastq.gz')])
    assert single_sample(str(reads / 'A_R1.fastq.gz'), name='other').name == 'other'


def test_single_sample_errors(tmp_path):
    touch(tmp_path, 'A_R1.fastq.gz', 'assembly.fasta')
    with pytest.raises(BrucellaIDError, match='Input file not found'):
        single_sample(str(tmp_path / 'missing.fastq'))
    with pytest.raises(BrucellaIDError, match='is a fasta file'):
        single_sample(str(tmp_path / 'assembly.fasta'))
    with pytest.raises(BrucellaIDError, match='the same file'):
        single_sample(str(tmp_path / 'A_R1.fastq.gz'), str(tmp_path / 'A_R1.fastq.gz'))
