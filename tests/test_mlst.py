"""The MLST allele calling and the sequence type lookup."""

import pytest

from brucellaid.errors import BrucellaIDError
from brucellaid.mlst import (
    NEW_PROFILE,
    UNKNOWN,
    AlleleCall,
    MlstResult,
    alleles_fasta,
    call_alleles,
    call_locus,
    counts_by_locus,
    load_profiles,
    lookup_st,
    parse_allele,
    parse_profiles,
    run_mlst,
)
from tests.conftest import read_fasta, requires_seal

HEADER = 'ST\tgap\taroA\n'


@pytest.mark.parametrize('name, expected', [
    ('gap_2', ('gap', '2')),
    ('int_hyp_12', ('int_hyp', '12')),
    ('gap', None),
    ('_2', None),
    ('gap_x', None),
])
def test_parse_allele(name, expected):
    assert parse_allele(name) == expected


def test_the_profile_table_of_the_package():
    loci, profiles = load_profiles()
    assert loci == ('gap', 'aroA', 'glk', 'dnaK', 'gyrB', 'trpE', 'cobQ', 'int_hyp', 'omp25')
    assert len(profiles) == 81
    assert profiles[('2', '1', '1', '2', '1', '3', '1', '1', '1')] == '1'
    assert load_profiles() == (loci, profiles)  # Read once


def test_every_allele_of_the_profile_table_is_in_the_allele_file():
    loci, profiles = load_profiles()
    available = set(read_fasta(alleles_fasta()))
    used = {'{}_{}'.format(locus, allele) for profile in profiles for locus, allele in zip(loci, profile, strict=True)}
    assert used <= available
    assert all(parse_allele(name) and parse_allele(name)[0] in loci for name in available)


def test_parse_profiles_ignores_comments_and_trailing_tabs():
    loci, profiles = parse_profiles(['# a comment\n', HEADER, '1\t2\t1\t\n', '\n', '2\t2\t2\n'])
    assert loci == ('gap', 'aroA')
    assert profiles == {('2', '1'): '1', ('2', '2'): '2'}


@pytest.mark.parametrize('lines, message', [
    (['sequence_type\tgap\n', '1\t2\n'], 'must start with a header'),
    (['ST\n'], 'must start with a header'),
    ([HEADER, '1\t2\n'], 'has 2 columns, and the header 3'),
    ([HEADER, '1\t2\t1\n', '2\t2\t1\n'], 'have the same alleles'),
    ([HEADER], 'holds no sequence type'),
])
def test_parse_profiles_errors(lines, message):
    with pytest.raises(BrucellaIDError, match=message):
        parse_profiles(lines)


def test_load_profiles_of_an_unreadable_file(tmp_path):
    with pytest.raises(BrucellaIDError, match='cannot be read'):
        load_profiles(str(tmp_path / 'nowhere.tsv'))


def test_counts_by_locus():
    assert counts_by_locus({'gap_1': 10, 'gap_2': 3, 'aroA_1': 5}, ('gap', 'aroA')) == {
        'gap': {'1': 10, '2': 3}, 'aroA': {'1': 5}}


@pytest.mark.parametrize('counts, message', [
    ({'gap_1': 1, 'rubbish': 0}, 'do not belong to a locus'),
    ({'gap_1': 1, 'glk_1': 2}, 'do not belong to a locus'),
    ({'gap_1': 1}, 'no allele for the locus'),
])
def test_counts_by_locus_errors(counts, message):
    with pytest.raises(BrucellaIDError, match=message):
        counts_by_locus(counts, ('gap', 'aroA'))


def test_call_locus_takes_the_allele_the_most_reads_match():
    call = call_locus('gap', {'1': 40, '2': 12, '3': 0})
    assert (call.allele, call.reads) == ('1', 40)
    assert (call.runner_up, call.runner_up_reads) == ('2', 12)
    assert call.note == ''
    assert call.name == 'gap_1'


def test_call_locus_reports_a_tie_rather_than_breaking_it():
    call = call_locus('gap', {'1': 20, '2': 20, '3': 1})
    assert call.allele is None
    assert call.name == ''
    assert 'tied at 20 read(s) between gap_1, gap_2' in call.note


def test_call_locus_below_the_minimum():
    assert 'no allele reached 5 read(s) (best: gap_2 with 4)' in call_locus('gap', {'1': 1, '2': 4}).note
    assert call_locus('gap', {'1': 0, '2': 0}, min_reads=1).note == 'no read matched this locus'


def test_call_locus_with_one_allele_in_the_scheme():
    call = call_locus('gap', {'1': 7})
    assert (call.allele, call.runner_up, call.runner_up_reads) == ('1', None, 0)


def test_call_alleles_follows_the_order_of_the_loci():
    calls = call_alleles({'gap_1': 9, 'aroA_2': 8, 'aroA_1': 1}, ('aroA', 'gap'))
    assert [call.locus for call in calls] == ['aroA', 'gap']


def test_lookup_st():
    profiles = {('2', '1'): '1'}
    known = [AlleleCall('gap', '2', 9), AlleleCall('aroA', '1', 9)]
    assert lookup_st(known, profiles) == ('1', '')
    unknown = [AlleleCall('gap', '7', 9), AlleleCall('aroA', '1', 9)]
    assert lookup_st(unknown, profiles) == (NEW_PROFILE, 'this allele profile is not in the profile table')
    missing = [AlleleCall('gap', None), AlleleCall('aroA', '1', 9)]
    assert lookup_st(missing, profiles) == (UNKNOWN, 'no allele at gap')


def test_mlst_result_properties():
    result = MlstResult(st='1', calls=[AlleleCall('gap', '2', 9), AlleleCall('aroA', None)])
    assert result.alleles == {'gap': '2', 'aroA': None}
    assert result.profile == '2--'
    assert not result.complete
    assert MlstResult().st == UNKNOWN
    assert not MlstResult().complete


@requires_seal
@pytest.mark.parametrize('st', ['1', '8'])
def test_run_mlst_types_the_reads_of_a_known_sequence_type(st, st_reads, mlst_loci, mlst_profiles):
    result = run_mlst([str(st_reads(st))], threads=2)
    assert result.st == st
    assert result.complete
    assert result.note == ''
    assert tuple(result.alleles[locus] for locus in mlst_loci) in mlst_profiles
    assert result.total_reads > 100
    assert result.assigned_reads == result.total_reads  # Every read of the input comes from an allele


@requires_seal
def test_run_mlst_of_reads_that_hold_no_allele(tmp_path):
    from tests.conftest import random_reads, random_sequence, write_fastq
    reads = write_fastq(tmp_path / 'other_R1.fastq.gz', random_reads(random_sequence(5000, seed=4), 200, seed=2))
    result = run_mlst([str(reads)], threads=2)
    assert result.st == UNKNOWN
    assert result.profile == '-' * 17
    assert 'no allele at gap' in result.note
