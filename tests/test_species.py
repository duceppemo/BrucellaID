"""The ST to species table, and the rule that decides the species from MASH and the ST together."""

import pytest

from brucellaid.mash import Identification, parse_dist
from brucellaid.mlst import UNKNOWN, MlstResult
from brucellaid.references import COLUMNS as REFERENCE_COLUMNS
from brucellaid.references import parse_references
from brucellaid.species import (
    AGREE,
    ARBITRATED,
    COLUMNS,
    CONFLICT,
    FROM_BOTH,
    FROM_MASH,
    FROM_ST,
    UNNAMED,
    UNRESOLVED,
    SpeciesError,
    StSpecies,
    call_species,
    candidates_of,
    closest_of,
    load_st_species,
    parse_st_species,
    st_evidence,
)

# B. canis and B. suis bv. 4 are the pair of the vSNP sketch MASH cannot separate (0.00024 apart); these
# fixtures are vSNP sketch IDs. B. abortus is far enough from both to be decided on its own. The default
# 432-genome set has the same canis/suis problem at 0.000264 and 13 other pairs besides; see species.py.
CANIS = 'NC_010103-NC_010104.fasta'
SUIS4 = 'B-REF-BS4-40.fasta'
ABORTUS = 'NC_006932-NC_006933.fasta'
REFERENCES = parse_references([
    '\t'.join(REFERENCE_COLUMNS) + '\n',
    '{}\tcanis\tBrucella canis\t\tATCC 23365\tNC_010103.1\n'.format(CANIS),
    '{}\tsuis4\tBrucella suis\t4\t40\t\n'.format(SUIS4),
    '{}\tab1\tBrucella abortus\t1\t9-941\tNC_006932.1\n'.format(ABORTUS),
])

TABLE = parse_st_species([
    '# a comment\n',
    '\t'.join(COLUMNS) + '\n',
    '20\tBrucella canis\t1.000\t70\t\t\n',
    '16\tBrucella suis\t1.000\t15\t2:8\t\n',
    '1\tBrucella abortus\t0.994\t466\t1:86;2:16\tBrucella melitensis:3\n',
    '17\tBrucella suis\t0.900\t20\t3:4\tBrucella melitensis:2\n',  # Too impure to be used
    '104\tBrucella suis\t1.000\t1\t2:1\t\n',  # Too few isolates to be used
    '50\t{}\t1.000\t16\t\t\n'.format(UNNAMED),  # No species named at all
])


def identification(*hits, min_hashes=150):
    """An Identification built from "<reference> <distance> <shared>/1000" triples, closest first."""
    lines = ''.join('{}\tq\t{}\t0\t{}/1000\n'.format(name, distance, shared) for name, distance, shared in hits)
    return Identification(parse_dist(lines), min_hashes=min_hashes)


def test_the_table_of_the_package_separates_canis_from_suis():
    """The reason the table exists: the ST decides the pair no sketch can."""
    table = load_st_species()
    canis = {entry.st for entry in table.values() if entry.species == 'Brucella canis'}
    suis = {entry.st for entry in table.values() if entry.species == 'Brucella suis'}
    assert canis and suis
    assert not canis & suis
    assert all(table[st].reliable() for st in canis)


def test_the_table_of_the_package_is_read():
    table = load_st_species()
    assert len(table) > 100
    assert table['20'].species == 'Brucella canis'
    assert table is load_st_species()  # Read once


def test_parse_st_species_reads_a_row():
    entry = TABLE['1']
    assert (entry.st, entry.species, entry.purity, entry.isolates) == ('1', 'Brucella abortus', 0.994, 466)
    assert entry.biovars == '1:86;2:16'
    assert entry.agreeing == 463
    assert 'ST 1 is Brucella abortus in 463 of 466 isolate(s)' in entry.evidence
    assert 'Brucella melitensis:3' in entry.evidence


def test_parse_st_species_allows_a_row_without_its_trailing_empty_columns():
    table = parse_st_species(['\t'.join(COLUMNS) + '\n', '20\tBrucella canis\t1.000\t70\n'])
    assert table['20'] == StSpecies('20', 'Brucella canis', 1.0, 70)


@pytest.mark.parametrize('lines, message', [
    (['ST\tSpecies\n', '20\tBrucella canis\n'], 'must start with the columns'),
    (['\t'.join(COLUMNS) + '\n', '20\tBrucella canis\t1.0\t70\t\t\textra\n'], 'more than the 6'),
    (['\t'.join(COLUMNS) + '\n', '\tBrucella canis\t1.0\t70\t\t\n'], 'no ST or no species'),
    (['\t'.join(COLUMNS) + '\n', '20\t\t1.0\t70\t\t\n'], 'no ST or no species'),
    (['\t'.join(COLUMNS) + '\n', '20\tBrucella canis\tmost\t70\t\t\n'], 'not a number'),
    (['\t'.join(COLUMNS) + '\n', '20\tBrucella canis\t1.4\t70\t\t\n'], 'outside 0-1'),
    (['\t'.join(COLUMNS) + '\n', '20\tBrucella canis\t1.0\t0\t\t\n'], 'outside 0-1 or no isolate'),
    (['\t'.join(COLUMNS) + '\n', '20\tBrucella canis\t1.0\t70\t\t\n', '20\tBrucella suis\t1.0\t5\t\t\n'],
     'appears twice'),
    (['\t'.join(COLUMNS) + '\n'], 'holds no ST'),
])
def test_parse_st_species_refuses_a_broken_table(lines, message):
    with pytest.raises(SpeciesError, match=message):
        parse_st_species(lines)


def test_load_st_species_reports_a_missing_file():
    with pytest.raises(SpeciesError, match='cannot be read'):
        load_st_species('/nonexistent/st_species.tsv')


@pytest.mark.parametrize('st, reason', [
    ('20', ''),
    ('16', ''),
    (UNKNOWN, 'has no ST'),
    ('new', 'no ST yet'),
    ('99', 'no isolate of ST 99'),
    ('17', 'below the 95% agreement'),
    ('104', 'rests on 1 isolate'),
    ('50', 'no species named'),
])
def test_st_evidence_says_when_an_st_cannot_be_used(st, reason):
    entry, why_not = st_evidence(MlstResult(st=st), TABLE)
    assert why_not == '' if not reason else reason in why_not
    assert (entry is not None) == (reason == '')


def test_st_evidence_without_a_table_or_without_the_mlst():
    assert st_evidence(MlstResult(st='20'), None) == (None, 'the ST to species table is not in use')
    assert st_evidence(None, TABLE) == (None, 'the MLST was not run')


def test_candidates_are_the_species_mash_leaves_within_reach():
    found = identification((CANIS, 0.000239, 990), (SUIS4, 0.000420, 980), (ABORTUS, 0.0012, 950))
    # B. suis bv. 4 is within CLOSE_MATCH_MARGIN of B. canis, B. abortus is not
    assert candidates_of(found, REFERENCES) == ('Brucella canis', 'Brucella suis')
    assert candidates_of(identification((ABORTUS, 0.0012, 950)), REFERENCES) == ('Brucella abortus',)
    assert closest_of(found, REFERENCES, 'Brucella suis').code == 'suis4'
    assert closest_of(found, REFERENCES, 'Brucella ovis') is None


def coarse_tie(count=50, min_hashes=1):
    """
    An identification where so many references tie that the margin is measuring nothing.

    mash derives the distance from the shared-hash count, so every reference sharing the same few hashes
    gets the same distance. The closest row is then first by an accident of ordering.
    """
    hits = [(CANIS, 0.2, 4), (SUIS4, 0.2, 4), (ABORTUS, 0.2, 4)]
    hits += [('filler{}.fasta'.format(i), 0.2, 4) for i in range(count)]
    return identification(*hits, min_hashes=min_hashes)


def test_an_undiscriminating_distance_is_never_reported_as_agreeing_with_the_st():
    """
    The bug this guards: withholding the close matches made the call read "MASH+ST (agree)".

    That is the verdict of a clean 950/1000 match. Here the distance separated nothing, so the agreement
    it would be claiming was never expressed.
    """
    found = coarse_tie()
    assert found.identified and found.undiscriminating == 52
    assert found.close_matches == []        # withheld, so the note can count instead of listing
    # A tie is ordered by name, so B-REF-BS4-40.fasta leads it and ST 16, B. suis, agrees with that
    call = call_species(found, REFERENCES, mlst=MlstResult(st='16'), table=TABLE)
    assert (call.species, call.source, call.verdict) == ('Brucella suis', FROM_BOTH, ARBITRATED)
    assert 'cannot separate 52 of the references' in call.notes[0]
    assert 'the ST agrees with the closest of them' in call.notes[0]


def test_an_undiscriminating_distance_gives_the_species_to_a_disagreeing_st():
    found = coarse_tie()
    call = call_species(found, REFERENCES, mlst=MlstResult(st='20'), table=TABLE)
    assert (call.species, call.source, call.verdict) == ('Brucella canis', FROM_ST, ARBITRATED)
    assert "the species is the ST's" in call.notes[0]


def test_an_undiscriminating_distance_with_no_st_is_unresolved():
    found = coarse_tie()
    call = call_species(found, REFERENCES, mlst=None, table=TABLE)
    assert (call.species, call.source, call.verdict) == ('Brucella suis', FROM_MASH, UNRESOLVED)
    assert 'the ST cannot choose because the MLST was not run' in call.notes[0]


def test_the_st_chooses_between_the_species_mash_cannot_separate():
    """The canis/suis case: MASH lands on canis, and ST 16 says B. suis."""
    found = identification((CANIS, 0.000239, 990), (SUIS4, 0.000420, 980))
    call = call_species(found, REFERENCES, mlst=MlstResult(st='16'), table=TABLE)
    assert (call.species, call.source, call.verdict) == ('Brucella suis', FROM_ST, ARBITRATED)
    assert call.candidates == ('Brucella canis', 'Brucella suis')
    # The vSNP code follows the species the ST chose, not the genome MASH happened to land on
    assert call.code == 'suis4'
    assert 'cannot separate Brucella canis and Brucella suis' in call.notes[0]
    assert 'ST 16 is Brucella suis in 15 of 15 isolate(s)' in call.notes[0]


def test_the_st_confirms_the_species_mash_could_not_separate():
    found = identification((CANIS, 0.000239, 990), (SUIS4, 0.000420, 980))
    call = call_species(found, REFERENCES, mlst=MlstResult(st='20'), table=TABLE)
    assert (call.species, call.source, call.verdict) == ('Brucella canis', FROM_BOTH, ARBITRATED)
    assert call.code == 'canis'
    assert 'agrees with the closest reference' in call.notes[0]


def test_an_unambiguous_species_the_st_agrees_with_is_reported_without_a_note():
    call = call_species(identification((ABORTUS, 0.0012, 950)), REFERENCES, mlst=MlstResult(st='1'), table=TABLE)
    assert (call.species, call.source, call.verdict) == ('Brucella abortus', FROM_BOTH, AGREE)
    assert (call.code, call.st_species, call.notes) == ('ab1', 'Brucella abortus', ())


def test_an_st_that_contradicts_an_unambiguous_species_is_reported_but_not_obeyed():
    """ST 20 is B. canis, and B. canis is nowhere near the closest reference: the MASH call stands."""
    call = call_species(identification((ABORTUS, 0.0012, 950)), REFERENCES, mlst=MlstResult(st='20'), table=TABLE)
    assert (call.species, call.source, call.verdict) == ('Brucella abortus', FROM_MASH, CONFLICT)
    assert call.code == 'ab1'
    assert 'the ST disagrees' in call.notes[0]
    assert 'check the sample' in call.notes[0]


def test_two_species_within_reach_and_an_st_that_cannot_choose():
    found = identification((CANIS, 0.000239, 990), (SUIS4, 0.000420, 980))
    call = call_species(found, REFERENCES, mlst=MlstResult(st=UNKNOWN), table=TABLE)
    assert (call.species, call.source, call.verdict) == ('Brucella canis', FROM_MASH, UNRESOLVED)
    assert 'the ST cannot choose because the sample has no ST' in call.notes[0]


def test_without_the_table_the_species_is_the_mash_one_and_nothing_is_said():
    found = identification((CANIS, 0.000239, 990), (SUIS4, 0.000420, 980))
    call = call_species(found, REFERENCES, mlst=MlstResult(st='16'), table=None)
    assert (call.species, call.source, call.verdict) == ('Brucella canis', FROM_MASH, UNRESOLVED)
    assert 'the ST to species table is not in use' in call.notes[0]


def test_a_sample_that_is_not_brucella_gets_no_species_from_its_st():
    found = identification((CANIS, 0.03, 100))
    call = call_species(found, REFERENCES, mlst=MlstResult(st='20'), table=TABLE)
    assert (call.species, call.source, call.verdict) == (UNKNOWN, '', '')
    assert 'not identified as Brucella' in call.notes[0]
    assert call_species(found, REFERENCES, mlst=MlstResult(st=UNKNOWN), table=TABLE).notes == ()


def test_a_reference_the_table_does_not_name_leaves_the_species_to_the_st():
    found = identification(('something_else.fasta', 0.0005, 900))
    call = call_species(found, REFERENCES, mlst=MlstResult(st='20'), table=TABLE)
    assert (call.species, call.source, call.verdict) == ('Brucella canis', FROM_ST, ARBITRATED)
    assert call.code == UNKNOWN  # No genome of that species in the sketch
    assert 'has no species in the reference table' in call.notes[0]
    # With nothing from either side there is no species to report
    assert call_species(found, REFERENCES, mlst=MlstResult(st='17'), table=TABLE).species == UNKNOWN


def test_without_mash_there_is_no_call():
    call = call_species(None, REFERENCES, mlst=MlstResult(st='20'), table=TABLE)
    assert (call.species, call.source, call.st) == (UNKNOWN, '', '20')


def test_the_thresholds_can_be_moved():
    found = identification((CANIS, 0.000239, 990), (SUIS4, 0.000420, 980))
    # ST 17 is B. suis in 18 of 20 isolates: too impure by default, usable at 0.9
    assert call_species(found, REFERENCES, mlst=MlstResult(st='17'), table=TABLE).verdict == UNRESOLVED
    call = call_species(found, REFERENCES, mlst=MlstResult(st='17'), table=TABLE, min_purity=0.9)
    assert (call.species, call.verdict) == ('Brucella suis', ARBITRATED)
    # ST 104 rests on one isolate
    assert call_species(found, REFERENCES, mlst=MlstResult(st='104'), table=TABLE,
                        min_isolates=1).species == 'Brucella suis'
