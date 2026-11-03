"""The reference table of the MASH sketch."""

import pytest

from brucellaid.errors import BrucellaIDError
from brucellaid.mash import sketch_references
from brucellaid.references import (
    COLUMNS,
    EXTENDED_REFERENCES_FILE,
    EXTENDED_SKETCH_FILE,
    Reference,
    compare_with_sketch,
    data_path,
    load_references,
    lookup,
    parse_references,
    reference_sketch,
)

HEADER = '\t'.join(COLUMNS)


def table(*rows):
    return ['# a comment\n', HEADER + '\n', '\n'] + [row + '\n' for row in rows]


def test_the_table_of_the_package():
    references = load_references()
    assert len(references) == 16
    assert load_references() is references  # Read once
    abortus = references['NC_006932-NC_006933.fasta']
    assert abortus.code == 'ab1'
    assert abortus.label == 'Brucella abortus bv. 1'
    assert abortus.description == 'Brucella abortus bv. 1 str. 9-941'
    assert abortus.accessions == 'NC_006932.1;NC_006933.1'
    assert {r.species for r in references.values()} == {
        'Brucella abortus', 'Brucella melitensis', 'Brucella suis', 'Brucella canis', 'Brucella ovis',
        'Brucella ceti', 'Brucella microti', 'Brucella neotomae'}
    # Every code is unique: the legacy xlsx report writes it as the species
    assert len({r.code for r in references.values()}) == len(references)


def test_a_reference_without_a_biovar_or_a_strain():
    reference = Reference('x.fasta', 'x', 'Brucella species')
    assert reference.label == 'Brucella species'
    assert reference.description == 'Brucella species'


def test_parse_references_fills_in_the_missing_trailing_columns():
    references = parse_references(table('a.fasta\ta\tBrucella a\t1\tstrain', 'b.fasta\tb\tBrucella b'))
    assert references['a.fasta'].accessions == ''
    assert references['b.fasta'] == Reference('b.fasta', 'b', 'Brucella b', '', '', '')


@pytest.mark.parametrize('lines, message', [
    (['reference\tcode\n', 'a\tb\n'], 'must start with the columns'),
    (table('a.fasta\ta\tBrucella a\t1\ts\tNC_1\textra'), 'more than the'),
    (table('a.fasta\t\tBrucella a'), 'no reference, code or species'),
    (table('a.fasta\ta\tBrucella a', 'a.fasta\tb\tBrucella b'), 'appears twice'),
    ([HEADER + '\n'], 'holds no reference'),
])
def test_parse_references_errors(lines, message):
    with pytest.raises(BrucellaIDError, match=message):
        parse_references(lines)


def test_load_references_of_an_unreadable_file(tmp_path):
    with pytest.raises(BrucellaIDError, match='cannot be read'):
        load_references(str(tmp_path / 'nowhere.tsv'))


def test_data_path_of_a_missing_file():
    assert data_path('references.tsv').endswith('data/references.tsv')
    with pytest.raises(BrucellaIDError, match='is missing from the installed'):
        data_path('not_shipped.tsv')


def test_reference_sketch():
    assert reference_sketch().endswith('data/references.msh')
    assert reference_sketch('/my/own.msh') == '/my/own.msh'


def test_lookup_by_name_and_by_path():
    references = parse_references(table('a.fasta\ta\tBrucella a'))
    assert lookup(references, 'a.fasta').code == 'a'
    assert lookup(references, '/data/refs/a.fasta').code == 'a'  # mash reports the path it sketched
    assert lookup(references, 'b.fasta') is None


def test_compare_with_sketch():
    references = parse_references(table('a.fasta\ta\tBrucella a', 'b.fasta\tb\tBrucella b'))
    assert compare_with_sketch(['/refs/a.fasta', 'b.fasta'], references) == []
    messages = compare_with_sketch(['a.fasta', 'c.fasta'], references)
    assert len(messages) == 2
    assert 'not in the reference table' in messages[0] and 'c.fasta' in messages[0]
    assert 'not in the MASH sketch' in messages[1] and 'b.fasta' in messages[1]


def test_the_extended_table_of_the_package():
    references = load_references(data_path(EXTENDED_REFERENCES_FILE))
    assert len(references) == 428
    assert all(r.code and r.species for r in references.values())
    species = {r.species for r in references.values()}
    # 29 of the 30 named species, and "Brucella sp." for the unnamed lineages. The thirtieth,
    # B. papionis, has no assembly in INSDC, so no sketch can hold it.
    assert len(species) == 30
    assert 'Brucella sp.' in species
    assert 'Brucella anthropi' in species  # The clade formerly known as Ochrobactrum
    abortus5 = references['GCA_000163115.1_ASM16311v1_genomic.fna.gz']
    assert abortus5.label == 'Brucella abortus bv. 5'
    assert abortus5.description == 'Brucella abortus bv. 5 str. B3196'
    assert abortus5.accessions == 'GCA_000163115.1'


def test_the_extended_sketch_matches_the_extended_table():
    ids = sketch_references(reference_sketch(data_path(EXTENDED_SKETCH_FILE)))
    assert len(ids) == 428
    references = load_references(data_path(EXTENDED_REFERENCES_FILE))
    assert compare_with_sketch(ids, references) == []
