#!/usr/bin/env python3
"""
Compare what BrucellaID called from the simulated reads with what the genomes they came from really are.

    python3 validation/check_results.py validation/results validation/data

The species comes from validation/genomes.tsv. The sequence type is not listed anywhere: it is read off each
genome by looking for an exact match of every MLST allele of the package, which is an answer that owes
nothing to the read based caller being checked.
"""
import argparse
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from brucellaid.mlst import alleles_fasta, load_profiles  # noqa: E402

COMPLEMENT = str.maketrans('ACGT', 'TGCA')
UNCHECKED = '-'


def read_table(path):
    """[{column: value}] of a tab separated file, ignoring comment lines."""
    with open(path) as f:
        rows = [line.rstrip('\n').split('\t') for line in f if line.strip() and not line.startswith('#')]
    return [dict(zip(rows[0], row, strict=True)) for row in rows[1:]]


def read_fasta(path):
    """{name: sequence} of a fasta file."""
    records, name = {}, None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                name = line[1:].split()[0]
                records[name] = []
            elif line and name:
                records[name].append(line.upper())
    return {name: ''.join(parts) for name, parts in records.items()}


def profile_of_genome(genome, alleles, loci):
    """
    The MLST profile of an assembly, by exact match: for each locus, the allele that occurs in the genome.

    :return: ({locus: allele or None}, {locus: note}) with a note for the loci that have no exact match, or
        more than one
    """
    sequence = ''.join(read_fasta(genome).values())
    both = sequence + '\n' + sequence.translate(COMPLEMENT)[::-1]
    profile, notes = {}, {}
    for locus in loci:
        found = [name.rsplit('_', 1)[1] for name, allele in alleles.items()
                 if name.rsplit('_', 1)[0] == locus and allele in both]
        profile[locus] = found[0] if len(found) == 1 else None
        if len(found) != 1:
            notes[locus] = 'no exact match' if not found else 'matches {}'.format(', '.join(found))
    return profile, notes


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument('results', help='Folder holding brucellaid.tsv')
    parser.add_argument('data', help='Folder holding the downloaded genomes')
    parser.add_argument('--genomes', default=str(pathlib.Path(__file__).resolve().parent / 'genomes.tsv'))
    args = parser.parse_args()

    loci, profiles = load_profiles()
    alleles = read_fasta(alleles_fasta())
    expected = {row['sample']: row for row in read_table(args.genomes)}
    called = {row['Sample']: row for row in read_table(pathlib.Path(args.results) / 'brucellaid.tsv')}

    header = ['sample', 'species', 'biovar of the genome', 'reference', 'hashes', 'ST called',
              'ST in the genome', 'species from', 'result']
    rows, failures = [], []
    for sample, row in expected.items():
        if sample not in called:
            failures.append('{}: no result'.format(sample))
            continue
        result = called[sample]
        problems = []
        # The biovar of genomes.tsv is checked through the reference: BrucellaID reports no biovar of
        # its own, and the closest reference being right is what the biovar column used to stand for.
        for column, key in (('Species', 'species'), ('ClosestReference', 'reference')):
            want = row[key]
            if want == UNCHECKED:
                continue
            if result[column] != want:
                problems.append('{} is "{}", expected "{}"'.format(column, result[column], want))
        if row['species'] == UNCHECKED and result['Species'] != 'ND':
            problems.append('identified as "{}", and it is not Brucella'.format(result['Species']))
        # These genomes are what they say they are, so the two kinds of evidence must not contradict each
        # other: a conflict here would mean the ST to species table and the reference table disagree.
        if 'conflict' in result['SpeciesSource']:
            problems.append('the ST contradicts the identification: {}'.format(result['Notes']))

        genome = pathlib.Path(args.data) / '{}.fasta'.format(sample)
        in_genome, notes = profile_of_genome(genome, alleles, loci)
        if all(in_genome[locus] for locus in loci):
            key = tuple(in_genome[locus] for locus in loci)
            st = profiles.get(key, 'new')
        elif len(notes) == len(loci):
            st = 'ND (no locus)'
        else:
            st = 'ND ({})'.format(', '.join(sorted(notes)))
        if row['species'] != UNCHECKED and not st.startswith('ND') and result['SequenceType'] != st:
            problems.append('ST is "{}", and the genome carries ST {}'.format(result['SequenceType'], st))

        rows.append([sample, result['Species'] or 'ND', row['biovar'], result['ClosestReference'],
                     result['SharedHashes'], result['SequenceType'], st, result['SpeciesSource'],
                     'ok' if not problems else 'FAILED'])
        failures += ['{}: {}'.format(sample, problem) for problem in problems]

    widths = [max(len(str(row[i])) for row in [header] + rows) for i in range(len(header))]
    for row in [header, ['-' * w for w in widths]] + rows:
        print('  '.join(str(cell).ljust(width) for cell, width in zip(row, widths, strict=True)).rstrip())
    print()
    if failures:
        print('{} problem(s):'.format(len(failures)))
        for failure in failures:
            print('  {}'.format(failure))
        return 1
    print('{} sample(s), all as expected.'.format(len(rows)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
