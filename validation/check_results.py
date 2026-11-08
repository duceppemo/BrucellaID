#!/usr/bin/env python3
"""
Compare what BrucellaID called from the reads with what the samples really are.

    python3 validation/check_results.py validation/results validation/data

The simulated samples are checked against validation/genomes.tsv, and their sequence type against the
genome the reads came from: it is read off the assembly by looking for an exact match of every MLST allele
of the package, which is an answer that owes nothing to the read based caller being checked. A sample
simulated at a fixed depth is also checked for what that depth demands: the k-mer copy threshold must have
been derived from the coverage, and the sketch must have stayed clean.

The real read sets are checked against validation/reads.tsv, and each one against the simulated sample
whose expected reference is the same genome: real reads and simulated reads of one reference must agree on
the species and on the reference.
"""
import argparse
import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from brucellaid.mlst import alleles_fasta, load_profiles  # noqa: E402

COMPLEMENT = str.maketrans('ACGT', 'TGCA')
UNCHECKED = '-'
# A simulated sample must keep at least this many of the 1000 sketch hashes: its reads hold nothing but the
# genome plus 0.2% substitutions, so a short sketch means the copy threshold let the error k-mers in.
# Measured on four genomes at 150x and two at 250x: the derived threshold keeps 1000/1000 on all six, and a
# fixed threshold of 2 keeps 902 at best (860-902 at 150x, 728-748 at 250x), so 950 tells the two apart
# with about fifty hashes of room on either side.
MIN_CLEAN_HASHES = 950
# The thresholds are this project's own arithmetic, kept apart from brucellaid.mash on purpose: if the
# shipped formula drifts, the validation must notice, not drift along. So these three are restated here
# rather than imported, and they are changed by hand when the rule is re-derived -- which is what caught
# this file the last time, when it still expected one copy per 20x and the rule had become one per 5x.
# The derivation is on the wiki, under "The k-mer copy threshold".
DEPTH_PER_COPY, COPY_FLOOR, COPY_CEILING = 5, 2, 48
# A real library at or above this estimated coverage must have raised the threshold above the floor. One
# copy per 5x clears a floor of 2 from about 13x, so any real library in reads.tsv is well past it.
ADAPTIVE_FROM = 15


def expected_copies(coverage):
    """The k-mer copy threshold a read set of this coverage must get: one more copy per 5x, within [2, 48]."""
    return min(COPY_CEILING, max(COPY_FLOOR, round(coverage / DEPTH_PER_COPY)))


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


def genome_st(genome, alleles, loci, profiles):
    """The sequence type carried by an assembly, read off the genome itself."""
    in_genome, notes = profile_of_genome(genome, alleles, loci)
    if all(in_genome[locus] for locus in loci):
        return profiles.get(tuple(in_genome[locus] for locus in loci), 'new')
    if len(notes) == len(loci):
        return 'ND (no locus)'
    return 'ND ({})'.format(', '.join(sorted(notes)))


def identification_of(record):
    """(shared hashes, min k-mer copies, estimated coverage) of one sample of brucellaid.json."""
    best = (record or {}).get('closest_reference') or {}
    return best.get('shared_hashes'), best.get('min_kmer_copies'), best.get('estimated_coverage')


def copies_cell(record):
    """"7 (150x)" for the table: the k-mer copy threshold, and the coverage it came from."""
    _, copies, coverage = identification_of(record)
    if copies is None:
        return '-'
    return '{} ({:.0f}x)'.format(copies, coverage) if coverage else str(copies)


def check_expected(row, result, problems):
    """The checks genomes.tsv and reads.tsv share: species, closest reference, and no conflict."""
    for column, key in (('Species', 'species'), ('ClosestReference', 'reference')):
        want = row.get(key, UNCHECKED)
        if want != UNCHECKED and result[column] != want:
            problems.append('{} is "{}", expected "{}"'.format(column, result[column], want))
    # These samples are what they say they are, so the two kinds of evidence must not contradict each
    # other: a conflict here would mean the ST to species table and the reference table disagree.
    if 'conflict' in result['SpeciesSource']:
        problems.append('the ST contradicts the identification: {}'.format(result['Notes']))


def check_depth(row, record, problems):
    """A sample simulated at a fixed depth: the copy threshold must fit the depth, the sketch stay clean."""
    depth = float(row['coverage'])
    shared, copies, estimated = identification_of(record)
    if estimated is None or abs(estimated - depth) > 0.2 * depth:
        problems.append('estimated coverage is {}, the reads hold {:.0f}x'.format(
            '{:.0f}x'.format(estimated) if estimated else 'unknown', depth))
    # Against the ESTIMATE, not the true depth. The tool only ever sees the estimate, and the check above
    # already tolerates that being 20% out -- which at one copy per 5x is six copies at 150x, so comparing
    # with expected_copies(depth) would fail a sample whose coverage it had just accepted. What is being
    # checked here is that the shipped formula was applied, so the formula's own input is what to apply it
    # to. (The tolerance of 1 absorbs the rounding of a fractional estimate.)
    elif abs(copies - expected_copies(estimated)) > 1:
        problems.append('the k-mer copy threshold is {}, expected {} at the {:.0f}x estimated'.format(
            copies, expected_copies(estimated), estimated))
    if shared is not None and shared < MIN_CLEAN_HASHES:
        problems.append('only {}/1000 hashes shared with the reference: the sketch is polluted '
                        'at {:.0f}x'.format(shared, depth))


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument('results', help='Folder holding brucellaid.tsv and brucellaid.json')
    parser.add_argument('data', help='Folder holding the downloaded genomes')
    parser.add_argument('--genomes', default=str(pathlib.Path(__file__).resolve().parent / 'genomes.tsv'))
    parser.add_argument('--reads', default=str(pathlib.Path(__file__).resolve().parent / 'reads.tsv'))
    parser.add_argument('--no-real', action='store_true',
                        help='The run held no real read sets: check the simulated samples only')
    args = parser.parse_args()

    loci, profiles = load_profiles()
    alleles = read_fasta(alleles_fasta())
    genomes = read_table(args.genomes)
    real = [] if args.no_real else read_table(args.reads)
    called = {row['Sample']: row for row in read_table(pathlib.Path(args.results) / 'brucellaid.tsv')}
    with open(pathlib.Path(args.results) / 'brucellaid.json') as f:
        records = {record['sample']: record for record in json.load(f)['samples']}

    header = ['sample', 'kind', 'species', 'reference', 'hashes', 'k-mer copies', 'ST called',
              'ST expected', 'species from', 'result']
    rows, failures = [], []

    def show(sample, kind, result, st, problems):
        rows.append([sample, kind, result['Species'] or 'ND', result['ClosestReference'],
                     result['SharedHashes'], copies_cell(records.get(sample)), result['SequenceType'], st,
                     result['SpeciesSource'], 'ok' if not problems else 'FAILED'])
        failures.extend('{}: {}'.format(sample, problem) for problem in problems)

    for row in genomes:
        sample = row['sample']
        if sample not in called:
            failures.append('{}: no result'.format(sample))
            continue
        result, problems = called[sample], []
        check_expected(row, result, problems)
        if row['species'] == UNCHECKED and result['Species'] != 'ND':
            problems.append('identified as "{}", and it is not Brucella'.format(result['Species']))
        if row['coverage'] != UNCHECKED:
            check_depth(row, records.get(sample), problems)
        st = genome_st(pathlib.Path(args.data) / '{}.fasta'.format(sample), alleles, loci, profiles)
        if row['species'] != UNCHECKED and not st.startswith('ND') and result['SequenceType'] != st:
            problems.append('ST is "{}", and the genome carries ST {}'.format(result['SequenceType'], st))
        show(sample, 'simulated {:.0f}x'.format(float(row['coverage'])) if row['coverage'] != UNCHECKED
             else 'simulated', result, st, problems)

    twins = {row['reference']: row['sample'] for row in genomes if row['reference'] != UNCHECKED}
    for row in real:
        sample = row['sample']
        if sample not in called:
            failures.append('{}: no result'.format(sample))
            continue
        result, problems = called[sample], []
        check_expected(row, result, problems)
        if row['st'] != UNCHECKED and result['SequenceType'] != row['st']:
            problems.append('ST is "{}", this library carries ST {}'.format(result['SequenceType'], row['st']))
        shared, copies, estimated = identification_of(records.get(sample))
        if estimated and estimated >= ADAPTIVE_FROM and copies is not None and copies <= COPY_FLOOR:
            problems.append('the k-mer copy threshold stayed at {} on an estimated {:.0f}x'.format(
                copies, estimated))
        twin = called.get(twins.get(row['reference'], ''))
        if twin:
            for column in ('Species', 'ClosestReference'):
                if result[column] != twin[column]:
                    problems.append('real reads say {} "{}", the simulated {} says "{}"'.format(
                        column, result[column], twins[row['reference']], twin[column]))
        show(sample, 'real', result, row['st'], problems)

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
