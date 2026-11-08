#!/usr/bin/env python3
"""
Derive LOOSE_MATCH_DISTANCE from scripts/loose_distance_calibration.tsv.

    python3 -I scripts/recalibrate_loose_distance.py

The rule, unchanged since 2.0: the threshold is the 95th percentile of the MASH distance of calls that
agree with the ENA submitter, rounded up to a round number, subject to the note staying rare. Past that
distance the closest reference is further away than a correctly identified library normally sits, which is
worth a note even though the species is usually still right.

The percentiles are weighted: the 4,020 identified runs of the public sweep were binned by the distance
they had under the previous k-mer copy threshold, a sample was drawn from each bin, and each row carries
the weight pop/n of its bin. See the header of the data file.

Printing the table this produces is how the figures quoted on the wiki's "Validation at scale" page, and
in the comment on LOOSE_MATCH_DISTANCE, are checked.
"""
import csv
import pathlib
import sys

DATA = pathlib.Path(__file__).resolve().parent / 'loose_distance_calibration.tsv'
MIN_HASHES = 150
CANDIDATES = (0.0005, 0.0007, 0.001, 0.002, 0.005)


def rows():
    with open(DATA) as f:
        for row in csv.DictReader((line for line in f if not line.startswith('#')), delimiter='\t'):
            if int(row['shared']) < MIN_HASHES:
                continue  # not identified: no species was reported, so there is no call to calibrate on
            yield (float(row['distance']), int(row['shared']), float(row['weight']),
                   row['agrees'] == 'yes', row['run'], float(row['distance_old']))


def weighted_quantile(pairs, q):
    """The q-quantile of [(value, weight)]."""
    pairs = sorted(pairs)
    target, cumulative = q * sum(w for _, w in pairs), 0.0
    for value, weight in pairs:
        cumulative += weight
        if cumulative >= target:
            return value
    return pairs[-1][0]


def main():
    points = list(rows())
    if not points:
        print('no data in {}'.format(DATA), file=sys.stderr)
        return 1
    population = sum(p[2] for p in points)
    agreeing = [p for p in points if p[3]]
    print('{} runs measured, standing for {:.0f} identified runs; {:.0f} of those agree with the '
          'submitter'.format(len(points), population, sum(p[2] for p in agreeing)))
    print()

    print('distance of the calls that agree with the submitter:')
    for q in (0.5, 0.9, 0.95, 0.99):
        print('  p{:<3.0f} {:.5f}'.format(q * 100, weighted_quantile([(p[0], p[2]) for p in agreeing], q)))
    print()

    print('{:<11} {:>11} {:>10} {:>11}'.format('threshold', 'flagged', 'of those,', 'measured'))
    print('{:<11} {:>11} {:>10} {:>11}'.format('', 'share', 'agreeing', 'runs past'))
    for threshold in CANDIDATES:
        past = [p for p in points if p[0] > threshold]
        weight = sum(p[2] for p in past)
        agree = sum(p[2] for p in past if p[3])
        print('{:<11} {:>10.1f}% {:>9.0f}% {:>11}'.format(
            threshold, 100 * weight / population, 100 * agree / weight if weight else 0, len(past)))
    print()

    # The point the note rests on, and the one that is NOT a restatement of the threshold: a loose match
    # does not predict a wrong species. mash derives the distance from the shared-hash count, so "their
    # shared hashes were low" says nothing the distance has not already said -- the agreement rate does.
    print('agreement with the submitter, by distance (this is what justifies the note not changing '
          'the call):')
    print('  {:<18} {:>9} {:>9}'.format('distance', 'runs', 'agree'))
    bands = [(0, 0.0005), (0.0005, 0.001), (0.001, 0.005), (0.005, 0.02), (0.02, 1.01)]
    for low, high in bands:
        band = [p for p in points if low <= p[0] < high]
        if not band:
            continue
        weight = sum(p[2] for p in band)
        print('  {:<18} {:>9.0f} {:>8.1f}%'.format(
            '{:g}-{:g}'.format(low, high), weight, 100 * sum(p[2] for p in band if p[3]) / weight))
    print()

    clean = [p for p in points if p[5] < 0.001]
    print('the clean bulk: the {:.0f} runs that were already within 0.001 under the old threshold are '
          'now all within {:.5f}'.format(sum(p[2] for p in clean), max(p[0] for p in clean)))
    print()

    # The same arithmetic on the old distances, as a check that the sampling is not skewing the answer.
    # It will not land exactly on the published 2.0 figures, which were measured over all 4,020 runs
    # rather than over a 207-run sample; how close it comes is the sampling error of that sample.
    old_p95 = weighted_quantile([(p[5], p[2]) for p in agreeing], 0.95)
    old_past = sum(p[2] for p in points if p[5] > 0.005)
    print('the same weighting applied to the pre-2.1 distances in this file gives p95 {:.5f} and {:.1f}% '
          'past 0.005,'.format(old_p95, 100 * old_past / population))
    print('against the 0.00464 and 4.9% measured over all 4,020 runs: that gap is this sample\'s '
          'sampling error.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
