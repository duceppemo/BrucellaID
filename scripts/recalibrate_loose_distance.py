#!/usr/bin/env python3
"""
Derive LOOSE_MATCH_DISTANCE from scripts/loose_distance_calibration.tsv.

    python3 -I scripts/recalibrate_loose_distance.py

The rule, unchanged since 2.0: the threshold is the 95th percentile of the MASH distance of calls that
agree with the ENA submitter, taken to a round number, subject to the note staying rare. Past that
distance the closest reference is further away than a correctly identified library normally sits, which is
worth a note even though the species is usually still right.

The data file holds every identified run of the public sweep, so these are population figures. Until the
sweep was re-run under the current k-mer copy threshold it held a weighted 207-run sample instead; that
sample put the 95th percentile at 0.00066 against the 0.00068 measured here, which is the one piece of
evidence that the sampling had been sound.
"""
import csv
import pathlib
import sys

DATA = pathlib.Path(__file__).resolve().parent / 'loose_distance_calibration.tsv'
CANDIDATES = (0.0005, 0.0007, 0.001, 0.002, 0.005)
BANDS = ((0, 0.0005), (0.0005, 0.001), (0.001, 0.005), (0.005, 0.02), (0.02, 1.01))


def rows():
    with open(DATA) as f:
        for row in csv.DictReader((line for line in f if not line.startswith('#')), delimiter='\t'):
            yield float(row['distance']), int(row['shared']), row['agrees'] == 'yes', row['run']


def quantile(values, q):
    values = sorted(values)
    return values[min(len(values) - 1, int(q * len(values)))]


def main():
    points = list(rows())
    if not points:
        print('no data in {}'.format(DATA), file=sys.stderr)
        return 1
    agreeing = [p for p in points if p[2]]
    print('{} identified run(s); {} agree with the submitter'.format(len(points), len(agreeing)))
    print()

    print('distance of the calls that agree with the submitter:')
    for q in (0.5, 0.9, 0.95, 0.99):
        print('  p{:<3.0f} {:.5f}'.format(q * 100, quantile([p[0] for p in agreeing], q)))
    print()

    print('{:<11} {:>10} {:>11} {:>10}'.format('threshold', 'flagged', 'of ident.', 'agreeing'))
    for threshold in CANDIDATES:
        past = [p for p in points if p[0] > threshold]
        agree = sum(1 for p in past if p[2])
        print('{:<11} {:>10} {:>10.1f}% {:>9.0f}%'.format(
            threshold, len(past), 100 * len(past) / len(points),
            100 * agree / len(past) if past else 0))
    print()

    # What the note rests on, and the one thing here that is NOT a restatement of the threshold. mash
    # derives the distance from the shared-hash count, so "their shared hashes were low" says nothing the
    # distance has not said already; how often a call at a given distance is still right does.
    print('agreement with the submitter, by distance:')
    print('  {:<18} {:>8} {:>9}'.format('distance', 'runs', 'agree'))
    for low, high in BANDS:
        band = [p for p in points if low <= p[0] < high]
        if not band:
            continue
        print('  {:<18} {:>8} {:>8.1f}%'.format(
            '{:g}-{:g}'.format(low, high), len(band),
            100 * sum(1 for p in band if p[2]) / len(band)))
    print()
    print('Agreement falls with distance rather than staying flat, so a loose match does carry some')
    print('signal about the species being wrong. It is weak -- the worst band is still mostly right --')
    print('which is why the note is reported and not acted on.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
