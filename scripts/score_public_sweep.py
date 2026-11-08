#!/usr/bin/env python3
"""
Derive the head-to-head figures of the at-scale validation from scripts/public_sweep_disagreements.tsv.

    python3 -I scripts/score_public_sweep.py

Every number this prints is quoted on the wiki's "Validation at scale" page, under "Every usable public
read set in ENA". It is shipped because three of those figures had drifted from the data -- they were
computed before the sweep's retry pass filled in its last runs, and nothing recomputed them afterwards.

The agreement rule is the sweep's own: a call matches when its binomial equals the submitter's, and a
submitter label of "Brucella sp." matches any Brucella, because such a submission names no species to be
right or wrong about. Those runs are therefore excluded from "which set is right".
"""
import collections
import csv
import pathlib
import sys

DATA = pathlib.Path(__file__).resolve().parent / 'public_sweep_disagreements.tsv'
NO_SPECIES = ('Brucella sp.', 'Brucella')


def rows():
    with open(DATA) as f:
        return list(csv.DictReader((line for line in f if not line.startswith('#')), delimiter='\t'))


def matches(call, submitted):
    """Whether a call agrees with the submitter label, by the sweep's rule."""
    if call == 'ND':
        return False
    return call == submitted or (submitted in NO_SPECIES and call.startswith('Brucella'))


def main():
    disagreements = rows()
    print('{} runs the two sets called differently'.format(len(disagreements)))
    print()

    vague = [r for r in disagreements if r['submitted'] in NO_SPECIES]
    judgeable = [r for r in disagreements if r['submitted'] not in NO_SPECIES]
    extended = [r for r in judgeable
                if matches(r['extended_call'], r['submitted']) and not matches(r['vsnp_call'], r['submitted'])]
    vsnp = [r for r in judgeable
            if matches(r['vsnp_call'], r['submitted']) and not matches(r['extended_call'], r['submitted'])]
    raw = [r for r in disagreements
           if matches(r['extended_call'], r['submitted']) and not matches(r['vsnp_call'], r['submitted'])]
    print('which set is right, over the {} with a species-level answer:'.format(len(judgeable)))
    print('  extended set    {}'.format(len(extended)))
    print('  vSNP set        {}'.format(len(vsnp)))
    print('  ({} before setting aside the {} of those wins whose submitter named no species;'
          ' {} of the {} runs here have such a label)'.format(
              len(raw), len(raw) - len(extended), len(vague), len(disagreements)))
    print()

    nd = [r for r in disagreements if r['vsnp_call'] == 'ND']
    named = [r for r in disagreements if r['vsnp_call'] != 'ND']
    corrected = [r for r in named if r['submitted'] not in NO_SPECIES
                 and matches(r['extended_call'], r['submitted'])
                 and not matches(r['vsnp_call'], r['submitted'])]
    print('how the {} split:'.format(len(disagreements)))
    print('  {} the vSNP set could not name at all, and the extended set names'.format(len(nd)))
    print('  {} both named, differently -- the extended set corrects {}'.format(len(named), len(corrected)))
    print()

    distances = sorted(float(r['vsnp_distance']) for r in nd if r['vsnp_distance'])
    if distances:
        print('  the ND runs sit at {:.2f}-{:.2f} from the nearest of vSNP\'s 16'.format(
            distances[0], distances[-1]))
    print('  by the species the extended set gives them:')
    for species, n in collections.Counter(r['extended_call'] for r in nd).most_common():
        print('    {:<34} {}'.format(species, n))
    print()

    # The marine case the page singles out, because it is the one where the thin set does not merely fail
    # to name a species but names a different one.
    for species in ('Brucella pinnipedialis',):
        group = [r for r in disagreements if r['submitted'] == species]
        if not group:
            continue
        print('runs submitted as {} that the two sets disagree on: {}'.format(species, len(group)))
        for label, key in (('vSNP', 'vsnp_call'), ('extended', 'extended_call')):
            counts = collections.Counter(r[key] for r in group)
            right = sum(n for call, n in counts.items() if matches(call, species))
            print('  {:<9} {} right; {}'.format(
                label, right, ', '.join('{} x{}'.format(c, n) for c, n in counts.most_common())))
    return 0


if __name__ == '__main__':
    sys.exit(main())
