#!/usr/bin/env python3
"""
Simulate paired-end reads from a fasta file, so that the validation knows what the answer must be.

    python3 validation/simulate_reads.py genome.fasta out_prefix --coverage 20

Fragments are drawn at random positions, on either strand, and bases are substituted at a fixed rate. One
seed per sample, so a validation run gives the same reads every time.
"""
import argparse
import gzip
import random
import sys

COMPLEMENT = str.maketrans('ACGTN', 'TGCAN')


def read_fasta(path):
    """[sequence] of a fasta file, in file order."""
    sequences, parts = [], None
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if parts:
                    sequences.append(''.join(parts))
                parts = []
            elif parts is not None and line.strip():
                parts.append(line.strip().upper())
    if parts:
        sequences.append(''.join(parts))
    return sequences


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument('fasta')
    parser.add_argument('prefix', help='Reads are written to <prefix>_R1.fastq.gz and <prefix>_R2.fastq.gz')
    parser.add_argument('--coverage', type=float, default=20)
    parser.add_argument('--length', type=int, default=150)
    parser.add_argument('--fragment', type=int, default=350)
    parser.add_argument('--error', type=float, default=0.002, help='Substitutions per base')
    parser.add_argument('--seed', type=int, default=1)
    args = parser.parse_args()

    sequences = [s for s in read_fasta(args.fasta) if len(s) > args.fragment]
    if not sequences:
        raise SystemExit('No sequence longer than {} bp in {}'.format(args.fragment, args.fasta))
    total = sum(len(s) for s in sequences)
    pairs = int(total * args.coverage / (2 * args.length))
    rng = random.Random(args.seed)
    weights = [len(s) for s in sequences]
    with gzip.open(args.prefix + '_R1.fastq.gz', 'wt') as r1, gzip.open(args.prefix + '_R2.fastq.gz', 'wt') as r2:
        for i in range(pairs):
            sequence = rng.choices(sequences, weights=weights)[0]
            start = rng.randrange(0, len(sequence) - args.fragment)
            fragment = sequence[start:start + args.fragment]
            if rng.random() < 0.5:
                fragment = fragment.translate(COMPLEMENT)[::-1]
            forward, reverse = fragment[:args.length], fragment[-args.length:].translate(COMPLEMENT)[::-1]
            for handle, read, mate in ((r1, forward, 1), (r2, reverse, 2)):
                read = ''.join(rng.choice('ACGT') if rng.random() < args.error else base for base in read)
                handle.write('@{}_{}/{}\n{}\n+\n{}\n'.format(args.prefix.rsplit('/', 1)[-1], i, mate, read,
                                                             'I' * len(read)))
    print('{}: {:,} pairs of {} bp, {:.0f}x of {:,} bp'.format(args.prefix, pairs, args.length, args.coverage,
                                                               total))
    return 0


if __name__ == '__main__':
    sys.exit(main())
