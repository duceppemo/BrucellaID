#!/usr/bin/env bash
# A run you can look at, on data this script makes itself: reads covering the nine MLST alleles of one
# sequence type, which BrucellaID has to type as that sequence type.
#
#   bash examples/run_example.sh          # a few seconds, nothing is downloaded
#   ST=8 bash examples/run_example.sh
#
# The sample holds no genomic reads, only the alleles, so the species comes back as ND: there is nothing for
# MASH to recognize. For the species identification, see the tutorial in the wiki, which runs on a real
# sample, or validation/run_validation.sh, which simulates reads from whole genomes.
set -euo pipefail

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
data="$here/data"
results="$here/results"
st=${ST:-1}

command -v brucellaid >/dev/null || { echo "brucellaid is not installed: pip install brucellaid" >&2; exit 1; }

mkdir -p "$data"
python3 - "$st" "$data/ST${st}_R1.fastq.gz" <<'PY'
"""Reads tiled over the alleles of one sequence type, taken from the data the package ships."""
import gzip
import sys

from brucellaid.mlst import alleles_fasta, load_profiles

st, output = sys.argv[1], sys.argv[2]
loci, profiles = load_profiles()
wanted = next((alleles for alleles, number in profiles.items() if number == st), None)
if wanted is None:
    raise SystemExit('ST {} is not in the profile table'.format(st))

sequences, name = {}, None
with open(alleles_fasta()) as f:
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            name = line[1:].split()[0]
            sequences[name] = []
        elif line:
            sequences[name].append(line)
sequences = {name: ''.join(parts) for name, parts in sequences.items()}

length, step, reads = 150, 20, []
for locus, allele in zip(loci, wanted, strict=True):
    sequence = sequences['{}_{}'.format(locus, allele)]
    starts = list(range(0, len(sequence) - length + 1, step))
    starts += [] if (len(sequence) - length) in starts else [len(sequence) - length]
    reads += [('{}_{}_{}'.format(locus, allele, start), sequence[start:start + length]) for start in starts]
with gzip.open(output, 'wt') as f:
    for read_name, sequence in reads:
        f.write('@{}\n{}\n+\n{}\n'.format(read_name, sequence, 'I' * len(sequence)))
print('{}: {} reads covering the {} alleles of ST {}'.format(output, len(reads), len(loci), st))
PY

rm -rf "$results"
brucellaid -r1 "$data/ST${st}_R1.fastq.gz" -o "$results" -s "ST${st}_example"

echo
echo "== $results/brucellaid_mlst.tsv"
grep -v '^#' "$results/brucellaid_mlst.tsv" | column -t -s $'\t'
echo
echo "The allele profile above is the one of ST ${st}, and the SequenceType column of"
echo "$results/brucellaid.tsv says so."
