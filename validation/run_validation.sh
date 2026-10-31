#!/usr/bin/env bash
# Download the genomes of genomes.tsv, simulate reads from each one, run BrucellaID and check what it called.
#
#   bash validation/run_validation.sh            # about 15 minutes and 1 GB in validation/data
#   COVERAGE=10 bash validation/run_validation.sh
#
# Nothing is written outside validation/data and validation/results, both of which are git-ignored.
set -euo pipefail

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
data="$here/data"
results="$here/results"
coverage=${COVERAGE:-20}
threads=${THREADS:-$( (nproc 2>/dev/null) || echo 4)}
efetch="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

command -v brucellaid >/dev/null || { echo "brucellaid is not installed: pip install -e ." >&2; exit 1; }
command -v mash >/dev/null || { echo "mash is not installed" >&2; exit 1; }
command -v seal.sh >/dev/null || { echo "seal.sh (BBTools) is not installed" >&2; exit 1; }

mkdir -p "$data/reads"
while IFS=$'\t' read -r sample accessions reference species biovar note; do
    case "$sample" in '#'*|'sample'|'') continue;; esac
    genome="$data/$sample.fasta"
    if [ ! -s "$genome" ]; then
        echo "Downloading $sample ($accessions)"
        # NCBI allows three requests per second without an API key
        sleep 1
        curl -fsS "$efetch?db=nuccore&id=${accessions//;/,}&rettype=fasta&retmode=text" -o "$genome.part"
        mv "$genome.part" "$genome"
    fi
    if [ ! -s "$data/reads/${sample}_R1.fastq.gz" ]; then
        python3 "$here/simulate_reads.py" "$genome" "$data/reads/$sample" --coverage "$coverage"
    fi
done < "$here/genomes.tsv"

rm -rf "$results"
brucellaid -i "$data/reads" -o "$results" -t "$threads" -p 2 > /dev/null
python3 "$here/check_results.py" "$results" "$data"
