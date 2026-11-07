#!/usr/bin/env bash
# Download the genomes of genomes.tsv, simulate reads from each one, download the real read sets of
# reads.tsv, run BrucellaID over all of it and check what it called. The run is pinned to the 16 vSNP
# genomes (--reference-set vsnp), which is what this suite is about; see the note on the run itself.
#
#   bash validation/run_validation.sh            # about 40 minutes and 1.5 GB in validation/data
#   COVERAGE=10 bash validation/run_validation.sh
#   REAL=0 bash validation/run_validation.sh     # leave the real read sets out: simulated reads only
#
# Nothing is written outside validation/data and validation/results, both of which are git-ignored.
set -euo pipefail

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
data="$here/data"
results="$here/results"
coverage=${COVERAGE:-20}
threads=${THREADS:-$( (nproc 2>/dev/null) || echo 4)}
real=${REAL:-1}
efetch="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
ena="https://www.ebi.ac.uk/ena/portal/api/filereport"

command -v brucellaid >/dev/null || { echo "brucellaid is not installed: pip install -e ." >&2; exit 1; }
command -v mash >/dev/null || { echo "mash is not installed" >&2; exit 1; }
command -v seal.sh >/dev/null || { echo "seal.sh (BBTools) is not installed" >&2; exit 1; }

mkdir -p "$data/reads"
declare -A fetched
while IFS=$'\t' read -r sample accessions reference species biovar depth note; do
    case "$sample" in '#'*|'sample'|'') continue;; esac
    genome="$data/$sample.fasta"
    if [ ! -s "$genome" ]; then
        if [ -n "${fetched[$accessions]:-}" ]; then
            cp "${fetched[$accessions]}" "$genome"   # the same genome, simulated at another depth
        else
            echo "Downloading $sample ($accessions)"
            # NCBI allows three requests per second without an API key
            sleep 1
            curl -fsS "$efetch?db=nuccore&id=${accessions//;/,}&rettype=fasta&retmode=text" -o "$genome.part"
            mv "$genome.part" "$genome"
        fi
    fi
    fetched[$accessions]="$genome"
    [ "$depth" = "-" ] && depth=$coverage
    if [ ! -s "$data/reads/${sample}_R1.fastq.gz" ]; then
        python3 "$here/simulate_reads.py" "$genome" "$data/reads/$sample" --coverage "$depth"
    fi
done < "$here/genomes.tsv"

if [ "$real" != 0 ]; then
    while IFS=$'\t' read -r run rest; do
        case "$run" in '#'*|'sample'|'') continue;; esac
        if [ ! -s "$data/reads/${run}_1.fastq.gz" ]; then
            echo "Downloading $run from the ENA"
            urls=$(curl -fsS "$ena?accession=$run&result=read_run&fields=fastq_ftp&format=tsv" \
                   | tail -1 | cut -f2 | tr ';' ' ')
            [ -n "$urls" ] || { echo "The ENA returned no fastq URL for $run" >&2; exit 1; }
            for url in $urls; do
                name=$(basename "$url")
                # The ENA throttles with transient 403s and cut transfers: retry, and resume the file
                curl -fsS --retry 5 --retry-all-errors --retry-delay 10 -C - \
                     "https://$url" -o "$data/reads/$name.part"
                mv "$data/reads/$name.part" "$data/reads/$name"
            done
        fi
    done < "$here/reads.tsv"
fi

rm -rf "$results"
# --reference-set vsnp, not the 2.0 default: this suite measures the 16 genomes read for read, and its
# expectations are their sketch IDs. The negative control depends on it too -- Brucella anthropi
# ATCC 49188 is here to come back ND as the nearest thing to the genus that is outside the vSNP set, and
# the extended set carries it as a reference and so names it correctly. What the default set does is
# measured at scale instead: wiki/Validation.md, "At scale" and "On real public read sets".
brucellaid --reference-set vsnp -i "$data/reads" -o "$results" -t "$threads" -p 2 > /dev/null
flags=(); [ "$real" = 0 ] && flags+=(--no-real)
python3 "$here/check_results.py" "$results" "$data" "${flags[@]}"
