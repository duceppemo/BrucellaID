# Validation

`run_validation.sh` downloads the genomes listed in `genomes.tsv` from NCBI, simulates paired-end reads
from each one, runs BrucellaID over the folder and checks what it called:

```bash
bash validation/run_validation.sh            # about 25 minutes, 1 GB in validation/data
COVERAGE=10 THREADS=8 bash validation/run_validation.sh
```

Nothing is written outside `validation/data` and `validation/results`, which are git-ignored.

| File | What it is |
| --- | --- |
| `genomes.tsv` | The genomes, their accessions, and the species the identification must return |
| `simulate_reads.py` | Paired-end reads from a fasta file: random fragments, 0.2% substitutions, one seed per sample |
| `check_results.py` | Compares the reports with `genomes.tsv`, and with the MLST profile read off each genome |
| `run_validation.sh` | The three steps above, in order |

The species comes from `genomes.tsv`. The **sequence type is not listed anywhere**: `check_results.py`
reads it off each genome by looking for an exact match of every MLST allele of the package, which is an
answer that owes nothing to the read-based caller it is checking. A genome that is not in the MASH sketch
is only expected to give the right species, and a genome that is not *Brucella* is expected not to be
identified at all.

The results are on the [Validation](https://github.com/duceppemo/BrucellaID/wiki/Validation) page of the
wiki.
