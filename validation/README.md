# Validation

`run_validation.sh` downloads the genomes listed in `genomes.tsv` from NCBI, simulates paired-end reads
from each one, downloads the real read sets listed in `reads.tsv` from the ENA, runs BrucellaID over the
folder and checks what it called:

```bash
bash validation/run_validation.sh            # about 40 minutes, 1.5 GB in validation/data
COVERAGE=10 THREADS=8 bash validation/run_validation.sh
REAL=0 bash validation/run_validation.sh     # simulated reads only, no ENA download
```

Nothing is written outside `validation/data` and `validation/results`, which are git-ignored.

| File | What it is |
| --- | --- |
| `genomes.tsv` | The genomes, their accessions, the species the identification must return, and the simulated depth |
| `reads.tsv` | Real read sets from the ENA: the species the submitter gave each library, and the answer it must keep giving |
| `simulate_reads.py` | Paired-end reads from a fasta file: random fragments, 0.2% substitutions, one seed per sample |
| `check_results.py` | Compares the reports with the two tables, and with the MLST profile read off each genome |
| `run_validation.sh` | The steps above, in order |

The species comes from `genomes.tsv`. The **sequence type is not listed anywhere**: `check_results.py`
reads it off each genome by looking for an exact match of every MLST allele of the package, which is an
answer that owes nothing to the read-based caller it is checking. A genome that is not in the MASH sketch
is only expected to give the right species, and a genome that is not *Brucella* is expected not to be
identified at all.

Two things exercise the depth-adaptive k-mer copy threshold, which a run at the default 20x never would:

* Two genomes are simulated deep no matter what `COVERAGE` asks for: one at **150x** and one, of another
  species, at **250x**. At those depths the threshold must rise (about one more required copy per 20x) or
  the error k-mers crowd the bottom-1000 sketch, and `check_results.py` checks all of it: the estimated
  coverage is within 20% of the truth, the threshold fits the depth, and the sample keeps at least 950 of
  its 1000 sketch hashes. The floor comes from a sweep of four genomes at 150x and two at 250x: the
  derived threshold keeps 1000/1000 on all six, the old fixed threshold of 2 keeps 860-902 at 150x and
  728-748 at 250x, with the distance inflated to 0.0025-0.0082.
* The **real read sets** run next to the simulated ones. Each is checked against the answer it is known to
  give, a library estimated at 60x or more must have raised the threshold above the floor of 2, and each
  one is compared with the simulated sample of the same reference genome: real reads and simulated reads
  of one genome must agree on the species and on the reference.

The results are on the [Validation](https://github.com/duceppemo/BrucellaID/wiki/Validation) page of the
wiki.
