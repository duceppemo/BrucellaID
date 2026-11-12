# Running the at-scale sweep again

`run_validation.sh` beside this file is the 22-sample suite, and it is the one to run. The at-scale
measurement is different: every usable public *Brucella* read set in ENA, 1.6 TB of reads, both reference
sets, with MLST. It does not fit in a repository and there is no script here that reproduces it. These are
the notes anyone repeating it will want, since the data and the drivers live outside the repo.

What it establishes, and the figures themselves, are on the
[Validation at scale](https://github.com/duceppemo/BrucellaID/wiki/Validation-at-scale) page.

## What is shipped, and what is not

Each derived constant ships the measurements behind it and a script that re-derives the published number:

| in the repo | derives |
| --- | --- |
| `scripts/kmer_threshold_sweep.tsv`, `make_kmer_figure.py` | the k-mer copy threshold, and the figure |
| `scripts/kmer_threshold_simulated.tsv` | why simulated reads cannot set it |
| `scripts/loose_distance_calibration.tsv`, `recalibrate_loose_distance.py` | `--loose-distance` |
| `scripts/public_sweep_disagreements.tsv`, `score_public_sweep.py` | the head-to-head figures |

Run either script with `python3 -I scripts/<name>.py`. If a number on the wiki disagrees with what they
print, the wiki is wrong: three figures were found stale that way, having been computed before a sweep's
retry pass finished and never recomputed.

The reads, the sweep output and the assemblies are not shipped. They were kept at
`/media/8tb/analyses/` — `kmer_rederive/reads` (1.5 TB), `brucellaid_211_resweep` (the 2.1.1 sweep of both
sets), `kmer_rederive/assemblies` and `assemblies_unpaired`.

## Four things that cost a day each

**ENA publishes some runs as three fastq files.** `_1` holds the unpaired reads and `_2`/`_3` are the
mates, with identical read counts; 18 of 4,044 *Brucella* runs are like this. A sweep that globs `*_1` and
assumes `_2` is its mate therefore mis-pairs exactly the runs that carry unpaired reads, feeding the
orphan file in as R1. The counts then differ, Seal refuses the pair, and the sample gets no sequence type
with nothing to say why — which looks like an archive fault and is not. Take the pair from the manifest's
file count. BrucellaID warns when a `_3` sits beside a pair it has formed, which is the only point at
which the layout is visible.

**`pgrep -f` and `pkill -f` match their own invocation**, and anything else whose command line merely
mentions the pattern. A wait loop that polls `pgrep -f fetch_all.py` will match the shell that launched it
if that command also mentioned the string, and wait forever. Scan `/proc` for the real process, skip your
own PID, and do not put the needle in the command line as a literal.

**Killing a child of a wrapper loop achieves nothing.** A driver that relaunches a worker in rounds will
simply start the next round. Kill the wrapper first, then the worker, then confirm neither is left.

**Judge a long job by its output, not its process list.** Count the files it should have produced and look
at the newest one's timestamp. A job whose log is still empty has probably died at the first step — a
missing tool on `PATH` is the usual reason, since `mash`, `seal.sh` and `reformat.sh` are in the conda
environment and not on the default `PATH`.
