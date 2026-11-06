# Usage

```
brucellaid (-r1 FILE | -i FOLDER) -o FOLDER [options]
```

## Input
One sample at a time, or a folder:

```bash
brucellaid -r1 sample_R1.fastq.gz -r2 sample_R2.fastq.gz -o results/   # Paired-end
brucellaid -r1 sample.fastq.gz -o results/                             # Single-end
brucellaid -i reads/ -o results/                                       # A folder of samples
```

| Option | What it does |
| --- | --- |
| `-r1`, `--r1 FILE` | Reads: single-end, or R1 of paired-end. fastq, gzipped or not |
| `-r2`, `--r2 FILE` | R2 reads, for paired-end fastq files |
| `-i`, `--input FOLDER` | Folder of fastq files, searched recursively |
| `-s`, `--sample NAME` | With `-r1`: the sample name used in the reports |

With `-i`, each single-end file and each R1/R2 pair is one sample. The sample name is the file name without
the extension and without the read suffix (`_R1`, `_1`, `_S12_L001_R1_001`). Two files that would give the
same sample name are an error rather than a guess, and fasta files are skipped with a warning: BrucellaID
reads fastq only (see the [FAQ](FAQ)).

## Output
| Option | What it does |
| --- | --- |
| `-o`, `--output FOLDER` | Folder for the reports. Created if needed |
| `--no-html` | Do not write `brucellaid_report.html`, the quality assurance report |
| `--pdf` | Also write that report as a PDF. Needs reportlab: `pip install 'brucellaid[pdf]'` |
| `--no-md5` | Do not checksum the input files for the reports (a second or two per gigabyte) |
| `--operator NAME` | The name to put on the report. Default: the user account |
| `--xlsx` | Also write `stat_alignment_summary_<date>.xlsx`, the report of the first version |

The main table is printed to standard output as well, so `brucellaid ... > results.tsv` works: the messages
go to standard error. See [Output files](Output-files).

Six files are written by default: the three tables, the JSON report, a MultiQC file and
`brucellaid_report.html`, a self-contained quality assurance report holding the closest reference
genomes, every allele call with its margin, the input files with their checksums and the checksum of
every reference data file.
Run `multiqc results/` on the output folder to fold the run into a MultiQC report.

## Identification and typing
| Option | Default | What it does |
| --- | --- | --- |
| `--min-hashes N` | 150 | Of the reference's 1000 min-hashes, how many the reads must share for a species to be reported |
| `--min-kmer-copies N` | from the coverage | Copies of a k-mer the reads must hold for MASH to sketch it. The default derives it from the number of bases in the files: one more copy per 20x, from 2 to 20 |
| `--min-reads N` | 5 | Reads behind an allele call for it to be kept |
| `--no-mlst` | off | Identify the species only. The ST then takes no part in the species call either |
| `--no-st-species` | off | Do not let the ST take part in the species call: report the closest reference's species, whatever else was just as close |
| `--min-st-purity F` | 0.95 | Share of the isolates of an ST that must carry the same species for that ST to take part |
| `--min-st-isolates N` | 3 | Isolates an ST must rest on to take part |
| `-k`, `--kmer-size N` | 25 | K-mer size used to assign the reads to the alleles |
| `-m`, `--mismatch N` | 0 | Mismatches allowed between an allele k-mer and a read |

`--mismatch` is 0 on purpose: the alleles of a locus differ by a few bases, and a k-mer allowed one mismatch
would match the called allele and its neighbours alike. Raising it makes the allele calls less, not more,
reliable.

The ST is only consulted where the MASH distance leaves two species within reach, which in the sketch of
the package means *B. canis* and *B. suis* bv. 4 and nothing else. `--no-st-species` restores the
behaviour of version 1.0.0, where the closest reference decided alone. See
[How it works](How-it-works) for what the `SpeciesSource` column reports.

Lowering `--min-hashes` reports a species for a sample that holds little *Brucella* (a mixed culture, a low
yield extraction); the closest reference is then still the closest of the 16, which is not the same as the
right one. The `brucellaid_mash.tsv` report shows what the call rests on.

`--min-kmer-copies` is almost never worth setting by hand: the right value depends on the coverage, which
is why it is derived from the input rather than fixed. Set it only to reproduce an older run
(`--min-kmer-copies 2` is what versions up to 1.1.0 always used) or when the sample is not a single
*Brucella* genome, so that coverage is not what the number of bases implies. See
[How it works](How-it-works) for why the value has to move with the depth.

## Reference data
The package ships all five files, and each can be replaced:

| Option | What it takes |
| --- | --- |
| `--reference-set NAME` | Which packaged reference set to identify against: `default` (the 16 genomes of vSNP, 8 species) or `extended` (431 curated genomes, the whole *Brucella* clade — see [How it works](How-it-works#the-extended-reference-set)) |
| `--sketch FILE` | A MASH sketch of reference genomes. Keep MASH's defaults (k = 21, 1000 hashes) |
| `--references FILE` | The table saying what each genome of the sketch is |
| `--alleles FILE` | A fasta file of MLST alleles, named `<locus>_<number>` |
| `--profiles FILE` | A profile table: a header of `ST` and the loci, then one row per sequence type |
| `--st-species FILE` | A table of the species of each ST: `ST`, `Species`, `Purity`, `Isolates`, `Biovars`, `OtherSpecies` |

`--alleles` and `--profiles` go together: every allele a profile uses must be in the fasta file, and every
locus of the fasta file must be a column of the table. The genomes of `--sketch` and the rows of
`--references` should name the same genomes; BrucellaID warns when they do not, and a sample closest to a
genome that is not in the table is reported without a species.

`--reference-set extended` names packaged files, so it cannot be combined with `--sketch` or
`--references`: pass your own files, or pick a packaged set, not both.

`scripts/make_mlst_data.py` rebuilds the MLST files from PubMLST, `scripts/make_st_species.py` rebuilds the
ST to species table from PubMLST's isolate database, and `scripts/make_reference_data.py` and
`scripts/make_extended_references.py` rebuild the two sketches from the accessions of their reference
tables.

## Performance and other options
| Option | Default | What it does |
| --- | --- | --- |
| `-t`, `--threads N` | all available | Threads given to MASH and to Seal |
| `-p`, `-j`, `--parallel`, `--jobs N` | 1 | With `-i`: samples analysed at the same time, sharing the threads |
| `--memory SIZE` | from the allele file | Java heap for Seal, per job, e.g. `2g` |
| `--verbose` | off | Debug messages, including the MASH and Seal command lines |
| `-v`, `--version` | | Print the version |

A sample takes a few seconds at 20x. Most of it is reading the reads twice, once for MASH and once for Seal,
so `-p` helps more than `-t` on a folder of samples.

## Exit codes
| Code | Meaning |
| --- | --- |
| 0 | Every sample was analysed |
| 1 | The run could not start, or at least one sample failed (the others are still reported) |
| 2 | The arguments are wrong |
| 130 | Interrupted |
