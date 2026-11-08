# Output files

A run writes six files to the output folder, and one more with `--pdf` or `--xlsx`. Every table is tab
separated, with the command that produced it on a first `#` line.

## `brucellaid.tsv` — one row per sample
Also printed to standard output.

| Column | What it holds |
| --- | --- |
| `Sample` | The sample name |
| `Species` | The species, or `ND`. Usually the closest reference genome's; see `SpeciesSource` |
| `SpeciesSource` | What the species rests on: `MASH`, `MASH+ST (agree)`, `MASH+ST (arbitrated)`, `ST (arbitrated)`, `MASH (unresolved)` or `MASH (conflict)`. See [How it works](How-it-works#the-species-call), and [Troubleshooting](Troubleshooting) for the last two |
| `ClosestReference` | The reference genome the reads are closest to, whichever species is reported |
| `MashDistance` | The MASH distance to it: roughly the fraction of positions that differ |
| `SharedHashes` | How many of the reference's 1000 min-hashes the reads share |
| `SequenceType` | The ST of the profile table, `new` for a profile it does not hold, or `ND` |
| `AlleleProfile` | The nine allele numbers in the order of the loci, `-` for a locus not called |
| `STSpecies` | The species the isolates of that ST carry, empty when the ST could not be used (the note says why) |
| `Reads` | Reads in the input |
| `Notes` | Close matches, a closest reference that is not close, a distance too coarse to separate the references at all, loci without a call, and why a sample was not identified |
| `Error` | Why the sample failed, when it did. The other samples are still reported |

```
Sample      Species           SpeciesSource    ClosestReference                           MashDistance  SharedHashes  SequenceType  AlleleProfile      STSpecies         Reads    Notes  Error
ERR3454415  Brucella abortus  MASH+ST (agree)  GCA_000742195.1_ASM74219v1_genomic.fna.gz  0.000144      994/1000      1             2-1-1-2-1-3-1-1-1  Brucella abortus  1175082
            nearly as close: 9 other Brucella abortus genome(s)
```

A *B. canis* sample, where the distance alone cannot choose and the ST does, with the notes it carries:

```
ERR2136546  Brucella canis  MASH+ST (arbitrated)  GCA_000530495.1_Oliveri_genomic.fna.gz  0.000119  995/1000  20  1-6-4-1-5-3-5-4-6  Brucella canis  625576
            nearly as close: GCA_000741995.1_ASM74199v1_genomic.fna.gz (Brucella suis, 986/1000 hashes),
            3 other Brucella canis genome(s); MASH cannot separate Brucella canis and Brucella suis:
            ST 20 is Brucella canis in 70 of 70 isolate(s), which agrees with the closest reference
```

Both rows are from the [Tutorial](Tutorial) run, against the 432-genome default set, which is why the
references are GCA accessions. With `--reference-set vsnp` the same samples name vSNP's own sketch IDs
instead (`NC_006932-NC_006933.fasta`, `NC_010103-NC_010104.fasta`).

And a sample whose closest reference is **not close**, which is the other note the column carries. This is
`SRR24420453` from the ENA sweep, a library that holds a minority *Brucella* population among a great deal
of something else:

```
SRR24420453  Brucella melitensis  MASH+ST (agree)  GCA_002191375.1_ASM219137v1_genomic.fna.gz  0.061326  160/1000  11  3-2-3-2-1-5-3-2-10  Brucella melitensis  53748612
             loose match: the closest reference is 0.06133 away and shares 160/1000 hashes, where a read
             set of a species the sketch holds is expected to match within 0.005 (--loose-distance). The
             reads may hold more than Brucella, this lineage may have no near relative in the sketch, or
             a minority Brucella population may be being filtered out with the sequencing errors (try
             --min-kmer-copies either side of the derived value and see whether the distance moves); the
             species named is still the closest one; nearly as close: 14 other Brucella melitensis
             genome(s)
```

The species is still reported, and the nine loci still typed it as ST 11, which is *B. melitensis* and
agrees — but 160 of 1000 hashes at 0.061 is nothing like the 989-to-996 a clean library of a species in the
sketch returns, and the note is what says so. The three causes it lists are the three there are; see
[The k-mer copy threshold](K-mer-threshold) for the third.

**`Status`, not `Notes`, is the field to alert on.** Against the 432-genome default set the `Notes` column
is populated on nearly every sample, because a sample's nearest neighbours are mostly other genomes of its
own species and "nearly as close" is therefore the normal case. In 1.x the column was usually empty, so a
rule of "a non-empty `Notes` means look at this one" fired rarely; it now fires on almost everything and
needs replacing with a check on `Status`.

**There is no biovar column.** BrucellaID reports the species and the sequence type, and does not claim a
biovar: see [How it works](How-it-works#no-biovar). `brucellaid_mash.tsv` gives the biovar of each
*reference genome*, which is a fact about the reference and not about the sample.

## `brucellaid_mlst.tsv` — one row per sample and locus
What each allele call rests on.

| Column | What it holds |
| --- | --- |
| `Sample` | The sample name |
| `Locus` | One of the nine loci: `gap`, `aroA`, `glk`, `dnaK`, `gyrB`, `trpE`, `cobQ`, `int_hyp`, `omp25` |
| `Allele` | The allele called, or `ND` |
| `Reads` | Reads whose best match among the alleles of the locus is that allele |
| `RunnerUp`, `RunnerUpReads` | The next allele, and its reads: the margin the call has |
| `Note` | Why no allele was called: too few reads, or a tie |

A call whose runner-up is close behind is still a call, but it is worth a look: see
[How it works](How-it-works#the-mlst-alleles), which gives the margins the validated calls had.

## `brucellaid_mash.tsv` — one row per sample and reference genome
Every distance, closest first. This is the file to read when a species call surprises you, and the one the
quality assurance report points at when it lists only the closest few: the *Brucella* genomes are close to
one another, and the second line is often not far behind the first.

| Column | What it holds |
| --- | --- |
| `Sample` | The sample name |
| `Reference` | The reference genome, named as the reference table names it |
| `Species`, `Biovar` | What that **reference genome** is. The biovar here is the reference's own, not a call about the sample |
| `MashDistance` | The distance from the reads to that reference: roughly the fraction of positions that differ |
| `PValue` | MASH's p-value for the shared hashes, which is 0 for anything in the genus |
| `SharedHashes`, `TotalHashes` | Hashes shared, and the sketch size they are out of (1000) |

## `brucellaid.json` — everything, for a pipeline
The command, the start and end times, who ran it and on what, the versions of BrucellaID, Python, MASH,
BBTools and Java, the parameters, the reference table, the MD5 of every reference data file, and per
sample: the species and what decided it, the closest reference with its accessions, the k-mer copy
threshold the sample was sketched with and the coverage it was derived from (`min_kmer_copies`,
`estimated_coverage`), the close matches, every allele call, the input files with their checksums, the
status and the error, if any.

```json
{
  "tool": "brucellaid",
  "samples": [
    {
      "sample": "ERR3454415",
      "species": "Brucella abortus",
      "species_code": "abortus",
      "identified": true,
      "species_call": {"species": "Brucella abortus", "source": "MASH+ST", "verdict": "agree",
                       "mash_candidates": ["Brucella abortus"], "sequence_type": "1",
                       "st_species": "Brucella abortus", "st_purity": 0.994, "st_isolates": 466, "notes": []},
      "closest_reference": {"reference": "GCA_000742195.1_ASM74219v1_genomic.fna.gz",
                            "description": "Brucella abortus str. 2308",
                            "mash_distance": 0.000143503, "shared_hashes": 994, "total_hashes": 1000,
                            "min_kmer_copies": 14, "estimated_coverage": 69.84226787878788},
      "mlst": {"sequence_type": "1", "profile": "2-1-1-2-1-3-1-1-1", "complete": true}
    }
  ]
}
```

## `brucellaid_report.html` — the quality assurance report
**[View an example report](https://duceppemo.github.io/BrucellaID/example_report.html)** — the four samples of the [Tutorial](Tutorial), with the machine it ran
on removed. It is served from [GitHub Pages](https://duceppemo.github.io/BrucellaID/) because GitHub shows an HTML file in a repository as
source rather than rendering it; [the file itself](https://github.com/duceppemo/BrucellaID/blob/main/docs/example_report.html) is in `docs/`.

Each sample's identification facts include `K-mer copies required`: the threshold MASH sketched with, and
the estimated coverage it was derived from.

One self-contained HTML file: no stylesheet, no script, no font and no image is fetched from anywhere, so it
renders the same on a machine with no network, in an email attachment or out of an archive in ten years. It
prints to clean pages as it stands, and `--no-html` skips it.

It is written for the question "can I show how this result was reached?", so it holds the evidence rather
than the conclusion:

1. **Summary**: the run, the operator, and the samples grouped by what they were called — same species and
   same ST together, the largest group first, the failed samples last — with a **Reviewed by** box for a
   name, a date and a signature.
2. **Notes, warnings and errors**: one row per note, errors first. A note on a sample whose status is `ok`
   is a remark, not a problem.
3. **Samples**, one block each, **one to a printed page**: the identification with the five closest
   reference genomes and any other near enough that the distance cannot separate them (the closest row
   highlighted, the near ones shaded) — `brucellaid_mash.tsv` has the remaining distances; how the
   species was decided, including what the ST contributed and, when it was left out, why; the nine allele
   calls with the allele each one had to beat and the margin it had; and the input files with their size,
   date and MD5.
4. **Definitions and methods**: what a MASH distance is, why there is no biovar, how Seal
   counts, what `new` and `ND` mean, how the species is decided when two are equally close, and what the
   report does *not* establish.
5. **Run information**: operator, user account, computer, operating system, working directory, start and
   end time, the exact command, every parameter, the versions of BrucellaID, Python, MASH, BBTools and
   Java, and **the MD5 of every reference data file the run measured against** — which is what lets two
   runs a year apart be shown to have used the same data, or not.

Each sample has a `Status`: `ok`, `warning` (not identified, no ST, a `new` profile, two species the
evidence could not separate, or the two kinds of evidence disagreeing) or `failed`.

`--no-md5` skips the checksums of the input files, which cost a second or two per gigabyte. `--operator
NAME` puts a name other than the user account on the report.

## `brucellaid_report.pdf` — the same report as a PDF (`--pdf`)
The same evidence, laid out for paper, with "Page 3 of 8" in the footer of every page so a loose sheet can
be placed. It needs reportlab, which is not installed by default:

```bash
pip install 'brucellaid[pdf]'     # or: conda install -c conda-forge reportlab
```
The HTML report holds the same information and needs nothing installed, so the PDF is for the places that
want a file to sign and file rather than a page to read. **[Example PDF report](https://duceppemo.github.io/BrucellaID/example_report.pdf)** — the same four
samples, 8 pages.

## `brucellaid_mqc.json` — for MultiQC
A [MultiQC custom content](https://docs.seqera.io/multiqc/custom_content) file. Run `multiqc` on the output
folder and the run appears as a **BrucellaID** section: one row per sample with the species,
`SpeciesSource`, `STSpecies`, sequence type, shared hashes, closest reference and status.

```bash
brucellaid -i reads/ -o results/
multiqc results/
```
JSON rather than TSV on purpose: a sequence type is a label, and MultiQC would right-align a column of them
and offer to average it.

## `stat_alignment_summary_<date>.xlsx` — the legacy report (`--xlsx`)
The report the first version of BrucellaID wrote, for whatever still reads it: the four columns
`time_stamp`, `sample_name`, `species` and `mlst_type`, in Courier New 10. The `species` column holds the
reference table's **code** rather than the species name, exactly as before. It needs the `xlsxwriter`
package.

**Which code depends on the reference set**, and that changed in 2.0. The default 432-genome set codes the
species epithet and the biovar where NCBI records one (`abortus`, `abortus5`, `anthropi`); vSNP's own codes
(`ab1`, `mel3`, `suis2`, ...) exist only in the 16-genome table. A pipeline that matches on `ab1` needs
`--reference-set vsnp`, and `--xlsx` logs a warning when it is writing epithets rather than those codes.
The same applies to `species_code` in `brucellaid.json`. Note that the two kinds of code identify different
things: a vSNP code names one reference genome and is unique across the 16, while an epithet names a species,
so the 432 genomes of the default set carry 37 distinct codes between them.
