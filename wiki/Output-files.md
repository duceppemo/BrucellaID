# Output files

A run writes six files to the output folder, and one more with `--pdf` or `--xlsx`. Every table is tab
separated, with the command that produced it on a first `#` line.

## `brucellaid.tsv` — one row per sample
Also printed to standard output.

| Column | What it holds |
| --- | --- |
| `Sample` | The sample name |
| `Species` | The species, or `ND`. Usually the closest reference genome's; see `SpeciesSource` |
| `Biovar` | The biovar of the closest reference **of that species**, when it has one |
| `SpeciesSource` | What the species rests on: `MASH`, `MASH+ST (agree)`, `MASH+ST (arbitrated)`, `ST (arbitrated)`, `MASH (unresolved)` or `MASH (conflict)`. See [How it works](How-it-works) |
| `ClosestReference` | The reference genome the reads are closest to, whichever species is reported |
| `MashDistance` | The MASH distance to it: roughly the fraction of positions that differ |
| `SharedHashes` | How many of the reference's 1000 min-hashes the reads share |
| `SequenceType` | The ST of the profile table, `new` for a profile it does not hold, or `ND` |
| `AlleleProfile` | The nine allele numbers in the order of the loci, `-` for a locus not called |
| `STSpecies` | The species the isolates of that ST carry, empty when the ST could not be used (the note says why) |
| `Reads` | Reads in the input |
| `Notes` | Close matches, loci without a call, and why a sample was not identified |
| `Error` | Why the sample failed, when it did. The other samples are still reported |

```
Sample  Species           Biovar  SpeciesSource    ClosestReference           MashDistance  SharedHashes  SequenceType  AlleleProfile      STSpecies         Reads   Notes  Error
ab1     Brucella abortus  1       MASH+ST (agree)  NC_006932-NC_006933.fasta  0.000024      999/1000      1             2-1-1-2-1-3-1-1-1  Brucella abortus  438192
```

A *B. canis* sample, where the distance alone cannot choose and the ST does, with the note it carries:

```
canis   Brucella canis            MASH+ST (arbitrated)  NC_010103-NC_010104.fasta  0.000000  1000/1000  20  1-6-4-1-5-3-5-4-6  Brucella canis  66251
        nearly as close: B-REF-BS4-40.fasta (Brucella suis bv. 4, 990/1000 hashes); MASH cannot separate
        Brucella canis and Brucella suis: ST 20 is Brucella canis in 70 of 70 isolate(s), which agrees with
        the closest reference
```

## `brucellaid_mlst.tsv` — one row per sample and locus
What each allele call rests on.

| Column | What it holds |
| --- | --- |
| `Allele` | The allele called, or `ND` |
| `Reads` | Reads whose best match among the alleles of the locus is that allele |
| `RunnerUp`, `RunnerUpReads` | The next allele, and its reads: the margin the call has |
| `Note` | Why no allele was called: too few reads, or a tie |

A call whose runner-up is close behind is still a call, but it is worth a look: see
[How it works](How-it-works).

## `brucellaid_mash.tsv` — one row per sample and reference genome
Every distance, closest first, with the species of each reference. This is the file to read when a species
call surprises you: the *Brucella* genomes are close to one another, and the second line is often not far
behind the first.

## `brucellaid.json` — everything, for a pipeline
The command, the start and end times, who ran it and on what, the versions of BrucellaID, Python, MASH,
BBTools and Java, the parameters, the reference table, the MD5 of every reference data file, and per
sample: the species and what decided it, the closest reference with its accessions, the close matches,
every allele call, the input files with their checksums, the status and the error, if any.

```json
{
  "tool": "brucellaid",
  "samples": [
    {
      "sample": "ab1",
      "species": "Brucella abortus",
      "biovar": "1",
      "species_code": "ab1",
      "identified": true,
      "species_call": {"species": "Brucella abortus", "source": "MASH+ST", "verdict": "agree",
                       "mash_candidates": ["Brucella abortus"], "sequence_type": "1",
                       "st_species": "Brucella abortus", "st_purity": 0.994, "st_isolates": 466, "notes": []},
      "closest_reference": {"reference": "NC_006932-NC_006933.fasta", "mash_distance": 2.4e-05, "shared_hashes": 999},
      "mlst": {"sequence_type": "1", "profile": "2-1-1-2-1-3-1-1-1", "complete": true}
    }
  ]
}
```

## `brucellaid_report.html` — the quality assurance report
**[View an example report](https://htmlpreview.github.io/?https://github.com/duceppemo/BrucellaID/blob/main/assets/example_report.html)** — the four samples of the [Tutorial](Tutorial), with the machine it ran
on removed ([the file itself](https://github.com/duceppemo/BrucellaID/blob/main/assets/example_report.html); GitHub shows HTML as source, hence the viewer).

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
3. **Samples**, one block each: the identification with **every** distance MASH measured and not only the
   winning one (the closest row highlighted, and the rows near enough that the distance does not separate
   them shaded); how the species was decided, including what the ST contributed and, when it was left out,
   why; the nine allele calls with the allele each one had to beat and the margin it had; and the input
   files with their size, date and MD5.
4. **Definitions and methods**: what a MASH distance is, why the biovar is not a biovar call, how Seal
   counts, what `new` and `ND` mean, how the species is decided when two are equally close, and what the
   report does *not* establish.
5. **Run information**: operator, user account, computer, operating system, working directory, start and
   end time, the exact command, every parameter, the versions of BrucellaID, Python, MASH, BBTools and
   Java, and **the MD5 of every reference data file the run measured against** — which is what lets two
   runs a year apart be shown to have used the same data, or not.

Each sample has a `Status`: `ok`, `warning` (not identified, no ST, a `new` profile, or the two kinds of
evidence disagreeing) or `failed`.

`--no-md5` skips the checksums of the input files, which cost a second or two per gigabyte. `--operator
NAME` puts a name other than the user account on the report.

## `brucellaid_report.pdf` — the same report as a PDF (`--pdf`)
The same evidence, laid out for paper, with "Page 3 of 9" in the footer of every page so a loose sheet can
be placed. It needs reportlab, which is not installed by default:

```bash
pip install 'brucellaid[pdf]'     # or: conda install -c conda-forge reportlab
```
The HTML report holds the same information and needs nothing installed, so the PDF is for the places that
want a file to sign and file rather than a page to read.

## `brucellaid_mqc.json` — for MultiQC
A [MultiQC custom content](https://docs.seqera.io/multiqc/custom_content) file. Run `multiqc` on the output
folder and the run appears as a **BrucellaID** section: one row per sample with the species, biovar,
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
vSNP species code (`ab1`, `mel3`, `suis2`, ...), not the species name, exactly as before. It needs the
`xlsxwriter` package.
