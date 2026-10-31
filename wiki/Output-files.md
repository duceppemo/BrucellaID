# Output files

A run writes four files to the output folder, and a fifth with `--xlsx`. Every table is tab separated, with
the command that produced it on a first `#` line.

## `brucellaid.tsv` — one row per sample
Also printed to standard output.

| Column | What it holds |
| --- | --- |
| `Sample` | The sample name |
| `Species` | The species of the closest reference genome, or `ND` |
| `Biovar` | Its biovar, when the reference has one |
| `ClosestReference` | The reference genome the reads are closest to |
| `MashDistance` | The MASH distance to it: roughly the fraction of positions that differ |
| `SharedHashes` | How many of the reference's 1000 min-hashes the reads share |
| `SequenceType` | The ST of the profile table, `new` for a profile it does not hold, or `ND` |
| `AlleleProfile` | The nine allele numbers in the order of the loci, `-` for a locus not called |
| `Reads` | Reads in the input |
| `Notes` | Close matches, loci without a call, and why a sample was not identified |
| `Error` | Why the sample failed, when it did. The other samples are still reported |

```
Sample  Species           Biovar  ClosestReference           MashDistance  SharedHashes  SequenceType  AlleleProfile      Reads   Notes  Error
ab1     Brucella abortus  1       NC_006932-NC_006933.fasta  0.000024      999/1000      1             2-1-1-2-1-3-1-1-1  438192
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
The command, the start and end times, the versions of BrucellaID, Python, MASH, BBTools and Java, the
parameters, the reference table, and per sample: the species, the closest reference with its accessions, the
close matches, every allele call and the error, if any.

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
      "closest_reference": {"reference": "NC_006932-NC_006933.fasta", "mash_distance": 2.4e-05, "shared_hashes": 999},
      "mlst": {"sequence_type": "1", "profile": "2-1-1-2-1-3-1-1-1", "complete": true}
    }
  ]
}
```

## `stat_alignment_summary_<date>.xlsx` — the legacy report (`--xlsx`)
The report the first version of BrucellaID wrote, for whatever still reads it: the four columns
`time_stamp`, `sample_name`, `species` and `mlst_type`, in Courier New 10. The `species` column holds the
vSNP species code (`ab1`, `mel3`, `suis2`, ...), not the species name, exactly as before. It needs the
`xlsxwriter` package.
