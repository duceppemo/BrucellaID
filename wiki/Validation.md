# Validation

BrucellaID is validated three ways, each answering a question the others cannot, and each on its own page.

| | What it measures | Headline |
| --- | --- | --- |
| **[Validation suite](Validation-suite)** | 22 samples whose answer is known, reproducible from the repository in about 40 minutes | 22 of 22 as expected |
| **[Validation at scale](Validation-at-scale)** | simulated reads of all 2,022 assemblies the reference set was curated from, then every usable public *Brucella* read set in ENA | 1,538 of 1,543 held-out samples against the labels as written, all 1,543 against what the sequences demonstrably are; 98.7% of 4,042 real runs |
| **[The k-mer copy threshold](K-mer-threshold)** | how deep a read set has to be before a k-mer is believed, measured rather than assumed | 81 sweeps over 63 real libraries, 5x to 1903x |

The first is the one to run yourself: `bash validation/run_validation.sh` downloads its own data, simulates
its own reads and checks every answer against a truth the tool had no part in producing. The second is
measured at a scale the repository's script is not meant for, and is reported rather than reproduced. The
third is a derivation rather than a test: it is why `--min-kmer-copies` defaults to what it does.

**What none of them establishes** is on each page, under its own heading. In short: no metagenomic or mixed
sample is tested anywhere, no biovar is tested because BrucellaID reports none
([How it works](How-it-works#no-biovar)), and the public sweep scores agreement with submitters rather than
truth — the 20 mislabelled read sets it found are a reminder that the denominator has errors in it too.
Every mislabelled record either validation turned up is listed with its evidence in
[`curation/mislabelled_records.md`](https://github.com/duceppemo/BrucellaID/blob/main/curation/mislabelled_records.md).
