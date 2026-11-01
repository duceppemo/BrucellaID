# Troubleshooting

## `mash: error while loading shared libraries: libgsl.so.25`
The bioconda builds of mash 2.3 are linked against `libgsl.so.25`, but their recipe asks for
`gsl >=2.7,<2.8.0a0`, which gsl 2.7.1 satisfies while providing `libgsl.so.27`. Pin gsl in the environment:

```bash
conda install -n brucellaid -c conda-forge "gsl=2.7.0"
```
`environment.yml` already pins it, and installing BrucellaID from bioconda pulls a working combination.

## `"mash" was not found in your PATH` / `"seal.sh" was not found in your PATH`
Neither tool is on PyPI, so `pip install brucellaid` does not bring them:

```bash
conda install -c conda-forge -c bioconda mash "gsl=2.7.0" bbmap
```
BrucellaID also looks next to the running Python, so `/path/to/env/bin/brucellaid` works without activating
the environment.

## Every sample comes back as `ND`
The reads share fewer than 150 of the reference's 1000 hashes. In order of likelihood:

1. **It is not *Brucella*.** `brucellaid_mash.tsv` shows how far the closest reference is. A distance above
   about 0.05 is another genus; *Brucella anthropi* (formerly *Ochrobactrum anthropi*), the nearest
   neighbour, sits there.
2. **There is too little of it**: a mixed sample, a low yield extraction, or a very shallow run. Check the
   `Reads` column. Below roughly 5x, both the sketch and the allele calls get thin.
3. **The reads are of poor quality.** MASH drops the k-mers it sees once, so a read set full of errors
   loses most of its k-mers. Trim the reads and run it again.
4. **A sample is not what you think.** Check the `brucellaid_mash.tsv` rows: a sample that is *Brucella*
   but far from all 16 references is worth knowing about.

## The sequence type is `ND`
At least one locus has no call. The `Notes` column of `brucellaid.tsv` and the `Note` column of
`brucellaid_mlst.tsv` say which locus and why:

* **`no read matched this locus`** — the locus is not in the read set at all. That is a sample that is not
  *Brucella*, or reads that do not cover it.
* **`no allele reached 5 read(s)`** — too shallow. Lower `--min-reads` to see what the call would have been,
  and treat it as provisional.
* **`tied at N read(s) between ...`** — two alleles are indistinguishable in these reads, usually because
  nothing covers the position that separates them. More depth, or longer reads, settles it.

## The sequence type is `new`
The nine alleles were called, but the combination is not in the profile table. Either the sample really
carries a new profile (submit it to [PubMLST](https://pubmlst.org/organisms/brucella-spp)), or the profile
table is older than the sample. Rebuild it from PubMLST:

```bash
python3 scripts/make_mlst_data.py
```
Check the per-locus margins in `brucellaid_mlst.tsv` first: a call whose runner-up is one or two reads
behind is the kind of call a novel allele produces.

## `SpeciesSource` says `MASH (conflict)`
The sequence type belongs to a species that was not among the ones the MASH distance left within reach, and
the two statements cannot both be right. The species reported is the MASH one, because it rests on a
thousand min-hashes rather than nine loci, but the sample is worth a look:

* **A mixture.** Two species in one tube give one species' hashes and another's alleles. Check the shared
  hashes in `brucellaid.tsv` and the per-locus margins in `brucellaid_mlst.tsv`: a mixture usually shows up
  as both.
* **A mislabelled isolate in the ST table.** The species of each ST comes from the submitters' labels in the
  PubMLST isolate database. `brucellaid/data/mlst_st_species.tsv` gives the counts, so an ST resting on a
  handful of isolates can be read for what it is; `--min-st-isolates` raises the bar.
* **A genuinely new combination**, for an ST seen so far in one species only.

`--no-st-species` reports the MASH species without consulting the ST at all.

## The biovar is not the one the laboratory reported
Read `Biovar` as the biovar of the closest reference genome, not as a biovar call. The *Brucella* biovar
references differ by very little, and the note on the sample says when another was nearly as close. Use
classical biotyping, or [brucetyper](https://github.com/duceppemo/brucetyper), for a real biovar.

## `Seal failed: java.lang.OutOfMemoryError`
Give Seal more memory with `--memory 4g`, or run fewer samples at a time (`-p`). The default heap is
computed from the size of the allele file and shared between the jobs, within 85% of the machine.

## `The temporary folder ... contains a space, a comma, "xmx" or "xms"`
`seal.sh` runs its command line through `eval`, and BBTools reads any argument holding `xmx` as a memory
setting. BrucellaID works around unusual paths by linking the inputs under neutral names in a temporary
folder, but it cannot do that if the temporary folder itself is unusable. Set `TMPDIR` to a plain path.

## A folder gives fewer samples than expected
* fasta files are skipped: BrucellaID reads fastq (see the [FAQ](FAQ)).
* Two files that give the same sample name are an error, not a guess. Rename or move one.
* An `_R2` file without its `_R1` is used as single-end reads, with a warning.
* Hidden folders (`.git`, `.snakemake`) and hidden files are skipped, and so is the output folder when it is
  inside the input folder.

## Nothing above
Run it again with `--verbose`: the log then holds the MASH and Seal command lines and their output. Open an
issue with that log, the command and the version
([bug report](https://github.com/duceppemo/BrucellaID/issues/new/choose)).
