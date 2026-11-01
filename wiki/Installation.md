# Installation

BrucellaID needs Python 3.10 or newer, [MASH](https://github.com/marbl/mash) and
[BBTools](https://sourceforge.net/projects/bbmap/) (which provides `seal.sh`, and Java with it). Neither MASH
nor BBTools is on PyPI, so they come from bioconda; BrucellaID itself needs nothing but the Python standard
library. Two output formats need a package each, and neither is required: `xlsxwriter` for the legacy xlsx
report (`--xlsx`) and `reportlab` for the PDF quality assurance report (`--pdf`). The HTML version of that
report holds the same information and needs nothing.

## pip, with the tools from conda (recommended)
```bash
conda create -n brucellaid -c conda-forge -c bioconda mash "gsl=2.7.0" bbmap python">=3.10"
conda activate brucellaid
pip install brucellaid
brucellaid --version
```
The `gsl=2.7.0` pin is not cosmetic: the bioconda builds of mash 2.3 are linked against `libgsl.so.25`, while
their recipe accepts gsl 2.7.1, which provides `libgsl.so.27`. Without the pin, mash may fail to start with
`error while loading shared libraries: libgsl.so.25`. See [Troubleshooting](Troubleshooting).

The optional output formats: `pip install 'brucellaid[xlsx]'`, `pip install 'brucellaid[pdf]'`, or
`pip install 'brucellaid[all]'` for both. With conda, `conda install -c conda-forge xlsxwriter reportlab`.

## conda, in one command — not yet
```bash
conda create -n brucellaid -c conda-forge -c bioconda brucellaid    # once the recipe is merged
```
This will bring MASH and BBTools with it, and is the shortest way in. The recipe is
[awaiting review on bioconda](https://github.com/bioconda/bioconda-recipes/pull/71040); until it is merged the command above finds nothing, so use pip as
above. The same goes for the ready-made container: biocontainers builds one from the bioconda package, so
there is none yet.

## From a clone
```bash
git clone https://github.com/duceppemo/BrucellaID
cd BrucellaID
conda env create -f environment.yml
conda activate brucellaid
pip install --no-deps -e .
pytest            # optional: 310 tests
```
`python -m brucellaid` works too, and takes the same arguments as the `brucellaid` command.

## Docker or Singularity
Once the bioconda recipe is merged, biocontainers builds an image from it automatically:
```bash
docker run --rm -v "$PWD":/data quay.io/biocontainers/brucellaid:<tag> brucellaid -i /data/reads -o /data/results
```
The tag will be listed on biocontainers.pro. Until then, build your own from the pip install above, or run
BrucellaID from a clone.

## What is installed
The package ships its reference data and the logo of the reports, about 350 kB in all:

| File | What it is |
| --- | --- |
| `references.msh` | The MASH sketch of the 16 *Brucella* reference genomes |
| `references.tsv` | What each genome is: species, biovar, strain, NCBI accessions |
| `mlst_alleles.fasta` | The 309 alleles of the nine MLST loci |
| `mlst_profiles.tsv` | The 122 sequence types |
| `mlst_st_species.tsv` | The species of each sequence type, from PubMLST's isolate database |
| `logo.svg`, `logo.png` | The logo the quality assurance reports carry |

`--sketch`, `--references`, `--alleles`, `--profiles` and `--st-species` take files of your own instead
(see [Usage](Usage)).
