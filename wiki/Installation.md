# Installation

BrucellaID needs Python 3.10 or newer, [MASH](https://github.com/marbl/mash) and
[BBTools](https://sourceforge.net/projects/bbmap/) (which provides `seal.sh`, and Java with it). Neither MASH
nor BBTools is on PyPI, so they come from bioconda; BrucellaID itself needs nothing but the Python standard
library, plus `xlsxwriter` if you want the legacy xlsx report.

## conda (recommended)
```bash
conda create -n brucellaid -c conda-forge -c bioconda brucellaid
conda activate brucellaid
brucellaid --version
```
This brings MASH and BBTools with it.

## pip, with the tools from conda
```bash
conda create -n brucellaid -c conda-forge -c bioconda mash "gsl=2.7.0" bbmap python">=3.10"
conda activate brucellaid
pip install brucellaid
```
The `gsl=2.7.0` pin is not cosmetic: the bioconda builds of mash 2.3 are linked against `libgsl.so.25`, while
their recipe accepts gsl 2.7.1, which provides `libgsl.so.27`. Without the pin, mash may fail to start with
`error while loading shared libraries: libgsl.so.25`. See [Troubleshooting](Troubleshooting).

Add the legacy xlsx report with `pip install brucellaid[xlsx]`.

## From a clone
```bash
git clone https://github.com/duceppemo/BrucellaID
cd BrucellaID
conda env create -f environment.yml
conda activate brucellaid
pip install --no-deps -e .
pytest            # optional: 194 tests
```
`python -m brucellaid` works too, and takes the same arguments as the `brucellaid` command.

## Docker or Singularity
The bioconda package gives a container for free:
```bash
docker run --rm -v "$PWD":/data quay.io/biocontainers/brucellaid:<tag> brucellaid -i /data/reads -o /data/results
```
The tag is listed on the [biocontainers page](https://biocontainers.pro/tools/brucellaid).

## What is installed
The package ships its reference data, about 290 kB in all:

| File | What it is |
| --- | --- |
| `references.msh` | The MASH sketch of the 16 *Brucella* reference genomes |
| `references.tsv` | What each genome is: species, biovar, strain, NCBI accessions |
| `mlst_alleles.fasta` | The 309 alleles of the nine MLST loci |
| `mlst_profiles.tsv` | The 122 sequence types |

`--sketch`, `--references`, `--alleles` and `--profiles` take files of your own instead
(see [Usage](Usage)).
