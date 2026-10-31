#!/usr/bin/env python3
"""
Rebuild the MASH sketch of the reference genomes from the accessions of the reference table.

    python3 scripts/make_reference_data.py --genomes /tmp/refs --out /tmp/references.msh
    python3 scripts/make_reference_data.py --genomes /tmp/refs \
        --fasta Bceti1Cudo.fasta=/data/Bceti1Cudo.fasta --fasta B-REF-BS4-40.fasta=/data/B-REF-BS4-40.fasta

Each genome is downloaded from NCBI into the folder given by --genomes, under exactly the name the
reference table gives it, and "mash sketch" is run from inside that folder: the ID MASH stores is the file
name it was given, and that ID is what the reference table keys on.

Two of the 16 genomes of the shipped sketch (Bceti1Cudo.fasta and B-REF-BS4-40.fasta) come from vSNP under
its own names rather than from an accession, so they cannot be downloaded. Pass them with --fasta, or the
script stops rather than quietly writing a sketch with 14 genomes in it.

The sketch keeps MASH's defaults, k = 21 and 1000 hashes: the threshold of 150 shared hashes is read
against that sketch size.
"""
import argparse
import pathlib
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.request

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from brucellaid.references import load_references  # noqa: E402

EFETCH = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
DATA = pathlib.Path(__file__).resolve().parent.parent / 'brucellaid' / 'data'
NUCCORE = ('NC_', 'NZ_', 'AM', 'CP', 'AE', 'KF')  # An assembly accession (GCF_...) is not in nuccore


def download(accessions, path):
    """The fasta of a list of NCBI nucleotide accessions, written to path."""
    url = '{}?db=nuccore&id={}&rettype=fasta&retmode=text'.format(EFETCH, ','.join(accessions))
    try:
        with urllib.request.urlopen(url, timeout=300) as response, open(str(path) + '.part', 'wb') as f:
            shutil.copyfileobj(response, f)
    except (urllib.error.URLError, TimeoutError) as e:
        raise SystemExit('Cannot download {}: {}'.format(', '.join(accessions), e)) from None
    pathlib.Path(str(path) + '.part').rename(path)
    with open(path) as f:
        headers = sum(1 for line in f if line.startswith('>'))
    if headers != len(accessions):
        raise SystemExit('{}: {} sequence(s) downloaded for {} accession(s)'.format(
            path.name, headers, len(accessions)))
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument('--genomes', required=True, help='Folder for the downloaded genomes (created if needed)')
    parser.add_argument('--out', default=str(DATA / 'references.msh'), help='Sketch to write')
    parser.add_argument('--references', default=None, help='A reference table of your own')
    parser.add_argument('--fasta', action='append', default=[], metavar='NAME=PATH',
                        help='A genome that has no accession, as it is named in the reference table')
    args = parser.parse_args()

    references = load_references(args.references)
    given = dict(entry.split('=', 1) for entry in args.fasta)
    unknown = sorted(set(given) - set(references))
    if unknown:
        raise SystemExit('--fasta names {} , which is not in the reference table'.format(', '.join(unknown)))

    folder = pathlib.Path(args.genomes)
    folder.mkdir(parents=True, exist_ok=True)
    missing, files = [], []
    for name, reference in references.items():
        path = folder / name
        accessions = [a for a in reference.accessions.split(';') if a]
        if name in given:
            shutil.copyfile(given[name], path)
            print('{}: copied from {}'.format(name, given[name]))
        elif not accessions or not all(a.startswith(NUCCORE) for a in accessions):
            missing.append(name)
            continue
        elif path.exists():
            print('{}: already downloaded'.format(name))
        else:
            time.sleep(1)  # NCBI allows three requests per second without an API key
            download(accessions, path)
            print('{}: downloaded {}'.format(name, ', '.join(accessions)))
        files.append(path.name)

    if missing:
        raise SystemExit('{} genome(s) of the reference table have no nucleotide accession and were not given '
                         'with --fasta: {}'.format(len(missing), ', '.join(missing)))

    mash = shutil.which('mash') or str(pathlib.Path(sys.executable).parent / 'mash')
    out = pathlib.Path(args.out).resolve()
    subprocess.run([mash, 'sketch', '-o', str(out.with_suffix('')), *sorted(files)],
                   cwd=folder, check=True, capture_output=True)
    info = subprocess.run([mash, 'info', '-t', str(out)], check=True, capture_output=True, text=True)
    ids = [line.split('\t')[2] for line in info.stdout.splitlines() if line.strip() and not line.startswith('#')]
    if sorted(ids) != sorted(references):
        raise SystemExit('The sketch holds {}, and the reference table {}'.format(sorted(ids), sorted(references)))
    print('{}: {} genome(s), {:,} bytes'.format(out, len(ids), out.stat().st_size))
    return 0


if __name__ == '__main__':
    sys.exit(main())
