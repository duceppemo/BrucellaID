#!/usr/bin/env python3
"""
Rebuild the extended MASH sketch from the accessions of the extended reference table.

    python3 scripts/make_extended_references.py --genomes /tmp/extended_refs

Each of the 428 genomes is downloaded from the NCBI assembly FTP into the folder given by --genomes,
under exactly the name the reference table gives it, and "mash sketch" is run from inside that folder:
the ID MASH stores is the file name it was given, and that ID is what the reference table keys on.
A genome already in the folder is not downloaded again, so an interrupted run resumes where it stopped.

How the 428 were chosen is scripts/extended_references_curation.tsv: the dereplicated representatives
of the 2,022 NCBI taxon-234 assemblies behind the mashID Brucella pack (figshare article 33970168,
version 4), minus the genomes whose NCBI label their own sequence contradicts, minus the unlabelled
twins of named species, with one cluster-consensus relabel.

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

from brucellaid.references import EXTENDED_REFERENCES_FILE, load_references  # noqa: E402

FTP = 'https://ftp.ncbi.nlm.nih.gov/genomes/all'
DATA = pathlib.Path(__file__).resolve().parent.parent / 'brucellaid' / 'data'


def url_of(reference):
    """
    The NCBI FTP URL of a reference, whose file name is <accession>_<assembly name>_genomic.fna.gz.

    GCA_000163115.1_ASM16311v1_genomic.fna.gz lives under .../GCA/000/163/115/GCA_000163115.1_ASM16311v1/.
    """
    stem = reference[:-len('_genomic.fna.gz')]
    accession = '_'.join(stem.split('_')[:2])
    digits = accession.split('_')[1].split('.')[0]
    return '/'.join([FTP, accession[:3], digits[0:3], digits[3:6], digits[6:9], stem, reference])


def download(reference, path, tries=3):
    """One genome, written to path, with a pause and another try when NCBI drops the connection."""
    for attempt in range(1, tries + 1):
        try:
            with urllib.request.urlopen(url_of(reference), timeout=300) as response, \
                    open(str(path) + '.part', 'wb') as f:
                shutil.copyfileobj(response, f)
            pathlib.Path(str(path) + '.part').rename(path)
            return
        except (urllib.error.URLError, TimeoutError) as e:
            if attempt == tries:
                raise SystemExit('Cannot download {}: {}'.format(reference, e)) from None
            time.sleep(10 * attempt)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument('--genomes', required=True, help='Folder for the downloaded genomes (created if needed)')
    parser.add_argument('--out', default=str(DATA / 'references_extended.msh'), help='Sketch to write')
    parser.add_argument('--references', default=str(DATA / EXTENDED_REFERENCES_FILE),
                        help='A reference table of your own')
    args = parser.parse_args()

    references = load_references(args.references)
    genomes = pathlib.Path(args.genomes)
    genomes.mkdir(parents=True, exist_ok=True)
    for number, name in enumerate(references, start=1):
        path = genomes / name
        if path.is_file():
            continue
        print('[{}/{}] {}'.format(number, len(references), name), flush=True)
        download(name, path)
        time.sleep(0.4)  # NCBI asks for no more than three requests a second

    names = sorted(references)
    list_file = genomes / 'genomes.txt'
    list_file.write_text('\n'.join(names) + '\n')
    out = pathlib.Path(args.out).resolve()
    subprocess.run(['mash', 'sketch', '-k', '21', '-s', '1000', '-o', str(out.with_suffix('')),
                    '-l', list_file.name], cwd=genomes, check=True)
    print('Wrote {}: {} genome(s)'.format(out, len(names)))


if __name__ == '__main__':
    main()
