#!/usr/bin/env python3
"""
Rebuild the MLST data of the package from PubMLST.

    python3 scripts/make_mlst_data.py              # writes brucellaid/data/mlst_alleles.fasta and mlst_profiles.tsv
    python3 scripts/make_mlst_data.py --check      # says what would change, and writes nothing

The source is the "MLST (9 locus)" scheme of the Brucella database at https://pubmlst.org, read through its
REST API. Without an API key the data is the one submitted up to a cut-off date, which the API states and this
script prints: that is enough for a release, but note the date in wiki/Changelog.md and in
brucellaid/data/README.md.

After running it, check the diff, run the tests, and look at a validation sample: a retired allele that is
still in a local copy of the database is exactly the kind of thing this script is here to remove (vSNP's copy
held a gyrB_23 that PubMLST has deleted, and it matched the B. abortus 9-941 genome better than the real one).
"""
import argparse
import json
import pathlib
import sys
import urllib.error
import urllib.request

API = 'https://rest.pubmlst.org/db/pubmlst_brucella_seqdef'
SCHEME = 1
DATA = pathlib.Path(__file__).resolve().parent.parent / 'brucellaid' / 'data'
WRAP = 70


def fetch(url):
    """The body of a PubMLST URL, as text."""
    request = urllib.request.Request(url, headers={'User-Agent': 'brucellaid make_mlst_data.py'})
    try:
        with urllib.request.urlopen(request, timeout=120) as response:
            return response.read().decode('utf-8')
    except (urllib.error.URLError, TimeoutError) as e:
        raise SystemExit('Cannot read {}: {}'.format(url, e)) from None


def read_fasta(text):
    """[(name, sequence)] of a fasta file, in file order."""
    records, name, parts = [], None, []
    for line in text.splitlines():
        line = line.strip()
        if line.startswith('>'):
            if name is not None:
                records.append((name, ''.join(parts)))
            name, parts = line[1:].split()[0], []
        elif line:
            parts.append(line.upper())
    if name is not None:
        records.append((name, ''.join(parts)))
    return records


def allele_number(name):
    return int(name.rsplit('_', 1)[1])


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument('--check', action='store_true', help='Say what would change, and write nothing')
    args = parser.parse_args()

    scheme = json.loads(fetch('{}/schemes/{}'.format(API, SCHEME)))
    print('{}: "{}", last updated {}'.format(API, scheme['description'], scheme.get('last_updated', 'unknown')))
    if scheme.get('message'):
        print('PubMLST says: {}'.format(scheme['message'].strip()))
    loci = [url.rsplit('/', 1)[1] for url in scheme['loci']]
    print('{} loci: {}'.format(len(loci), ', '.join(loci)))

    alleles = []
    for locus in loci:
        downloaded = read_fasta(fetch('{}/loci/{}/alleles_fasta'.format(API, locus)))
        records = sorted(downloaded, key=lambda record: allele_number(record[0]))
        print('  {}: {} alleles'.format(locus, len(records)))
        alleles += records

    profiles = fetch('{}/schemes/{}/profiles_csv'.format(API, SCHEME)).replace('\r\n', '\n').strip('\n')
    rows = [line.split('\t') for line in profiles.split('\n') if line.strip()]
    header = rows[0]
    if header != ['ST'] + loci:
        raise SystemExit('The profile table has columns {}, and the scheme has loci {}'.format(header, loci))
    print('{} sequence types'.format(len(rows) - 1))

    # Every allele a profile uses must be one of the alleles downloaded, or a sample could never be typed
    available = {name for name, _ in alleles}
    used = {'{}_{}'.format(locus, allele) for row in rows[1:]
            for locus, allele in zip(header[1:], row[1:], strict=True)}
    missing = sorted(used - available)
    if missing:
        raise SystemExit('{} allele(s) used by a profile were not downloaded: {}'.format(len(missing), missing))

    fasta = ''.join('>{}\n{}\n'.format(name, '\n'.join(sequence[i:i + WRAP] for i in range(0, len(sequence), WRAP)))
                    for name, sequence in alleles)
    table = ''.join('\t'.join(row) + '\n' for row in rows)
    for path, content in ((DATA / 'mlst_alleles.fasta', fasta), (DATA / 'mlst_profiles.tsv', table)):
        old = path.read_text() if path.exists() else ''
        if old == content:
            print('{}: unchanged'.format(path.name))
        elif args.check:
            print('{}: would change ({:,} -> {:,} bytes)'.format(path.name, len(old), len(content)))
        else:
            path.write_text(content)
            print('{}: written ({:,} -> {:,} bytes)'.format(path.name, len(old), len(content)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
