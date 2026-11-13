#!/usr/bin/env python3
"""
Build the ST to species table of the package from the PubMLST Brucella isolate database.

    python3 scripts/make_st_species.py             # writes brucellaid/data/mlst_st_species.tsv
    python3 scripts/make_st_species.py --check     # says what would change, and writes nothing
    python3 scripts/make_st_species.py --report    # prints which species the STs separate

The "MLST (9 locus)" profile table of PubMLST (scripts/make_mlst_data.py) says which alleles make an ST, but
not which species carries it: the scheme has no species field. The isolate database does. Every isolate record
carries the species and the biovar the submitter gave it, and the ST the curators assigned, so counting the
isolates of each ST gives the species that ST is seen in, and how consistently.

That is what lets the ST arbitrate an identification MASH cannot make. The classical Brucella species are
0.02 to 0.4% apart: in the 429-genome reference set BrucellaID identifies against by default, the nearest
B. canis and B. suis references are 0.000264 apart, and in vSNP's 16 the pair is 0.00024 apart. Either way
that is inside the variation within either species, so no sketch separates them. The ST does,
where an ST belongs to one species only -- which is exactly what the "purity" column of this table reports, and
what --report measures.

The species labels are the submitters', and are the same labels that turned out to be wrong for some NCBI
assemblies, so purity is reported rather than assumed: an ST whose isolates disagree is a fact about the data,
and the caller should refuse to arbitrate with it instead of taking a majority of three.

Without an API key the data is the one submitted up to a cut-off date, which the API states and this script
prints. Note it in wiki/Changelog.md and in brucellaid/data/README.md.

The isolate records are fetched one by one -- the REST API has no bulk export -- so they are cached on disk and
a second run costs nothing. Delete the cache folder to fetch them again.
"""
import argparse
import collections
import json
import pathlib
import sys
import time
import urllib.error
import urllib.request

API = 'https://rest.pubmlst.org/db/pubmlst_brucella_isolates'
SCHEME_DESCRIPTION = 'MLST (9 locus)'
DATA = pathlib.Path(__file__).resolve().parent.parent / 'brucellaid' / 'data'
TABLE = 'mlst_st_species.tsv'
HEADER = ['ST', 'Species', 'Purity', 'Isolates', 'Biovars', 'OtherSpecies']
# PubMLST rate limits the API: five requests at a time is answered with "429 Too Many Requests" within a
# second. The records are read one at a time with a pause between them, which is slow (half an hour for the
# whole database) but only has to be done once, because they are cached.
DELAY = 0.25
RETRIES = 5


def fetch(url, delay=DELAY, retries=RETRIES):
    """
    The body of a PubMLST URL, as text.

    A 429 is the server asking us to slow down, so it is waited out with a doubling pause rather than treated
    as a failure. Anything else is reported as it is.
    """
    request = urllib.request.Request(url, headers={'User-Agent': 'brucellaid make_st_species.py'})
    wait = max(delay, 1.0)
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                return response.read().decode('utf-8')
        except urllib.error.HTTPError as e:
            if e.code != 429 or attempt == retries - 1:
                raise SystemExit('Cannot read {}: {}'.format(url, e)) from None
            print('  {}: waiting {:.0f}s'.format(e, wait), file=sys.stderr)
            time.sleep(wait)
            wait *= 2
        except (urllib.error.URLError, TimeoutError) as e:
            raise SystemExit('Cannot read {}: {}'.format(url, e)) from None
    raise SystemExit('Cannot read {}: too many retries'.format(url))


def isolate_urls():
    """Every isolate URL of the database, and the message the API attaches to the listing."""
    listing = json.loads(fetch('{}/isolates?return_all=1'.format(API)))
    return listing['isolates'], listing.get('message', '')


def read_isolate(url, cache, delay=DELAY):
    """One isolate record, from the cache folder if it is already there."""
    path = cache / '{}.json'.format(url.rstrip('/').rsplit('/', 1)[-1])
    if path.is_file():
        try:
            return json.loads(path.read_text()), True
        except ValueError:
            pass  # A truncated file from an interrupted run: fetch it again.
    record = json.loads(fetch(url, delay=delay))
    path.write_text(json.dumps(record))
    time.sleep(delay)
    return record, False


def read_isolates(urls, cache, delay=DELAY):
    """Every isolate record, fetching the ones the cache does not have."""
    cache.mkdir(parents=True, exist_ok=True)
    records, fetched = [], 0
    for done, url in enumerate(urls, start=1):
        record, cached = read_isolate(url, cache, delay=delay)
        records.append(record)
        fetched += not cached
        if done % 250 == 0 or done == len(urls):
            print('  {}/{} isolates ({} fetched, {} from the cache)'.format(
                done, len(urls), fetched, done - fetched), file=sys.stderr)
    return records


def species_name(label):
    """
    "B. abortus" written out in full, as the reference table writes it.

    PubMLST abbreviates the genus and sometimes adds a biovar, a biotype or a note to the species field; only
    the species is kept here, because that is what the table is keyed on.
    """
    label = ' '.join(label.split())
    if label.startswith('B. '):
        label = 'Brucella ' + label[3:]
    words = label.split()
    if len(words) >= 2 and words[0] == 'Brucella':
        return 'Brucella {}'.format(words[1].rstrip(','))
    return label


def designated_st(record):
    """The ST of the nine-locus scheme, or None when the isolate has none."""
    for scheme in record.get('schemes', []):
        if scheme.get('description') == SCHEME_DESCRIPTION:
            st = scheme.get('fields', {}).get('ST')
            return str(st) if st is not None else None
    return None


def observations(records):
    """
    [(ST, species, biovar)] of the isolates that have both an ST and a species.

    An isolate with no ST (its alleles were never completed) or no species says nothing about the mapping, and
    is counted as skipped by the caller.
    """
    seen = []
    for record in records:
        st = designated_st(record)
        provenance = record.get('provenance', {})
        label = provenance.get('Brucella_species', '')
        if st is None or not label:
            continue
        seen.append((st, species_name(label), str(provenance.get('biovar', '') or '')))
    return seen


def tabulate(seen):
    """
    One row per ST: the species of its isolates, in HEADER order.

    "Purity" is the share of the isolates of that ST that carry the commonest species. 1.000 means every
    isolate of the ST agrees; anything less is a disagreement the caller has to decide what to do with.
    """
    by_st = collections.defaultdict(list)
    for st, species, biovar in seen:
        by_st[st].append((species, biovar))

    rows = []
    for st in sorted(by_st, key=lambda s: (len(s), s)):
        isolates = by_st[st]
        counts = collections.Counter(species for species, _ in isolates)
        # Commonest first, then alphabetically, so a tie is broken the same way on every machine.
        ranked = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
        species, top = ranked[0]
        biovars = collections.Counter(biovar for s, biovar in isolates if s == species and biovar)
        rows.append([
            st,
            species,
            '{:.3f}'.format(top / len(isolates)),
            str(len(isolates)),
            ';'.join('{}:{}'.format(b, n) for b, n in sorted(biovars.items())),
            ';'.join('{}:{}'.format(s, n) for s, n in ranked[1:]),
        ])
    return rows


def write_table(rows, message, path):
    """The table, with the header comment that says where it came from."""
    lines = [
        '# The species of each ST of the nine-locus Brucella MLST scheme, counted over the isolates of the',
        '# PubMLST Brucella isolate database. Built by scripts/make_st_species.py.',
        '# Species: the commonest species among the isolates of that ST. Purity: the share of them that carry',
        '# it, so 1.000 is an ST seen in one species only. Biovars: the biovars recorded for that species, with',
        '# their counts. OtherSpecies: the rest, with their counts.',
    ]
    if message:
        lines.append('# PubMLST: {}'.format(' '.join(message.split())))
    lines.append('\t'.join(HEADER))
    lines.extend('\t'.join(row) for row in rows)
    text = '\n'.join(lines) + '\n'
    if path is None:
        return text
    path.write_text(text)
    return text


def report(rows, seen, skipped):
    """What the table says about the species the STs can and cannot separate."""
    print('\n{} isolates with an ST and a species, {} without; {} STs'.format(len(seen), skipped, len(rows)))

    pure = [r for r in rows if r[2] == '1.000']
    print('{} of {} STs are seen in one species only'.format(len(pure), len(rows)))

    # An ST shared by two species is an ST that cannot arbitrate between them. Reported with the isolate
    # counts, because one isolate of a second species is more likely a wrong label than a real sharing.
    shared = [r for r in rows if r[2] != '1.000']
    if shared:
        print('\nSTs seen in more than one species:')
        for st, species, purity, n, _, other in shared:
            print('  ST {:>4}  {:>3} isolates  {} {} ({})  + {}'.format(st, n, purity, species, '', other))

    # The pairs the identification has to tell apart: for each pair of species, whether any ST holds both.
    by_species = collections.defaultdict(set)
    for st, species, _ in seen:
        by_species[species].add(st)
    names = sorted(by_species)
    print('\nSTs per species:')
    for name in names:
        print('  {:<28} {:>3} STs, {:>4} isolates'.format(
            name, len(by_species[name]), sum(1 for _, s, _ in seen if s == name)))

    print('\nSpecies pairs that share at least one ST:')
    overlaps = 0
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            common = by_species[a] & by_species[b]
            if common:
                overlaps += 1
                print('  {} <-> {}: ST {}'.format(a, b, ', '.join(sorted(common, key=lambda s: (len(s), s)))))
    if not overlaps:
        print('  none: every ST belongs to one species')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1],
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--check', action='store_true', help='say what would change, and write nothing')
    parser.add_argument('--report', action='store_true', help='print which species the STs separate')
    parser.add_argument('--cache', metavar='DIR', default=str(pathlib.Path('.pubmlst_isolates')),
                        help='where the isolate records are kept between runs (default: %(default)s)')
    parser.add_argument('--delay', type=float, default=DELAY, metavar='SECONDS',
                        help='pause between requests, to stay under PubMLST\'s rate limit '
                             '(default: %(default)s)')
    args = parser.parse_args(argv)

    print('Reading the isolate list of {}'.format(API), file=sys.stderr)
    urls, message = isolate_urls()
    if message:
        print('PubMLST: {}'.format(' '.join(message.split())), file=sys.stderr)
    records = read_isolates(urls, pathlib.Path(args.cache), delay=args.delay)

    seen = observations(records)
    rows = tabulate(seen)
    text = write_table(rows, message, None)

    path = DATA / TABLE
    if args.check:
        old = path.read_text() if path.is_file() else ''
        print('{}: {}'.format(path, 'unchanged' if old == text else 'would change'))
    else:
        write_table(rows, message, path)
        print('{}: {} STs'.format(path, len(rows)))

    if args.report:
        report(rows, seen, len(records) - len(seen))
    return 0


if __name__ == '__main__':
    sys.exit(main())
