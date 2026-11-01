#!/usr/bin/env python3
"""
Make the example quality assurance report of the wiki from a real run, with the machine it ran on removed.

    python3 scripts/make_example_report.py RUN_FOLDER -o docs/example_report.html

A quality assurance report names the computer it was produced on, the account that produced it and the full
path of every file it touched. That is the point of it, and it is exactly what should not be published in a
repository. This script takes a real report and blanks those values -- and only those -- so that the example
shows what the report holds without publishing anyone's machine. The result is served by GitHub Pages
out of docs/, because GitHub shows an HTML file in a repository as source rather than rendering it.

It reads the values to remove out of brucellaid.json rather than guessing at them with a pattern, so it
removes what that run actually recorded and cannot quietly miss a path in a corner of the report. Each
removed value is replaced by a marker that says what used to be there, and a banner is added at the top of
the page saying that the environment was removed, so the example is never mistaken for a complete report.
"""
import argparse
import html
import json
import pathlib
import sys

REPORT = 'brucellaid_report.html'
JSON = 'brucellaid.json'
BANNER = (
    '<div class="panel noprint" style="border-color:var(--accent);background:var(--accent-soft)">'
    '<b>This is an example report.</b> It is a real run of the four public read sets of the '
    '<a href="https://github.com/duceppemo/BrucellaID/wiki/Tutorial">tutorial</a>, with the computer it ran '
    'on removed: the user account, the host name, the operating system and every directory path read '
    '<i>[removed]</i> below, and the input files are named without the folder they were in. Everything '
    'else is as BrucellaID wrote it. A report of your own names them all, which is what makes it worth '
    'signing.</div>')


def replacements(report):
    """
    {text to remove: what to put in its place}, read from the JSON report of the same run.

    Longest first, so that a working directory is removed before a path that merely starts with it.
    """
    pairs = {}
    # One short marker for all of them: the label to the left of the value already says what is gone,
    # and the footer reads "[removed]@[removed]" rather than running off the line.
    for key in ('working_directory', 'host', 'user', 'system'):
        value = report.get(key)
        if value:
            pairs[value] = '[removed]'
    # A reference data file keeps its name, which is what a reader needs, and loses the folder it sits in
    for item in report.get('reference_data', []):
        path = item.get('path', '')
        if path:
            pairs[path] = '.../brucellaid/data/{}'.format(pathlib.PurePath(path).name)
    # An input file keeps its name too: the report shows the same name in the table just above the path
    for sample in report.get('samples', []):
        for item in sample.get('input_files', []):
            path = item.get('path', '')
            if path:
                pairs[path] = pathlib.PurePath(path).name
        for path in sample.get('files', []):
            pairs.setdefault(path, pathlib.PurePath(path).name)
    return dict(sorted(pairs.items(), key=lambda pair: -len(pair[0])))


def redact(markup, pairs):
    """
    Remove each value from the page, in the plain and the HTML escaped form.

    The report escapes everything it writes, so a path holding an "&" appears escaped in the file; replacing
    both forms means a value cannot survive because of the character it happens to contain.
    """
    removed = {}
    for value, marker in pairs.items():
        for form in (html.escape(value, quote=True), value):
            if form in markup:
                markup = markup.replace(form, html.escape(marker, quote=True))
                removed[value] = removed.get(value, 0) + 1
    return markup, removed


def add_banner(markup):
    """Put the banner under the subtitle, where a reader cannot miss it."""
    anchor = '</p>\n<nav class="toc">'
    if anchor not in markup:
        raise SystemExit('The report does not look like one of brucellaid: no subtitle before the contents.')
    return markup.replace(anchor, '</p>\n{}\n<nav class="toc">'.format(BANNER), 1)


def check(markup, pairs):
    """Make sure nothing that had to go is still in the page."""
    left = [value for value in pairs
            if value in markup or html.escape(value, quote=True) in markup]
    if left:
        raise SystemExit('Still in the page after redacting: {}'.format(', '.join(left)))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument('run', help='the output folder of a run, holding {} and {}'.format(REPORT, JSON))
    parser.add_argument('-o', '--out', default='docs/example_report.html', metavar='FILE',
                        help='where to write the example (default: %(default)s)')
    args = parser.parse_args(argv)

    folder = pathlib.Path(args.run)
    report_path, json_path = folder / REPORT, folder / JSON
    for path in (report_path, json_path):
        if not path.is_file():
            raise SystemExit('{} is not there: run brucellaid with its reports first.'.format(path))

    report = json.loads(json_path.read_text())
    pairs = replacements(report)
    markup, removed = redact(report_path.read_text(encoding='utf-8'), pairs)
    check(markup, pairs)
    markup = add_banner(markup)

    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(markup, encoding='utf-8')
    print('{}: {:,} bytes, {} value(s) removed'.format(out, len(markup), len(removed)), file=sys.stderr)
    for value in pairs:
        print('  {} -> {}'.format(value, pairs[value]), file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main())
