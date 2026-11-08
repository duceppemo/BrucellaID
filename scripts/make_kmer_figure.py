#!/usr/bin/env python3
"""
Draw the figure behind the k-mer copy threshold: assets/kmer_threshold.svg and kmer_threshold_dark.svg.

    python3 scripts/make_kmer_figure.py          # reads scripts/kmer_threshold_sweep.tsv

Two panels, because the rule rests on two facts that no single panel shows.

Left: what each threshold recovers, as a share of the best that library's own sweep reached. A shallow
library peaks early and then collapses -- at 53x, requiring 64 copies leaves 1% of the sketch -- while a
deep one climbs slowly and only arrives at 40. So no single constant can serve both, which is the whole
argument for deriving the threshold from the depth.

Right: the plateau of each of the 81 curves, meaning every threshold within 1% of that curve's best, drawn
as a vertical span at its depth. A rule has to land inside that span using only an estimate of the depth.
The two rules are drawn over it: one more copy per 20x capped at 20, which runs along the bottom and falls
out of the span for most of the range, and one more copy per 5x capped at 48, which stays inside it.

No matplotlib: the repository depends on the standard library, and assets/make_logo.py already writes its
SVG by hand. Text is left as <text> rather than outlines because, unlike the logo, this figure is only ever
read on a page where a missing font degrades to a different font rather than to nothing.
"""
import pathlib
import statistics
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
SWEEP = ROOT / 'scripts' / 'kmer_threshold_sweep.tsv'
TOLERANCE = 0.99          # inside the plateau means within 1% of that curve's best recovery
FLOOR, SLOPE, CEILING = 2, 5, 48      # the rule this figure justifies
OLD_FLOOR, OLD_SLOPE, OLD_CEILING = 2, 20, 20

# The validated palette. The depth ramp is one hue stepped light to dark (ordinal: magnitude with an
# order), and each mode gets its own steps rather than a flip of the other's, because the step nearest
# the surface disappears into it. The two rules are categorical slots 1 and 2.
LIGHT = {
    'surface': '#fcfcfb', 'ink': '#0b0b0b', 'muted': '#52514e', 'grid': '#e5e4e0',
    'ramp': ['#86b6ef', '#5598e7', '#2a78d6', '#1c5cab', '#0d366b'],
    'band': '#c8c7c0', 'new': '#2a78d6', 'old': '#eb6834',
}
DARK = {
    'surface': '#1a1a19', 'ink': '#f0efec', 'muted': '#c3c2b7', 'grid': '#333330',
    'ramp': ['#9ec5f4', '#6da7ec', '#3987e5', '#256abf', '#184f95'],
    'band': '#4a4945', 'new': '#3987e5', 'old': '#d95926',
}
BANDS = [(45, 60, '53x'), (95, 120, '105x'), (190, 260, '230x'), (400, 560, '500x'), (1200, 2000, '1500x')]

W, H = 1060, 430
PAD_L, PAD_R, PAD_T, PAD_B = 58, 152, 54, 58
GAP = 70
PANEL_W = (W - PAD_L - PAD_R - GAP) / 2


def rule(depth, floor, slope, ceiling):
    return int(min(ceiling, max(floor, round(depth / slope))))


# The figure is of the real-library sweep only. The simulated measurements live in their own file
# (scripts/kmer_threshold_simulated.tsv) because they set nothing, and they are swept to a different
# maximum threshold, so mixing them in would both change the curve count and censor the spans. This set is
# stated rather than inferred so that adding a source to the data file cannot silently redraw the figure.
REAL_SOURCES = ('subsample', 'ena')


def read_sweep(path):
    """[(library, depth, {copies: shared})] -- one entry per library per depth, real libraries only."""
    curves = {}
    for line in path.read_text().splitlines():
        if line.startswith('#') or line.startswith('source\t'):
            continue
        source, library, depth, copies, shared, _total = line.split('\t')
        if source not in REAL_SOURCES:
            continue
        key = (source, library, round(float(depth)))
        curves.setdefault(key, [float(depth), {}])[1][int(copies)] = int(shared)
    out = []
    for (_source, library, _d), (depth, curve) in curves.items():
        if len(curve) >= 6 and max(curve.values()) > 0:
            out.append((library, depth, curve))
    return sorted(out, key=lambda c: c[1])


def plateau(curve):
    peak = max(curve.values())
    inside = sorted(m for m, s in curve.items() if s >= peak * TOLERANCE)
    return inside[0], inside[-1]


def esc(text):
    return str(text).replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')


class Axes:
    """A panel: user space in, pixels out, on a log x and a linear y."""

    def __init__(self, x0, y0, w, h, xlo, xhi, ylo, yhi):
        import math
        self.x0, self.y0, self.w, self.h = x0, y0, w, h
        self.lxlo, self.lxhi = math.log10(xlo), math.log10(xhi)
        self.ylo, self.yhi = ylo, yhi
        self._log10 = math.log10

    def x(self, v):
        return self.x0 + self.w * (self._log10(v) - self.lxlo) / (self.lxhi - self.lxlo)

    def y(self, v):
        return self.y0 + self.h * (1 - (v - self.ylo) / (self.yhi - self.ylo))


def frame(parts, ax, c, xticks, yticks, xlabel, ylabel, xfmt=str, yfmt=str):
    """Recessive grid and axis labels: the data carries the ink, the frame stays out of the way."""
    for v in yticks:
        y = ax.y(v)
        parts.append('<line x1="{:.1f}" y1="{:.1f}" x2="{:.1f}" y2="{:.1f}" stroke="{}" stroke-width="1"/>'
                     .format(ax.x0, y, ax.x0 + ax.w, y, c['grid']))
        parts.append('<text x="{:.1f}" y="{:.1f}" fill="{}" font-size="11" text-anchor="end" '
                     'dominant-baseline="middle">{}</text>'.format(ax.x0 - 8, y, c['muted'], esc(yfmt(v))))
    for v in xticks:
        x = ax.x(v)
        parts.append('<text x="{:.1f}" y="{:.1f}" fill="{}" font-size="11" text-anchor="middle">{}</text>'
                     .format(x, ax.y0 + ax.h + 17, c['muted'], esc(xfmt(v))))
    parts.append('<text x="{:.1f}" y="{:.1f}" fill="{}" font-size="12" text-anchor="middle">{}</text>'
                 .format(ax.x0 + ax.w / 2, ax.y0 + ax.h + 40, c['muted'], esc(xlabel)))
    parts.append('<text transform="translate({:.1f} {:.1f}) rotate(-90)" fill="{}" font-size="12" '
                 'text-anchor="middle">{}</text>'
                 .format(ax.x0 - 40, ax.y0 + ax.h / 2, c['muted'], esc(ylabel)))


def panel_recovery(parts, curves, ax, c):
    """Median recovery against threshold, one line per depth band, each labelled at its own end."""
    frame(parts, ax, c, [1, 2, 4, 8, 16, 32, 64], [0, 25, 50, 75, 100],
          'k-mer copies required (--min-kmer-copies)', 'of the best this library reached',
          yfmt=lambda v: '{:.0f}%'.format(v))
    ends = []
    for (lo_d, hi_d, label), colour in zip(BANDS, c['ramp'], strict=True):
        sel = [cv for _lib, d, cv in curves if lo_d <= d < hi_d]
        if not sel:
            continue
        pts = []
        for m in sorted({m for cv in sel for m in cv}):
            vals = [100 * cv[m] / max(cv.values()) for cv in sel if m in cv]
            if len(vals) >= max(2, len(sel) // 2):
                pts.append((ax.x(m), ax.y(statistics.median(vals))))
        if len(pts) < 2:
            continue
        path = 'M' + ' L'.join('{:.1f} {:.1f}'.format(x, y) for x, y in pts)
        parts.append('<path d="{}" fill="none" stroke="{}" stroke-width="2" stroke-linecap="round" '
                     'stroke-linejoin="round"/>'.format(path, colour))
        ends.append([pts[-1][0], pts[-1][1], label, colour])

    # Direct labels, so identity never rests on colour alone -- but three of the bands finish at the same
    # 100% and their labels would sit on top of one another, so they are pushed apart and joined to their
    # own line by a leader.
    ends_raw = [list(e) for e in ends]
    ends.sort(key=lambda e: e[1])
    for i in range(1, len(ends)):
        if ends[i][1] - ends[i - 1][1] < 14:
            ends[i][1] = ends[i - 1][1] + 14
    for x, ly, label, colour in ends:
        y = next(e[1] for e in ends_raw if e[2] == label)
        parts.append('<circle cx="{:.1f}" cy="{:.1f}" r="3.5" fill="{}" stroke="{}" stroke-width="2"/>'
                     .format(x, y, colour, c['surface']))
        if abs(ly - y) > 1:
            parts.append('<path d="M{:.1f} {:.1f} L{:.1f} {:.1f}" fill="none" stroke="{}" '
                         'stroke-width="1" opacity="0.55"/>'.format(x + 4, y, x + 13, ly, colour))
        parts.append('<text x="{:.1f}" y="{:.1f}" fill="{}" font-size="11" dominant-baseline="middle">'
                     '{}</text>'.format(x + 16, ly, c['ink'], esc(label)))


def panel_plateau(parts, curves, ax, c):
    """Every curve's plateau as a span at its depth, with the two rules drawn over it."""
    frame(parts, ax, c, [5, 10, 50, 100, 500, 2000], [0, 8, 16, 24, 32, 40, 48, 56, 64],
          'sequencing depth', 'k-mer copies required',
          xfmt=lambda v: '{:.0f}x'.format(v))
    for _lib, depth, curve in curves:
        lo, hi = plateau(curve)
        x = ax.x(depth)
        parts.append('<line x1="{:.1f}" y1="{:.1f}" x2="{:.1f}" y2="{:.1f}" stroke="{}" stroke-width="3" '
                     'stroke-linecap="round" opacity="0.75"/>'
                     .format(x, ax.y(lo), x, ax.y(hi), c['band']))
    for label, floor, slope, ceiling, colour, dash in (
            ('one per 20x, cap 20', OLD_FLOOR, OLD_SLOPE, OLD_CEILING, c['old'], '5 4'),
            ('one per 5x, cap 48', FLOOR, SLOPE, CEILING, c['new'], '')):
        pts = []
        d = 5.0
        while d <= 2000:
            pts.append((ax.x(d), ax.y(rule(d, floor, slope, ceiling))))
            d *= 1.06
        path = 'M' + ' L'.join('{:.1f} {:.1f}'.format(x, y) for x, y in pts)
        parts.append('<path d="{}" fill="none" stroke="{}" stroke-width="2" stroke-linejoin="round"{}/>'
                     .format(path, colour, ' stroke-dasharray="{}"'.format(dash) if dash else ''))
        lx, ly = pts[-1]
        parts.append('<text x="{:.1f}" y="{:.1f}" fill="{}" font-size="11" dominant-baseline="middle">'
                     '{}</text>'.format(lx + 7, ly, c['ink'], esc(label)))


def build(curves, c):
    parts = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {} {}" width="{}" height="{}" '
             'role="img" aria-label="How much of the reference sketch each k-mer copy threshold recovers, '
             'and where the plateau sits against sequencing depth">'.format(W, H, W, H),
             '<rect width="{}" height="{}" fill="{}"/>'.format(W, H, c['surface'])]
    plateau_h = H - PAD_T - PAD_B
    left = Axes(PAD_L, PAD_T, PANEL_W, plateau_h, 1, 64, 0, 104)
    right = Axes(PAD_L + PANEL_W + GAP, PAD_T, PANEL_W, plateau_h, 5, 2000, 0, 68)
    parts.append('<text x="{}" y="26" fill="{}" font-size="15" font-weight="600">What the k-mer copy '
                 'threshold recovers, and what it has to hit</text>'.format(PAD_L - 40, c['ink']))
    parts.append('<text x="{}" y="44" fill="{}" font-size="11.5">81 sweeps over 63 real Brucella '
                 'libraries, 5x to 1903x. Left: median recovery per depth band. Right: every sweep’s '
                 'plateau, and the two rules.</text>'.format(PAD_L - 40, c['muted']))
    panel_recovery(parts, curves, left, c)
    panel_plateau(parts, curves, right, c)
    parts.append('</svg>')
    return '\n'.join(parts) + '\n'


def main():
    if not SWEEP.exists():
        raise SystemExit('{} is missing: it holds the measurements this figure draws.'.format(SWEEP))
    curves = read_sweep(SWEEP)
    if len(curves) < 20:
        raise SystemExit('only {} usable curves in {}'.format(len(curves), SWEEP))
    for name, colours in (('kmer_threshold.svg', LIGHT), ('kmer_threshold_dark.svg', DARK)):
        path = ROOT / 'assets' / name
        path.write_text(build(curves, colours), encoding='utf-8')
        print('{}: {:,} bytes'.format(path, path.stat().st_size), file=sys.stderr)
    print('{} curves from {} libraries, {:.0f}x to {:.0f}x'.format(
        len(curves), len({c[0] for c in curves}), curves[0][1], curves[-1][1]), file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main())
