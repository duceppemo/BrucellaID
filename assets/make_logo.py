#!/usr/bin/env python3
"""
Generate the logo: assets/logo.svg and assets/logo_dark.svg for the README and the wiki, and
brucellaid/data/logo.svg and logo.png for the reports BrucellaID writes.

The text is written as outlines rather than as <text>, so the logo renders identically everywhere: a viewer
without DejaVu Sans would otherwise fall back to a wider font and the tagline would run past the canvas. The
canvas is then trimmed to the ink, so the logo centres properly wherever it is placed.

    python3 assets/make_logo.py        # needs fonttools, cairosvg, pillow and the DejaVu fonts

The copy in the package is the same drawing with its four colours written as CSS variables, so that one
file serves the light and the dark report instead of two; a viewer that knows nothing of the variables,
such as a browser opening the file on its own, gets the light colours from the fallbacks. The png is for
the PDF report, which reportlab draws without needing anything that can read an svg.

Edit the layout below and run it again; the generated files are what the repository ships.
"""
import math
import pathlib
import sys

from fontTools.misc.transform import Transform
from fontTools.pens.svgPathPen import SVGPathPen
from fontTools.pens.transformPen import TransformPen
from fontTools.ttLib import TTFont

REGULAR = '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'
BOLD = '/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf'

# Colours: the sample and "ID" share one plum, so the name says what the icon shows
INK = '#2f2a3d'        # "Brucella"
ACCENT = '#7d4a9e'     # the sample at the centre, the reference it is closest to, and "ID"
REFERENCE = '#9aa5ad'  # the other reference genomes of the sketch
MUTED = '#5d6a72'      # the tagline
DARK = {INK: '#f0edf5', ACCENT: '#b185d8', REFERENCE: '#6d787f', MUTED: '#a9b6bd'}
# The same four colours as CSS variables, for the copy the reports inline
THEMED = {INK: 'var(--logo-ink, {})'.format(INK), ACCENT: 'var(--logo-accent, {})'.format(ACCENT),
          REFERENCE: 'var(--logo-reference, {})'.format(REFERENCE),
          MUTED: 'var(--logo-muted, {})'.format(MUTED)}
PNG_WIDTH = 1200  # Wide enough for the PDF report at 300 dpi

TAGLINE = 'Brucella species identification and MLST from reads'
MARGIN = 10  # Uniform margin around the ink, in user units

# The icon: the 16 reference genomes of the sketch around the sample, with the closest one highlighted
CENTRE = (62, 66)
RADIUS = 46
REFERENCES = 16
CLOSEST = 2  # Which of the 16 the sample is closest to, counting clockwise from the top


def short(value):
    """One decimal is plenty at these sizes, and it halves the file."""
    return '{:.1f}'.format(value).rstrip('0').rstrip('.') or '0'


def runs_to_paths(runs, x, baseline, size, font_path, tracking=0.0):
    """Each (text, colour) run as a <path>, laid out from x along the baseline."""
    font = TTFont(font_path)
    glyphs, cmap, upem = font.getGlyphSet(), font.getBestCmap(), font['head'].unitsPerEm
    scale = size / upem
    paths = []
    for text, colour in runs:
        pen = SVGPathPen(glyphs, ntos=short)
        for character in text:
            glyph = glyphs[cmap[ord(character)]]
            glyph.draw(TransformPen(pen, Transform(scale, 0, 0, -scale, x, baseline)))
            x += glyph.width * scale + tracking
        paths.append('  <path fill="{}" d="{}"/>'.format(colour, pen.getCommands()))
    return paths


def icon():
    """The sample at the centre, every reference genome around it, and the distance to the closest one."""
    points = []
    for i in range(REFERENCES):
        angle = -math.pi / 2 + 2 * math.pi * i / REFERENCES
        points.append((CENTRE[0] + RADIUS * math.cos(angle), CENTRE[1] + RADIUS * math.sin(angle)))
    closest = points[CLOSEST]
    others = [point for i, point in enumerate(points) if i != CLOSEST]
    return [
        '  <!-- The reference genomes of the MASH sketch, around the sample -->',
        '  <g fill="{}">'.format(REFERENCE),
        *['    <circle cx="{}" cy="{}" r="5"/>'.format(short(x), short(y)) for x, y in others],
        '  </g>',
        '  <!-- The distance to the closest reference, and that reference -->',
        '  <path d="M{} {}L{} {}" stroke="{}" stroke-width="4" stroke-linecap="round"/>'.format(
            short(CENTRE[0]), short(CENTRE[1]), short(closest[0]), short(closest[1]), ACCENT),
        '  <circle cx="{}" cy="{}" r="8" fill="{}"/>'.format(short(closest[0]), short(closest[1]), ACCENT),
        '  <!-- The sample -->',
        '  <circle cx="{}" cy="{}" r="11" fill="{}"/>'.format(short(CENTRE[0]), short(CENTRE[1]), INK),
    ]


def build(view_box):
    wordmark = runs_to_paths([('Brucella', INK), ('ID', ACCENT)], 130, 66, 52, BOLD, tracking=-1)
    tagline = runs_to_paths([(TAGLINE, MUTED)], 132, 100, 18, REGULAR)
    width, height = view_box[2], view_box[3]
    return '\n'.join([
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="{} {} {} {}" width="{:.0f}" height="{:.0f}" role="img"'
        .format(*('{:.1f}'.format(v) for v in view_box), width, height),
        '     aria-label="BrucellaID: {}">'.format(TAGLINE),
        '  <title>BrucellaID</title>',
        '  <!-- Generated by assets/make_logo.py: the text is outlines, so no font is needed to render it -->',
        *icon(), *wordmark, *tagline,
        '</svg>', ''])


def ink_box(svg, scale=4):
    """The bounding box of everything drawn, in user units, by rendering the svg and looking at the pixels."""
    import io

    import cairosvg
    from PIL import Image, ImageChops
    png = cairosvg.svg2png(bytestring=svg.encode(), scale=scale, background_color='white')
    image = Image.open(io.BytesIO(png)).convert('RGB')
    box = ImageChops.difference(image, Image.new('RGB', image.size, 'white')).getbbox()
    offset_x, offset_y = float(svg.split('viewBox="')[1].split()[0]), float(svg.split('viewBox="')[1].split()[1])
    return [box[0] / scale + offset_x, box[1] / scale + offset_y,
            box[2] / scale + offset_x, box[3] / scale + offset_y]


def write_png(svg, path, width=PNG_WIDTH):
    """The light logo as a transparent png, for the PDF report."""
    import cairosvg
    cairosvg.svg2png(bytestring=svg.encode(), write_to=str(path), output_width=width)


def main():
    # First pass on a canvas large enough for anything, then trim it to the ink
    box = ink_box(build([0, 0, 900, 200]))
    view_box = [box[0] - MARGIN, box[1] - MARGIN, box[2] - box[0] + 2 * MARGIN, box[3] - box[1] + 2 * MARGIN]
    svg = build(view_box)

    assets = pathlib.Path(__file__).resolve().parent
    (assets / 'logo.svg').write_text(svg)
    dark = svg
    for light, replacement in DARK.items():
        dark = dark.replace(light, replacement)
    (assets / 'logo_dark.svg').write_text(dark)

    # The copies the reports use, inside the package
    data = assets.parent / 'brucellaid' / 'data'
    themed = svg
    for light, replacement in THEMED.items():
        themed = themed.replace('"{}"'.format(light), '"{}"'.format(replacement))
    (data / 'logo.svg').write_text(themed)
    write_png(svg, data / 'logo.png')

    print('viewBox {:.1f} {:.1f} {:.1f} {:.1f}'.format(*view_box))
    for path in (assets / 'logo.svg', assets / 'logo_dark.svg', data / 'logo.svg', data / 'logo.png'):
        print('{}: {:,} bytes'.format(path.relative_to(assets.parent), path.stat().st_size))
    return 0


if __name__ == '__main__':
    sys.exit(main())
