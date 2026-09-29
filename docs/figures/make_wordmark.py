#!/usr/bin/env python3
"""Write the FlutterWing wordmark assets as path-only SVGs (no font needed to view them).

    uv run --with fonttools --with uharfbuzz docs/figures/make_wordmark.py [--out docs/figures]

The mark is the V-g crossing: the blue damping trace dipping through the dashed zero line, the
hollow coral circle at the flutter crossing, drawn on a 64 x 64 grid. Writes wordmark.svg (mark +
"FlutterWing"), mark.svg (the mark alone) and social_preview.svg (1280 x 640, for the repository's
social preview). wordmark.svg and mark.svg are opaque white cards like every other asset, so they
stay readable on GitHub's dark theme; the social preview sits on paper. social_preview.png is
rasterised afterwards with resvg, the rasteriser of make_coupling_diagram.py: `uv run --with resvg-py
python -c "import resvg_py; open('docs/figures/social_preview.png','wb').write(bytes(
resvg_py.svg_to_bytes(svg_path='docs/figures/social_preview.svg', width=1280, height=640)))"`
(Quick Look's qlmanage crops non-square SVGs).
Text is set in Latin Modern Roman, regular weight throughout (the OpenType Computer Modern, GUST
Font License; the identity uses no bold): the fonts are
fetched once from CTAN (mirrors.ctan.org/fonts/lm.zip) into a cache directory, shaped with
HarfBuzz and converted to outlines with fontTools. Colours are the identity tokens of
util/fw_style.m: ink #1D1D1F, blue #0066CC, coral #E35336, paper #F5F5F7, grey #D2D2D7,
muted #6E6E73.
"""
import argparse, io, os, sys, tempfile, urllib.request, zipfile
import uharfbuzz as hb
from fontTools.ttLib import TTFont
from fontTools.pens.svgPathPen import SVGPathPen
from fontTools.pens.transformPen import TransformPen

INK, BLUE, CORAL, PAPER, GREY, MUTED, WHITE = '#1D1D1F', '#0066CC', '#E35336', '#F5F5F7', '#D2D2D7', '#6E6E73', '#FFFFFF'
FACES = {'bold': 'lmroman10-bold.otf', 'regular': 'lmroman10-regular.otf', 'italic': 'lmroman10-italic.otf',
         'mono': 'lmmono10-regular.otf'}
LM_URL = 'https://mirrors.ctan.org/fonts/lm.zip'


def fonts_dir():
    d = os.path.join(tempfile.gettempdir(), 'flutterwing-latin-modern')
    if all(os.path.exists(os.path.join(d, f)) for f in FACES.values()):
        return d
    os.makedirs(d, exist_ok=True)
    print('fetching Latin Modern from CTAN ...', file=sys.stderr)
    data = urllib.request.urlopen(LM_URL, timeout=120).read()
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        for f in FACES.values():
            with z.open('lm/fonts/opentype/public/lm/' + f) as src, open(os.path.join(d, f), 'wb') as dst:
                dst.write(src.read())
    return d


class Face:
    def __init__(self, path):
        self.tt = TTFont(path); self.gs = self.tt.getGlyphSet(); self.order = self.tt.getGlyphOrder()
        self.upem = self.tt['head'].unitsPerEm
        self.hb = hb.Font(hb.Face(hb.Blob.from_file_path(path)))
        self.cap = self.tt['OS/2'].sCapHeight / self.upem
        self.asc = self.tt['hhea'].ascent / self.upem; self.desc = -self.tt['hhea'].descent / self.upem

    def run(self, text, size, x=0.0, y=0.0):
        """Outline path data of `text` at `size` px with the baseline at (x, y); returns (d, advance)."""
        buf = hb.Buffer(); buf.add_str(text); buf.guess_segment_properties()
        hb.shape(self.hb, buf, {'kern': True, 'liga': True})
        s = size / self.upem; parts = []; cx = x
        for info, pos in zip(buf.glyph_infos, buf.glyph_positions):
            pen = SVGPathPen(self.gs, ntos=lambda v: f'{v:.2f}')
            self.gs[self.order[info.codepoint]].draw(TransformPen(pen, (s, 0, 0, -s, cx + pos.x_offset * s, y - pos.y_offset * s)))
            d = pen.getCommands()
            if d: parts.append(d)
            cx += pos.x_advance * s
        return ' '.join(parts), cx - x


def mark_svg(px, x=0, y=0):
    """The mark on a 64 x 64 design grid, placed at (x, y) with side length px: the V-g crossing."""
    s = px / 64
    body = (f'<line x1="5" y1="40" x2="59" y2="40" stroke="{INK}" stroke-width="2.2" stroke-dasharray="4 3.2" stroke-linecap="round"/>'
            f'<path d="M5 15 C 22 16, 30 24, 37 40 S 50 59, 59 61" fill="none" stroke="{BLUE}" stroke-width="4.6" stroke-linecap="round"/>'
            f'<circle cx="37" cy="40" r="6" fill="{WHITE}" stroke="{CORAL}" stroke-width="3.8"/>')
    return f'<g transform="translate({x:.2f} {y:.2f}) scale({s:.4f})">{body}</g>'


def lockup(bold, cap_px, x=0, y=0):
    """Mark + wordmark with the cap height cap_px and the baseline at y; returns (svg, width, top, bottom)."""
    size = cap_px / bold.cap
    mark_px = 1.28 * cap_px; gap = 0.38 * cap_px
    # the mark's visual centre sits on the cap-height middle; its 64-grid content spans roughly y = 12..62
    my = y - cap_px / 2 - mark_px * 37 / 64
    d, adv = bold.run('FlutterWing', size, x + mark_px + gap, y)
    svg = mark_svg(mark_px, x, my) + f'<path fill="{INK}" d="{d}"/>'
    return svg, mark_px + gap + adv, min(my, y - cap_px), max(my + mark_px, y + bold.desc * size)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--out', default=os.path.dirname(os.path.abspath(__file__)))
    a = ap.parse_args()
    fd = fonts_dir(); italic = Face(os.path.join(fd, FACES['italic'])); regular = Face(os.path.join(fd, FACES['regular'])); bold = regular   # the wordmark is set in the regular weight
    hdr = f'<!-- FlutterWing identity asset, generated by docs/figures/make_wordmark.py. Text: Latin Modern Roman outlines (GUST Font License). -->'

    # mark alone, an opaque white card
    with open(os.path.join(a.out, 'mark.svg'), 'w') as f:
        f.write(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 64 64" width="64" height="64">{hdr}<rect width="64" height="64" fill="{WHITE}"/>{mark_svg(64)}</svg>\n')

    # wordmark lockup, cap height 40 px, an opaque white card
    cap = 40; svg, w, top, bottom = lockup(bold, cap, 0, 0)
    pad = 10; vb = f'{-pad} {top - pad:.2f} {w + 2 * pad:.2f} {bottom - top + 2 * pad:.2f}'
    card = f'<rect x="{-pad}" y="{top - pad:.2f}" width="{w + 2 * pad:.2f}" height="{bottom - top + 2 * pad:.2f}" fill="{WHITE}"/>'
    with open(os.path.join(a.out, 'wordmark.svg'), 'w') as f:
        f.write(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{vb}" width="{w + 2 * pad:.0f}" height="{bottom - top + 2 * pad:.0f}" role="img" aria-label="FlutterWing">{hdr}{card}{svg}</svg>\n')

    # social preview 1280 x 640
    W, H = 1280, 640
    parts = [f'<rect width="{W}" height="{H}" fill="{PAPER}"/>']
    lk, lw, _, _ = lockup(bold, 78, 84, 232); parts.append(lk)
    y = 300
    for line in ('A minimal, self-contained MATLAB aeroservoelastic', 'modeling framework for ASE control research'):
        d, _ = italic.run(line, 31, 84, y); parts.append(f'<path fill="{INK}" d="{d}"/>'); y += 40
    d, _ = regular.run('MATLAB · GPL-3.0 · AIAA SciTech 2026 benchmark model', 22, 84, 584); parts.append(f'<path fill="{MUTED}" d="{d}"/>')
    # a faint V-g crossing on the right: the damping curve dips through zero
    pts = []; import math
    for v in range(90, 191):
        g = 4.2 + 1.2 * math.sin((v - 90) / 60) - 0.0032 * max(v - 128, 0) ** 2 * (1 + (v - 128) / 200)
        pts.append((880 + (v - 90) / 100 * 360, 330 - g * 9))
    vc = next(v for v in range(90, 191) if 4.2 + 1.2 * math.sin((v - 90) / 60) - 0.0032 * max(v - 128, 0) ** 2 * (1 + (v - 128) / 200) < 0)
    poly = ' '.join(f'{x:.1f},{y:.1f}' for x, y in pts if 0 < y < H)
    parts.append(f'<line x1="850" y1="330" x2="1250" y2="330" stroke="{INK}" stroke-opacity="0.35" stroke-width="2" stroke-dasharray="10 8"/>')
    parts.append(f'<polyline points="{poly}" fill="none" stroke="{BLUE}" stroke-opacity="0.45" stroke-width="4" stroke-linejoin="round"/>')
    parts.append(f'<circle cx="{880 + (vc - 90) / 100 * 360:.1f}" cy="330" r="13" fill="{PAPER}" stroke="{CORAL}" stroke-width="6"/>')
    with open(os.path.join(a.out, 'social_preview.svg'), 'w') as f:
        f.write(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}">{hdr}{"".join(parts)}</svg>\n')
    print(f'wrote mark.svg, wordmark.svg ({w + 2 * pad:.0f} x {bottom - top + 2 * pad:.0f} px), social_preview.svg in {a.out}')


if __name__ == '__main__':
    main()
