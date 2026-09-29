#!/usr/bin/env python3
"""Draw the DLM-FEM coupling schematic of tutorial chapter 4: docs/figures/DLM_FEM_Coupling.svg and .png.

    uv run --with fonttools --with uharfbuzz --with resvg-py docs/figures/make_coupling_diagram.py [--out DIR]

One beam element (nodes 1 and 2 with the six DOFs d_1..d_6: bending displacement, bending slope and
torsion angle at each node) and one DLM panel in the wing plane ahead of the beam, with its
quarter-chord line (bound vortex points BVP_1, BVP_2, load point l) and three-quarter-chord line
(JP_1, JP_2, control point j). In blue, what the coupling exchanges: the pressure C_p at l becomes
the force F_z and the moment M_y (arm x_l) at the beam station under l and reaches d_1..d_6 through
the shape functions at y_l (build_S_ele); the nodal motion moves j and sets the downwash w_j
(build_DReDIm_ele). Rotations (the slopes d_2, d_5 about x, the twists d_3, d_6 and the moment M_y
about y) are arrows with two heads along their axis, right-hand rule.
Structural frame, x upstream, y spanwise, z down, right-handed, seen by the wing_scene camera (aft of
the wing, above, towards the root: azimuth 37, elevation 26, orthographic), physical up on top.
Text is set from Latin Modern outlines through make_wordmark.py, so the SVG needs no font; the PNG
is the 14 x 9.8 cm canvas at 300 dpi, rasterised by resvg. Colours are the fw_style tokens, on white.
"""
import argparse, math, os, sys
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_wordmark import Face, fonts_dir, FACES, INK, BLUE, GREY, MUTED, WHITE

PT = 72 / 2.54                                  # pt per cm: 1 canvas unit = 1 pt, so sizes are points
W_CM, H_CM, DPI = 14.0, 9.8, 300
W, H = W_CM * PT, H_CM * PT
AZ, EL = 37.0, 26.0                             # the wing_scene view
LW_HAIR, LW_THIN, LW, LW_ARROW, LW_BEAM = 0.5, 0.8, 1.2, 1.8, 2.2
FS_SMALL, FS = 8, 9
DASH, DOT = '3 2', '0.5 1.5'


def camera(az=AZ, el=EL):
    """Screen right and up unit vectors in structural coordinates for MATLAB's view(az, el)."""
    a, e = math.radians(az), math.radians(el)
    r_a = (math.cos(a), math.sin(a), 0.0)                                   # aero frame, z up
    u_a = (-math.sin(a) * math.sin(e), math.cos(a) * math.sin(e), math.cos(e))
    return (-r_a[0], r_a[1], -r_a[2]), (-u_a[0], u_a[1], -u_a[2])         # structural = (-x_a, y_a, -z_a)


R, U = camera()
def dot(a, b): return sum(x * y for x, y in zip(a, b))
def add(*vs): return tuple(sum(c) for c in zip(*vs))
def mul(a, k): return tuple(x * k for x in a)
def proj(p): return (dot(R, p), dot(U, p))


class Tex:
    """Mini TeX: $..$ sets letters in italic, _x and _{xy} are subscripts, \\ldots and \\rightarrow glyphs."""
    SUB, SUB_DY = 0.7, 0.22

    def __init__(self):
        fd = fonts_dir()
        self.rm = Face(os.path.join(fd, FACES['regular']))
        self.it = Face(os.path.join(fd, FACES['italic']))

    def runs(self, s):
        out, math_, i = [], False, 0

        def emit(face, ch, sc=1.0, dy=0.0):
            if out and out[-1][0] is face and out[-1][2] == sc and out[-1][3] == dy:
                out[-1] = (face, out[-1][1] + ch, sc, dy)
            else:
                out.append((face, ch, sc, dy))
        while i < len(s):
            ch = s[i]
            if ch == '$':
                math_ = not math_; i += 1
            elif ch == '\\':
                for name, rep in (('ldots', '…'), ('rightarrow', '→'), ('infty', '∞')):
                    if s.startswith(name, i + 1):
                        emit(self.rm, rep); i += 1 + len(name); break
                else:
                    raise ValueError('unknown macro in ' + s)
            elif ch == '_' and math_:
                if s[i + 1] == '{':
                    j = s.index('}', i); grp, i = s[i + 2:j], j + 1
                else:
                    grp, i = s[i + 1], i + 2
                for g in grp:
                    emit(self.it if g.isalpha() else self.rm, g, self.SUB, self.SUB_DY)
            else:
                emit(self.it if (math_ and ch.isalpha()) else self.rm, ch); i += 1
        return out

    def width(self, s, size):
        return sum(face.run(t, size * sc)[1] for face, t, sc, dy in self.runs(s))

    def svg(self, s, size, x, y, color=INK, anchor='start'):
        x -= {'start': 0.0, 'middle': 0.5, 'end': 1.0}[anchor] * self.width(s, size)
        paths = []
        for face, t, sc, dy in self.runs(s):
            d, adv = face.run(t, size * sc, x, y + dy * size)
            if d:
                paths.append(f'<path fill="{color}" d="{d}"/>')
            x += adv
        return ''.join(paths)


class Canvas:
    def __init__(self, tex, scale, origin):
        self.tex, self.s, self.o, self.parts = tex, scale, origin, []

    def pt(self, p):                              # 3-D structural point -> canvas point (pt, y down)
        sx, sy = proj(p)
        return (self.o[0] + self.s * sx, self.o[1] - self.s * sy)

    @staticmethod
    def _pts(pts):
        return ' '.join(f'{x:.2f},{y:.2f}' for x, y in pts)

    def poly(self, pts, color, w, dash=None, fill='none', close=False):
        tag = 'polygon' if close else 'polyline'
        extra = f' stroke-dasharray="{dash}"' if dash else ''
        self.parts.append(f'<{tag} points="{self._pts(pts)}" fill="{fill}" stroke="{color}" stroke-width="{w}"'
                          f' stroke-linejoin="round" stroke-linecap="round"{extra}/>')

    def line3(self, p, q, color, w, dash=None):
        self.poly([self.pt(p), self.pt(q)], color, w, dash)

    def circle3(self, p, r, fill, stroke=None, w=0.7):
        x, y = self.pt(p)
        st = f' stroke="{stroke}" stroke-width="{w}"' if stroke else ''
        self.parts.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{r}" fill="{fill}"{st}/>')

    def arrow_pts(self, base, tip, color, w, heads=1):
        """Straight arrow with flat heads at its tip: one for a vector, two for a rotation (right-hand rule)."""
        hl = max(4.2, 3.2 * w); hw = 0.4 * hl; step = 0.85 * hl
        tx, ty = tip[0] - base[0], tip[1] - base[1]; n = math.hypot(tx, ty); tx, ty = tx / n, ty / n
        end = (tip[0] - tx * (hl + step * (heads - 1)), tip[1] - ty * (hl + step * (heads - 1)))
        self.poly([base, end], color, w)
        for k in range(heads):
            t = (tip[0] - tx * step * k, tip[1] - ty * step * k)
            b = (t[0] - tx * hl, t[1] - ty * hl)
            self.parts.append(f'<polygon points="{self._pts([t, (b[0] - ty * hw, b[1] + tx * hw), (b[0] + ty * hw, b[1] - tx * hw)])}" fill="{color}"/>')

    def arrow3(self, base, d, L, color, w, heads=1):
        self.arrow_pts(self.pt(base), self.pt(add(base, mul(d, L))), color, w, heads)

    def label(self, s, p, dx=0.0, dy=0.0, size=FS, color=INK, anchor='start'):
        x, y = self.pt(p)
        self.parts.append(self.tex.svg(s, size, x + dx, y + dy, color, anchor))


def draw(cv):
    X, Y, Z = (1, 0, 0), (0, 1, 0), (0, 0, 1)
    UP = (0, 0, -1)
    y1, y2, yl = 2.0, 6.0, 4.0                     # element nodes and the load station along the beam
    xa, xb, ya, yb = 1.6, 2.8, 3.0, 5.0            # panel in the wing plane: TE-side edge, LE edge, inboard, outboard
    xl, xj = xb - 0.25 * (xb - xa), xb - 0.75 * (xb - xa)
    O, N1, N2, B = (0, 0, 0), (0, y1, 0), (0, y2, 0), (0, yl, 0)
    P1, P2, P3, P4 = (xb, ya, 0), (xa, ya, 0), (xa, yb, 0), (xb, yb, 0)   # build_PaPs corner order
    BVP1, BVP2, Lp = (xl, ya, 0), (xl, yb, 0), (xl, yl, 0)
    JP1, JP2, J = (xj, ya, 0), (xj, yb, 0), (xj, yl, 0)
    xd = -0.45                                     # the y_l dimension line, downstream of the beam
    off = 0.16                                     # the rotation arrows along the beam float this far above it
    L_rot, L_up, L_dn = 0.8, 0.9, 0.9

    # ---- DLM panel, farthest from the camera, drawn first ----
    cv.poly([cv.pt(p) for p in (P1, P2, P3, P4)], INK, LW_THIN, fill=GREY, close=True)
    cv.line3(BVP1, BVP2, INK, LW_THIN, DASH)
    cv.line3(JP1, JP2, MUTED, LW_THIN, DASH)
    cv.line3(B, Lp, MUTED, LW_HAIR, DOT)                                # the arm x_l, in the wing plane
    cv.label('$x_l$', add(mul(B, 0.5), mul(Lp, 0.5)), 3, -4, FS, MUTED, 'start')
    for p in (BVP1, BVP2, Lp, JP1, JP2, J):
        cv.circle3(p, 1.9, WHITE, INK, 0.7)
    cv.label('$BVP_1$', BVP1, -4, 3, FS_SMALL, INK, 'end')
    cv.label('$BVP_2$', BVP2, 4, 3, FS_SMALL, INK, 'start')
    cv.label('$JP_1$', JP1, -4, 3, FS_SMALL, INK, 'end')
    cv.label('$JP_2$', JP2, 4, 3, FS_SMALL, INK, 'start')
    cv.label('$l$', Lp, -5, 3, FS, INK, 'end')
    cv.label('$j$', J, 5, 3, FS, INK, 'start')
    cv.arrow3(Lp, UP, 1.0, BLUE, LW_ARROW)
    cv.label('$C_p$', add(Lp, mul(UP, 1.0)), 5, 3, FS, BLUE, 'start')
    cv.arrow3(J, Z, 0.8, BLUE, LW_ARROW)
    cv.label('$w_j$', add(J, mul(Z, 0.8)), 5, 3, FS, BLUE, 'start')

    # ---- wing plane: dimension y_l, the beam beyond the element ----
    cv.line3((xd, 0, 0), (xd, yl, 0), MUTED, LW_HAIR)
    for yy in (0, yl):
        cv.line3((xd + 0.12, yy, 0), (xd - 0.12, yy, 0), MUTED, LW_HAIR)
    cv.label('$y_l$', (xd, 0.5 * yl, 0), 4, 9, FS, MUTED, 'start')
    cv.line3((0, 0.95, 0), N1, MUTED, LW_HAIR, DASH)
    cv.line3(N2, (0, y2 + 0.8, 0), MUTED, LW_HAIR, DASH)

    # ---- beam element, nodes and DOFs ----
    cv.line3(N1, N2, INK, LW_BEAM)
    for N, k in ((N1, 0), (N2, 3)):
        cv.circle3(N, 2.3, INK)
        cv.arrow3(N, Z, L_dn, INK, LW_THIN)                              # d1, d4: bending displacement u_z
        cv.arrow3(N, X, L_rot, INK, LW_THIN, heads=2)                    # d2, d5: bending slope, rotation about x
        Na = add(N, mul(UP, off))
        cv.arrow3(Na, Y, L_rot, INK, LW_THIN, heads=2)                   # d3, d6: torsion, rotation about y
        cv.label(f'$d_{k + 1}$', add(N, mul(Z, L_dn)), 0, 10, FS, INK, 'middle')
        cv.label(f'$d_{k + 2}$', add(N, mul(X, L_rot)), -4, 2, FS, INK, 'end')
        cv.label(f'$d_{k + 3}$', add(Na, mul(Y, L_rot)), 2, -5, FS, INK, 'start')
    cv.label('1', N1, 7, 12, FS, INK, 'start')
    cv.label('2', N2, 7, 12, FS, INK, 'start')

    # ---- what the coupling exchanges at the beam ----
    cv.arrow3(B, Z, 1.1, BLUE, LW_ARROW)
    cv.label('$F_z$', add(B, mul(Z, 1.1)), 5, 3, FS, BLUE, 'start')
    Ba = add(B, mul(UP, off))
    cv.arrow3(Ba, Y, L_up, BLUE, LW_ARROW, heads=2)
    cv.label('$M_y$', add(Ba, mul(Y, L_up)), 3, -5, FS, BLUE, 'start')

    # ---- frame triad at the origin ----
    for d, L, s, dx, dy, anc in ((X, 0.9, '$x$', -4, 3, 'end'), (Y, 0.8, '$y$', -3, -6, 'end'), (Z, 0.7, '$z$', 0, 10, 'middle')):
        cv.arrow3(O, d, L, INK, LW_THIN)
        cv.label(s, add(O, mul(d, L)), dx, dy, FS, INK, anc)

    return [O, (0.9, 0, 0), (0, 0, 0.7), (0, y2 + 0.8, 0), P1, P2, P3, P4, add(Lp, mul(UP, 1.0)),
            add(N1, mul(Z, L_dn)), add(N2, mul(Z, L_dn)), add(B, mul(Z, 1.1)), (xd, 0, 0),
            add(N1, mul(X, L_rot)), add(N2, mul(X, L_rot))]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'figures'))
    a = ap.parse_args()
    tex = Tex()
    # fit the scene: bounding box of its key points, with room for the labels
    probe = Canvas(tex, 1.0, (0.0, 0.0)); pts = [proj(p) for p in draw(probe)]
    x0, x1 = min(p[0] for p in pts), max(p[0] for p in pts)
    y0, y1 = min(p[1] for p in pts), max(p[1] for p in pts)
    ml, mr, mt, mb = 26, 40, 22, 30
    s = min((W - ml - mr) / (x1 - x0), (H - mt - mb) / (y1 - y0))
    ox = ml + 0.5 * (W - ml - mr - s * (x1 - x0)) - s * x0
    oy = mt + 0.5 * (H - mt - mb - s * (y1 - y0)) + s * y1
    cv = Canvas(tex, s, (ox, oy)); draw(cv)
    hdr = '<!-- DLM-FEM coupling schematic of the FlutterWing tutorial, generated by docs/figures/make_coupling_diagram.py. Text: Latin Modern Roman outlines (GUST Font License). -->'
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W:.2f} {H:.2f}" width="{W:.2f}" height="{H:.2f}">{hdr}'
           f'<rect width="{W:.2f}" height="{H:.2f}" fill="{WHITE}"/>{"".join(cv.parts)}</svg>\n')
    os.makedirs(a.out, exist_ok=True)
    svg_file, png_file = os.path.join(a.out, 'DLM_FEM_Coupling.svg'), os.path.join(a.out, 'DLM_FEM_Coupling.png')
    with open(svg_file, 'w') as f:
        f.write(svg)
    import resvg_py
    wpx, hpx = round(W_CM / 2.54 * DPI), round(H_CM / 2.54 * DPI)
    png = resvg_py.svg_to_bytes(svg_string=svg, width=wpx, height=hpx)
    with open(png_file, 'wb') as f:
        f.write(bytes(png))
    print(f'wrote {svg_file} and {png_file} ({wpx} x {hpx} px, scale {s:.1f} pt/m)')


if __name__ == '__main__':
    main()
