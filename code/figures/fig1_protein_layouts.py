"""Assemble the Figure 1 candidates (7.0 in wide, 300 dpi PNG + PDF), in the style of the approved Figure 2
(code/figures/fig2_v3_*.py: Liberation Sans 7 pt, residue labels 7.5 pt, subsite numbers 8 pt bold black on the
rings with a white halo, panel letters 9 pt bold, panel titles 7.5 pt; sugar gold #F4B942, catalytic magenta
#B0308F, faint grey protein).

Inputs: ChimeraX renders and page coordinates (work/whole_<pdb>[_opaque].{png,json}, fig1_protein_cx.py), ring
centroids in the family frame (work/centroids_<fam>.tsv, fig1_protein_frame.py), RDKit depictions
(fig1_protein_chem.py, called here at the pixel size of their axes).

Layout A (proteins first): a GH13 6JQB | b GH31 3WEO (cartoon under a transparent surface); c GH13 frame |
  d GH31 frame; e chemistry (maltohexaose over -3..+3 with the scissile bond; glucose unit with R at O2/O3/O6 and
  the four substituents in one row).
Layout B (question first): a chemistry (substituents in a block at right); b GH13 6JQB opaque surface | c GH13
  frame; d GH31 3WEO opaque surface | e GH31 frame.
Run: ~/miniconda3/bin/python code/figures/fig1_protein_layouts.py [A] [B]
"""
import csv
import json
import math
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import gemmi
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import Patch, Polygon
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fig1_protein_common as C
import fig1_protein_frame as F
import fig1_protein_chem as K

plt.rcParams.update({'font.family': C.FONT, 'font.size': C.FS, 'axes.linewidth': 0.6,
                     'pdf.fonttype': 42, 'ps.fonttype': 42})
FS = C.FS
W_IN = 7.0
DPI = 300
CHEM_DPI = 600
HALO = [pe.withStroke(linewidth=2.0, foreground='white')]       # fig2_v3 label halo
TITLES = {'6JQB': 'GH13 maltotetraose-forming amylase, PDB 6JQB (1.10 Å)',
          '3WEO': 'GH31 α-glucosidase, PDB 3WEO (1.45 Å)',
          'GH13': 'GH13: sugar units of 79 entries in the 1AMY frame',
          'GH31': 'GH31: sugar units of 18 entries in the 2QMJ frame'}
DOM_NAME = {'6JQB': [('A', 'A'), ('B', 'B'), ('C', 'C')],
            '3WEO': [('N', 'N'), ('A', 'A'), ('C', 'C'), ('C2', 'C′')]}
DOMAINS = {'6JQB': [('A', 1, 113), ('B', 114, 167), ('A', 168, 361), ('C', 362, 418)],
           '3WEO': [('N', 38, 313), ('A', 314, 655), ('C', 656, 756), ('C2', 757, 909)]}


def letter(fig, x, y, s):
    fig.text(x, y, s, fontsize=C.FS_LET, fontweight='bold', va='top', ha='left')


def title(fig, x, y, s):
    fig.text(x, y, s, fontsize=C.FS_TITLE, va='top', ha='left')


# ------------------------------------------------------------------ label placement
def place(items, obstacles, bounds, rmin, rmax, avoid_boxes=()):
    """items: dict(anchor=(x,y), text, w, h, prefer=(ux,uy) or None) in data units.
    obstacles: (N,2) points with radius r_obs (data units).  Greedy, longest-constrained first."""
    angs = np.radians(np.arange(0, 360, 10))
    radii = np.linspace(rmin, rmax, 6)
    chosen = []
    for it in items:
        best = None
        ax_, ay_ = it['anchor']
        for r in radii:
            for t in angs:
                u = np.array([math.cos(t), math.sin(t)])
                cx = ax_ + u[0] * (r + it['w'] / 2 * abs(u[0]))
                cy = ay_ + u[1] * (r + it['h'] / 2 * abs(u[1]))
                b = (cx - it['w'] / 2, cy - it['h'] / 2, cx + it['w'] / 2, cy + it['h'] / 2)
                cost = r * 0.3
                if b[0] < bounds[0] or b[1] < bounds[1] or b[2] > bounds[2] or b[3] > bounds[3]:
                    cost += 1e6
                pad = 0.4
                dx = np.maximum(np.maximum(b[0] - pad - obstacles[:, 0], obstacles[:, 0] - b[2] - pad), 0)
                dy = np.maximum(np.maximum(b[1] - pad - obstacles[:, 1], obstacles[:, 1] - b[3] - pad), 0)
                cost += 30 * (np.hypot(dx, dy) < obstacles[:, 2]).sum()
                for cb in list(avoid_boxes) + [c['box'] for c in chosen]:
                    w = min(b[2], cb[2]) - max(b[0], cb[0]) + 0.6
                    h = min(b[3], cb[3]) - max(b[1], cb[1]) + 0.6
                    if w > 0 and h > 0:
                        cost += 500 + 50 * w * h
                for c in chosen:          # leader crossing
                    if _cross((ax_, ay_), (cx, cy), c['anchor'], c['xy']):
                        cost += 200
                if it.get('prefer') is not None:
                    cost -= it.get('pw', 3) * float(u @ np.asarray(it['prefer']))
                if best is None or cost < best[0]:
                    best = (cost, (cx, cy), b)
        chosen.append(dict(anchor=(ax_, ay_), xy=best[1], box=best[2], it=it))
    return chosen


def _cross(p1, p2, p3, p4):
    def ccw(a, b, c):
        return (c[1] - a[1]) * (b[0] - a[0]) > (b[1] - a[1]) * (c[0] - a[0])
    return ccw(p1, p3, p4) != ccw(p2, p3, p4) and ccw(p1, p2, p3) != ccw(p1, p2, p4)


def text_size(ax, s, fs=FS, bold=True):
    """width, height of a label in data units of ax (equal aspect assumed)"""
    fig = ax.figure
    t = ax.text(0, 0, s, fontsize=fs, fontweight='bold' if bold else 'normal')
    bb = t.get_window_extent(renderer=fig.canvas.get_renderer())
    t.remove()
    inv = ax.transData.inverted()
    (x0, y0), (x1, y1) = inv.transform([(bb.x0, bb.y0), (bb.x1, bb.y1)])
    return abs(x1 - x0), abs(y1 - y0)


def scalebar(ax, x, y, length, label):
    ax.plot([x, x + length], [y, y], color='black', lw=1.0, solid_capstyle='butt', clip_on=False)
    ax.text(x + length + 0.6, y, label, ha='left', va='center', fontsize=FS, clip_on=False)


# ------------------------------------------------------------------ panels
def domain_anchors(pdb, meta, units_xy):
    """page position for each domain letter: median of the domain's Calpha on the page, pushed away from the
    sugar until 12 A clear of every sugar atom (coordinates of work/whole_<pdb>_rot.pdb)."""
    st = gemmi.read_structure(str(C.WORK / f'whole_{pdb}_rot.pdb'))
    ca = {}
    for ch in st[0]:
        for r in ch:
            a = r.find_atom('CA', '*')
            if a is not None and r.het_flag == 'A':
                ca[r.seqid.num] = (a.pos.x, a.pos.y)
    out = {}
    sc = units_xy.mean(0)
    for dn, a, b in DOMAINS[pdb]:
        out.setdefault(dn, []).extend([ca[i] for i in range(a, b + 1) if i in ca])
    res = {}
    for dn, pts in out.items():
        pts = np.array(pts)
        p = np.median(pts, 0)
        u = p - sc; u = u / (np.linalg.norm(u) + 1e-9)
        for _ in range(60):
            if np.min(np.linalg.norm(units_xy - p, axis=1)) >= 12:
                break
            p = p + u
        res[dn] = p
    return res


def whole_panel(ax, pdb, style=''):
    meta = json.load(open(C.WORK / f'whole_{pdb}{style}.json'))
    hw, hh = meta['half']
    ax.imshow(plt.imread(meta['png']), extent=[-hw, hw, -hh, hh], interpolation='lanczos')
    ax.set_xlim(-hw, hw); ax.set_ylim(-hh, hh); ax.set_aspect('equal'); ax.set_axis_off()
    shown = sorted([u for u in meta['units'] if u['sub'] in C.SUBSITES], key=lambda u: C.SUBSITES.index(u['sub']))
    sug = np.array([[x, y] for u in meta['units'] for x, y, _ in u['atoms']])
    cen = np.array([u['xyz'][:2] for u in shown])
    d = cen[-1] - cen[0]; nrm = np.array([-d[1], d[0]]) / np.linalg.norm(d)
    # subsite numbers: on the ring or just beside it (fig2_v3: "at the ring centre or beside it")
    obs = np.array([[x, y, 0.35] for x, y in sug])
    items = []
    for k, u in enumerate(shown):
        w, h = text_size(ax, C.sub_label(u['sub']), fs=C.FS_SUB)
        items.append(dict(anchor=tuple(u['xyz'][:2]), text=C.sub_label(u['sub']), w=w, h=h,
                          prefer=tuple(nrm if k % 2 == 0 else -nrm), pw=2))
    ch = place(items, obs, (-hw, -hh, hw, hh), 0.0, 3.5)
    for c in ch:
        if np.hypot(c['xy'][0] - c['anchor'][0], c['xy'][1] - c['anchor'][1]) > 2.6:
            ax.plot([c['anchor'][0], c['xy'][0]], [c['anchor'][1], c['xy'][1]], color='black', lw=0.5, zorder=3)
        ax.text(*c['xy'], c['it']['text'], color='black', fontsize=C.FS_SUB, fontweight='bold', ha='center',
                va='center', path_effects=HALO, zorder=4)
    for dn, p in domain_anchors(pdb, meta, sug).items():
        name = dict(DOM_NAME[pdb])[dn]
        ax.text(p[0], p[1], name, fontsize=C.FS_RES, color='0.30', fontweight='bold', ha='center', va='center',
                path_effects=HALO, zorder=4)
    scalebar(ax, -hw + 1, -hh + 1.5, 20, '20 Å')
    return meta


def frame_panel(ax, fam):
    rows, X = F.load_page(fam)
    subs = [r['subsite'] for r in rows]
    lo = X[:, :2].min(0) - 5; hi = X[:, :2].max(0) + 5
    zmid = float(np.median(X[:, 2]))
    xs, ys, Fd, cat = F.reference_outline(fam, zmid, ((lo[0], hi[0]), (lo[1], hi[1])))
    ax.contourf(xs, ys, Fd, levels=[0.5, 10], colors=['#EFEFEF'], zorder=0)
    ax.contour(xs, ys, Fd, levels=[0.5], colors=['#B4B4B4'], linewidths=0.5, zorder=1)
    ax.scatter(X[:, 0], X[:, 1], s=7, color=C.SUGAR_C, alpha=0.9, lw=0.25, edgecolor=C.SUGAR_EDGE, zorder=3)
    ax.scatter(cat[:, 0], cat[:, 1], marker='+', s=36, color=C.CAT_C, lw=1.1, zorder=4)
    ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_aspect('equal'); ax.set_axis_off()
    cnt = Counter(subs)
    obs = np.array([[x, y, 0.45] for x, y in X[:, :2]] + [[x, y, 0.8] for x, y in cat[:, :2]])
    meds = {s: np.median(X[np.array([x == s for x in subs]), :2], 0) for s in C.SUBSITES if s in subs}
    ks = list(meds)
    dvec = meds[ks[-1]] - meds[ks[0]]; nrm = np.array([-dvec[1], dvec[0]]) / np.linalg.norm(dvec)
    items = []
    for s in C.SUBSITES:
        m = np.array([x == s for x in subs])
        if not m.any():
            continue
        side = nrm if ks.index(s) % 2 == 0 else -nrm      # alternate sides of the chain of clusters
        txt = f'{C.sub_label(s)}\n{cnt[s]}'
        w1, h1 = text_size(ax, C.sub_label(s), fs=C.FS_SUB)
        w2, h2 = text_size(ax, str(cnt[s]), fs=FS, bold=False)
        items.append(dict(anchor=tuple(meds[s]), text=s, w=max(w1, w2), h=h1 + h2, n=cnt[s],
                          prefer=tuple(side), pw=8))
    ch = place(items, obs, (lo[0], lo[1], hi[0], hi[1]), 1.5, 5.0)
    for c in ch:
        x, y = c['xy']; it = c['it']
        b = c['box']; ax0, ay0 = c['anchor']        # leader from the cluster median to the label box
        tx, ty = np.clip(ax0, b[0], b[2]), np.clip(ay0, b[1], b[3])
        ax.plot([ax0, tx], [ay0, ty], color='0.25', lw=0.5, zorder=4.5, solid_capstyle='round')
        h1 = it['h'] * 0.55
        ax.text(x, y + it['h'] / 2 - h1 / 2, C.sub_label(it['text']), fontsize=C.FS_SUB, fontweight='bold',
                ha='center', va='center', path_effects=HALO, zorder=5)
        ax.text(x, y - it['h'] / 2 + (it['h'] - h1) / 2, str(it['n']), fontsize=FS, ha='center', va='center',
                path_effects=HALO, zorder=5)
    scalebar(ax, lo[0] + 1, lo[1] + 1.2, 5, '5 Å')
    return dict(n=len(rows), counts=dict(cnt), entries=len({r['pdb'] for r in rows}))


def chem_chain(ax, tag):
    """maltohexaose -3..+3 on axes ax (axes size sets the image size)"""
    fig = ax.figure
    bb = ax.get_position()
    w_in, h_in = bb.width * fig.get_figwidth(), bb.height * fig.get_figheight()
    w, h = int(w_in * CHEM_DPI), int(h_in * CHEM_DPI)
    fpx = FS / 72 * CHEM_DPI
    png = C.WORK / f'chem_chain_{tag}.png'
    info = K.chain_panel(w, int(h * 0.84), fpx, str(png), bond_px=0.6 / 72 * CHEM_DPI)
    ax.set_xlim(0, w); ax.set_ylim(h, 0); ax.set_axis_off()
    hh = int(h * 0.84)
    ax.imshow(plt.imread(png), extent=[0, w, hh, 0], zorder=1, interpolation='antialiased')
    for s, u in info['units'].items():     # ring fill in the sugar gold of the structure panels
        ax.add_patch(Polygon(u['ring'], closed=True, facecolor=C.SUGAR_C, alpha=0.40, lw=0, zorder=2))
    ax.set_ylim(h, 0)
    ybot = max(y for (_, _), (x, y) in info['atoms'].items()) + fpx * 2.4   # below the lowest atom and its note
    for s, u in info['units'].items():
        ax.text(u['cen'][0], ybot, C.sub_label(s), color='black', fontsize=C.FS_SUB,
                fontweight='bold', ha='center', va='center')
    # scissile bond C1(-1)-O4(+1): short bar across the bond and a scissors mark
    (x1, y1), (x2, y2) = info['cut']
    mx, my = (x1 + x2) / 2, (y1 + y2) / 2
    v = np.array([x2 - x1, y2 - y1]); v /= np.linalg.norm(v); n = np.array([-v[1], v[0]])
    L = 0.09 * h
    ax.plot([mx - n[0] * L, mx + n[0] * L], [my - n[1] * L, my + n[1] * L], color='black', lw=1.0, zorder=2)
    ax.text(mx + n[0] * L * 1.25, my + n[1] * L * 1.25, '✂', fontsize=FS + 2, ha='center', va='center',
            rotation=0, zorder=2, fontfamily='DejaVu Sans')       # Liberation Sans has no scissors glyph
    return info


def chem_substituents(fig, rect, tag, grid='2x2'):
    """glucose unit with R at O2/O3/O6 and the four substituents inside rect (figure fractions).
    grid '2x2': glucose left, substituents 2 x 2; grid 'row': glucose and the four substituents in one row."""
    x0, y0, w0, h0 = rect
    fpx = FS / 72 * CHEM_DPI
    bond = 0.6 / 72 * CHEM_DPI
    name_h = (FS + 3) / 72 / fig.get_figheight()          # height of a name line, figure fraction

    def img(axrect, fn, *args):
        a = fig.add_axes(axrect)
        bb = a.get_position()
        W, H = int(bb.width * fig.get_figwidth() * CHEM_DPI), int(bb.height * fig.get_figheight() * CHEM_DPI)
        png = C.WORK / f'{fn.__name__}_{tag}_{len(fig.axes)}.png'
        lead = [x for x in args if x is not None]
        out = fn(*lead, W, H, fpx, str(png), bond)
        a.imshow(plt.imread(png), extent=[0, W, H, 0], interpolation='antialiased')
        a.set_xlim(0, W); a.set_ylim(H, 0); a.set_axis_off()
        return a, out

    if grid == 'row':
        wts = [1.15, 1.0, 0.8, 1.75, 1.0]
        ws = [w0 * v / sum(wts) for v in wts]
        xs = np.concatenate([[x0], x0 + np.cumsum(ws)[:-1]])
        ga, opos = img([xs[0], y0, ws[0], h0], K.glucose_R_panel, None)
        cells = [(xs[k + 1], y0 + name_h, ws[k + 1], h0 - name_h) for k in range(4)]
    else:
        gw = w0 * 0.28
        ga, opos = img([x0, y0 + h0 * 0.08, gw, h0 * 0.84], K.glucose_R_panel, None)
        sw, sh = (w0 - gw) / 3, h0 / 2
        top = [0, 1, 3]                                   # hydroxypropyl, acetyl, carboxymethyl on top
        cells = [None] * 4
        for j, k in enumerate(top):
            cells[k] = (x0 + gw + j * sw, y0 + sh + name_h, sw, sh - name_h)
        cells[2] = (x0 + gw, y0 + name_h, w0 - gw, sh - name_h)   # octenyl succinyl across
    for nm in ('O2', 'O3', 'O6'):
        p = np.array(opos[nm]); r = np.array(opos['R' + nm[1:]])
        u = p - r; u = np.array([u[1], -u[0]]) / np.linalg.norm(u)     # perpendicular to the O-R bond
        q = p + u * fpx * 1.1
        ga.text(q[0], q[1], nm[1], fontsize=FS, color='0.25', ha='center', va='center', style='italic')
    for k, (name, smi) in enumerate(K.SUBST):
        cx, cy, cw, ch = cells[k]
        img([cx, cy, cw, ch], K.substituent, smi)
        fig.text(cx + cw / 2, cy - name_h * 0.45, name, ha='center', va='center', fontsize=FS)



def legend(fig, x, y, ncol):
    h = [Line2D([], [], color=C.SUGAR_C, lw=3.2, label='sugar unit, subsite −3 to +3 (number on the ring)'),
         Line2D([], [], color=C.CTX_C, lw=2.0, label='sugar unit outside −3 to +3'),
         Line2D([], [], color=C.CAT_C, lw=2.2, label='catalytic residues'),
         Line2D([], [], color=C.CAT_C, marker='+', ls='', ms=6, mew=1.1, label='catalytic carboxylates, frame reference'),
         Patch(facecolor='#EFEFEF', edgecolor='#B4B4B4', lw=0.5, label='frame reference, ±5 Å section')]
    fig.legend(handles=h, loc='center', bbox_to_anchor=(x, y), ncol=ncol, frameon=False, fontsize=FS,
               handlelength=1.8, columnspacing=1.1, handletextpad=0.5)


# ------------------------------------------------------------------ layouts
def layout_A():
    H = 8.25
    fig = plt.figure(figsize=(W_IN, H))
    fy = lambda y: y / H
    t = H - 0.04
    letter(fig, 0.0, fy(t), 'a'); title(fig, 0.035, fy(t - 0.02), TITLES['6JQB'])
    letter(fig, 0.51, fy(t), 'b'); title(fig, 0.545, fy(t - 0.02), TITLES['3WEO'])
    r1 = 2.55
    a = fig.add_axes([0.0, fy(t - 0.2 - r1), 0.49, fy(r1)])
    b = fig.add_axes([0.51, fy(t - 0.2 - r1), 0.49, fy(r1)])
    whole_panel(a, '6JQB'); whole_panel(b, '3WEO')
    t2 = t - 0.2 - r1 - 0.06
    letter(fig, 0.0, fy(t2), 'c'); title(fig, 0.035, fy(t2 - 0.02), TITLES['GH13'])
    letter(fig, 0.51, fy(t2), 'd'); title(fig, 0.545, fy(t2 - 0.02), TITLES['GH31'])
    r2 = 2.05
    c = fig.add_axes([0.0, fy(t2 - 0.2 - r2), 0.49, fy(r2)])
    d = fig.add_axes([0.51, fy(t2 - 0.2 - r2), 0.49, fy(r2)])
    sc = frame_panel(c, 'GH13'); sd = frame_panel(d, 'GH31')
    legend(fig, 0.5, fy(t2 - 0.2 - r2 - 0.2), 3)
    t3 = t2 - 0.2 - r2 - 0.48
    letter(fig, 0.0, fy(t3), 'e')
    r3 = 1.55
    e = fig.add_axes([0.03, fy(t3 - r3), 0.94, fy(r3)])
    chem_chain(e, 'A')
    r4 = 1.2
    t4 = t3 - r3 - 0.05
    chem_substituents(fig, (0.02, fy(t4 - r4), 0.96, fy(r4)), 'A', grid='row')
    save(fig, 'fig1_protein_A')
    return sc, sd


def layout_B():
    H = 8.6
    fig = plt.figure(figsize=(W_IN, H))
    fy = lambda y: y / H
    t = H - 0.04
    r1 = 1.85
    letter(fig, 0.0, fy(t), 'a')
    a = fig.add_axes([0.0, fy(t - r1 + 0.1), 0.52, fy(r1 - 0.15)])
    chem_chain(a, 'B')
    chem_substituents(fig, (0.52, fy(t - r1), 0.455, fy(r1)), 'B', grid='2x2')
    rr = 2.75
    top = t - r1 - 0.12
    for k, (pdb, fam) in enumerate((('6JQB', 'GH13'), ('3WEO', 'GH31'))):
        y = top - k * (rr + 0.3)
        l1, l2 = ('b', 'c') if k == 0 else ('d', 'e')
        letter(fig, 0.0, fy(y), l1); title(fig, 0.035, fy(y - 0.02), TITLES[pdb])
        letter(fig, 0.55, fy(y), l2); title(fig, 0.585, fy(y - 0.02), TITLES[fam])
        a1 = fig.add_axes([0.0, fy(y - 0.2 - rr), 0.53, fy(rr)])
        a2 = fig.add_axes([0.55, fy(y - 0.2 - rr + 0.25), 0.45, fy(rr - 0.5)])
        whole_panel(a1, pdb, style='_opaque')
        frame_panel(a2, fam)
    legend(fig, 0.5, fy(top - 2 * rr - 0.3 - 0.35), 3)
    save(fig, 'fig1_protein_B')


def save(fig, name):
    # tight bounding box trims empty margins only; the panels span the 7.0 in width
    fig.savefig(C.OUT / f'{name}.png', dpi=DPI, facecolor='white', bbox_inches='tight', pad_inches=0.03)
    fig.savefig(C.OUT / f'{name}.pdf', dpi=DPI, facecolor='white', bbox_inches='tight', pad_inches=0.03)
    plt.close(fig)
    print('wrote', C.OUT / f'{name}.png')


if __name__ == '__main__':
    which = sys.argv[1:] or ['A', 'B']
    if 'A' in which:
        layout_A()
    if 'B' in which:
        layout_B()
