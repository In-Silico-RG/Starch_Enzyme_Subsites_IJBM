"""Assemble the Figure 1 candidates (7.0 in wide, 300 dpi PNG + PDF) from the rendered structure panels
(fig1_protein_struct.py), the frame centroids (fig1_protein_frame.py) and the RDKit chemistry
(fig1_protein_chem.py).

Layout A: proteins first.  a GH13 6JQB | b GH31 3WEO ; c GH13 frame | d GH31 frame | e catalytic centre 6JQB ;
          f chemistry (maltohexaose at -3..+3 with the scissile bond; glucose unit with R at O2/O3/O6 and the
          four substituents).
Layout B: question first.  a chemistry ; b GH13 6JQB | c GH13 frame | d 6JQB catalytic centre ;
          e GH31 3WEO | f GH31 frame | g 3WEO catalytic centre.
Labels are placed by a greedy search over candidate positions around each anchor (page coordinates in A of the
orthoscopic render), penalising overlap with sugar atoms, other labels and the panel edge.
Run: ~/miniconda3/bin/python code/figures/fig1_protein_layouts.py
"""
import csv
import json
import math
import sys
from collections import Counter
from pathlib import Path

import numpy as np
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

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7, 'axes.linewidth': 0.5,
                     'pdf.fonttype': 42, 'ps.fonttype': 42})
FS = 7
LET = 9
W_IN = 7.0
DPI = 300
CHEM_DPI = 600
HALO = [pe.withStroke(linewidth=2.2, foreground='white')]
THREE = {'ASP': 'Asp', 'GLU': 'Glu'}
DOM_NAME = {'6JQB': [('A', 'A'), ('B', 'B'), ('C', 'C')],
            '3WEO': [('N', 'N'), ('A', 'A'), ('C', 'C prox.'), ('C2', 'C dist.')]}
# reference numbering (frame reference) of the catalytic residues shown
REFNUM = {('6JQB', 193): 'D179', ('6JQB', 219): 'E204', ('6JQB', 294): 'D289',
          ('3WEO', 469): 'D443', ('3WEO', 568): 'D542'}


def letter(ax, s, dx=-0.0, dy=1.0):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=LET, fontweight='bold', va='top', ha='left')


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
                    cost -= 3 * float(u @ np.asarray(it['prefer']))
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
    fig.canvas.draw()
    bb = t.get_window_extent()
    t.remove()
    inv = ax.transData.inverted()
    (x0, y0), (x1, y1) = inv.transform([(bb.x0, bb.y0), (bb.x1, bb.y1)])
    return abs(x1 - x0), abs(y1 - y0)


def scalebar(ax, x, y, length, label, color='black'):
    ax.plot([x, x + length], [y, y], color=color, lw=1.0, solid_capstyle='butt')
    ax.text(x + length / 2, y + 0.8, label, ha='center', va='bottom', fontsize=FS)


# ------------------------------------------------------------------ panels
def whole_panel(ax, pdb, legend_loc='lower right', style=''):
    meta = json.load(open(C.WORK / f'whole_{pdb}{style}.json'))
    hw, hh = meta['half']
    img = plt.imread(meta['png'])
    ax.imshow(img, extent=[-hw, hw, -hh, hh], interpolation='lanczos')
    ax.set_xlim(-hw, hw); ax.set_ylim(-hh, hh); ax.set_aspect('equal'); ax.set_axis_off()
    shown = [u for u in meta['units'] if u['sub'] in C.SUBSITES]
    shown.sort(key=lambda u: C.SUBSITES.index(u['sub']))
    obs = np.array([[x, y, 0.9] for u in meta['units'] for x, y, _ in u['atoms']] +
                   [[x, y, 0.9] for c in meta['cats'] for x, y, _ in c['atoms']])
    cen = np.array([u['xyz'][:2] for u in shown])
    d = cen[-1] - cen[0]; nrm = np.array([-d[1], d[0]]) / np.linalg.norm(d)
    items = []
    for k, u in enumerate(shown):
        w, h = text_size(ax, C.sub_label(u['sub']))
        items.append(dict(anchor=tuple(u['xyz'][:2]), text=C.sub_label(u['sub']), w=w, h=h, sub=u['sub'],
                          prefer=tuple(nrm if k % 2 == 0 else -nrm)))
    ch = place(items, obs, (-hw, -hh, hw, hh), 3.0, 9.0)
    for c in ch:
        col = C.SUB_COL[c['it']['sub']]
        ax.plot([c['anchor'][0], c['xy'][0]], [c['anchor'][1], c['xy'][1]], color=col, lw=0.6, zorder=3,
                path_effects=[pe.withStroke(linewidth=1.6, foreground='white')])
        ax.text(*c['xy'], c['it']['text'], color=col, fontsize=FS, fontweight='bold', ha='center', va='center',
                path_effects=HALO, zorder=4)
    handles = [Patch(facecolor=C.DOM_COL[k], edgecolor='0.35', lw=0.4, label=n) for k, n in DOM_NAME[pdb]]
    handles.append(Line2D([], [], color=C.CAT_C, lw=1.6, label='catalytic'))
    leg = ax.legend(handles=handles, loc=legend_loc, fontsize=FS, frameon=False, handlelength=1.0,
                    handleheight=0.8, borderaxespad=0.1, labelspacing=0.25, handletextpad=0.4)
    scalebar(ax, -hw + 2, -hh + 2, 20, '20 Å')
    return meta


def close_panel(ax, pdb):
    meta = json.load(open(C.WORK / f'close_{pdb}.json'))
    hw, hh = meta['half']
    ax.imshow(plt.imread(meta['png']), extent=[-hw, hw, -hh, hh], interpolation='lanczos')
    ax.set_xlim(-hw, hw); ax.set_ylim(-hh, hh); ax.set_aspect('equal'); ax.set_axis_off()
    obs = np.array([[x, y, 0.75] for u in meta['units'] for x, y, _ in u['atoms']] +
                   [[x, y, 0.75] for c in meta['cats'] for x, y, _ in c['atoms']])
    items = []
    for c in meta['cats']:
        s = f"{THREE[c['name']]}{c['seq']}"
        w, h = text_size(ax, s)
        # anchor: carboxylate C (CG of Asp, CD of Glu)
        a = c['X']['CG' if c['name'] == 'ASP' else 'CD']
        items.append(dict(anchor=tuple(a[:2]), text=s, w=w, h=h, col=C.CAT_C, prefer=None))
    for u in meta['units']:
        s = C.sub_label(u['sub'])
        w, h = text_size(ax, s)
        items.append(dict(anchor=tuple(u['xyz'][:2]), text=s, w=w, h=h, col=C.SUB_COL[u['sub']], prefer=None,
                          inside=True))
    ch = place(items, obs, (-hw, -hh, hw, hh), 0.8, 4.0)
    for c in ch:
        it = c['it']
        ax.plot([c['anchor'][0], c['xy'][0]], [c['anchor'][1], c['xy'][1]], color=it['col'], lw=0.5, zorder=3)
        ax.text(*c['xy'], it['text'], color=it['col'], fontsize=FS, fontweight='bold', ha='center', va='center',
                path_effects=HALO, zorder=4)
    scalebar(ax, -hw + 0.5, -hh + 0.4, 5, '5 Å')
    return meta


def frame_panel(ax, fam, legend=True, lim=None):
    rows, X = F.load_page(fam)
    subs = [r['subsite'] for r in rows]
    if lim is None:
        lo = X[:, :2].min(0) - 5; hi = X[:, :2].max(0) + 5
    else:
        lo, hi = np.array(lim[0]), np.array(lim[1])
    zmid = float(np.median(X[:, 2]))
    xs, ys, Fd, cat = F.reference_outline(fam, zmid, ((lo[0], hi[0]), (lo[1], hi[1])))
    ax.contourf(xs, ys, Fd, levels=[0.5, 10], colors=['#ECECEC'], zorder=0)
    ax.contour(xs, ys, Fd, levels=[0.5], colors=['#A8A8A8'], linewidths=0.5, zorder=1)
    cnt = Counter(subs)
    for s in C.SUBSITES:
        m = np.array([x == s for x in subs])
        if m.any():
            ax.scatter(X[m, 0], X[m, 1], s=5, color=C.SUB_COL[s], alpha=0.8, lw=0.15, edgecolor='white',
                       zorder=3, label=f'{C.sub_label(s)}  {cnt[s]}')
    ax.scatter(cat[:, 0], cat[:, 1], marker='+', s=28, color='black', lw=0.8, zorder=4)
    ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_aspect('equal'); ax.set_axis_off()
    if legend:
        h, l = ax.get_legend_handles_labels()
        h.append(Line2D([], [], marker='+', ls='', color='black', ms=5, mew=0.8))
        l.append('cat.')
        ax.legend(h, l, loc='best', fontsize=FS, frameon=False, handletextpad=0.1, borderaxespad=0.1,
                  labelspacing=0.2, markerscale=1.6, handlelength=1.0)
    scalebar(ax, lo[0] + 1, lo[1] + 1, 5, '5 Å')
    return dict(n=len(rows), counts=dict(cnt), entries=len({r['pdb'] for r in rows}))


def chem_chain(ax, tag):
    """maltohexaose -3..+3 on axes ax (axes size sets the image size)"""
    fig = ax.figure
    bb = ax.get_position()
    w_in, h_in = bb.width * fig.get_figwidth(), bb.height * fig.get_figheight()
    w, h = int(w_in * CHEM_DPI), int(h_in * CHEM_DPI)
    fpx = FS / 72 * CHEM_DPI
    png = C.WORK / f'chem_chain_{tag}.png'
    info = K.chain_panel(w, int(h * 0.80), fpx * 0.92, str(png), bond_px=0.6 / 72 * CHEM_DPI)
    ax.set_xlim(0, w); ax.set_ylim(h, 0); ax.set_axis_off()
    hh = int(h * 0.80)
    for s, u in info['units'].items():
        ax.add_patch(Polygon(u['ring'], closed=True, facecolor=C.SUB_COL[s], alpha=0.35, lw=0, zorder=0))
    ax.imshow(plt.imread(png), extent=[0, w, hh, 0], zorder=1, interpolation='lanczos')
    ax.set_ylim(h, 0)
    for s, u in info['units'].items():
        ax.text(u['cen'][0], hh + (h - hh) * 0.45, C.sub_label(s), color=C.SUB_COL[s], fontsize=FS,
                fontweight='bold', ha='center', va='center')
    # scissile bond C1(-1)-O4(+1): short bar across the bond and a scissors mark
    (x1, y1), (x2, y2) = info['cut']
    mx, my = (x1 + x2) / 2, (y1 + y2) / 2
    v = np.array([x2 - x1, y2 - y1]); v /= np.linalg.norm(v); n = np.array([-v[1], v[0]])
    L = 0.09 * h
    ax.plot([mx - n[0] * L, mx + n[0] * L], [my - n[1] * L, my + n[1] * L], color='black', lw=1.0, zorder=2)
    ax.text(mx + n[0] * L * 1.25, my + n[1] * L * 1.25, '✂', fontsize=FS + 2, ha='center', va='center',
            rotation=0, zorder=2)
    # O2, O3, O6 named on unit -2 (k = 1)
    ring_c = np.array(info['units']['-2']['cen'])
    for nm in ('O2', 'O3', 'O6'):
        p = np.array(info['atoms'][(1, nm)])
        u = p - ring_c; u /= np.linalg.norm(u)
        q = p + u * fpx * 1.7
        ax.text(q[0], q[1], nm, fontsize=FS - 0.5, style='italic', color='0.25', ha='center', va='center',
                zorder=3, path_effects=HALO)
    return info


def chem_substituents(fig, rect, tag):
    """glucose unit with R at O2/O3/O6 and the four substituents in a 1 + 2x2 arrangement inside rect"""
    x0, y0, w0, h0 = rect
    fpx = FS / 72 * CHEM_DPI
    bond = 0.6 / 72 * CHEM_DPI
    # glucose unit
    gw = w0 * 0.34
    ax = fig.add_axes([x0, y0 + h0 * 0.08, gw, h0 * 0.84])
    bb = ax.get_position()
    W, H = int(bb.width * fig.get_figwidth() * CHEM_DPI), int(bb.height * fig.get_figheight() * CHEM_DPI)
    png = C.WORK / f'chem_glcR_{tag}.png'
    opos = K.glucose_R_panel(W, H, fpx * 0.92, str(png), bond)
    ax.imshow(plt.imread(png), extent=[0, W, H, 0], interpolation='lanczos')
    ax.set_xlim(0, W); ax.set_ylim(H, 0); ax.set_axis_off()
    for nm in ('O2', 'O3', 'O6'):
        p = np.array(opos[nm]); r = np.array(opos['R' + nm[1:]])
        u = p - r; u = np.array([u[1], -u[0]]) / np.linalg.norm(u)     # perpendicular to the O-R bond
        q = p + u * fpx * 1.1
        ax.text(q[0], q[1], nm[1], fontsize=FS - 1, color='0.25', ha='center', va='center', style='italic')
    axes = [ax]
    # four substituents
    sw, sh = (w0 - gw) / 2, h0 / 2
    for k, (name, smi) in enumerate(K.SUBST):
        i, j = divmod(k, 2)
        a = fig.add_axes([x0 + gw + j * sw, y0 + h0 - (i + 1) * sh + sh * 0.22, sw, sh * 0.78])
        bb = a.get_position()
        W, H = int(bb.width * fig.get_figwidth() * CHEM_DPI), int(bb.height * fig.get_figheight() * CHEM_DPI)
        png = C.WORK / f'chem_sub{k}_{tag}.png'
        K.substituent(smi, W, H, fpx * 0.92, str(png), bond)
        a.imshow(plt.imread(png), extent=[0, W, H, 0], interpolation='lanczos')
        a.set_xlim(0, W); a.set_ylim(H, 0); a.set_axis_off()
        fig.text(x0 + gw + (j + 0.5) * sw, y0 + h0 - (i + 1) * sh + sh * 0.12, name, ha='center', va='center',
                 fontsize=FS)
        axes.append(a)
    return axes


# ------------------------------------------------------------------ layouts
def layout_A():
    H = 8.3
    fig = plt.figure(figsize=(W_IN, H))
    fy = lambda y: y / H
    # row 1: whole proteins (y from top)
    r1 = 2.75
    a = fig.add_axes([0.0, fy(H - 0.05 - r1), 0.5, fy(r1)])
    b = fig.add_axes([0.5, fy(H - 0.05 - r1), 0.5, fy(r1)])
    whole_panel(a, '6JQB'); letter(a, 'a')
    whole_panel(b, '3WEO'); letter(b, 'b')
    r2 = 2.35
    top2 = H - 0.15 - r1
    c = fig.add_axes([0.0, fy(top2 - r2), 0.36, fy(r2)])
    d = fig.add_axes([0.36, fy(top2 - r2), 0.30, fy(r2)])
    e = fig.add_axes([0.67, fy(top2 - r2), 0.33, fy(r2)])
    sc = frame_panel(c, 'GH13'); letter(c, 'c')
    sd = frame_panel(d, 'GH31'); letter(d, 'd')
    close_panel(e, '6JQB'); letter(e, 'e')
    r3 = 1.45
    top3 = top2 - r2 - 0.12
    f = fig.add_axes([0.0, fy(top3 - r3), 1.0, fy(r3)])
    chem_chain(f, 'A'); letter(f, 'f', dy=1.0)
    r4 = 1.45
    top4 = top3 - r3 - 0.05
    chem_substituents(fig, (0.04, fy(top4 - r4), 0.94, fy(r4)), 'A')
    save(fig, 'fig1_protein_A')
    return sc, sd


def layout_B():
    H = 7.9
    fig = plt.figure(figsize=(W_IN, H))
    fy = lambda y: y / H
    r1 = 1.30
    a = fig.add_axes([0.0, fy(H - 0.05 - r1), 0.62, fy(r1)])
    chem_chain(a, 'B'); letter(a, 'a')
    chem_substituents(fig, (0.63, fy(H - 0.05 - r1 - 0.05), 0.37, fy(r1 + 0.05)), 'B')
    rr = 3.1
    top = H - 0.2 - r1
    for k, (pdb, fam) in enumerate((('6JQB', 'GH13'), ('3WEO', 'GH31'))):
        y = top - (k + 1) * rr - k * 0.08
        a1 = fig.add_axes([0.0, fy(y), 0.44, fy(rr)])
        a2 = fig.add_axes([0.45, fy(y + 0.25), 0.28, fy(rr - 0.5)])
        a3 = fig.add_axes([0.74, fy(y + 0.25), 0.26, fy(rr - 0.5)])
        whole_panel(a1, pdb, style='_opaque')
        frame_panel(a2, fam)
        close_panel(a3, pdb)
        for ax, l in zip((a1, a2, a3), 'bcd' if k == 0 else 'efg'):
            letter(ax, l)
    save(fig, 'fig1_protein_B')


def save(fig, name):
    fig.savefig(C.OUT / f'{name}.png', dpi=DPI, facecolor='white')
    fig.savefig(C.OUT / f'{name}.pdf', dpi=DPI, facecolor='white')
    plt.close(fig)
    print('wrote', C.OUT / f'{name}.png')


if __name__ == '__main__':
    which = sys.argv[1:] or ['A', 'B']
    if 'A' in which:
        layout_A()
    if 'B' in which:
        layout_B()
