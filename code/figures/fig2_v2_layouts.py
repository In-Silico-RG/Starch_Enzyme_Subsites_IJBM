"""Figure 2 (v2) candidates A, B, C for the short draft (two-column width 7.0 in).

Inputs
  data/figures/fig2_v2_work/per_protein_counts.tsv, contact_rows.tsv  (code/figures/fig2_v2_data.py)
  data/census/cif/{3BC9,1G94,3WEO}.cif                                (deposited coordinates)
Structure panels: code/figures/fig2_v2_struct.py (PyMOL ray-trace + label placement here).
Data panel: residue types per position, one vote per protein (UniProt accession), engineered mutants excluded.

Run:  ~/miniconda3/bin/python code/figures/fig2_v2_data.py
      ~/miniconda3/bin/python code/figures/fig2_v2_layouts.py [A B C]
Writes figures/candidates/fig2_v2_{A,B,C}.{png,pdf} (PNG at 300 dpi).
"""
import csv
import sys
import math
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, FancyBboxPatch

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fig2_v2_struct as S

P = S.P
OUT = P / 'figures/candidates'
OUT.mkdir(parents=True, exist_ok=True)
WORK = S.WORK
COUNTS = WORK / 'per_protein_counts.tsv'

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7, 'pdf.fonttype': 42, 'ps.fonttype': 42,
                     'axes.linewidth': 0.6})
THREE = {'ALA': 'Ala', 'ARG': 'Arg', 'ASN': 'Asn', 'ASP': 'Asp', 'CYS': 'Cys', 'GLN': 'Gln', 'GLU': 'Glu',
         'GLY': 'Gly', 'HIS': 'His', 'ILE': 'Ile', 'LEU': 'Leu', 'LYS': 'Lys', 'MET': 'Met', 'PHE': 'Phe',
         'PRO': 'Pro', 'SER': 'Ser', 'THR': 'Thr', 'TRP': 'Trp', 'TYR': 'Tyr', 'VAL': 'Val'}
ONE = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
       'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
       'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}
ORDER = {
    'GH13': [('H', ['D179', 'D289', 'H288', 'R177', 'H92', 'E204', 'K182']), ('Ar', ['Y51']),
             ('Al', ['F180', 'W206'])],
    'GH31': [('H', ['D327', 'H600', 'D443', 'D542']), ('Ar', ['W441', 'W539', 'F575']),
             ('Al', ['I328', 'I364'])],
}
CLASSNAME = {'H': 'H', 'Ar': 'Ar', 'Al': 'Al'}
MINUS = '−'


def sub_txt(s):
    return s.replace('-', MINUS)


# ------------------------------------------------------------------ data panel
def read_counts():
    d = {}
    for r in csv.DictReader(open(COUNTS), delimiter='\t'):
        votes = []
        for part in r['votes'].split(', '):
            t, v = part.rsplit(' ', 1)
            votes.append((t, float(v)))
        d[(r['family'], r['position'])] = dict(sub=r['subsite'], cls=r['class'], n=int(r['n_proteins']),
                                               votes=votes, nent=int(r['n_entries_native']))
    return d


def tint(c, f):
    r = np.array(to_rgb(c))
    return tuple(r + (1 - r) * f)


def family_groups(fam):
    return [(cls, [(fam, p) for p in poss]) for cls, poss in ORDER[fam]]


def draw_bars(ax, groups_in, data, fs=7, show_title=True, title=None, bar_h=0.70, gap=0.5, axis_label=True,
              fam_prefix=False, pitch=0.20, top_frac=None):
    """Stacked horizontal bars of residue-type fractions per position (one vote per protein).
    groups_in: [(class, [(family, position), ...]), ...]. Row label: residue of the reference structure and
    subsite. Bar segments: residue types, majority first, tints of the class colour. Right: n = number of
    proteins, then the types too narrow to label, with their protein counts.
    The axes height is set so that one row is `pitch` inches, top edge kept at top_frac (figure fraction)."""
    rows, groups, y = [], [], 0.0
    for cls, poss in groups_in:
        y0 = y
        for fp in poss:
            rows.append((y, fp, cls)); y += 1
        groups.append((cls, y0, y - 1))
        y += gap
    ymax = y - gap - 1
    top = -0.5 - (1.0 if show_title else 0)
    bottom = ymax + (1.7 if axis_label else 0.6)
    fig = ax.figure
    pos = ax.get_position()
    Hfig = fig.get_size_inches()[1]
    h = (bottom - top) * pitch / Hfig
    t = pos.y1 if top_frac is None else top_frac
    ax.set_position((pos.x0, t - h, pos.width, h))
    ax.set_xlim(0, 1)
    ax.set_ylim(bottom, top)
    ax.axis('off')
    fig = ax.figure
    rnd = fig.canvas.get_renderer()
    bb = ax.get_window_extent(rnd)
    x_per_px = 1.0 / bb.width

    def tw(txt, size=fs, weight='normal'):
        t = ax.text(0, 0, txt, fontsize=size, fontweight=weight)
        w = t.get_window_extent(rnd).width; t.remove(); return w * x_per_px

    tints = [0.0, 0.45, 0.65, 0.80, 0.90]
    lab_w = 0
    for y, (fam, pos), cls in rows:
        d = data[(fam, pos)]
        base = S.COL[cls]
        x = 0.0
        notes = []
        for k, (t, v) in enumerate(d['votes']):
            w = v / d['n']
            col = tint(base, tints[min(k, len(tints) - 1)])
            ax.add_patch(Rectangle((x, y - bar_h / 2), w, bar_h, fc=col, ec='white', lw=0.7))
            tc = 'white' if (k == 0 and cls == 'H') else 'black'
            for txt in (THREE[t], ONE[t]):
                if tw(txt) + 0.012 < w:
                    ax.text(x + w / 2, y, txt, ha='center', va='center', fontsize=fs, color=tc)
                    break
            else:
                notes.append(f'{THREE[t]} {v:g}')
            x += w
        ax.text(1.02, y, f'{d["n"]}', ha='left', va='center', fontsize=fs, color='0.1')
        if notes:
            ax.text(1.02 + tw('00') + 0.03, y, ', '.join(notes), ha='left', va='center', fontsize=fs,
                    color='0.25')
        three = THREE[{'D': 'ASP', 'H': 'HIS', 'R': 'ARG', 'E': 'GLU', 'K': 'LYS', 'Y': 'TYR', 'F': 'PHE',
                       'W': 'TRP', 'I': 'ILE'}[pos[0]]]
        lab = f'{three}{pos[1:]} ({sub_txt(d["sub"])})'
        if fam_prefix:
            lab = f'{fam}  ' + lab
        ax.text(-0.025, y, lab, ha='right', va='center', fontsize=fs)
        lab_w = max(lab_w, tw(lab))
    xg = -0.025 - lab_w - 0.03
    for cls, y0, y1 in groups:
        ax.add_patch(Rectangle((xg - 0.012, y0 - bar_h / 2), 0.012, (y1 - y0) + bar_h, fc=S.COL[cls], ec='none',
                               clip_on=False))
        ax.text(xg - 0.03, (y0 + y1) / 2, cls, ha='right', va='center', fontsize=fs + 0.5,
                fontweight='bold', color=S.COL[cls])
    if show_title:
        ax.text(0.0, top + 0.35, title, ha='left', va='center', fontsize=fs + 0.5, fontweight='bold')
        ax.text(1.02, top + 0.35, 'n', ha='left', va='center', fontsize=fs, style='italic', color='0.1')
    yb = ymax + 0.5
    ax.plot([0, 1], [yb, yb], color='0.35', lw=0.6, clip_on=False)
    for xv, t in ((0, '0'), (0.5, '0.5'), (1, '1')):
        ax.plot([xv, xv], [yb, yb + 0.15], color='0.35', lw=0.6, clip_on=False)
        ax.text(xv, yb + 0.25, t, ha='center', va='top', fontsize=fs, clip_on=False)
    if axis_label:
        ax.text(0.5, yb + 1.0, 'fraction of proteins', ha='center', va='top', fontsize=fs, clip_on=False)
    return xg - 0.03 - tw('Ar', fs + 0.5, 'bold')


# ------------------------------------------------------------------ structure panel
def bonds_mid(Xs, names):
    pts = []
    names = list(Xs)
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a, b = Xs[names[i]], Xs[names[j]]
            if np.linalg.norm(a - b) < 1.9:
                pts += [a * 0.75 + b * 0.25, (a + b) / 2, a * 0.25 + b * 0.75]
    return pts


def structure_panel(ax, info, png, half, px, fs=7.5, fs_o=7, margin=3.0, olabels=True, crop=None):
    img = plt.imread(png)
    ax.imshow(img, extent=(0, px, px, 0), interpolation='lanczos')
    sc = px / (2 * half)
    tp = lambda x: S.to_px(S.rot(info, x)[:2], half, px)
    pts = S.shown_points(info)
    P2 = np.array([S.to_px(p[:2], half, px) for p in pts])
    if crop is None:
        x0, y0 = P2.min(0) - margin * sc; x1, y1 = P2.max(0) + margin * sc
    else:
        x0, y0, x1, y1 = crop
    ax.set_xlim(x0, x1); ax.set_ylim(y1, y0); ax.axis('off')
    return (x0, y0, x1, y1), sc, tp


def add_labels(ax, info, half, px, bounds, fs=7.5, fs_o=7, ring_labels=True, olabels=True, stack_note=True):
    fig = ax.figure
    rnd = fig.canvas.get_renderer()
    sc = px / (2 * half)
    tp = lambda x: S.to_px(S.rot(info, x)[:2], half, px)
    bb = ax.get_window_extent(rnd)
    d_per_disp = (bounds[2] - bounds[0]) / bb.width

    def measure(txt, size, weight='normal'):
        t = ax.text(0, 0, txt, fontsize=size, fontweight=weight)
        e = t.get_window_extent(rnd); t.remove()
        return e.width * d_per_disp + 6 * d_per_disp, e.height * d_per_disp + 2 * d_per_disp

    # obstacles: atoms and bond points of everything drawn
    obs = []
    for k, rd in enumerate(info['residues']):
        for x in list(rd['X'].values()) + bonds_mid(rd['X'], None):
            p = tp(x); obs.append([p[0], p[1], 0.55 * sc, 100 + k])
    for k, u in enumerate(info['units']):
        for x in list(u['X'].values()) + bonds_mid(u['X'], None):
            p = tp(x); obs.append([p[0], p[1], 0.6 * sc, k])
    obs = np.array(obs)
    segs = []
    for hb in info['hbonds']:
        a = tp(hb['res']['X'][hb['a']]); b = tp(hb['unit']['X'][hb['b']])
        segs.append((a, b))
    sugar_c = np.mean([tp(x) for u in info['units'] for x in u['X'].values()], 0)
    items = []
    # residue labels
    for rd in info['residues']:
        sc_atoms = {n: x for n, x in rd['X'].items() if n != 'CA'} or rd['X']
        far = max(sc_atoms.values(), key=lambda x: np.linalg.norm(tp(x) - sugar_c))
        anc = tp(far)
        w, h = measure(rd['label'], fs, 'bold')
        v = anc - sugar_c; v = v / (np.linalg.norm(v) + 1e-9)
        # candidate attachment atoms: the outer half of the side chain (by distance from the sugar on the page)
        dist = {n: np.linalg.norm(tp(x) - sugar_c) for n, x in sc_atoms.items()}
        med = np.median(list(dist.values()))
        anchors = [tp(x) for n, x in sc_atoms.items() if dist[n] >= med] or [anc]
        items.append(dict(kind='res', text=rd['label'], color=S.COL[rd['cls']], anchor=anc, anchors=anchors,
                          w=w, h=h, radii=[0.4, 0.9, 1.5, 2.2, 3.0, 3.8], prefer=v, priority=2, size=fs,
                          weight='bold',
                          group=100 + info['residues'].index(rd)))
    # ring labels
    if ring_labels:
        for k, u in enumerate(info['units']):
            c = np.mean([tp(u['X'][n]) for n in u['ring']], 0)
            w, h = measure(sub_txt(u['sub']), fs, 'bold')
            items.append(dict(kind='ring', text=sub_txt(u['sub']), color='black', anchor=c, w=w, h=h,
                              radii=[1.1, 1.5, 2.0, 2.6], prefer=None, priority=3, size=fs, weight='bold',
                              center_ok=False, own=k))
    # sugar oxygen labels (hydrogen-bonded ones)
    if olabels:
        seen = set()
        for hb in info['hbonds']:
            u, b = hb['unit'], hb['b']
            if (u['sub'], b) in seen:
                continue
            seen.add((u['sub'], b))
            name = u['names'].get(b, b)
            if name.startswith('N'):
                name = 'N'
            anc = tp(u['X'][b])
            w, h = measure(name, fs_o)
            items.append(dict(kind='O', text=name, color=S.OLAB, anchor=anc, w=w, h=h,
                              radii=[0.35, 0.7, 1.1, 1.6], prefer=None, priority=1, size=fs_o, weight='normal',
                              weightc=1.0))
    boxes, total = S.place_labels(items, obs, bounds, sc, segs=segs)
    for it, b in zip(items, boxes):
        cx, cy = (b[0] + b[2]) / 2, (b[1] + b[3]) / 2
        st = S._attach(b, np.array(it['anchor']))
        lead = np.linalg.norm(st - it['anchor'])
        if it['kind'] == 'res' and lead > 0.35 * sc:
            ax.plot([st[0], it['anchor'][0]], [st[1], it['anchor'][1]], color=it['color'], lw=0.6, zorder=4,
                    solid_capstyle='round')
        if it['kind'] == 'O' and lead > 0.45 * sc:
            ax.plot([st[0], it['anchor'][0]], [st[1], it['anchor'][1]], color=S.OLAB, lw=0.45, zorder=4)
        ax.text(cx, cy, it['text'], ha='center', va='center', fontsize=it['size'], fontweight=it['weight'],
                color=it['color'], zorder=5,
                bbox=dict(boxstyle='round,pad=0.08', fc='white', ec='none', alpha=0.8))
    return items, boxes, total


def check_overlaps(items, boxes):
    bad = []
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            if S._overlap(boxes[i], boxes[j]) > 0:
                bad.append((items[i]['text'], items[j]['text']))
    return bad


def panel_title(ax, letter, text, fs=8):
    ax.text(0.0, 1.0, letter, transform=ax.transAxes, fontsize=fs + 1.5, fontweight='bold', va='top', ha='left')
    ax.text(0.07, 1.0, text, transform=ax.transAxes, fontsize=fs, va='top', ha='left')


def legend_handles():
    return [Line2D([], [], color=S.COL['H'], lw=3, label='H'),
            Line2D([], [], color=S.COL['Ar'], lw=3, label='Ar'),
            Line2D([], [], color=S.COL['Al'], lw=3, label='Al'),
            Line2D([], [], color='0.15', lw=1, ls=(0, (1.5, 1.5)), label='N/O–N/O ≤ 3.5 Å')]


# ------------------------------------------------------------------ stacking inset
INSET_TILT = 62


def stack_view(info):
    """View direction that shows the stacked aromatic ring and the sugar ring from the side."""
    srd, c_ar, c_su = info['stack']
    n_ar = S.plane(np.array([srd['X'][n] for n in S.AROM[srd['resname']]]))[1]
    sub = info['spec']['stack'][1]
    n_su = info['rings'][sub][1]
    if n_ar @ n_su < 0:
        n_ar = -n_ar
    v = c_ar - c_su
    d = np.linalg.norm(v)
    best, bs = None, -1
    for ez in S.fib_dirs(6000):
        gap = np.linalg.norm(v - (v @ ez) * ez) / d          # fraction of the centroid distance seen on page
        if gap < 0.75:
            continue
        sc = min(abs(n_ar @ ez), abs(n_su @ ez))             # both rings seen as rings, not as lines
        if sc > bs:
            bs, best = sc, ez
    return best


def inset_render(pdb, tag):
    base = S.prepare(pdb)
    ez = stack_view(base)
    spos, ssub = base['spec']['stack']
    info = S.prepare(pdb, positions=[spos], view_override=ez)
    info['units'] = [u for u in info['units'] if u['sub'] == ssub]
    png = WORK / f'{pdb}_stack_{tag}.png'
    half = S.render(info, png, px=1400, slab=8, cartoon=False, draw_hbonds=False, sticks=0.2)
    return info, png, half


def draw_inset(ax, info, png, half, px=1400, fs=7):
    img = plt.imread(png)
    ax.imshow(img, extent=(0, px, px, 0), interpolation='lanczos')
    tp = lambda x: S.to_px(S.rot(info, x)[:2], half, px)
    srd = info['residues'][0]
    u = info['units'][0]
    pts = np.array([tp(x) for x in list(srd['X'].values()) + list(u['X'].values())])
    sc = px / (2 * half)
    x0, y0 = pts.min(0) - 1.0 * sc; x1, y1 = pts.max(0) + 1.0 * sc
    c_ar = tp(np.mean([srd['X'][n] for n in S.AROM[srd['resname']]], 0))
    c_su = tp(np.mean([u['X'][n] for n in u['ring']], 0))
    ax.plot([c_ar[0], c_su[0]], [c_ar[1], c_su[1]], color='0.1', lw=0.8, ls=(0, (1.5, 1.2)), zorder=4)
    for c in (c_ar, c_su):
        ax.plot(*c, 'o', ms=2.2, color='0.1', zorder=5)
    d, ang = S.stack_geometry(info)
    y0 -= 2.6 * sc   # room for the two-line inset title
    ax.set_xlim(x0, x1); ax.set_ylim(y1, y0)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_color('0.55'); s.set_linewidth(0.6)
    return d, ang


def fit_rect(cell, aspect, fig):
    """Largest rect with given data aspect (h/w in display) inside cell (figure fraction), centred."""
    W, H = fig.get_size_inches()
    cx, cy, cw, ch = cell
    wi, hi = cw * W, ch * H
    if hi / wi > aspect:
        h = wi * aspect; return (cx, cy + ch - h / H, cw, h / H)
    w = hi / aspect
    return (cx + (cw - w / W) / 2, cy, w / W, ch)


# ------------------------------------------------------------------ composition helpers
PX = 2400


def build_structure(pdb, tag):
    info = S.prepare(pdb)
    png = WORK / f'{pdb}_{tag}.png'
    half = S.render(info, png, px=PX)
    return info, png, half


def place_structure(fig, cell, info, png, half, letter, title, inset=None, fs=7.5, fs_o=7, margin=2.6,
                    inset_frac=0.34, inset_corner=None):
    sc = PX / (2 * half)
    pts = S.shown_points(info)
    P2 = np.array([S.to_px(p[:2], half, PX) for p in pts])
    x0, y0 = P2.min(0) - margin * sc; x1, y1 = P2.max(0) + margin * sc
    title_h = 0.22  # inches reserved for the title line
    W, H = fig.get_size_inches()
    cell_img = (cell[0], cell[1], cell[2], cell[3] - title_h / H)
    rect = fit_rect(cell_img, (y1 - y0) / (x1 - x0), fig)
    ax = fig.add_axes(rect)
    ax.imshow(plt.imread(png), extent=(0, PX, PX, 0), interpolation='lanczos')
    ax.set_xlim(x0, x1); ax.set_ylim(y1, y0); ax.axis('off')
    bounds = (x0, y0, x1, y1)
    extra_obs = []
    if inset is not None:
        iinfo, ipng, ihalf = inset
        # corner with fewest drawn atoms
        wpx, hpx = (x1 - x0) * inset_frac, (y1 - y0) * inset_frac
        corners = {'tl': (x0, y0), 'tr': (x1 - wpx, y0), 'bl': (x0, y1 - hpx), 'br': (x1 - wpx, y1 - hpx)}
        def n_in(c):
            return int(((P2[:, 0] > c[0] - 0.5 * sc) & (P2[:, 0] < c[0] + wpx + 0.5 * sc) &
                        (P2[:, 1] > c[1] - 0.5 * sc) & (P2[:, 1] < c[1] + hpx + 0.5 * sc)).sum())
        key = inset_corner or min(corners, key=lambda k: n_in(corners[k]))
        cxy = corners[key]
        # figure-fraction rect of the inset
        def to_fig(x, y):
            d = ax.transData.transform((x, y)); return fig.transFigure.inverted().transform(d)
        a = to_fig(cxy[0], cxy[1] + hpx); b = to_fig(cxy[0] + wpx, cxy[1])
        iax = fig.add_axes((a[0], a[1], b[0] - a[0], b[1] - a[1]))
        d, ang = draw_inset(iax, iinfo, ipng, ihalf)
        # match inset aspect by expanding limits
        xl, yl = iax.get_xlim(), iax.get_ylim()
        wd, hd = xl[1] - xl[0], yl[0] - yl[1]
        target = (b[1] - a[1]) * H / ((b[0] - a[0]) * W)
        if hd / wd < target:
            e = (wd * target - hd) / 2; iax.set_ylim(yl[0] + e, yl[1] - e)
        else:
            e = (hd / target - wd) / 2; iax.set_xlim(xl[0] - e, xl[1] + e)
        srd = iinfo['residues'][0]
        iax.text(0.04, 0.97, f"{srd['label']} / {sub_txt(iinfo['units'][0]['sub'])} ring\n{d:.1f} \u00c5, {ang:.0f}\u00b0", transform=iax.transAxes, linespacing=1.1,
                 fontsize=7, va='top', ha='left', color=S.COL[srd['cls']], fontweight='bold')
        extra_obs.append((cxy[0] - 0.3 * sc, cxy[1] - 0.3 * sc, cxy[0] + wpx + 0.3 * sc, cxy[1] + hpx + 0.3 * sc))
        inset_info = (d, ang, key)
    else:
        inset_info = None
    items, boxes, total = add_labels_ex(ax, info, half, PX, bounds, fs=fs, fs_o=fs_o, extra_obs=extra_obs)
    tx = fig.transFigure.inverted().transform(ax.transAxes.transform((0, 1)))
    fig.text(rect[0], rect[1] + rect[3] + 0.05 / H, letter, fontsize=10, fontweight='bold', va='bottom', ha='left')
    fig.text(rect[0] + 0.17 / W, rect[1] + rect[3] + 0.05 / H, title, fontsize=7.5, va='bottom', ha='left')
    bad = check_overlaps(items, boxes)
    print(f'  {info["pdb"]}: label cost {total:.0f}, overlapping label pairs: {bad}')
    return ax, inset_info


def add_labels_ex(ax, info, half, px, bounds, fs, fs_o, extra_obs):
    orig = S.place_labels

    def wrapped(items, obstacles, bnds, scale, segs=(), rounds=4):
        return orig(items, obstacles, bnds, scale, segs=segs, rounds=rounds, forbid=extra_obs)
    S.place_labels = wrapped
    try:
        return add_labels(ax, info, half, px, bounds, fs=fs, fs_o=fs_o)
    finally:
        S.place_labels = orig


def caption_line(fig, y, fs=7):
    fig.legend(handles=legend_handles(), loc='lower center', bbox_to_anchor=(0.5, y), ncol=4, fontsize=fs,
               frameon=False, handlelength=2.0, columnspacing=1.6)


# ------------------------------------------------------------------ candidates
def _title(pdb):
    e = S.ENTRIES[pdb]
    return f'{e["family"]}, PDB {pdb} ({e["resolution"]:.2f} \u00c5)'


def candidate_A(gh13='3BC9', tag='A'):
    """Two structure panels side by side (a GH13, b GH31, stacking insets); data panel c below, one block per
    family."""
    data = read_counts()
    a = build_structure(gh13, tag); b = build_structure('3WEO', tag)
    ia = inset_render(gh13, tag); ib = inset_render('3WEO', tag)
    W, H = 7.0, 7.3
    fig = plt.figure(figsize=(W, H))
    top = 3.75 / H
    axa, _ = place_structure(fig, (0.005, 1 - top, 0.49, top - 0.005), *a, 'a', _title(gh13), inset=ia)
    axb, _ = place_structure(fig, (0.505, 1 - top, 0.49, top - 0.005), *b, 'b', _title('3WEO'), inset=ib)
    ylo = min(axa.get_position().y0, axb.get_position().y0)
    caption_line(fig, ylo - 0.045)
    yt = ylo - 0.065
    ax13 = fig.add_axes((0.17, 0.0, 0.25, 0.1))
    ax31 = fig.add_axes((0.71, 0.0, 0.235, 0.1))
    draw_bars(ax13, family_groups('GH13'), data, title='GH13', top_frac=yt)
    draw_bars(ax31, family_groups('GH31'), data, title='GH31', top_frac=yt)
    fig.text(0.005, yt, 'c', fontsize=10, fontweight='bold', va='top')
    save(fig, tag)


def candidate_B():
    """Structures stacked in the left column (larger), data panel in the right column."""
    data = read_counts()
    a = build_structure('3BC9', 'B'); b = build_structure('3WEO', 'B')
    ia = inset_render('3BC9', 'B'); ib = inset_render('3WEO', 'B')
    W, H = 7.0, 7.0
    fig = plt.figure(figsize=(W, H))
    lw = 3.85 / W
    place_structure(fig, (0.003, 0.5, lw, 0.497), *a, 'a', _title('3BC9'), inset=ia, inset_frac=0.33)
    place_structure(fig, (0.003, 0.008, lw, 0.487), *b, 'b', _title('3WEO'), inset=ib, inset_frac=0.33)
    x0 = lw + 0.165
    ax13 = fig.add_axes((x0, 0.0, 0.165, 0.1))
    ax31 = fig.add_axes((x0, 0.0, 0.165, 0.1))
    yt = 0.985
    draw_bars(ax13, family_groups('GH13'), data, title='GH13', top_frac=yt - 0.01, axis_label=False)
    y2 = ax13.get_position().y0 - 0.025
    draw_bars(ax31, family_groups('GH31'), data, title='GH31', top_frac=y2)
    fig.text(lw + 0.012, yt, 'c', fontsize=10, fontweight='bold', va='top')
    fig.legend(handles=legend_handles(), loc='upper left', bbox_to_anchor=(lw + 0.03, ax31.get_position().y0 - 0.02), ncol=2,
               fontsize=7, frameon=False, handlelength=2.0, columnspacing=1.2)
    save(fig, 'B')


def candidate_C():
    """GH13 entry 1G94 (every H position hydrogen-bonds the sugar in this entry); data panel arranged by
    contact class across both families (H left; Ar and Al right), which is the form of the claim."""
    data = read_counts()
    a = build_structure('1G94', 'C'); b = build_structure('3WEO', 'C')
    ia = inset_render('1G94', 'C'); ib = inset_render('3WEO', 'C')
    W, H = 7.0, 7.4
    fig = plt.figure(figsize=(W, H))
    top = 3.75 / H
    axa, _ = place_structure(fig, (0.005, 1 - top, 0.49, top - 0.005), *a, 'a', _title('1G94'), inset=ia)
    axb, _ = place_structure(fig, (0.505, 1 - top, 0.49, top - 0.005), *b, 'b', _title('3WEO'), inset=ib)
    ylo = min(axa.get_position().y0, axb.get_position().y0)
    caption_line(fig, ylo - 0.045)
    yt = ylo - 0.065
    hgroups = [('H', [('GH13', p) for p in ORDER['GH13'][0][1]] + [('GH31', p) for p in ORDER['GH31'][0][1]])]
    rgroups = [('Ar', [('GH13', 'Y51')] + [('GH31', p) for p in ORDER['GH31'][1][1]]),
               ('Al', [('GH13', p) for p in ORDER['GH13'][2][1]] + [('GH31', p) for p in ORDER['GH31'][2][1]])]
    axh = fig.add_axes((0.215, 0.0, 0.23, 0.1))
    axr = fig.add_axes((0.70, 0.0, 0.19, 0.1))
    draw_bars(axh, hgroups, data, title='side-chain hydrogen bond', fam_prefix=True, top_frac=yt, pitch=0.19)
    draw_bars(axr, rgroups, data, title='apolar contact', fam_prefix=True, top_frac=yt, pitch=0.19)
    fig.text(0.005, yt, 'c', fontsize=10, fontweight='bold', va='top')
    save(fig, 'C')


def save(fig, tag):
    for ext in ('png', 'pdf'):
        fig.savefig(OUT / f'fig2_v2_{tag}.{ext}', dpi=300, facecolor='white', bbox_inches='tight',
                    pad_inches=0.04)
    from PIL import Image
    im = Image.open(OUT / f'fig2_v2_{tag}.png')
    print('wrote', OUT / f'fig2_v2_{tag}.png', im.size, f'= {im.size[0] / 300:.2f} x {im.size[1] / 300:.2f} in')
    plt.close(fig)


if __name__ == '__main__':
    which = sys.argv[1:] or ['A', 'B', 'C']
    for w in which:
        globals()[f'candidate_{w}']()
