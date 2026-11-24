"""Common code for the descriptive structure figures (candidates A-D), in the style of Figure 2 v3
(code/figures/fig2_v3_*.py; Recipe A of figures/references/REFERENCES.md).

Everything is drawn in the family frame of data/calc/03_frames (GH13: 1AMY chain A; GH31: 2QMJ chain A):
x_frame = R x + t with R, t of superposition.tsv for the chain (family, pdb, chain). The viewing directions are
those of the approved fig2_v3 panels (chosen on 1G94 and 3WEO), carried into the frame by the rotation of those
two chains, so every candidate shows the site from the same side as fig2_v3.

Colours, ChimeraX helpers (orthographic camera set from Python, camera matrix, scene -> pixel projection) and
the reference-numbering helpers are imported from fig2_v3_common; the ChimeraX runner and the label placer are
re-implemented here so that nothing is written into the fig2_v3 folder.
"""
import csv
import subprocess
import sys
from pathlib import Path

import gemmi
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fig2_v3_common as C            # noqa: E402

P = C.P
CIF = C.CIF
D44 = C.D44
D45 = P / 'data/analysis/45_hydroxyl_map'
FR = P / 'data/calc/03_frames'
OUT = P / 'figures/candidates/descriptive'
WORK = OUT / 'work'
WORK.mkdir(parents=True, exist_ok=True)

FONT = 'Liberation Sans'
plt.rcParams.update({'font.family': FONT, 'font.size': 7, 'pdf.fonttype': 42, 'axes.linewidth': 0.6})
REF = {'GH13': ('1AMY', 'A'), 'GH31': ('2QMJ', 'A')}
# fig2_v3 viewing directions (in the coordinates of 1G94 and 3WEO), from fig2_v3_struct.VIEW
VIEW_ENTRY = {'GH13': ('1G94', 'A', [-0.2821, 0.6736, -0.6832]),
              'GH31': ('3WEO', 'A', [-0.7836, 0.5235, 0.3346])}
TEXT_COL = {'#E69F00': '#B87800', '#CC79A7': '#B5508E'}
HB_COL = '#4A4A4A'
BACKBONE = {'N', 'C', 'O', 'OXT'}


# ---------------------------------------------------------------- frames
_MATS = None


def mats():
    global _MATS
    if _MATS is None:
        _MATS = {}
        for r in csv.DictReader(open(FR / 'superposition.tsv'), delimiter='\t'):
            if not r['status'].startswith('accepted'):
                continue
            v = [float(x) for x in r['R11..R33_t1..t3'].split()]
            _MATS[(r['family'], r['pdb'], r['chain'])] = (np.array(v[:9]).reshape(3, 3), np.array(v[9:12]))
    return _MATS


def to_frame(fam, pdb, chain, X):
    R, t = mats()[(fam, pdb, chain)]
    return np.atleast_2d(X) @ R.T + t


def cx_matrix(fam, pdb, chain):
    """ChimeraX 'view matrix models' argument (3x4 row-major) placing the model in the frame."""
    R, t = mats()[(fam, pdb, chain)]
    M = np.column_stack([R, t])
    return ','.join(f'{v:.5f}' for v in M.ravel())


def view_dir(fam):
    pdb, ch, v = VIEW_ENTRY[fam]
    R, _ = mats()[(fam, pdb, ch)]
    d = R @ np.array(v, float)
    return d / np.linalg.norm(d)


# ---------------------------------------------------------------- tables
def read_tsv(path):
    return list(csv.DictReader(open(path), delimiter='\t'))


def cell_states(tag='all_glucoselike'):
    out = {}
    for r in read_tsv(D45 / 'cell_summary.tsv'):
        if r['set'] == tag:
            out[(r['family'], r['subsite'], r['label'])] = r
    return out


def position_identity():
    return {(r['family'], r['position']): r for r in read_tsv(D45 / 'position_identity.tsv')}


def is_fixed(fam, pos, ident, nmin=3):
    """45's criterion for a lineage-invariant residue: n_cl50 >= NMIN (3) and m50 = 1.0 (INVM)."""
    r = ident.get((fam, pos))
    return bool(r) and float(r['n_cl50']) >= nmin and float(r['m50']) >= 1.0


def units_labelled():
    return read_tsv(D45 / 'units_labelled.tsv')


def unit_labels(row):
    """label -> atom name from the 'labels' column of units_labelled (e.g. 'O5=C7B(C)')."""
    out = {}
    for tok in row['labels'].split():
        lab, rest = tok.split('=')
        out[lab] = rest.split('(')[0]
    return out


def structure(pdb):
    st = gemmi.read_structure(str(CIF / f'{pdb}.cif'))
    st.remove_hydrogens()
    return st


def residue_atoms(model, chain, seq, name=None, keep=None):
    r = C.find(model, chain, seq, name)
    return {a.name: np.array(a.pos.tolist()) for a in C.first_alt_heavy(r) if keep is None or a.name in keep}, r


# ---------------------------------------------------------------- ChimeraX
def run_chimerax(lines, stem, timeout=1800):
    py = WORK / f'{stem}_cx.py'
    body = ['from chimerax.core.commands import run']
    for ln in lines:
        body.append(ln[3:].strip() if ln.startswith('PY:') else f'run(session, {ln!r})')
    py.write_text('\n'.join(body) + '\n')
    r = subprocess.run(['chimerax', '--offscreen', '--nogui', '--exit', str(py)], capture_output=True, text=True,
                       stdin=subprocess.DEVNULL, timeout=timeout)
    txt = r.stdout + r.stderr
    if r.returncode != 0 or 'Traceback' in txt or 'ERROR' in txt:
        print(txt[-3000:])
    return py


def camera(B, c, W, px):
    return [f'windowsize {px} {px}', *C.ortho_cmds(W), C.camera_cmd(B, c)]


def field(B, pts, margin, aspect):
    v = pts @ B.T
    lo, hi = v.min(0), v.max(0)
    mid = (lo + hi) / 2
    ext = hi - lo
    W = max(ext[0] + 2 * margin, (ext[1] + 2 * margin) / aspect)
    c = mid @ B
    z = (pts - c) @ B[2]
    return c, W, float(z.min()), float(z.max())


def bild_cylinders(path, segs, rgba, r):
    """segs: list of (p, q); rgba: list of '#rrggbb' or (hex, transparency)."""
    L = []
    last = None
    for (p, q), col in zip(segs, rgba):
        hexc, tr = (col, 0.0) if isinstance(col, str) else col
        if (hexc, tr) != last:
            rr, gg, bb = (int(hexc[i:i + 2], 16) / 255 for i in (1, 3, 5))
            L.append(f'.color {rr:.3f} {gg:.3f} {bb:.3f}')
            L.append(f'.transparency {tr:.2f}')
            last = (hexc, tr)
        L.append('.cylinder ' + ' '.join(f'{v:.3f}' for v in (*p, *q)) + f' {r:.3f}')
    Path(path).write_text('\n'.join(L) + '\n')


# ---------------------------------------------------------------- labels
def text_size_px(txt, fs, pt, weight='normal'):
    fp = FontProperties(family=FONT, weight=weight, size=fs)
    e = TextPath((0, 0), txt, size=fs, prop=fp).get_extents()
    return e.width * pt + 4, fs * 0.95 * pt + 2


def segment_samples(X, bond_cut=1.9):
    """Atom positions plus three samples on each bond (distance < bond_cut)."""
    names = list(X) if isinstance(X, dict) else range(len(X))
    arr = np.array([X[n] for n in names])
    out = list(arr)
    for i in range(len(arr)):
        for j in range(i + 1, len(arr)):
            if np.linalg.norm(arr[i] - arr[j]) < bond_cut:
                for t in (0.25, 0.5, 0.75):
                    out.append(arr[i] * (1 - t) + arr[j] * t)
    return out


def place_labels(items, OP, OW, OR, crop, scale, forbid=(), npass=5, leader_min=1.0):
    """Greedy label placement (the algorithm of fig2_v3_struct.place_labels, generalised).

    items: dicts with txt, ancs (k x 2 px), col, fs, weight, kind ('res' or 'sub'), w, h (px).
    OP, OW, OR: obstacle points (px), weights, radii (px). crop: page box in px. scale: px per A.
    forbid: boxes (x0, y0, x1, y1) no label may touch (thumbnail, legend)."""
    x0, y0, x1, y1 = crop
    angles = np.radians(np.arange(0, 360, 15))

    def cands(it):
        out = []
        rs = it.get('rs') or ([0.0, 0.4, 0.8, 1.3] if it['kind'] == 'sub' else [0.25, 0.6, 1.0, 1.6, 2.3, 3.0])
        for anc in it['ancs']:
            for r in rs:
                for a in (angles if r > 0 else [0.0]):
                    d = np.array([np.cos(a), -np.sin(a)])
                    if r > 0:
                        t = min(it['w'] / 2 / max(abs(d[0]), 1e-6), it['h'] / 2 / max(abs(d[1]), 1e-6))
                        ctr = anc + d * (r * scale + t)
                    else:
                        ctr = anc.copy()
                    out.append((ctr, r, anc))
        return out

    def box_of(ctr, it):
        return (ctr[0] - it['w'] / 2, ctr[1] - it['h'] / 2, ctr[0] + it['w'] / 2, ctr[1] + it['h'] / 2)

    def overlap(b1, b2):
        dx = min(b1[2], b2[2]) - max(b1[0], b2[0])
        dy = min(b1[3], b2[3]) - max(b1[1], b2[1])
        return max(dx, 0) * max(dy, 0)

    def cost(it, ctr, r, anc, placed):
        b = box_of(ctr, it)
        pad = 3
        if b[0] < x0 + pad or b[2] > x1 - pad or b[1] < y0 + pad or b[3] > y1 - pad:
            return 1e6
        for fb in forbid:
            if overlap(b, (fb[0] - 4, fb[1] - 4, fb[2] + 4, fb[3] + 4)) > 0:
                return 1e6
        inside = ((OP[:, 0] + OR > b[0]) & (OP[:, 0] - OR < b[2]) & (OP[:, 1] + OR > b[1]) & (OP[:, 1] - OR < b[3]))
        s = 10 * OW[inside].sum()
        for ob in placed:
            if overlap(b, ob) > 0:
                s += 1e5
            s += 0.02 * overlap((b[0] - 8, b[1] - 8, b[2] + 8, b[3] + 8), ob)
        s += 12.0 * r * (1 if it['kind'] == 'res' else 2)
        if it['kind'] == 'res' and r >= leader_min:
            e = np.array([np.clip(anc[0], b[0], b[2]), np.clip(anc[1], b[1], b[3])])
            v = e - anc
            L = np.linalg.norm(v)
            if L > 1:
                t = np.clip(((OP - anc) @ v) / L ** 2, 0, 1)
                dist = np.linalg.norm(anc + t[:, None] * v - OP, axis=1)
                hit = (dist < OR * 0.8) & (t > 0.15) & (t < 0.98)
                s += 6 * OW[hit].sum()
                tt = np.linspace(0, 1, 12)[:, None]
                pts = anc + tt * v
                for ob in placed:
                    if np.any((pts[:, 0] > ob[0]) & (pts[:, 0] < ob[2]) & (pts[:, 1] > ob[1]) & (pts[:, 1] < ob[3])):
                        s += 1e4
        return s

    sol = {}
    order = sorted(range(len(items)), key=lambda k: (items[k]['kind'] != 'sub', k))
    for _ in range(npass):
        for k in order:
            placed = [box_of(sol[j][0], items[j]) for j in sol if j != k]
            sol[k] = min(cands(items[k]), key=lambda cr: cost(items[k], cr[0], cr[1], cr[2], placed))
    res = []
    for k, it in enumerate(items):
        ctr, r, anc = sol[k]
        placed = [box_of(sol[j][0], items[j]) for j in sol if j != k]
        res.append(dict(it, ctr=ctr, r=r, anc=anc, cost=cost(it, ctr, r, anc, placed), box=box_of(ctr, it)))
    return res


def draw_labels(ax, labs, leader_min=1.0):
    import matplotlib.patheffects as pe
    halo = [pe.withStroke(linewidth=2.0, foreground='white')]
    for L in labs:
        ax.text(L['ctr'][0], L['ctr'][1], L['txt'], color=L['col'], fontsize=L['fs'], weight=L['weight'],
                ha='center', va='center', zorder=6, path_effects=halo, style=L.get('style', 'normal'))
        if L['kind'] == 'res' and L['r'] >= leader_min:
            b = L['box']
            tx = np.clip(L['anc'][0], b[0], b[2]); ty = np.clip(L['anc'][1], b[1], b[3])
            ax.plot([L['anc'][0], tx], [L['anc'][1], ty], color=L['col'], lw=0.5, zorder=5, solid_capstyle='round')
    bad = [L['txt'] for L in labs if L['cost'] > 1e4]
    if bad:
        print('    labels with overlap/edge cost:', bad)


def ref_label(resname, pos):
    """Present residue with reference number (fig2_v3 convention); unmapped residues by author number."""
    return C.label_text(resname, pos, None)
