"""Shared code for the Figure 2 rendering-style variants (one structure panel, GH13 entry 1G94).

Every atom drawn comes from the deposited mmCIF (data/census/cif/1G94.cif, model 1, first altloc, no H),
selected exactly as in fig2_v2_struct.py: the protein residues mapped to the figure positions
(data/figures/fig2_v2_work/contact_rows.tsv) and the sugar units -1, +1, +2. Hydrogen bonds: side-chain
N/O of a shown residue within 3.5 A of a sugar O/N, closest partner per side-chain atom (fig2_v2_struct.prepare).

What this module adds to fig2_v2_struct.py:
 * declutter_view(): a search over viewing directions that counts on-page crossings of covalent bonds of
   different groups (residues, sugar units) and close approaches of atoms of different groups, keeping the
   ring-tilt and stacking constraints of the earlier search (stricter: -1 ring |n.z| 0.55-0.95, +1 ring >= 0.50,
   the aromatic ring of Tyr51 not edge-on, |n.z| >= 0.40); the -1 -> +1 vector is then turned to point right.
 * a slight perspective (field of view FOV) used by every 3D renderer; the label overlay projects atoms with the
   same perspective, so the labels sit on the rendered atoms.
 * compose(): one 3.4 in panel (matplotlib) with the rendered image, the labels (greedy placement of
   fig2_v2_struct.place_labels), and a one-line class legend; written at 300 dpi and 900 dpi.

Run with ~/miniconda3/bin/python (gemmi, matplotlib, numpy). Renderers: PyMOL (pymol -cq) and ChimeraX
(chimerax --offscreen --nogui).
"""
import math
import sys
import subprocess
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fig2_v2_struct as S

P = S.P
PDB = '1G94'
OUT = P / 'figures/candidates/styles'
OUT.mkdir(parents=True, exist_ok=True)
WORK = P / 'data/figures/fig2_style_work'
WORK.mkdir(parents=True, exist_ok=True)

COL = S.COL                      # Okabe-Ito: H blue, Ar orange, Al reddish purple
SUGAR_C = '#D4D4D4'              # sugar carbons: neutral light grey
O_COL = '#E0242B'
N_COL = '#2F4FD8'
OLAB = '#8A1414'
PX = 2400                        # render size (square)
FOV = 16.0                       # degrees; slight perspective
PANEL_W = 3.4                    # inches

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7, 'pdf.fonttype': 42, 'axes.linewidth': 0.6})


# ------------------------------------------------------------------ geometry
def _bonds(X, el=None):
    names = list(X)
    out = []
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            if np.linalg.norm(X[names[i]] - X[names[j]]) < 1.9:
                out.append((names[i], names[j]))
    return out


def groups(info):
    """List of (group id, {name: xyz}, bonds) for residues and sugar units; inter-unit bonds go to group -1."""
    G = []
    for k, rd in enumerate(info['residues']):
        G.append((100 + k, rd['X'], _bonds(rd['X'])))
    for k, u in enumerate(info['units']):
        G.append((k, u['X'], _bonds(u['X'])))
    return G


def bond_segments(info):
    segs, gid = [], []
    for g, X, bl in groups(info):
        for a, b in bl:
            segs.append((X[a], X[b])); gid.append(g)
    # glycosidic / N-linked bonds between units
    us = info['units']
    for i in range(len(us)):
        for j in range(i + 1, len(us)):
            for a, x in us[i]['X'].items():
                for b, y in us[j]['X'].items():
                    if np.linalg.norm(x - y) < 1.9:
                        segs.append((x, y)); gid.append(-1)
    return np.array(segs), np.array(gid)


def _crossings(A, B, C, D):
    """Vectorised 2D segment intersection count between segments AB (n) and CD (n)."""
    def orient(p, q, r):
        return np.sign((q[..., 0] - p[..., 0]) * (r[..., 1] - p[..., 1]) - (q[..., 1] - p[..., 1]) * (r[..., 0] - p[..., 0]))
    o1 = orient(A, B, C); o2 = orient(A, B, D); o3 = orient(C, D, A); o4 = orient(C, D, B)
    return (o1 * o2 < 0) & (o3 * o4 < 0)


def declutter_view(info, n=12000):
    """Search viewing directions; returns the best basis (rows ex, ey, ez) and a report dict."""
    segs, gid = bond_segments(info)
    iu, ju = np.triu_indices(len(segs), 1)
    keep = (gid[iu] != gid[ju])
    iu, ju = iu[keep], ju[keep]
    center = info['center']
    rings = info['rings']
    c_m1, n_m1 = rings['-1']; c_p1, n_p1 = rings['+1']
    srd, sc_ar, sc_su = info['stack']
    sn_su = rings[info['spec']['stack'][1]][1]
    n_ar = S.plane(np.array([srd['X'][n] for n in S.AROM[srd['resname']]]))[1]
    anchors = np.array([rings[s][0] for s in rings] + [rd['centroid'] for rd in info['residues']])
    atoms, agid = [], []
    for g, X, _ in groups(info):
        for x in X.values():
            atoms.append(x); agid.append(g)
    atoms = np.array(atoms); agid = np.array(agid)
    diffg = agid[:, None] != agid[None, :]

    def score(ez, report=False):
        B = S.basis(ez)
        Sg = (segs - center) @ B.T
        cr = int(_crossings(Sg[iu, 0, :2], Sg[iu, 1, :2], Sg[ju, 0, :2], Sg[ju, 1, :2]).sum())
        Q = (atoms - center) @ B.T
        dd = np.linalg.norm(Q[:, None, :2] - Q[None, :, :2], axis=2)
        close = int(((dd < 0.9) & diffg).sum() // 2)
        A = (anchors - center) @ B.T
        D = np.linalg.norm(A[:, None, :2] - A[None, :, :2], axis=2) + 1e3 * np.eye(len(A))
        gap = np.linalg.norm(((sc_ar - sc_su) @ B.T)[:2])
        side = 1 - abs(sn_su @ ez)
        tilt = abs(n_m1 @ ez); tilt1 = abs(n_p1 @ ez)
        s = (-1.0 * cr - 0.5 * close + 1.0 * min(D.min(), 3.0) + 0.6 * min(gap, 4.5) + 0.5 * side
             + 0.3 * tilt + 0.3 * tilt1)
        if report:
            return dict(crossings=cr, close=close, min_anchor=float(D.min()), stack_gap=float(gap), tilt=float(tilt), tilt1=float(tilt1), arom=float(abs(n_ar @ ez)),
                        score=float(s))
        return s

    old = np.asarray(info['spec']['view'], float); old /= np.linalg.norm(old)
    best, bs = old, -1e9
    for ez in S.fib_dirs(n):
        tilt = abs(n_m1 @ ez); tilt1 = abs(n_p1 @ ez)
        # rings readable: -1 and +1 rings tilted (neither edge-on nor flat), aromatic ring not edge-on
        if not (0.55 <= tilt <= 0.95 and tilt1 >= 0.50 and abs(n_ar @ ez) >= 0.40):
            continue
        s = score(ez)
        if s > bs:
            best, bs = ez, s
    rep = dict(old=score(old, True), new=score(best, True), ez=best.tolist())
    B = S.basis(best)
    v = (rings['+1'][0] - c_m1) @ B.T
    ang = math.atan2(v[1], v[0])
    Rz = np.array([[math.cos(-ang), -math.sin(-ang), 0], [math.sin(-ang), math.cos(-ang), 0], [0, 0, 1]])
    return Rz @ B, rep


_CACHE = WORK / 'view_1G94.npy'


def setup(declutter=True):
    """prepare() + decluttered view (cached) + rotated PDB; returns info with 'half', 'dist', 'pdbfile'."""
    info = S.prepare(PDB)
    if declutter:
        if _CACHE.exists():
            B = np.load(_CACHE)
        else:
            B, rep = declutter_view(info)
            np.save(_CACHE, B)
            print('declutter: old view', rep['old'], '\n           new view', rep['new'], '\n           ez', rep['ez'])
        info['B'] = B
    pdbfile = WORK / f'{PDB}_style_rot.pdb'
    S.write_rotated(info, pdbfile)            # rigid rotation only (before the perspective is switched on)
    info['pdbfile'] = pdbfile
    pts = np.array([S.rot(info, x) for rd in info['residues'] for x in rd['X'].values()] +
                   [S.rot(info, x) for u in info['units'] for x in u['X'].values()])
    # centre the field on the shown atoms
    c = (pts.max(0) + pts.min(0)) / 2
    c[2] = 0
    info['center'] = info['center'] + c @ info['B']
    S.write_rotated(info, pdbfile)
    pts = pts - c
    half = float(np.abs(pts[:, :2]).max() + 2.6)
    info['half'] = half
    info['dist'] = half / math.tan(math.radians(FOV / 2))
    info['zrange'] = (float(pts[:, 2].min()), float(pts[:, 2].max()))
    return info


def proj(info, x):
    """Rotated coordinates with the perspective scale applied to x, y (z kept)."""
    v = (np.asarray(x) - info['center']) @ info['B'].T
    f = 1.0 if info.get('ortho') else info['dist'] / (info['dist'] - v[2])
    return np.array([v[0] * f, v[1] * f, v[2]])


def tp(info, x, px=PX):
    return S.to_px(proj(info, x)[:2], info['half'], px)


# ------------------------------------------------------------------ selections
def res_sel_pymol(rd):
    return f'(polymer.protein and chain {rd["chain"]} and resi {rd["seq"]} and not name N+C+O+OXT)'


def unit_sel_pymol(u):
    return f'(chain {u["chain"]} and resi {u["seq"]} and resn {u["comp"]} and name {"+".join(u["X"])})'


def res_sel_cx(rd):
    return f'/{rd["chain"]}:{rd["seq"]}@{",".join(rd["X"])}'


def unit_sel_cx(u):
    return f'/{u["chain"]}:{u["seq"]}@{",".join(u["X"])}'


def hbond_pairs(info):
    out = []
    for hb in info['hbonds']:
        rd, u = hb['res'], hb['unit']
        out.append(((rd['chain'], rd['seq'], hb['a']), (u['chain'], u['seq'], u['comp'], hb['b']), hb['d']))
    return out


def hexrgb(h):
    return S.hexcol(h)


def run_pymol(lines, stem):
    pml = WORK / f'{stem}.pml'
    pml.write_text('\n'.join(lines) + '\n')
    subprocess.run(['pymol', '-cq', str(pml)], check=True, stdout=subprocess.DEVNULL)


def run_chimerax(lines, stem):
    cxc = WORK / f'{stem}.cxc'
    cxc.write_text('\n'.join(lines + ['exit']) + '\n')
    r = subprocess.run(['chimerax', '--offscreen', '--nogui', '--silent', str(cxc)], capture_output=True, text=True)
    if r.returncode != 0 or 'Error' in r.stdout + r.stderr:
        print(r.stdout[-3000:], r.stderr[-3000:])


def pymol_camera(info, slab=None):
    d = info['dist']
    z0, z1 = info['zrange']
    front, back = d - z1 - 3.0, d - z0 + 3.0
    if slab is not None:
        front, back = d - slab, d + slab
    return [f'set field_of_view, {FOV}', 'set orthoscopic, 0',
            f'set_view (1,0,0, 0,1,0, 0,0,1, 0,0,{-d:.3f}, 0,0,0, {front:.3f},{back:.3f}, {-FOV})']


def chimerax_camera(info, png, px=PX, supersample=3):
    d = info['dist']
    z0, z1 = info['zrange']
    cam = (f'camera ortho fieldWidth {2 * info["half"]:.3f}' if info.get('ortho') else f'camera fieldOfView {FOV}')
    return [f'windowsize {px} {px}', cam,
            f'view matrix camera 1,0,0,0,0,1,0,0,0,0,1,{d:.3f}',
            f'clip near {-(z1 + 3.0):.2f} far {-(z0 - 6.0):.2f} position 0,0,0',
            f'save {png} width {px} height {px} supersample {supersample} transparentBackground false']


# ------------------------------------------------------------------ panel composition
def _bonds_mid(X):
    pts = []
    for a, b in _bonds(X):
        pts += [X[a] * 0.75 + X[b] * 0.25, (X[a] + X[b]) / 2, X[a] * 0.25 + X[b] * 0.75]
    return pts


def sub_txt(s):
    return s.replace('-', '−')


def add_labels(ax, info, bounds, fs=7, fs_o=6.5, olabels=True, label_box=True, res_weight='bold',
               ring_color='black', px=PX):
    fig = ax.figure
    rnd = fig.canvas.get_renderer()
    sc = px / (2 * info['half'])
    T = lambda x: tp(info, x, px)
    bb = ax.get_window_extent(rnd)
    dpd = (bounds[2] - bounds[0]) / bb.width

    def measure(txt, size, weight='normal'):
        t = ax.text(0, 0, txt, fontsize=size, fontweight=weight)
        e = t.get_window_extent(rnd); t.remove()
        return e.width * dpd + 6 * dpd, e.height * dpd + 2 * dpd

    obs = []
    for k, rd in enumerate(info['residues']):
        for x in list(rd['X'].values()) + _bonds_mid(rd['X']):
            p = T(x); obs.append([p[0], p[1], 0.55 * sc, 100 + k])
    for k, u in enumerate(info['units']):
        for x in list(u['X'].values()) + _bonds_mid(u['X']):
            p = T(x); obs.append([p[0], p[1], 0.6 * sc, k])
    obs = np.array(obs)
    segs = [(T(hb['res']['X'][hb['a']]), T(hb['unit']['X'][hb['b']])) for hb in info['hbonds']]
    sugar_c = np.mean([T(x) for u in info['units'] for x in u['X'].values()], 0)
    items = []
    for k, rd in enumerate(info['residues']):
        sca = {n: x for n, x in rd['X'].items() if n != 'CA'} or rd['X']
        dist = {n: np.linalg.norm(T(x) - sugar_c) for n, x in sca.items()}
        far = max(sca, key=lambda n: dist[n])
        anc = T(sca[far])
        med = np.median(list(dist.values()))
        anchors = [T(x) for n, x in sca.items() if dist[n] >= med] or [anc]
        w, h = measure(rd['label'], fs, res_weight)
        v = anc - sugar_c; v = v / (np.linalg.norm(v) + 1e-9)
        items.append(dict(kind='res', text=rd['label'], color=COL[rd['cls']], anchor=anc, anchors=anchors, w=w, h=h,
                          radii=[0.4, 0.9, 1.5, 2.2, 3.0, 3.8], prefer=v, priority=2, size=fs, weight=res_weight,
                          group=100 + k))
    for k, u in enumerate(info['units']):
        c = np.mean([T(u['X'][n]) for n in u['ring']], 0)
        w, h = measure(sub_txt(u['sub']), fs, 'bold')
        items.append(dict(kind='ring', text=sub_txt(u['sub']), color=ring_color, anchor=c, w=w, h=h,
                          radii=[1.1, 1.5, 2.0, 2.6], prefer=None, priority=3, size=fs, weight='bold'))
    if olabels:
        seen = set()
        for hb in info['hbonds']:
            u, b = hb['unit'], hb['b']
            if (u['sub'], b) in seen:
                continue
            seen.add((u['sub'], b))
            name = u['names'].get(b, b)
            name = 'N' if name.startswith('N') else name
            w, h = measure(name, fs_o)
            items.append(dict(kind='O', text=name, color=OLAB, anchor=T(u['X'][b]), w=w, h=h,
                              radii=[0.35, 0.7, 1.1, 1.6], prefer=None, priority=1, size=fs_o, weight='normal'))
    boxes, total = S.place_labels(items, obs, bounds, sc, segs=segs)
    for it, b in zip(items, boxes):
        cx, cy = (b[0] + b[2]) / 2, (b[1] + b[3]) / 2
        st = S._attach(b, np.array(it['anchor']))
        lead = np.linalg.norm(st - it['anchor'])
        if it['kind'] == 'res' and lead > 0.35 * sc:
            ax.plot([st[0], it['anchor'][0]], [st[1], it['anchor'][1]], color=it['color'], lw=0.5, zorder=4,
                    solid_capstyle='round')
        if it['kind'] == 'O' and lead > 0.45 * sc:
            ax.plot([st[0], it['anchor'][0]], [st[1], it['anchor'][1]], color=OLAB, lw=0.4, zorder=4)
        kw = dict(bbox=dict(boxstyle='round,pad=0.06', fc='white', ec='none', alpha=0.85)) if label_box else {}
        ax.text(cx, cy, it['text'], ha='center', va='center', fontsize=it['size'], fontweight=it['weight'],
                color=it['color'], zorder=5, **kw)
    bad = [(items[i]['text'], items[j]['text']) for i in range(len(boxes)) for j in range(i + 1, len(boxes))
           if S._overlap(boxes[i], boxes[j]) > 0]
    return total, bad


def legend(fig, y, hb_style='dash', fs=6.5):
    h = [Line2D([], [], color=COL['H'], lw=3, label='H'), Line2D([], [], color=COL['Ar'], lw=3, label='Ar'),
         Line2D([], [], color=COL['Al'], lw=3, label='Al'),
         Line2D([], [], color=SUGAR_C, lw=3, label='sugar C', path_effects=None),
         Line2D([], [], color='0.15', lw=0.9, ls=(0, (1.6, 1.4)), label='N/O–N/O ≤ 3.5 Å')]
    fig.legend(handles=h, loc='center', bbox_to_anchor=(0.5, y), ncol=5, frameon=False, fontsize=fs,
               handlelength=1.6, columnspacing=1.0, handletextpad=0.4)


def compose(info, png, stem, margin=2.2, olabels=True, label_box=True, res_weight='bold', fs=7, fs_o=6.5,
            px=PX, extra_artists=None):
    """Panel 3.4 in wide: rendered image cropped to the shown atoms + margin, labels, legend."""
    sc = px / (2 * info['half'])
    pts = np.array([tp(info, x, px) for rd in info['residues'] for x in rd['X'].values()] +
                   [tp(info, x, px) for u in info['units'] for x in u['X'].values()])
    x0, y0 = pts.min(0) - margin * sc; x1, y1 = pts.max(0) + margin * sc
    x0, y0 = max(x0, 0), max(y0, 0); x1, y1 = min(x1, px), min(y1, px)
    leg_h = 0.20
    H = PANEL_W * (y1 - y0) / (x1 - x0) + leg_h
    fig = plt.figure(figsize=(PANEL_W, H))
    ax = fig.add_axes((0, leg_h / H, 1, 1 - leg_h / H))
    ax.imshow(plt.imread(png), extent=(0, px, px, 0), interpolation='lanczos')
    ax.set_xlim(x0, x1); ax.set_ylim(y1, y0); ax.axis('off')
    if extra_artists:
        extra_artists(ax, sc)
    total, bad = add_labels(ax, info, (x0, y0, x1, y1), fs=fs, fs_o=fs_o, olabels=olabels, label_box=label_box,
                            res_weight=res_weight, px=px)
    legend(fig, leg_h / H / 2)
    print(f'  {stem}: label cost {total:.0f}, overlapping label pairs {bad}')
    fig.savefig(OUT / f'{stem}.png', dpi=300)
    fig.savefig(OUT / f'{stem}_hires.png', dpi=900)
    plt.close(fig)
    return OUT / f'{stem}.png'


# ------------------------------------------------------------------ ChimeraX helpers
def cx_pbfile(info, path, color='#303030', radius=0.035, dashes=12):
    lines = [f'; radius = {radius}', f'; dashes = {dashes}', f'; color = {color}']
    for (rc, rs, ra), (uc, us, ucomp, ub), d in hbond_pairs(info):
        lines.append(f'/{rc}:{rs}@{ra} /{uc}:{us}@{ub}')
    Path(path).write_text('\n'.join(lines) + '\n')
    return path


def cx_atoms(info, stick=0.14, sugar_stick=0.10, sugar_ball=0.22, sugar_style='ball', hb=True):
    """ChimeraX commands that show and colour the residues, sugar units and hydrogen bonds."""
    L = [f'open {info["pdbfile"]}', 'hide atoms', 'hide cartoons', 'set bgColor white', 'style stick']
    for rd in info['residues']:
        s = res_sel_cx(rd)
        L += [f'show {s} atoms', f'color {s} & C {COL[rd["cls"]]}', f'size {s} stickRadius {stick}']
    for u in info['units']:
        s = unit_sel_cx(u)
        L += [f'show {s} atoms', f'style {s} {sugar_style}', f'color {s} & C {SUGAR_C}',
              f'size {s} stickRadius {sugar_stick} ballScale {sugar_ball}']
    L += [f'color #1 & O {O_COL}', f'color #1 & N {N_COL}']
    if hb:
        pb = cx_pbfile(info, WORK / 'hbonds_1G94.pb')
        L += [f'open {pb}']
    return L
