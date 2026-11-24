"""Figure 2 zoomed structure panels: the -1 and +1 sugar units and five key residues per family, larger.

Data (deposited coordinates only): data/census/cif/{1G94,3WEO}.cif, model 1, first altloc, no H; residues
mapped to the figure positions in data/figures/fig2_v2_work/contact_rows.tsv; sugar units -1/+1 as in
fig2_v2_struct.ENTRIES. Hydrogen bonds: side-chain N/O of a shown residue within 3.5 A of an O/N of the -1 or +1
unit, closest partner per side-chain atom (fig2_v2_struct.prepare).

Steps per family
 1. Residues: GH13 1G94 Asp179, Glu204, Asp289, Tyr51, Ala (Phe180 position) and one of His288 / Arg177;
    GH31 3WEO Asp443, Asp542, Phe575, Ile328 and one of Asp327 / His600. The alternative giving the cleaner
    best view (score below) is kept.
 2. View search over 20 000 directions (Fibonacci sphere). Kept only: -1 ring normal at |n.z| 0.45-0.95, +1 ring
    |n.z| >= 0.35, aromatic ring |n.z| >= 0.45 (not edge-on). Score: -3 x on-page crossings of bonds/H-bond
    dashes of different groups, -1 x atom pairs of different groups closer than 1.4 A on the page, +1.5 x the
    smallest such distance (capped 2 A), -1 x non-bonded sugar atom pairs closer
    than 0.9 A on the page (sugar drawn over itself), -2 x residue atoms/bond midpoints falling inside a
    sugar ring outline (grown by 0.5 A) on the page, +1 x the smallest distance between label anchors (capped 4 A).
    The -1 -> +1 vector is turned to point right.
 3. PyMOL ray trace, slight perspective (16 deg): side chains from CB (Ala CB as a small sphere), thin sticks
    (0.11 A), class-coloured carbons, N/O standard colours; sugar thin ball-and-stick, grey carbons; soft light,
    thin grey outlines; no cartoon.
 4. Labels (residues, ring positions) by the greedy placement of fig2_v2_struct.place_labels, with leaders.

Writes figures/candidates/styles/zoom_GH13.png, zoom_GH31.png (3.45 in, 300 dpi) and
zoom_pair.png (7 in, 300 dpi).
Run: ~/miniconda3/bin/python code/figures/fig2_zoom.py
"""
import math
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fig2_v2_struct as S
import fig2_style_common as C

OUT = C.OUT
WORK = C.WORK
PX = 2400
FOV = 16.0
SUGAR_C = '#BEBEBE'
O_COL = '#E0242B'
N_COL = '#2F4FD8'

PANELS = {
    # GH13 view: picked by eye from the 8 best-scoring distinct directions of the search (montage of the 2D
    # projections); the top-scoring one hid the -1 ring behind Asp289. His288 kept (with Arg177 every candidate
    # view put Arg177 across the -1 ring or the ring edge-on).
    'GH13': dict(pdb='1G94', fixed=['D179', 'E204', 'D289', 'Y51', 'F180'], alt=['H288'],
                 view={'H288': [0.432, -0.668, 0.606]}),
    'GH31': dict(pdb='3WEO', fixed=['D443', 'D542', 'F575', 'I328'], alt=['D327', 'H600']),
}
SUBS = ('-1', '+1')
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8, 'pdf.fonttype': 42})


def trim(info):
    info['units'] = [u for u in info['units'] if u['sub'] in SUBS]
    info['hbonds'] = [h for h in info['hbonds'] if h['unit']['sub'] in SUBS]
    return info


def view_search(info, n=20000, override=None):
    groups = []
    for k, rd in enumerate(info['residues']):
        X = {a: x for a, x in rd['X'].items() if a != 'CA'}
        groups.append((100 + k, X))
    for k, u in enumerate(info['units']):
        groups.append((k, u['X']))
    segs, gid = [], []
    for g, X in groups:
        for a, b in C._bonds(X):
            segs.append((X[a], X[b])); gid.append(g)
    hbp = set()
    for h in info['hbonds']:
        segs.append((h['res']['X'][h['a']], h['unit']['X'][h['b']])); gid.append(-5 - len(segs))
        hbp.add((id(h['res']), h['a'], id(h['unit']), h['b']))
    segs = np.array(segs); gid = np.array(gid)
    iu, ju = np.triu_indices(len(segs), 1)
    # skip pairs that share an atom (H-bond dash touching its own atoms)
    share = np.array([np.allclose(segs[i, 0], segs[j, 0]) or np.allclose(segs[i, 0], segs[j, 1]) or
                      np.allclose(segs[i, 1], segs[j, 0]) or np.allclose(segs[i, 1], segs[j, 1])
                      for i, j in zip(iu, ju)])
    keep = (gid[iu] != gid[ju]) & ~share
    # sugar units -1/+1 are covalently linked: treat them as one group for crossings
    sug = (gid[iu] >= 0) & (gid[iu] < 100) & (gid[ju] >= 0) & (gid[ju] < 100)
    keep &= ~sug
    iu, ju = iu[keep], ju[keep]
    atoms, agid = [], []
    for g, X in groups:
        for x in X.values():
            atoms.append(x); agid.append(g if g >= 100 else 0)
    atoms = np.array(atoms); agid = np.array(agid)
    diffg = agid[:, None] != agid[None, :]
    rings = info['rings']
    n_m1 = rings['-1'][1]; n_p1 = rings['+1'][1]
    srd = info['stack'][0]
    n_ar = S.plane(np.array([srd['X'][a] for a in S.AROM[srd['resname']]]))[1]
    anchors = np.array([rings[s][0] for s in SUBS] + [rd['centroid'] for rd in info['residues']])
    center = np.vstack([atoms, anchors]).mean(0)
    # sugar self-overlap: non-bonded sugar atom pairs (both units) that fall on top of each other on the page
    sX = np.array([x for u in info['units'] for x in u['X'].values()])
    sd3 = np.linalg.norm(sX[:, None] - sX[None], axis=2)
    nb = np.triu((sd3 > 1.9), 1)
    from matplotlib.path import Path as MPath
    ringX = [np.array([u['X'][a] for a in u['ring']]) for u in info['units']]
    # residue atoms + bond midpoints (side chains from CB)
    rX = []
    for g, X in groups:
        if g >= 100:
            rX += list(X.values()) + [(X[a] + X[b]) / 2 for a, b in C._bonds(X)]
    rX = np.array(rX)
    best, bs, brep = None, -1e9, None
    allc = []
    for ez in S.fib_dirs(n):
        t1, t2, ta = abs(n_m1 @ ez), abs(n_p1 @ ez), abs(n_ar @ ez)
        if not (0.45 <= t1 <= 0.95 and t2 >= 0.35 and ta >= 0.45):
            continue
        B = S.basis(ez)
        Sg = (segs - center) @ B.T
        cr = int(C._crossings(Sg[iu, 0, :2], Sg[iu, 1, :2], Sg[ju, 0, :2], Sg[ju, 1, :2]).sum())
        Q = (atoms - center) @ B.T
        dd = np.linalg.norm(Q[:, None, :2] - Q[None, :, :2], axis=2)
        dd = np.where(diffg, dd, 99)
        close = int((dd < 1.4).sum() // 2)
        dmin = float(dd.min())
        A = (anchors - center) @ B.T
        D = np.linalg.norm(A[:, None, :2] - A[None, :, :2], axis=2) + 1e3 * np.eye(len(A))
        Sx = (sX - center) @ B.T
        sd = np.linalg.norm(Sx[:, None, :2] - Sx[None, :, :2], axis=2)
        selfov = int(((sd < 0.9) & nb).sum())
        Rr = (rX - center) @ B.T
        ringhit = 0
        for R in ringX:
            Rp = (R - center) @ B.T
            # ring polygon grown by 0.5 A so atoms touching the ring edge also count
            cen = Rp[:, :2].mean(0)
            poly = cen + (Rp[:, :2] - cen) * (1 + 0.5 / np.linalg.norm(Rp[:, :2] - cen, axis=1, keepdims=True))
            ringhit += int(MPath(poly).contains_points(Rr[:, :2]).sum())
        s = (-3 * cr - 1.0 * close - 1.0 * selfov - 2.0 * ringhit + 1.5 * min(dmin, 2.0) + 1.0 * min(D.min(), 4.0)
             + 0.3 * t1 + 0.2 * ta)
        allc.append((s, ez))
        if s > bs:
            best, bs = ez, s
            brep = dict(crossings=cr, close=close, dmin=round(dmin, 2), anchor_min=round(float(D.min()), 2),
                        selfov=selfov, ringhit=ringhit, tilt_m1=round(float(t1), 2), tilt_p1=round(t2, 2), arom=round(ta, 2), score=round(s, 2))
    if override is not None:
        best = np.asarray(override, float); best /= np.linalg.norm(best)
    info['allc'] = allc
    B = S.basis(best)
    v = (rings['+1'][0] - rings['-1'][0]) @ B.T
    ang = math.atan2(v[1], v[0])
    Rz = np.array([[math.cos(-ang), -math.sin(-ang), 0], [math.sin(-ang), math.cos(-ang), 0], [0, 0, 1]])
    return Rz @ B, center, bs, brep, best


def setup(fam):
    spec = PANELS[fam]
    cands = []
    for a in spec['alt']:
        info = trim(S.prepare(spec['pdb'], positions=spec['fixed'] + [a]))
        B, center, sc, rep, ez = view_search(info, override=spec.get('view', {}).get(a))
        print(f'{fam} with {a}: {rep}')
        cands.append((sc, a, info, B, center, rep, ez))
    sc, a, info, B, center, rep, ez = max(cands, key=lambda c: c[0])
    print(f'{fam}: kept {a}; view direction {np.round(ez, 4).tolist()}')
    info['B'], info['center'] = B, center
    info['choice'], info['viewrep'] = a, rep
    pts = np.array([S.rot(info, x) for rd in info['residues'] for n, x in rd['X'].items() if n != 'CA'] +
                   [S.rot(info, x) for u in info['units'] for x in u['X'].values()])
    c = (pts.max(0) + pts.min(0)) / 2; c[2] = 0
    info['center'] = center + c @ B
    pts = pts - c
    info['half'] = float(np.abs(pts[:, :2]).max() + 3.0)
    info['dist'] = info['half'] / math.tan(math.radians(FOV / 2))
    info['zrange'] = (float(pts[:, 2].min()), float(pts[:, 2].max()))
    info['pdbfile'] = WORK / f'zoom_{spec["pdb"]}_rot.pdb'
    S.write_rotated(info, info['pdbfile'])
    for rd in info['residues']:          # side chains from CB: CA is not drawn and not used for labels
        if len(rd['X']) > 2:
            rd['X'] = {k: v for k, v in rd['X'].items() if k != 'CA'}
        else:
            rd['X'] = {k: v for k, v in rd['X'].items() if k != 'CA'}
    return info


def render(info, png, sil=2.0):
    """ChimeraX offscreen: lighting soft (ambient occlusion, 64 shadow directions), thin silhouettes."""
    L = [f'open {info["pdbfile"]}', 'hide atoms', 'hide cartoons', 'set bgColor white', 'style stick']
    for rd in info['residues']:
        s = f'/{rd["chain"]}:{rd["seq"]}@{",".join(rd["X"])}'
        L += [f'show {s} atoms', f'color {s} & C {C.COL[rd["cls"]]}', f'size {s} stickRadius 0.12']
        if rd['resname'] == 'ALA':
            L += [f'style {s} ball', f'size {s} ballScale 0.30']
    for u in info['units']:
        s = C.unit_sel_cx(u)
        L += [f'show {s} atoms', f'style {s} ball', f'color {s} & C {SUGAR_C}',
              f'size {s} stickRadius 0.085 ballScale 0.14']
    L += [f'color #1 & O {O_COL}', f'color #1 & N {N_COL}']
    pb = WORK / f'zoom_{info["pdb"]}_hb.pb'
    pbl = ['; radius = 0.03', '; dashes = 12', '; color = #3a3a3a']
    for hb in info['hbonds']:
        rd, u = hb['res'], hb['unit']
        pbl.append(f'/{rd["chain"]}:{rd["seq"]}@{hb["a"]} /{u["chain"]}:{u["seq"]}@{hb["b"]}')
    pb.write_text('\n'.join(pbl) + '\n')
    L += [f'open {pb}', 'lighting soft', 'lighting intensity 0.45 ambientIntensity 1.3', 'lighting depthCue false', 'material dull',
          f'graphics silhouettes true width {sil} color #303030 depthJump 0.01', 'graphics quality 20']
    d = info['dist']
    L += [f'windowsize {PX} {PX}', f'camera fieldOfView {FOV}',
          f'view matrix camera 1,0,0,0,0,1,0,0,0,0,1,{d:.3f}', 'clip off',
          f'save {png} width {PX} height {PX} supersample 3 transparentBackground false']
    C.run_chimerax(L, Path(png).stem)


def render_pymol(info, png):
    L = ['reinitialize', f'load {info["pdbfile"]}, m', 'remove hydro', 'hide everything', 'bg_color white',
         'set ray_opaque_background, 1', 'set antialias, 2',
         'set light_count, 4', 'set ambient, 0.50', 'set direct, 0.30', 'set reflect, 0.45', 'set specular, 0.18',
         'set spec_power, 90', 'set spec_count, 1', 'set ray_shadows, 0', 'set depth_cue, 0',
         'set ray_trace_mode, 1', 'set ray_trace_gain, 0.10', 'set ray_trace_color, grey45',
         'set stick_radius, 0.11', 'set stick_quality, 40', 'set sphere_quality, 4', 'set valence, 0',
         'set dash_gap, 0.17', 'set dash_length, 0.13', 'set dash_radius, 0.028', 'set dash_color, grey25',
         'set dash_round_ends, 1',
         f'set_color sugc, {S.hexcol(SUGAR_C)}', f'set_color ocol, {S.hexcol(O_COL)}',
         f'set_color ncol, {S.hexcol(N_COL)}']
    for u in info['units']:
        s = C.unit_sel_pymol(u)
        L += [f'show sticks, {s}', f'show spheres, {s}', f'set sphere_scale, 0.15, {s}',
              f'set stick_radius, 0.075, {s}', f'color sugc, {s} and elem C', f'color ocol, {s} and elem O',
              f'color ncol, {s} and elem N']
    for rd in info['residues']:
        s = f'(polymer.protein and chain {rd["chain"]} and resi {rd["seq"]} and not name N+C+O+OXT+CA)'
        L += [f'set_color cls_{rd["cls"]}, {S.hexcol(C.COL[rd["cls"]])}', f'show sticks, {s}',
              f'color cls_{rd["cls"]}, {s} and elem C', f'color ocol, {s} and elem O', f'color ncol, {s} and elem N+S']
        if rd['resname'] == 'ALA':
            L += [f'show spheres, {s}', f'set sphere_scale, 0.30, {s}']
    for k, hb in enumerate(info['hbonds']):
        rd, u = hb['res'], hb['unit']
        L.append(f'distance hb{k}, polymer.protein and chain {rd["chain"]} and resi {rd["seq"]} and name {hb["a"]}, '
                 f'chain {u["chain"]} and resi {u["seq"]} and resn {u["comp"]} and name {hb["b"]}')
    d = info['dist']; z0, z1 = info['zrange']
    L += ['hide labels', f'set field_of_view, {FOV}', 'set orthoscopic, 0',
          f'set_view (1,0,0, 0,1,0, 0,0,1, 0,0,{-d:.3f}, 0,0,0, {d - z1 - 3:.3f},{d - z0 + 3:.3f}, {-FOV})',
          f'ray {PX},{PX}', f'png {png}, dpi=300']
    C.run_pymol(L, Path(png).stem)


def draw(ax, info, png, fs=8.5, margin=2.4):
    sc = PX / (2 * info['half'])
    pts = np.array([C.tp(info, x) for rd in info['residues'] for x in rd['X'].values()] +
                   [C.tp(info, x) for u in info['units'] for x in u['X'].values()])
    x0, y0 = pts.min(0) - margin * sc; x1, y1 = pts.max(0) + margin * sc
    ax.imshow(plt.imread(png), extent=(0, PX, PX, 0), interpolation='lanczos')
    ax.set_xlim(x0, x1); ax.set_ylim(y1, y0); ax.axis('off')
    return (x0, y0, x1, y1)


def labels(ax, info, bounds, fs=8.5):
    total, bad = C.add_labels(ax, info, bounds, fs=fs, olabels=False, px=PX)
    return total, bad


def title_of(fam, info):
    return f'{fam}, PDB {info["pdb"]} ({info["spec"]["resolution"]:.2f} Å)'


def legend(fig, y, fs=7.5):
    h = [Line2D([], [], color=C.COL['H'], lw=3, label='H'), Line2D([], [], color=C.COL['Ar'], lw=3, label='Ar'),
         Line2D([], [], color=C.COL['Al'], lw=3, label='Al'),
         Line2D([], [], color='0.2', lw=0.8, ls=(0, (1.5, 1.5)), label='N/O–N/O ≤ 3.5 Å')]
    fig.legend(handles=h, loc='center', bbox_to_anchor=(0.5, y), ncol=4, frameon=False, fontsize=fs,
               handlelength=1.8, columnspacing=1.4)


def single(fam, info, png, W=3.45):
    fig = plt.figure(figsize=(W, W))
    ax = fig.add_axes((0, 0.09, 1, 0.84))
    b = draw(ax, info, png)
    # keep the axes aspect equal to the crop
    asp = (b[3] - b[1]) / (b[2] - b[0])
    H = W * 0.84 * asp / 0.84 + 0.55
    fig.set_size_inches(W, H)
    ax.set_position((0, 0.3 / H, 1, (H - 0.55) / H))
    t, bad = labels(ax, info, b)
    fig.text(0.02, 1 - 0.05 / H, title_of(fam, info), va='top', fontsize=8)
    legend(fig, 0.15 / H)
    print(f'  zoom_{fam}: label cost {t:.0f}, overlapping labels {bad}')
    fig.savefig(OUT / f'zoom_{fam}.png', dpi=300)
    fig.savefig(OUT / f'zoom_{fam}_hires.png', dpi=900)
    plt.close(fig)


def pair(infos, pngs, W=7.0):
    crops = []
    for fam in ('GH13', 'GH31'):
        info = infos[fam]
        sc = PX / (2 * info['half'])
        pts = np.array([C.tp(info, x) for rd in info['residues'] for x in rd['X'].values()] +
                       [C.tp(info, x) for u in info['units'] for x in u['X'].values()])
        x0, y0 = pts.min(0) - 2.4 * sc; x1, y1 = pts.max(0) + 2.4 * sc
        crops.append((y1 - y0) / (x1 - x0))
    pw = (W - 0.15) / 2
    ph = pw * max(crops)
    H = ph + 0.30 + 0.30
    fig = plt.figure(figsize=(W, H))
    for k, fam in enumerate(('GH13', 'GH31')):
        info = infos[fam]
        ax = fig.add_axes(((0.0 + k * (pw + 0.15)) / W, 0.30 / H, pw / W, ph / H))
        b = draw(ax, info, pngs[fam])
        # centre the crop in the cell (equal aspect)
        cx, cy = (b[0] + b[2]) / 2, (b[1] + b[3]) / 2
        w = b[2] - b[0]; h = w * ph / pw
        if h < b[3] - b[1]:
            h = b[3] - b[1]; w = h * pw / ph
        bb = (cx - w / 2, cy - h / 2, cx + w / 2, cy + h / 2)
        ax.set_xlim(bb[0], bb[2]); ax.set_ylim(bb[3], bb[1])
        t, bad = labels(ax, info, bb)
        print(f'  pair {fam}: label cost {t:.0f}, overlapping labels {bad}')
        x = (k * (pw + 0.15)) / W
        fig.text(x + 0.01, 1 - 0.04 / H, 'ab'[k], fontsize=10, fontweight='bold', va='top')
        fig.text(x + 0.20 / W + 0.01, 1 - 0.06 / H, title_of(fam, info), fontsize=8, va='top')
    legend(fig, 0.14 / H)
    fig.savefig(OUT / 'zoom_pair.png', dpi=300)
    fig.savefig(OUT / 'zoom_pair_hires.png', dpi=600)
    plt.close(fig)


if __name__ == '__main__':
    infos, pngs = {}, {}
    for fam in ('GH13', 'GH31'):
        info = setup(fam)
        info['pdb'] = PANELS[fam]['pdb']
        png = WORK / f'zoom_{fam}_render.png'
        render(info, png)
        infos[fam], pngs[fam] = info, png
        single(fam, info, png)
    pair(infos, pngs)
