"""Figure 2 v3, panels a (GH13, PDB 1G94) and b (GH31, PDB 3WEO): active-site zooms in the style of the
published examples 01 (Nat Commun 2017 Fig. 3) and 03 (JACS 2022 Fig. 3), Recipe A of
figures/references/REFERENCES.md.

Data: deposited coordinates only (data/census/cif/{1G94,3WEO}.cif; model 1; ChimeraX reads the mmCIF
directly, gemmi reads the same file for the geometry used by the view search and the label placement).
Residues and sugar units: fig2_v3_common.ENTRIES (reference numbering via 44_rule_test/posmap.tsv).

Steps per panel
 1. View: direction chosen by eye among the 12 best of the search in fig2_v3_views.py (VIEW below; both are
    the solvent-side member of an antipodal pair); the -1 -> +2 vector is put along the page x axis.
 2. ChimeraX 1.11 offscreen, orthographic camera (set from Python: 'camera ortho' fails offscreen) on the same
    axes as the projection used for the labels, two layers:
    (i) cartoon only, #DCDCDC, soft lighting, depth cue, thin grey silhouettes, clipped from 1.5 A in front of
        the nearest shown atom to 10 A behind the panel centre; faded to 30 % strength in the figure;
    (ii) side chains (CA onwards; CB onwards for His288 of GH13, whose CA would be joined to that of
        Asp289 by a ChimeraX chain-trace bond) as sticks r = 0.16 A, carbons by contact class (H #0072B2, Ar #E69F00,
        Al #CC79A7), catalytic pair carbons magenta #B0308F; sugar units -1, +1, +2 as sticks r = 0.22 A,
        gold carbons #F4B942; O #E03C31, N #2E6DB4; lighting soft (ambient occlusion, multiShadow MS),
        black silhouettes; transparent background. 2000 x 2000 px, supersample SS.
    H-bonds (fig2_v3_common: side-chain N/O to sugar O/N <= 3.5 A, closest partner per atom) are drawn by
    matplotlib as thin dark-grey dashes between the two layers.
 3. Thumbnail: the whole chain from the same direction (surface #D5D9DE, soft lighting, silhouettes), the
    oligosaccharide as gold spheres; a black box marks the zoomed field; placed where it covers no atoms.
 4. Labels drawn by matplotlib on the exact orthographic projection: greedy placement over candidate positions
    around every side-chain atom, cost = drawn atoms/bonds under the label (sugar weighted most), overlap with
    other labels, thumbnail or panel edge, distance, and leader lines crossing bonds or labels; leader lines for
    labels more than 1 A from their atom. Subsite numbers placed on or beside their rings.

Writes figures/candidates/fig2_v3/fig2_v3_a.{png,pdf}, fig2_v3_b.{png,pdf} (3.42 in wide,
300 dpi, with legend) and the renders in fig2_v3/work/.
Run: ~/miniconda3/bin/python code/figures/fig2_v3_struct.py [--rerender]
"""
import sys
import json
import hashlib
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Rectangle

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fig2_v3_common as C

FONT = 'Liberation Sans'
plt.rcParams.update({'font.family': FONT, 'font.size': 7, 'pdf.fonttype': 42, 'axes.linewidth': 0.6})
PX = 2000
TPX = 900
MS = 256         # multiShadow directions (64 for drafts)
SS = 3
SIL = 3.0
SIL_CART = 1.5
CLIP_FAR = 10.0
CART_ALPHA = 0.70
PANEL_W = 3.42          # inches
ASPECT = 0.90           # panel height / width
MARGIN = 3.2            # A around the shown atoms
FS_RES = 7.5
FS_SUB = 8.0
LEADER_MIN = 1.0      # A; labels farther than this from their atom get a leader line
TEXT_COL = {'#E69F00': '#B87800', '#CC79A7': '#B5508E'}   # darker text for orange / pink on white

# viewing directions picked by eye from work/views_<fam>.png (fig2_v3_views.py)
VIEW = {
    'GH13': [-0.2821, 0.6736, -0.6832],      # candidate 0 (seen from the solvent side)
    'GH31': [-0.7836, 0.5235, 0.3346],       # candidate 1 (seen from the solvent side)
}
THUMB_POS = {'GH13': None, 'GH31': None}     # None = chosen by least clutter
THUMB = 0.27            # thumbnail size, fraction of the panel width
HB_COL = '#4A4A4A'


def cx_spec_res(info):
    """One ChimeraX spec '/A:n1,n2,...' (a space-separated list would bind a following '&' to its last item)."""
    ch = {rd['chain'] for rd in info['residues']}
    assert len(ch) == 1
    return f'/{ch.pop()}:' + ','.join(str(rd['seq']) for rd in info['residues'])


def setup(fam):
    info = C.prepare(fam)
    rings = info['rings']
    right = rings['+2'][0] - rings['-1'][0]
    B = C.basis(np.array(VIEW[fam], float), right)
    pts = np.array([x for rd in info['residues'] for n, x in rd['X'].items()] + list(info['sugarX'].values()))
    v = pts @ B.T
    lo, hi = v.min(0), v.max(0)
    mid = (lo + hi) / 2
    ext = hi - lo
    W = max(ext[0] + 2 * MARGIN, (ext[1] + 2 * MARGIN) / ASPECT)
    c = mid @ B                      # scene point at the page centre (depth = middle of the shown atoms)
    zmax = float((pts - c) @ B[2] @ np.ones(1) if False else ((pts - c) @ B[2]).max())
    zmin = float(((pts - c) @ B[2]).min())
    info.update(B=B, c=c, W=W, zmax=zmax, zmin=zmin)
    # where is the protein relative to the site (positive = between viewer and site)
    pz = (info['ca'] - c) @ B[2]
    print(f'{fam}: field {W:.1f} A, depth of shown atoms {zmin:.1f}..{zmax:.1f} A; protein CA in front of the '
          f'site: {np.mean(pz > zmax):.2f} of CA')
    return info


def write_pb(info, path):
    L = ['; radius = 0.045', '; dashes = 10', '; color = #555555']
    for h in info['hbonds']:
        rd, u = h['res'], h['unit']
        L.append(f'/{rd["chain"]}:{rd["seq"]}@{h["a"]} /{u["chain"]}:{u["seq"]}@{h["b"]}')
    path.write_text('\n'.join(L) + '\n')


def render_zoom(info, png_cart, png_stk):
    """Two layers from the same camera: (1) the cartoon alone, opaque pale grey, soft light, no outlines;
    (2) the sticks and H-bond dashes, soft light (ambient occlusion among the sticks), black silhouettes,
    transparent background. The figure lays (2) over a faded (1)."""
    spec = info['spec']
    pb = C.WORK / f'{info["fam"]}_hb.pb'
    write_pb(info, pb)
    uch = {u['chain'] for u in info['units']}
    assert len(uch) == 1
    sug = f'/{uch.pop()}:' + ','.join(sorted({str(u['seq']) for u in info['units']}))
    pc = spec['protein_chain']
    cam = [f'windowsize {PX} {PX}', *C.ortho_cmds(info['W']), C.camera_cmd(info['B'], info['c'])]
    clip = (f'clip near {-(info["zmax"] + 1.5):.2f} far {CLIP_FAR:.2f} position '
            f'{",".join(f"{v:.3f}" for v in info["c"])}')
    head = [f'open {C.CIF / (info["pdb"] + ".cif")}', 'delete solvent', 'delete H', 'set bgColor white',
            'graphics quality 6', 'hide atoms', 'hide cartoons']
    L1 = head + [f'cartoon /{pc}', f'cartoon style /{pc} width 1.6 thickness 0.35',
                 f'color /{pc} {C.CARTOON} target c',
                 f'lighting soft multiShadow {MS}', 'lighting depthCue true', 'lighting depthCueStart 0.2',
                 'lighting depthCueEnd 1.0', 'material dull',
                 f'graphics silhouettes true width {SIL_CART} color #8C8C8C depthJump 0.03', *cam, clip,
                 f'save {png_cart} width {PX} height {PX} supersample {SS} transparentBackground false']
    C.run_chimerax(L1, Path(png_cart).stem)
    L2 = head + [f'name site {cx_spec_res(info)}', f'name lig {sug}',
                 'show site atoms', 'hide site & @N,C,O,OXT', 'hide pbonds', 'show lig atoms', 'style site stick',
                 'style lig stick',
                 'size site stickRadius 0.16', 'size lig stickRadius 0.22']
    # ChimeraX draws a chain-trace bond between the CA atoms of consecutive residues whose backbone is hidden;
    # for such pairs (GH13 His288-Asp289) the CA of the first residue is hidden too.
    seqs = {rd['seq'] for rd in info['residues']}
    for rd in info['residues']:
        if rd['seq'] + 1 in seqs:
            L2.append(f'hide /{rd["chain"]}:{rd["seq"]}@CA')
    for rd in info['residues']:
        col = C.CAT_USED if rd['role'] else C.COL[rd['cls']]
        L2.append(f'color /{rd["chain"]}:{rd["seq"]} & C {col} target a')
    L2 += [f'color lig & C {C.SUGAR_C} target a',
           f'color (site | lig) & O {C.O_COL} target a', f'color (site | lig) & N {C.N_COL} target a',
           f'color site & S #E6C229 target a', 'hide pbonds',
           f'lighting soft multiShadow {MS}', 'lighting depthCue false', 'material dull',
           f'graphics silhouettes true width {SIL} color black depthJump 0.01', *cam, 'clip off',
           f'save {png_stk} width {PX} height {PX} supersample {SS} transparentBackground true']
    C.run_chimerax(L2, Path(png_stk).stem)


def thumb_geometry(info):
    B = info['B']
    v = info['ca'] @ B.T
    lo, hi = v[:, :2].min(0), v[:, :2].max(0)
    mid = np.r_[(lo + hi) / 2, (info['c'] @ B.T)[2]]
    Wt = max(hi - lo) + 9.0
    return mid @ B, Wt


def render_thumb(info, png):
    spec = info['spec']
    ct, Wt = thumb_geometry(info)
    pc = spec['protein_chain']
    sug_chains = sorted({u['chain'] for u in info['units']} | {ch for ch, _ in spec['context']})
    sc = ','.join(sug_chains)
    L = [f'open {C.CIF / (info["pdb"] + ".cif")}', 'delete solvent', 'delete H', 'set bgColor white',
         'hide atoms', 'hide cartoons', f'surface /{pc} enclose /{pc}', f'color /{pc} #D5D9DE target s',
         f'show /{sc} atoms', f'style /{sc} sphere', f'color /{sc} {C.SUGAR_C} target a',
         f'color /{sc} & O #E8833A target a',
         'lighting soft multiShadow 128', 'lighting depthCue false', 'material dull',
         'graphics silhouettes true width 3 color black depthJump 0.01',
         f'windowsize {TPX} {TPX}', *C.ortho_cmds(Wt), C.camera_cmd(info['B'], ct),
         'clip off',
         f'save {png} width {TPX} height {TPX} supersample {SS} transparentBackground true']
    C.run_chimerax(L, Path(png).stem)


def renders(info, rerender=False):
    key = json.dumps(dict(v=VIEW[info['fam']], W=round(info['W'], 3), m=MARGIN,
                          res=[r['pos'] for r in info['residues']], ver=8, sil=SIL, cf=CLIP_FAR, ms=MS, ss=SS), sort_keys=True)
    h = hashlib.md5(key.encode()).hexdigest()[:8]
    zp = (C.WORK / f'zoomcart_{info["fam"]}_{h}.png', C.WORK / f'zoomstk_{info["fam"]}_{h}.png')
    tp = C.WORK / f'thumb_{info["fam"]}_{h}.png'
    if rerender or not zp[1].exists():
        render_zoom(info, *zp)
    if rerender or not tp.exists():
        render_thumb(info, tp)
    return zp, tp


# ---------------------------------------------------------------- labels
def text_size_pt(txt, fs, weight='normal'):
    fp = FontProperties(family=FONT, weight=weight, size=fs)
    e = TextPath((0, 0), txt, size=fs, prop=fp).get_extents()
    return e.width, fs * 0.95


def obstacles(info, px):
    """Points (pixel coords) with weights: atoms and bond samples of everything drawn."""
    B, c, W = info['B'], info['c'], info['W']
    P, Wt, R = [], [], []

    def add(X, w, rad, bond_cut=1.9):
        names = list(X)
        arr = np.array([X[n] for n in names])
        for x in arr:
            P.append(x); Wt.append(w); R.append(rad)
        for i in range(len(arr)):
            for j in range(i + 1, len(arr)):
                if np.linalg.norm(arr[i] - arr[j]) < bond_cut:
                    for t in (0.25, 0.5, 0.75):
                        P.append(arr[i] * (1 - t) + arr[j] * t); Wt.append(w); R.append(rad)
    for rd in info['residues']:
        add(rd['X'], 3.0, 0.45)
    add({k: v for k, v in info['sugarX'].items()}, 6.0, 0.55)
    for h in info['hbonds']:
        a, b = h['res']['X'][h['a']], h['unit']['X'][h['b']]
        for t in np.linspace(0.2, 0.8, 4):
            P.append(a * (1 - t) + b * t); Wt.append(1.5); R.append(0.25)
    Q = C.to_px(np.array(P), B, c, W, px)
    scale = px / W
    return Q[:, :2], np.array(Wt), np.array(R) * scale


def place_labels(info, ax, crop, inset_box, px):
    """Greedy placement of residue and subsite labels.

    Residue labels: candidate boxes around every side-chain atom of the residue (CB onwards) at gaps of
    0.2-3.0 A in 24 directions. Subsite numbers: at the ring centre or beside it. Cost of a box: drawn atoms and
    bond samples under it (sugar x6, residues x3, H-bond dashes x1.5), any overlap with another label (1e5),
    the thumbnail or the panel edge (1e6), 12 per A of gap, and drawn bonds crossed by the leader line."""
    B, c, W = info['B'], info['c'], info['W']
    scale = px / W                                   # px per A
    x0, y0, x1, y1 = crop
    ppi = (x1 - x0) / PANEL_W                        # image px per inch on the page
    pt = ppi / 72.0                                  # image px per point
    OP, OW, OR = obstacles(info, px)
    items = []
    for rd in info['residues']:
        col = C.CAT_USED if rd['role'] else C.COL[rd['cls']]
        col = TEXT_COL.get(col, col)
        names = [n for n in rd['X'] if n != 'CA']
        anc = C.to_px(np.array([rd['X'][n] for n in names]), B, c, W, px)[:, :2]
        items.append(dict(kind='res', txt=rd['label'], ancs=anc, col=col, fs=FS_RES, weight='normal', own=None))
    for u in info['units']:
        cen = C.to_px(info['rings'][u['sub']][0], B, c, W, px)[0, :2]
        ringpx = C.to_px(np.array([u['X'][a] for a in u['ring']]), B, c, W, px)[:, :2]
        items.append(dict(kind='sub', txt=u['sub'].replace('-', '−'), ancs=cen[None], col='black', fs=FS_SUB,
                          weight='bold', own=ringpx))
    for it in items:
        w, h = text_size_pt(it['txt'], it['fs'], it['weight'])
        it['w'], it['h'] = w * pt + 4, h * pt + 2
    angles = np.radians(np.arange(0, 360, 15))

    def cands(it):
        out = []
        rs = [0.0, 0.4, 0.8, 1.3] if it['kind'] == 'sub' else [0.25, 0.6, 1.0, 1.6, 2.3, 3.0]
        for anc in it['ancs']:
            for r in rs:
                for a in (angles if r > 0 else [0.0]):
                    d = np.array([np.cos(a), -np.sin(a)])
                    if r > 0:
                        # box edge at distance r from the anchor along d
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

    def leader_end(b, anc):
        return np.array([np.clip(anc[0], b[0], b[2]), np.clip(anc[1], b[1], b[3])])

    def cost(it, ctr, r, anc, placed):
        b = box_of(ctr, it)
        pad = 3
        if b[0] < x0 + pad or b[2] > x1 - pad or b[1] < y0 + pad or b[3] > y1 - pad:
            return 1e6
        if inset_box is not None and overlap(b, (inset_box[0] - 4, inset_box[1] - 4, inset_box[2] + 4,
                                                 inset_box[3] + 4)) > 0:
            return 1e6
        inside = ((OP[:, 0] + OR > b[0]) & (OP[:, 0] - OR < b[2]) & (OP[:, 1] + OR > b[1]) & (OP[:, 1] - OR < b[3]))
        s = 10 * OW[inside].sum()
        for ob in placed:
            if overlap(b, ob) > 0:
                s += 1e5
            s += 0.02 * overlap((b[0] - 8, b[1] - 8, b[2] + 8, b[3] + 8), ob)
        s += 12.0 * r * (1 if it['kind'] == 'res' else 2)
        if it['kind'] == 'res' and r >= LEADER_MIN:
            e = leader_end(b, anc)
            v = e - anc
            L = np.linalg.norm(v)
            if L > 1:
                t = np.clip(((OP - anc) @ v) / L ** 2, 0, 1)
                dist = np.linalg.norm(anc + t[:, None] * v - OP, axis=1)
                hit = (dist < OR * 0.8) & (t > 0.15) & (t < 0.98)
                s += 6 * OW[hit].sum()
                for ob in placed:                     # leader through another label
                    tt = np.linspace(0, 1, 12)[:, None]
                    pts = anc + tt * v
                    if np.any((pts[:, 0] > ob[0]) & (pts[:, 0] < ob[2]) & (pts[:, 1] > ob[1]) & (pts[:, 1] < ob[3])):
                        s += 1e4
        return s

    sol = {}
    order = sorted(range(len(items)), key=lambda k: (items[k]['kind'] != 'sub', k))
    for _pass in range(5):
        for k in order:
            placed = [box_of(sol[j][0], items[j]) for j in sol if j != k]
            sol[k] = min(cands(items[k]), key=lambda cr: cost(items[k], cr[0], cr[1], cr[2], placed))
    res, total = [], 0
    for k, it in enumerate(items):
        ctr, r, anc = sol[k]
        placed = [box_of(sol[j][0], items[j]) for j in sol if j != k]
        cs = cost(it, ctr, r, anc, placed)
        total += cs
        res.append(dict(it, ctr=ctr, r=r, anc=anc, cost=cs, box=box_of(ctr, it)))
    return res, total


def thumb_placement(info, crop, tp, px):
    """Crop the thumbnail render to its content and choose the panel position (4 corners, 4 edge middles)
    that covers the fewest drawn atoms. Returns the image, its page box and the mapping scene -> page."""
    timg = plt.imread(tp)
    a = timg[..., 3] > 0.02
    rows, cols = np.where(a.any(1))[0], np.where(a.any(0))[0]
    v0, v1, u0, u1 = rows[0], rows[-1] + 1, cols[0], cols[-1] + 1
    sub = timg[v0:v1, u0:u1]
    x0, y0, x1, y1 = crop
    f = THUMB * (x1 - x0) / max(u1 - u0, v1 - v0)
    tw, th = (u1 - u0) * f, (v1 - v0) * f
    OP, OW, OR = obstacles(info, px)
    m = 0.01 * (x1 - x0)
    pos = {'tl': (x0 + m, y0 + m), 'tr': (x1 - tw - m, y0 + m), 'bl': (x0 + m, y1 - th - m),
           'br': (x1 - tw - m, y1 - th - m), 'tc': ((x0 + x1 - tw) / 2, y0 + m), 'bc': ((x0 + x1 - tw) / 2, y1 - th - m),
           'lc': (x0 + m, (y0 + y1 - th) / 2), 'rc': (x1 - tw - m, (y0 + y1 - th) / 2)}
    if THUMB_POS.get(info['fam']):
        pos = {THUMB_POS[info['fam']]: pos[THUMB_POS[info['fam']]]}
    best = None
    for name, (bx, by) in pos.items():
        b = (bx, by, bx + tw, by + th)
        inside = ((OP[:, 0] + OR > b[0]) & (OP[:, 0] - OR < b[2]) & (OP[:, 1] + OR > b[1]) & (OP[:, 1] - OR < b[3]))
        s = OW[inside].sum() + (0 if name in ('tl', 'tr', 'bl', 'br') else 5)
        if best is None or s < best[0]:
            best = (s, name, b)
    s, name, b = best
    ct, Wt = thumb_geometry(info)

    def scene_to_page(p):
        q = C.to_px(p, info['B'], ct, Wt, TPX)[:, :2]
        return np.column_stack([b[0] + (q[:, 0] - u0) * f, b[1] + (q[:, 1] - v0) * f])
    return sub, b, name, s, scene_to_page


def draw_panel(info, zp, tp, ax):
    px = PX
    B, c, W = info['B'], info['c'], info['W']
    cart = plt.imread(zp[0])[..., :3]
    stk = plt.imread(zp[1])
    cart = 1 - CART_ALPHA * (1 - cart)              # fade the cartoon towards white
    hpx = px * ASPECT
    crop = (0, (px - hpx) / 2, px, (px + hpx) / 2)
    ax.imshow(cart, extent=(0, px, px, 0), interpolation='lanczos', zorder=0)
    # hydrogen bonds: thin dashes under the sticks
    for h in info['hbonds']:
        q = C.to_px(np.array([h['res']['X'][h['a']], h['unit']['X'][h['b']]]), B, c, W, px)
        ax.plot(q[:, 0], q[:, 1], color=HB_COL, lw=0.75, ls=(0, (2.2, 1.6)), zorder=0.5, dash_capstyle='round')
    ax.imshow(stk, extent=(0, px, px, 0), interpolation='lanczos', zorder=1)
    ax.set_xlim(crop[0], crop[2]); ax.set_ylim(crop[3], crop[1]); ax.axis('off')
    # thumbnail with the zoomed field boxed
    sub, ib, corner, tcost, s2p = thumb_placement(info, crop, tp, px)
    ax.imshow(sub, extent=(ib[0], ib[2], ib[3], ib[1]), interpolation='lanczos', zorder=3)
    corners_px = np.array([[crop[0], crop[1]], [crop[2], crop[1]], [crop[2], crop[3]], [crop[0], crop[3]]])
    sc = np.column_stack([(corners_px[:, 0] - px / 2) * W / px, -(corners_px[:, 1] - px / 2) * W / px,
                          np.zeros(4)]) @ B + c
    q = s2p(sc)
    ax.add_patch(plt.Polygon(q, closed=True, fill=False, lw=0.9, ec='black', zorder=5))
    # labels
    labs, total = place_labels(info, ax, crop, ib, px)
    halo = [pe.withStroke(linewidth=2.0, foreground='white')]
    for L in labs:
        ax.text(L['ctr'][0], L['ctr'][1], L['txt'], color=L['col'], fontsize=L['fs'], weight=L['weight'],
                ha='center', va='center', zorder=6, path_effects=halo)
        if L['kind'] == 'res' and L['r'] >= LEADER_MIN:
            b = L['box']
            tx = np.clip(L['anc'][0], b[0], b[2]); ty = np.clip(L['anc'][1], b[1], b[3])
            ax.plot([L['anc'][0], tx], [L['anc'][1], ty], color=L['col'], lw=0.5, zorder=5,
                    solid_capstyle='round')
    bad = [L['txt'] for L in labs if L['cost'] > 1e4]
    print(f'  {info["fam"]}: thumbnail {corner} (cost {tcost:.0f}), label cost {total:.0f}, bad labels {bad}')
    for L in labs:
        print(f'    {L["txt"]:<22s} r={L["r"]:.1f} cost={L["cost"]:.0f}')
    return crop


def scale_bar(ax, info, crop):
    pass


LEG_H = 0.30            # inches for the legend row
TITLE_H = 0.22


def panel_into(sf, info, zp, tp, letter, Hin):
    """Draw one structure panel into a (sub)figure of height Hin inches (title row + image)."""
    ax = sf.add_axes((0, 0, 1, (Hin - TITLE_H) / Hin))
    draw_panel(info, zp, tp, ax)
    sf.text(0.005, 1 - 0.02 / Hin, letter, fontsize=9, weight='bold', va='top')
    sf.text(0.065, 1 - 0.04 / Hin, title_of(info), fontsize=7.5, va='top')


def struct_legend(sf, x, y, ncol=5, fs=7):
    from matplotlib.lines import Line2D
    h = [Line2D([], [], color=C.SUGAR_C, lw=3.2, label='sugar, subsites −1 to +2'),
         Line2D([], [], color=C.COL['H'], lw=2.2, label='H'),
         Line2D([], [], color=C.COL['Ar'], lw=2.2, label='Ar'),
         Line2D([], [], color=C.CAT_USED, lw=2.2, label='catalytic pair'),
         Line2D([], [], color=HB_COL, lw=0.75, ls=(0, (2.2, 1.6)), label='H-bond ≤ 3.5 Å')]
    sf.legend(handles=h, loc='center', bbox_to_anchor=(x, y), ncol=ncol, frameon=False, fontsize=fs,
              handlelength=1.8, columnspacing=1.1, handletextpad=0.5)


def panel_figure(info, zp, tp, stem, letter):
    Hp = PANEL_W * ASPECT + TITLE_H
    H = Hp + LEG_H + 0.12
    fig = plt.figure(figsize=(PANEL_W, H))
    sf_top, sf_leg = fig.subfigures(2, 1, height_ratios=[Hp, H - Hp])
    panel_into(sf_top, info, zp, tp, letter, Hp)
    struct_legend(sf_leg, 0.5, 0.5, ncol=3, fs=6.5)
    fig.savefig(C.OUT / f'{stem}.png', dpi=300)
    fig.savefig(C.OUT / f'{stem}.pdf', dpi=300)
    plt.close(fig)


def title_of(info):
    s = info['spec']
    return f'{info["fam"]} {s["short"]}, PDB {info["pdb"]} ({s["resolution"]:.2f} Å)'


def main():
    rer = '--rerender' in sys.argv
    out = {}
    for fam, letter in (('GH13', 'a'), ('GH31', 'b')):
        info = setup(fam)
        zp, tp = renders(info, rer)
        panel_figure(info, zp, tp, f'fig2_v3_{letter}', letter)
        out[fam] = (info, zp, tp)
    return out


if __name__ == '__main__':
    main()
