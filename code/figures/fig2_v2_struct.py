"""Structure panels of Figure 2 (v2): conserved subsite positions on deposited coordinates.

Every atom drawn comes from the deposited mmCIF (data/census/cif, model 1, first altloc, no H).
Steps for one entry:
 1. Residues: the protein residues mapped to the figure positions in this entry
    (data/figures/fig2_v2_work/contact_rows.tsv, written by fig2_v2_data.py); sugar units -1, +1 (+2)
    from data/analysis/33_subsites/units.tsv (RSCC >= 0.80).
 2. Hydrogen bonds drawn: side-chain N/O of a shown residue within 3.5 A of an O/N of a shown sugar unit
    (the criterion of code/analysis/a40_contact_types.py). Nothing else is drawn as a dashed line.
 3. View: a search over 30 000 directions scores (i) the smallest on-page distance between residue label
    anchors and ring centroids, (ii) the number of sugar atoms hidden by residue atoms in front, (iii) the
    on-page gap between the stacked aromatic ring and the sugar ring it stacks on (so the stack is seen from the
    side, not face-on), (iv) the tilt of the -1 ring (neither edge-on nor flat). The -1 -> +1 vector is then
    turned to point right.
 4. PyMOL (command line, ray-traced, orthoscopic) draws the whole chain as a transparent cartoon clipped to a
    slab around the site, the residues as sticks (side chain + CA) and the sugar units as sticks.
 5. Labels are placed by matplotlib on the exact orthoscopic projection with a greedy search over candidate
    positions that penalises overlap with atoms, dashes, other labels and crossing leader lines.
"""
import math
import csv
import subprocess
from pathlib import Path

import gemmi
import numpy as np

P = Path(__file__).resolve().parents[2]
CIF = P / 'data/census/cif'
WORK = P / 'data/figures/fig2_v2_work'
WORK.mkdir(parents=True, exist_ok=True)
ROWS = WORK / 'contact_rows.tsv'

# Okabe-Ito colour-blind-safe palette
COL = {'H': '#0072B2', 'Ar': '#E69F00', 'Al': '#CC79A7'}
SUGAR_C = '#C8C8C8'
O_COL = '#D62020'
N_COL = '#2B3BD6'
OLAB = '#7A0F0F'
HB = 3.5
BACKBONE = {'N', 'C', 'O', 'OXT'}
ONE3 = {'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys', 'Q': 'Gln', 'E': 'Glu', 'G': 'Gly', 'H': 'His',
        'I': 'Ile', 'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro', 'S': 'Ser', 'T': 'Thr', 'W': 'Trp',
        'Y': 'Tyr', 'V': 'Val'}
AROM = {'TYR': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'], 'PHE': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
        'TRP': ['CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2']}

CLASS = {'D179': 'H', 'D289': 'H', 'H288': 'H', 'R177': 'H', 'H92': 'H', 'E204': 'H', 'K182': 'H', 'Y51': 'Ar',
         'F180': 'Al', 'W206': 'Al', 'D327': 'H', 'H600': 'H', 'D443': 'H', 'D542': 'H',
         'W441': 'Ar', 'W539': 'Ar', 'F575': 'Ar', 'I328': 'Al', 'I364': 'Al'}

# Entries. view: the viewing direction (deposited frame) chosen from the 8 best-scoring distinct directions of
# the search in prepare() (montage data/figures/fig2_v2_work/views_<pdb>.png); remove it to re-run the search.
# units: (subsite, chain, seq, comp, ring atoms, {atom name: label (glucose-equivalent name)})
ENTRIES = {
    '3BC9': dict(family='GH13', ref='1AMY', resolution=1.35, view=[-0.9373, 0.33679, 0.0897], stack=('Y51', '-1'),
                 units=[('-1', 'A', 806, 'ACI', ['C1', 'C2', 'C3', 'C4', 'C5', 'C7'],
                         {'O2': 'O2', 'O3': 'O3', 'O4': 'O4', 'O6': 'O6', 'N1': 'N4'}),
                        ('+1', 'D', 3, 'GLD', ['C1', 'C2', 'C3', 'C4', 'C5', 'O5'],
                         {'O2': 'O2', 'O3': 'O3'}),
                        ('+2', 'D', 2, 'GLC', ['C1', 'C2', 'C3', 'C4', 'C5', 'O5'],
                         {'O2': 'O2', 'O3': 'O3', 'O4': 'O4', 'O6': 'O6'})],
                 context=[('C', 1, 'GLC'), ('D', 1, 'BGC')]),
    '1G94': dict(family='GH13', ref='1AMY', resolution=1.74, view=[-0.64175, 0.7157, 0.27557], stack=('Y51', '-1'),
                 units=[('-1', 'C', 3, 'DAF', ['C1H', 'C2H', 'C3H', 'C4H', 'C5H', 'C6H'],
                         {'O2H': 'O2', 'O3H': 'O3', 'O4H': 'O4', 'O7H': 'O6', 'N4': 'N4'}),
                        ('+1', 'C', 3, 'DAF', ['C1', 'C2', 'C3', 'C4', 'C5', 'O5'],
                         {'O2': 'O2', 'O3': 'O3'}),
                        ('+2', 'C', 2, 'GLC', ['C1', 'C2', 'C3', 'C4', 'C5', 'O5'],
                         {'O2': 'O2', 'O3': 'O3', 'O4': 'O4', 'O6': 'O6'})],
                 context=[('B', 1, 'GLC'), ('C', 1, 'GLC')]),
    '3WEO': dict(family='GH31', ref='2QMJ', resolution=1.45, view=[-0.99883, 0.0441, -0.0199], stack=('F575', '-1'),
                 units=[('-1', 'D', 7, 'AC1', ['C1B', 'C2B', 'C3B', 'C4A', 'C5B', 'C7B'],
                         {'O2B': 'O2', 'O3B': 'O3', 'O4': 'O4', 'O6B': 'O6', 'N4A': 'N4'}),
                        ('+1', 'D', 7, 'AC1', ['C1', 'C2', 'C3', 'C4', 'C5', 'O5'],
                         {'O2': 'O2', 'O3': 'O3'})],
                 context=[('D', 6, 'GLC')]),
}
# atoms belonging to each unit when a residue holds two rings (acarbose-type AC1 / DAF)
UNIT_ATOMS = {
    ('3WEO', '-1'): ['C1B', 'C2B', 'O2B', 'C3B', 'O3B', 'C4A', 'O4', 'C5B', 'C7B', 'C6B', 'O6B', 'N4A'],
    ('3WEO', '+1'): ['C1', 'O2', 'C2', 'C3', 'O3', 'C4', 'C5', 'O5', 'C6'],
    ('1G94', '-1'): ['C1H', 'C2H', 'O2H', 'C3H', 'O3H', 'C4H', 'O4H', 'C5H', 'C6H', 'C7H', 'O7H', 'N4'],
    ('1G94', '+1'): ['C1', 'C2', 'O2', 'C3', 'O3', 'C4', 'C5', 'O5', 'C6'],
}


def first_alt_heavy(res):
    out, seen = [], set()
    for a in res:
        if a.element.name == 'H' or a.name in seen:
            continue
        seen.add(a.name)
        out.append(a)
    return out


def find(model, chain, seq, name=None):
    for ch in model:
        if ch.name == chain:
            for r in ch:
                if r.seqid.num == int(seq) and (name is None or r.name == name):
                    return r
    raise KeyError((chain, seq, name))


def label_text(resname, ref):
    present = resname.capitalize()
    refname = ONE3[ref[0]]
    return f'{present}{ref[1:]}' if present == refname else f'{present} ({refname}{ref[1:]})'


def fib_dirs(n=30000):
    i = np.arange(n) + 0.5
    ph = np.arccos(1 - 2 * i / n)
    th = np.pi * (1 + 5 ** 0.5) * i
    return np.column_stack([np.cos(th) * np.sin(ph), np.sin(th) * np.sin(ph), np.cos(ph)])


def basis(ez):
    a = np.array([0, 1, 0]) if abs(ez[1]) < 0.9 else np.array([1, 0, 0])
    ex = np.cross(a, ez); ex /= np.linalg.norm(ex)
    ey = np.cross(ez, ex)
    return np.array([ex, ey, ez])


def plane(X):
    c = X.mean(0)
    return c, np.linalg.svd(X - c)[2][2]


def load(pdb):
    st = gemmi.read_structure(str(CIF / f'{pdb}.cif'))
    st.remove_hydrogens()
    return st


def prepare(pdb, positions=None, view_override=None, seed_score=None):
    """Collect atoms, hydrogen bonds and a viewing basis for one entry."""
    spec = ENTRIES[pdb]
    st = load(pdb)
    m = st[0]
    rows = {}
    for r in csv.DictReader(open(ROWS), delimiter='\t'):
        if r['pdb'] == pdb and (positions is None or r['position'] in positions):
            rows.setdefault(r['position'], r)
    residues = []
    for pos, r in rows.items():
        pr = find(m, r['chain'], r['seq'], r['resname'])
        atoms = [a for a in first_alt_heavy(pr) if a.name not in ('N', 'C', 'O', 'OXT')]
        residues.append(dict(pos=pos, cls=CLASS[pos], chain=r['chain'], seq=int(r['seq']), resname=r['resname'],
                             kind=r['kind'], label=label_text(r['resname'], pos),
                             X={a.name: np.array(a.pos.tolist()) for a in atoms},
                             el={a.name: a.element.name for a in atoms}))
    units = []
    for sub, ch, seq, comp, ring, names in spec['units']:
        res = find(m, ch, seq, comp)
        keep = UNIT_ATOMS.get((pdb, sub))
        atoms = [a for a in first_alt_heavy(res) if keep is None or a.name in keep]
        units.append(dict(sub=sub, chain=ch, seq=seq, comp=comp, ring=ring, names=names,
                          X={a.name: np.array(a.pos.tolist()) for a in atoms},
                          el={a.name: a.element.name for a in atoms}))
    # hydrogen bonds: residue side-chain N/O to sugar O/N <= 3.5 A
    # all pairs are kept in info['contacts']; the drawn set keeps, for each side-chain N/O, its closest sugar O/N
    contacts = []
    for rd in residues:
        for an, x in rd['X'].items():
            if an in BACKBONE or rd['el'][an] not in ('N', 'O'):
                continue
            for u in units:
                for bn, y in u['X'].items():
                    if u['el'][bn] in ('N', 'O'):
                        d = float(np.linalg.norm(x - y))
                        if d <= HB:
                            contacts.append(dict(res=rd, a=an, unit=u, b=bn, d=d))
    hbonds = []
    for key in {(id(c['res']), c['a']) for c in contacts}:
        cs = [c for c in contacts if (id(c['res']), c['a']) == key]
        hbonds.append(min(cs, key=lambda c: c['d']))
    hbonds.sort(key=lambda c: (c['res']['pos'], c['a']))
    # geometry for the view search
    rings = {u['sub']: plane(np.array([u['X'][n] for n in u['ring']])) for u in units}
    ring_pts = np.vstack([np.array([u['X'][n] for n in u['ring']]) for u in units])
    sugar_pts = np.vstack([np.array(list(u['X'].values())) for u in units])
    for rd in residues:
        sc = [x for n, x in rd['X'].items() if n != 'CA']
        rd['centroid'] = np.mean(sc if sc else list(rd['X'].values()), 0)
    anchors = np.array([rings[s][0] for s in rings] + [rd['centroid'] for rd in residues])
    resatoms = np.vstack([np.array(list(rd['X'].values())) for rd in residues])
    center = np.vstack([sugar_pts, anchors]).mean(0)
    spos, ssub = spec['stack']
    srd = next(rd for rd in residues if rd['pos'] == spos)
    sc_ar, sn_ar = plane(np.array([srd['X'][n] for n in AROM[srd['resname']]]))
    sc_su, sn_su = rings[ssub]
    c_m1, n_m1 = rings['-1']
    res_id = np.concatenate([[k] * len(rd['X']) for k, rd in enumerate(residues)])
    c_p1, n_p1 = rings['+1']
    best, bestscore = None, -1e9
    allc = []
    if view_override is None and spec.get('view') is not None:
        view_override = spec['view']
    dirs = fib_dirs() if view_override is None else [np.asarray(view_override) / np.linalg.norm(view_override)]
    for ez in dirs:
        tilt = abs(n_m1 @ ez); tilt1 = abs(n_p1 @ ez)
        if view_override is None and not (0.40 <= tilt <= 0.95 and tilt1 >= 0.35):
            continue
        B = basis(ez)
        A = (anchors - center) @ B.T
        D = np.linalg.norm(A[:, None, :2] - A[None, :, :2], axis=2) + 1e3 * np.eye(len(A))
        S = (sugar_pts - center) @ B.T
        Q = (resatoms - center) @ B.T
        dd = np.linalg.norm(Q[:, None, :2] - S[None, :, :2], axis=2)
        occl = int((dd < 1.1).any(axis=0).sum())          # sugar atoms covered by residue atoms, any depth
        rr = np.linalg.norm(Q[:, None, :2] - Q[None, :, :2], axis=2)
        cross = int(((rr < 0.9) & (res_id[:, None] != res_id[None, :])).sum() // 2)
        gap = np.linalg.norm(((sc_ar - sc_su) @ B.T)[:2])   # on-page separation of the stacked rings
        side = 1 - abs(sn_su @ ez)                            # stacked pair seen from the side
        score = (1.0 * min(D.min(), 3.0) + 0.25 * min(np.sort(D.ravel())[len(A)], 5) - 0.6 * occl - 0.45 * cross
                 + 0.8 * min(gap, 4.5) + 0.5 * side + 0.3 * tilt + 0.3 * tilt1)
        cand = dict(ez=ez, mind=D.min(), occl=occl, gap=gap, tilt=tilt, cross=cross, score=score)
        allc.append(cand)
        if score > bestscore:
            bestscore, best = score, cand
    alts = []
    for c in sorted(allc, key=lambda c: -c['score']):
        if all(np.degrees(np.arccos(np.clip(c['ez'] @ a['ez'], -1, 1))) > 30 for a in alts):
            alts.append(c)
        if len(alts) >= 8:
            break
    B = basis(best['ez'])
    # in-plane rotation: -1 -> +1 points right
    v = (rings['+1'][0] - c_m1) @ B.T
    ang = math.atan2(v[1], v[0])
    Rz = np.array([[math.cos(-ang), -math.sin(-ang), 0], [math.sin(-ang), math.cos(-ang), 0], [0, 0, 1]])
    B = Rz @ B
    print(f"{pdb}: view min anchor distance {best['mind']:.2f} A, sugar atoms covered {best['occl']}, residue overlaps {best['cross']}, "
          f"stack on-page gap {best['gap']:.2f} A, -1 ring tilt |n.z| {best['tilt']:.2f}")
    return dict(pdb=pdb, spec=spec, st=st, residues=residues, units=units, hbonds=hbonds, contacts=contacts, B=B, center=center,
                rings=rings, stack=(srd, sc_ar, sc_su), alts=alts)


def rot(info, x):
    return (np.asarray(x) - info['center']) @ info['B'].T


def write_rotated(info, path):
    st = info['st']
    out = st.clone()
    for mdl in out:
        for ch in mdl:
            for r in ch:
                for a in r:
                    v = rot(info, a.pos.tolist())
                    a.pos = gemmi.Position(*v)
    while len(out) > 1:
        del out[1]
    out.setup_entities()
    out.write_pdb(str(path))


def hexcol(h):
    h = h.lstrip('#')
    return [int(h[i:i + 2], 16) / 255 for i in (0, 2, 4)]


def render(info, png, px=2400, fov=20.0, half=None, slab=11.0, cartoon=True, only=None, sticks=0.17,
           draw_hbonds=True, context=False):
    """Ray-trace with PyMOL. Returns the half-width (A) of the square field so labels can be projected."""
    pdb = info['pdb']
    pdbfile = WORK / f'{pdb}_rot.pdb'
    write_rotated(info, pdbfile)
    residues = info['residues'] if only is None else [rd for rd in info['residues'] if rd['pos'] in only]
    units = info['units'] if only is None else [u for u in info['units'] if u['sub'] in only]
    shown = shown_points(info, residues, units)
    if half is None:
        half = np.abs(shown[:, :2]).max() + 2.5
    dist = half / math.tan(math.radians(fov / 2))
    L = ['reinitialize', f'load {pdbfile}, m', 'remove hydro', 'hide everything', 'bg_color white',
         'set ray_opaque_background, 1', 'set orthoscopic, 1', f'set field_of_view, {fov}',
         'set ray_shadows, 0', 'set antialias, 2', 'set ambient, 0.5', 'set specular, 0.2',
         'set spec_power, 60', 'set direct, 0.55', 'set reflect, 0.3', 'set depth_cue, 0',
         f'set stick_radius, {sticks}', 'set stick_quality, 24', 'set sphere_quality, 3',
         'set ray_trace_mode, 1', 'set ray_trace_gain, 0.4', 'set ray_trace_color, grey30',
         'set dash_gap, 0.22', 'set dash_length, 0.18', 'set dash_radius, 0.055', 'set dash_color, grey15',
         'set dash_round_ends, 1', 'set cartoon_transparency, 0.72', 'set cartoon_color, grey75',
         'set cartoon_loop_radius, 0.12', 'set cartoon_oval_length, 1.0', 'set cartoon_oval_width, 0.18',
         'set cartoon_rect_length, 1.0', 'set cartoon_rect_width, 0.25', 'set cartoon_fancy_helices, 1',
         f'set_color sugc, {hexcol(SUGAR_C)}', f'set_color ocol, {hexcol(O_COL)}',
         f'set_color ncol, {hexcol(N_COL)}']
    if cartoon:
        L += ['show cartoon, polymer.protein']
    for sub_ch, sub_seq, sub_comp in (info['spec'].get('context', []) if context else []):
        s = f'chain {sub_ch} and resi {sub_seq} and resn {sub_comp}'
        L += [f'show sticks, {s}', f'set stick_radius, 0.10, {s}', f'color grey60, {s}',
              f'set stick_transparency, 0.0, {s}']
    for u in units:
        names = '+'.join(u['X'])
        s = f'chain {u["chain"]} and resi {u["seq"]} and resn {u["comp"]} and name {names}'
        L += [f'show sticks, {s}', f'set stick_radius, 0.21, {s}', f'color sugc, {s} and elem C',
              f'color ocol, {s} and elem O', f'color ncol, {s} and elem N']
    for rd in residues:
        s = f'polymer.protein and chain {rd["chain"]} and resi {rd["seq"]} and not name N+C+O+OXT'
        cname = f'cls_{rd["cls"]}'
        L += [f'set_color {cname}, {hexcol(COL[rd["cls"]])}', f'show sticks, {s}',
              f'color {cname}, {s} and elem C', f'color ocol, {s} and elem O', f'color ncol, {s} and elem N+S']
    if draw_hbonds:
        k = 0
        for hb in info['hbonds']:
            if hb['res'] not in residues or hb['unit'] not in units:
                continue
            rd, u = hb['res'], hb['unit']
            L.append(f'distance hb{k}, polymer.protein and chain {rd["chain"]} and resi {rd["seq"]} and name {hb["a"]}, '
                     f'chain {u["chain"]} and resi {u["seq"]} and resn {u["comp"]} and name {hb["b"]}')
            k += 1
    L += ['hide labels', 'set label_size, 0']
    L += [f'set_view (1,0,0, 0,1,0, 0,0,1, 0,0,{-dist:.3f}, 0,0,0, {dist - slab:.3f},{dist + slab:.3f}, 1)',
          f'ray {px},{px}', f'png {png}, dpi=300']
    pml = WORK / (Path(png).stem + '.pml')
    text = '\n'.join(L) + '\n'
    if Path(png).exists() and pml.exists() and pml.read_text() == text:
        return half                      # identical scene already rendered
    pml.write_text(text)
    subprocess.run(['pymol', '-cq', str(pml)], check=True, stdout=subprocess.DEVNULL)
    return half


def shown_points(info, residues=None, units=None):
    residues = info['residues'] if residues is None else residues
    units = info['units'] if units is None else units
    pts = [rot(info, x) for rd in residues for x in rd['X'].values()]
    pts += [rot(info, x) for u in units for x in u['X'].values()]
    return np.array(pts)


def to_px(xy, half, px):
    s = px / (2 * half)
    return np.array([px / 2 + xy[0] * s, px / 2 - xy[1] * s])


# ---------------------------------------------------------------- label placement
def _overlap(a, b):
    """Overlap area of two boxes (x0, y0, x1, y1)."""
    w = min(a[2], b[2]) - max(a[0], b[0]); h = min(a[3], b[3]) - max(a[1], b[1])
    return w * h if (w > 0 and h > 0) else 0.0


def _seg_cross(p1, p2, p3, p4):
    def ccw(a, b, c):
        return (c[1] - a[1]) * (b[0] - a[0]) > (b[1] - a[1]) * (c[0] - a[0])
    return ccw(p1, p3, p4) != ccw(p2, p3, p4) and ccw(p1, p2, p3) != ccw(p1, p2, p4)


def _seg_box(p, q, box):
    x0, y0, x1, y1 = box
    edges = [((x0, y0), (x1, y0)), ((x1, y0), (x1, y1)), ((x1, y1), (x0, y1)), ((x0, y1), (x0, y0))]
    return any(_seg_cross(p, q, e[0], e[1]) for e in edges)


def _attach(box, anchor):
    """Point on the box boundary nearest the anchor (leader start)."""
    x = min(max(anchor[0], box[0]), box[2]); y = min(max(anchor[1], box[1]), box[3])
    return np.array([x, y])


def place_labels(items, obstacles, bounds, scale, segs=(), rounds=4, forbid=()):
    """items: dicts with anchor (px), w, h (px), radii (list, A), prefer (unit vector or None), weight.
    obstacles: (N, 3) array x, y, r (px).  bounds: (x0, y0, x1, y1).  scale: px per A.
    Returns the list of chosen boxes."""
    angles = np.radians(np.arange(0, 360, 12))
    cands = []
    for it in items:
        cs = []
        for a in it.get('anchors', [it['anchor']]):   # a label may attach to any of several atoms
            a = np.asarray(a, float)
            for r in it['radii']:
                for t in angles:
                    u = np.array([math.cos(t), -math.sin(t)])
                    # place box so that its nearest point is r from the anchor along u
                    off = u * (r * scale) + u * np.array([it['w'] / 2, it['h'] / 2]) * np.abs(u)
                    cs.append(np.concatenate([a + off, a]))
        cands.append(cs)

    def box_of(it, c):
        return (c[0] - it['w'] / 2, c[1] - it['h'] / 2, c[0] + it['w'] / 2, c[1] + it['h'] / 2)

    def cost(i, cc, chosen):
        it = items[i]
        c, anc = cc[:2], cc[2:]
        b = box_of(it, c)
        pad = 3
        bp = (b[0] - pad, b[1] - pad, b[2] + pad, b[3] + pad)
        x = 0.0
        if b[0] < bounds[0] or b[1] < bounds[1] or b[2] > bounds[2] or b[3] > bounds[3]:
            x += 1e6
        for fb in forbid:                      # e.g. the inset frame
            if _overlap(bp, fb) > 0:
                x += 1e6
        # atoms
        ob = obstacles
        dx = np.maximum(np.maximum(bp[0] - ob[:, 0], ob[:, 0] - bp[2]), 0)
        dy = np.maximum(np.maximum(bp[1] - ob[:, 1], ob[:, 1] - bp[3]), 0)
        inside = np.hypot(dx, dy) < ob[:, 2]
        x += 400 * inside.sum()
        # segments (H-bond dashes)
        for p, q in segs:
            if _seg_box(p, q, bp) or (bp[0] <= p[0] <= bp[2] and bp[1] <= p[1] <= bp[3]):
                x += 600
        st = _attach(b, anc)
        lead = np.linalg.norm(st - anc)
        x += 0.6 * lead / scale * 10
        if lead > 0.6 * scale:   # leader drawn over atoms of other groups
            tt = np.linspace(0, 1, 12)[2:-2]
            Lp = st[None, :] * (1 - tt[:, None]) + anc[None, :] * tt[:, None]
            dd = np.linalg.norm(ob[None, :, :2] - Lp[:, None, :], axis=2)
            hit = (dd < 0.5 * ob[None, :, 2]) & (ob[None, :, 3] != it.get('group', -99))
            x += 150 * hit.any(axis=0).sum()
        if it.get('prefer') is not None and lead > 1:
            x -= 4 * float(np.dot((c - anc) / (np.linalg.norm(c - anc) + 1e-9), it['prefer']))
        for j, cj in chosen.items():
            if j == i:
                continue
            bj = box_of(items[j], cj[:2]); aj = cj[2:]
            x += 5e4 * _overlap(bp, (bj[0] - pad, bj[1] - pad, bj[2] + pad, bj[3] + pad)) / (it['w'] * it['h'])
            if _overlap(bp, bj) > 0:
                x += 1e5
            sj = _attach(bj, aj)
            if np.linalg.norm(st - anc) > 2 and np.linalg.norm(sj - aj) > 2:
                if _seg_cross(st, anc, sj, aj):
                    x += 3000
            if np.linalg.norm(st - anc) > 2 and _seg_box(st, anc, bj):
                x += 3000
            if np.linalg.norm(sj - aj) > 2 and _seg_box(sj, aj, b):
                x += 3000
        return x * it.get('cw', 1.0)

    order = sorted(range(len(items)), key=lambda i: -items[i].get('priority', 0))
    chosen = {}
    for i in order:
        chosen[i] = min(cands[i], key=lambda c: cost(i, c, chosen))
    for _ in range(rounds):
        changed = False
        for i in order:
            best = min(cands[i], key=lambda c: cost(i, c, chosen))
            if not np.allclose(best, chosen[i]):
                chosen[i] = best; changed = True
        if not changed:
            break
    total = sum(cost(i, chosen[i], chosen) for i in chosen)
    import os
    if os.environ.get("FIG2_DEBUG"):
        for i in chosen:
            print("     cost", items[i]["text"], round(cost(i, chosen[i], chosen)))
    for i in chosen:
        items[i]['anchor'] = chosen[i][2:]
    return [box_of(items[i], chosen[i][:2]) for i in range(len(items))], total


def stack_geometry(info):
    srd, c_ar, c_su = info['stack']
    n_ar = plane(np.array([srd['X'][n] for n in AROM[srd['resname']]]))[1]
    sub = info['spec']['stack'][1]
    n_su = info['rings'][sub][1]
    d = float(np.linalg.norm(c_ar - c_su))
    ang = float(np.degrees(np.arccos(abs(n_ar @ n_su))))
    return d, ang
