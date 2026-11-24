"""Common definitions for Figure 2 v3 (active sites of GH13 and GH31, style of refs 01 and 03).

Data: deposited mmCIF files in data/census/cif (model 1, first altloc, no H). Residue numbering in the
labels is the reference numbering (GH13 1AMY, GH31 2QMJ); the mapping to the author numbering of each entry is
data/analysis/44_rule_test/posmap.tsv (checked below: residue name and author number must match).
Contact classes (H, Ar, Al) are class_cluster50 of data/analysis/44_rule_test/C_identity_T20.tsv.

Hydrogen bonds drawn: side-chain N/O of a shown residue within 3.5 A of an O/N of a shown sugar unit
(subsites -1, +1, +2), the closest partner per side-chain atom (the criterion of fig2_v2_struct.prepare).
"""
import csv
import math
import subprocess
from pathlib import Path

import gemmi
import numpy as np

P = Path(__file__).resolve().parents[2]
CIF = P / 'data/census/cif'
D44 = P / 'data/analysis/44_rule_test'
OUT = P / 'figures/candidates/fig2_v3'
WORK = OUT / 'work'
WORK.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------- colours (Okabe-Ito classes)
COL = {'H': '#0072B2', 'Ar': '#E69F00', 'Al': '#CC79A7'}
CAT = '#C76BB5'          # catalytic pair carbons (as requested); see CAT_USED
CAT_USED = '#B0308F'     # deeper magenta actually used: #C76BB5 is too close to the Al pink #CC79A7
SUGAR_C = '#F4B942'
CTX_C = '#C8C8C8'        # carbons of the rest of the oligosaccharide (thin)
O_COL = '#E03C31'
N_COL = '#2E6DB4'
CARTOON = '#DCDCDC'
HB = 3.5
BACKBONE = {'N', 'C', 'O', 'OXT'}
ONE3 = {'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys', 'Q': 'Gln', 'E': 'Glu', 'G': 'Gly', 'H': 'His',
        'I': 'Ile', 'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro', 'S': 'Ser', 'T': 'Thr', 'W': 'Trp',
        'Y': 'Tyr', 'V': 'Val'}
AROM = {'TYR': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'], 'PHE': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
        'TRP': ['CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'], 'HIS': ['CG', 'ND1', 'CD2', 'CE1', 'NE2']}

# ---------------------------------------------------------------- entries
# residues: (reference position, role); role 'nuc' or 'ab' for the catalytic pair
# units: (subsite, chain, seq, comp, ring atoms, atoms of the unit (None = whole residue))
# context: other residues of the same oligosaccharide drawn thin and grey
ENTRIES = {
    'GH13': dict(pdb='1G94', ref='1AMY', resolution=1.74, organism='Pseudoalteromonas haloplanktis α-amylase', short='α-amylase',
                 protein_chain='A',
                 residues=[('D179', 'nuc'), ('E204', 'ab'), ('D289', None), ('H288', None), ('H92', None),
                           ('R177', None), ('Y51', None), ('K182', None)],
                 units=[('-1', 'C', 3, 'DAF', ['C1H', 'C2H', 'C3H', 'C4H', 'C5H', 'C6H'],
                         ['C1H', 'C2H', 'O2H', 'C3H', 'O3H', 'C4H', 'O4H', 'C5H', 'C6H', 'C7H', 'O7H', 'N4']),
                        ('+1', 'C', 3, 'DAF', ['C1', 'C2', 'C3', 'C4', 'C5', 'O5'],
                         ['C1', 'C2', 'O2', 'C3', 'O3', 'C4', 'C5', 'O5', 'C6']),
                        ('+2', 'C', 2, 'GLC', ['C1', 'C2', 'C3', 'C4', 'C5', 'O5'], None)],
                 context=[('C', 1), ('B', 1), ('B', 2)]),
    'GH31': dict(pdb='3WEO', ref='2QMJ', resolution=1.45, organism='sugar beet α-glucosidase', short='α-glucosidase',
                 protein_chain='A',
                 residues=[('D443', 'nuc'), ('D542', 'ab'), ('D327', None), ('H600', None), ('R526', None),
                           ('D203', None), ('W441', None), ('F575', None)],
                 units=[('-1', 'D', 7, 'AC1', ['C1B', 'C2B', 'C3B', 'C4A', 'C5B', 'C7B'],
                         ['C1B', 'C2B', 'O2B', 'C3B', 'O3B', 'C4A', 'O4', 'C5B', 'C7B', 'C6B', 'O6B', 'N4A']),
                        ('+1', 'D', 7, 'AC1', ['C1', 'C2', 'C3', 'C4', 'C5', 'O5'],
                         ['C1', 'O2', 'C2', 'C3', 'O3', 'C4', 'C5', 'O5', 'C6']),
                        ('+2', 'D', 6, 'GLC', ['C1', 'C2', 'C3', 'C4', 'C5', 'O5'], None)],
                 context=[('D', 5), ('D', 4), ('D', 3), ('D', 2), ('D', 1)]),
}


def posmap():
    out = {}
    for r in csv.DictReader(open(D44 / 'posmap.tsv'), delimiter='\t'):
        out[(r['pdb'], r['position'])] = r
    return out


def classes():
    out = {}
    for r in csv.DictReader(open(D44 / 'C_identity_T20.tsv'), delimiter='\t'):
        out[(r['family'], r['position'])] = r
    return out


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


def plane(X):
    c = X.mean(0)
    return c, np.linalg.svd(X - c)[2][2]


def bonds(X, cut=1.9):
    names = list(X)
    out = []
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            if np.linalg.norm(X[names[i]] - X[names[j]]) < cut:
                out.append((names[i], names[j]))
    return out


def label_text(resname, pos, role):
    present = resname.capitalize()
    refname = ONE3[pos[0]]
    t = f'{present}{pos[1:]}' if present == refname else f'{present} ({refname}{pos[1:]})'
    if role == 'nuc':
        t += ' (nuc.)'
    elif role == 'ab':
        t += ' (a./b.)'
    return t


def prepare(fam):
    spec = ENTRIES[fam]
    pdb = spec['pdb']
    st = gemmi.read_structure(str(CIF / f'{pdb}.cif'))
    st.remove_hydrogens()
    m = st[0]
    pm, cl = posmap(), classes()
    residues = []
    for pos, role in spec['residues']:
        r = pm[(pdb, pos)]
        pr = find(m, r['chain'], r['seq'], r['resname'])
        atoms = [a for a in first_alt_heavy(pr) if a.name not in BACKBONE]
        cls = cl[(fam, pos)]['class_cluster50']
        residues.append(dict(pos=pos, role=role, cls=cls, chain=r['chain'], seq=int(r['seq']), resname=r['resname'],
                             label=label_text(r['resname'], pos, role),
                             X={a.name: np.array(a.pos.tolist()) for a in atoms},
                             el={a.name: a.element.name for a in atoms}))
    units = []
    for sub, ch, seq, comp, ring, keep in spec['units']:
        res = find(m, ch, seq, comp)
        atoms = [a for a in first_alt_heavy(res) if keep is None or a.name in keep]
        units.append(dict(sub=sub, chain=ch, seq=seq, comp=comp, ring=ring,
                          X={a.name: np.array(a.pos.tolist()) for a in atoms},
                          el={a.name: a.element.name for a in atoms}))
    # all atoms of the residues holding the units (drawn), and context residues
    sugar_res = sorted({(u['chain'], u['seq'], u['comp']) for u in units})
    sugarX = {}
    for ch, seq, comp in sugar_res:
        for a in first_alt_heavy(find(m, ch, seq, comp)):
            sugarX[(ch, seq, a.name)] = np.array(a.pos.tolist())
    ctxX = {}
    for ch, seq in spec['context']:
        for a in first_alt_heavy(find(m, ch, seq)):
            ctxX[(ch, seq, a.name)] = np.array(a.pos.tolist())
    # hydrogen bonds
    contacts = []
    for rd in residues:
        for an, x in rd['X'].items():
            if an == 'CA' or rd['el'][an] not in ('N', 'O'):
                continue
            for u in units:
                for bn, y in u['X'].items():
                    if u['el'][bn] in ('N', 'O'):
                        d = float(np.linalg.norm(x - y))
                        if d <= HB:
                            contacts.append(dict(res=rd, a=an, unit=u, b=bn, d=d))
    hbonds = []
    for key in sorted({(rd_i, c['a']) for rd_i, c in [(residues.index(c['res']), c) for c in contacts]}):
        cs = [c for c in contacts if residues.index(c['res']) == key[0] and c['a'] == key[1]]
        hbonds.append(min(cs, key=lambda c: c['d']))
    rings = {u['sub']: plane(np.array([u['X'][n] for n in u['ring']])) for u in units}
    # whole protein CA trace for the thumbnail extent
    ca = np.array([r['CA'][0].pos.tolist() for r in m[spec['protein_chain']] if r.find_atom('CA', '*')])
    return dict(fam=fam, spec=spec, pdb=pdb, residues=residues, units=units, sugarX=sugarX, ctxX=ctxX,
                hbonds=hbonds, rings=rings, ca=ca)


# ---------------------------------------------------------------- views
def fib_dirs(n):
    i = np.arange(n) + 0.5
    ph = np.arccos(1 - 2 * i / n)
    th = np.pi * (1 + 5 ** 0.5) * i
    return np.column_stack([np.cos(th) * np.sin(ph), np.sin(th) * np.sin(ph), np.cos(ph)])


def basis(ez, right=None):
    """Rows ex (screen right), ey (screen up), ez (toward the viewer). If `right` is given, ex is the projection
    of that vector on the screen plane."""
    ez = ez / np.linalg.norm(ez)
    if right is None:
        a = np.array([0, 1, 0]) if abs(ez[1]) < 0.9 else np.array([1, 0, 0])
        ex = np.cross(a, ez)
    else:
        ex = right - (right @ ez) * ez
    ex /= np.linalg.norm(ex)
    ey = np.cross(ez, ex)
    return np.array([ex, ey, ez])


def _crossings(A, B, C, D):
    def orient(p, q, r):
        return np.sign((q[..., 0] - p[..., 0]) * (r[..., 1] - p[..., 1]) - (q[..., 1] - p[..., 1]) * (r[..., 0] - p[..., 0]))
    o1 = orient(A, B, C); o2 = orient(A, B, D); o3 = orient(C, D, A); o4 = orient(C, D, B)
    return (o1 * o2 < 0) & (o3 * o4 < 0)


def geometry(info):
    """Segments with group ids, atoms with group ids, ring polygons, label anchors."""
    segs, gid = [], []
    for k, rd in enumerate(info['residues']):
        X = rd['X']
        for a, b in bonds(X):
            segs.append((X[a], X[b])); gid.append(100 + k)
    SX = info['sugarX']
    keys = list(SX)
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            if np.linalg.norm(SX[keys[i]] - SX[keys[j]]) < 1.9:
                segs.append((SX[keys[i]], SX[keys[j]])); gid.append(0)
    for k, h in enumerate(info['hbonds']):
        segs.append((h['res']['X'][h['a']], h['unit']['X'][h['b']])); gid.append(-10 - k)
    atoms, agid = [], []
    for k, rd in enumerate(info['residues']):
        for x in rd['X'].values():
            atoms.append(x); agid.append(100 + k)
    for x in SX.values():
        atoms.append(x); agid.append(0)
    return np.array(segs), np.array(gid), np.array(atoms), np.array(agid)


def tip(rd):
    """Label anchor of a residue: mean of the side-chain atoms beyond CB (CB for Ala)."""
    far = [x for n, x in rd['X'].items() if n not in ('CA', 'CB')]
    return np.mean(far if far else [rd['X']['CB']], 0)


def view_search(info, n=30000, top=12, dirs=None, sep=20, w_rr=0.0):
    from matplotlib.path import Path as MPath
    segs, gid, atoms, agid = geometry(info)
    iu, ju = np.triu_indices(len(segs), 1)
    share = np.zeros(len(iu), bool)
    for t in range(2):
        for s in range(2):
            share |= np.all(np.isclose(segs[iu, t], segs[ju, s]), axis=1)
    keep = (gid[iu] != gid[ju]) & ~share & ~((gid[iu] == 0) & (gid[ju] == 0))
    iu, ju = iu[keep], ju[keep]
    diffg = agid[:, None] != agid[None, :]
    rings = info['rings']
    ringX = [np.array([u['X'][a] for a in u['ring']]) for u in info['units']]
    arom = [(rd, plane(np.array([rd['X'][a] for a in AROM[rd['resname']]]))[1]) for rd in info['residues']
            if rd['resname'] in AROM and rd['resname'] != 'HIS']
    resX = np.array([x for rd in info['residues'] for n, x in rd['X'].items() if n != 'CA'] +
                    [(rd['X'][a] + rd['X'][b]) / 2 for rd in info['residues'] for a, b in bonds(rd['X'])])
    anchors = np.array([rings[s][0] for s in rings] + [tip(rd) for rd in info['residues']])
    center = atoms.mean(0)
    right = rings['+2'][0] - rings['-1'][0]
    out = []
    for ez in (fib_dirs(n) if dirs is None else dirs):
        tz = {s: abs(rings[s][1] @ ez) for s in rings}
        if not (0.45 <= tz['-1'] <= 0.97 and tz['+1'] >= 0.35 and tz['+2'] >= 0.35):
            continue
        ta = min([abs(nn @ ez) for _, nn in arom] or [1.0])
        if ta < 0.35:
            continue
        B = basis(ez, right)
        Sg = (segs - center) @ B.T
        cr = int(_crossings(Sg[iu, 0, :2], Sg[iu, 1, :2], Sg[ju, 0, :2], Sg[ju, 1, :2]).sum())
        Q = (atoms - center) @ B.T
        dd = np.linalg.norm(Q[:, None, :2] - Q[None, :, :2], axis=2)
        dd = np.where(diffg, dd, 99)
        close = int((dd < 1.3).sum() // 2)
        rr = (agid[:, None] >= 100) & (agid[None, :] >= 100)
        close_rr = int(((dd < 2.0) & rr).sum() // 2)
        Rr = (resX - center) @ B.T
        ringhit = 0
        for R in ringX:
            Rp = (R - center) @ B.T
            cen = Rp[:, :2].mean(0)
            poly = cen + (Rp[:, :2] - cen) * (1 + 0.6 / np.linalg.norm(Rp[:, :2] - cen, axis=1, keepdims=True))
            ringhit += int(MPath(poly).contains_points(Rr[:, :2]).sum())
        A = (anchors - center) @ B.T
        D = np.linalg.norm(A[:, None, :2] - A[None, :, :2], axis=2) + 1e3 * np.eye(len(A))
        ext = Q[:, :2].max(0) - Q[:, :2].min(0)
        aspect = ext[0] / ext[1]
        s = (-3 * cr - 1.0 * close - 2.5 * ringhit + 1.0 * min(D.min(), 4.0) + 0.5 * tz['-1'] + 0.4 * ta
             - 1.0 * abs(math.log(aspect / 1.15)) - w_rr * close_rr)
        out.append((s, ez, dict(cross=cr, close=close, close_rr=close_rr, ringhit=ringhit, anchor_min=round(float(D.min()), 2),
                                tilt_m1=round(float(tz['-1']), 2), arom=round(float(ta), 2),
                                aspect=round(float(aspect), 2), score=round(float(s), 2))))
    out.sort(key=lambda t: -t[0])
    picked = []
    for s, ez, rep in out:
        if all(np.degrees(np.arccos(np.clip(ez @ p[1], -1, 1))) > sep for p in picked):
            picked.append((s, ez, rep))
        if len(picked) >= top:
            break
    return picked, center, right


# ---------------------------------------------------------------- ChimeraX
def run_chimerax(lines, stem, timeout=1200):
    """Run ChimeraX commands offscreen from a Python script (lines starting with 'PY:' are raw Python).
    stdin is /dev/null so that a command error cannot leave ChimeraX waiting at its prompt."""
    py = WORK / f'{stem}_cx.py'
    body = ['from chimerax.core.commands import run']
    for ln in lines:
        if ln.startswith('PY:'):
            body.append(ln[3:].strip())
        else:
            body.append(f'run(session, {ln!r})')
    py.write_text('\n'.join(body) + '\n')
    r = subprocess.run(['chimerax', '--offscreen', '--nogui', '--exit', str(py)], capture_output=True, text=True,
                       stdin=subprocess.DEVNULL, timeout=timeout)
    txt = r.stdout + r.stderr
    if r.returncode != 0 or 'Traceback' in txt or 'ERROR' in txt:
        print(txt[-3000:])
    return py


def ortho_cmds(W):
    # 'camera ortho' fails offscreen in ChimeraX 1.11 (OffScreenRenderingContext has no 'stereo'); set it directly
    return ['PY: from chimerax.graphics import OrthographicCamera',
            'PY: _c = OrthographicCamera(); _c.position = session.main_view.camera.position',
            f'PY: _c.field_width = {W:.4f}; session.main_view.camera = _c']


def camera_cmd(B, c, dist=200.0):
    """ChimeraX camera matrix (camera-to-scene, 3x4 row-major) for screen axes B rows and look-at point c."""
    ex, ey, ez = B
    pos = c + dist * ez
    M = np.column_stack([ex, ey, ez, pos])
    return 'view matrix camera ' + ','.join(f'{v:.5f}' for v in M.ravel())


def to_px(p, B, c, W, px):
    """Pixel coordinates of scene point(s) p in an orthographic render of field width W, px x px pixels."""
    v = (np.atleast_2d(p) - c) @ B.T
    u = px / 2 + v[:, 0] / W * px
    w = px / 2 - v[:, 1] / W * px
    return np.column_stack([u, w, v[:, 2]])
