"""Figure 2 of the short draft: the three contact classes at subsites -1/+1, drawn from deposited coordinates.

Panel (a) GH13, panel (b) GH31. For each family the entry is the best-resolution X-ray entry of
data/analysis/40_contact_types/contact_types.tsv (units there already have RSCC >= 0.80) that contacts
every listed reference position, wild type (RCSB title checked: no 'mutant'):
  GH13: 6SAO, 1.20 A, alpha-amylase of Thamnidium elegans, acarbose-derived pseudo-oligosaccharide (AC1 B3 spans -1/+1)
  GH31: 3WEO, 1.45 A, sugar beet alpha-glucosidase, acarviosyl-maltohexaose (AC1 D7 spans -1/+1)
Selection done by the query documented in the docstring of select_entries() below (printed when run).

Method:
 1. Read the mmCIF (model 1, first altloc, no H) with gemmi; take the AC1 residue holding the -1 and +1 units
    (units.tsv: ring 1 = valienamine cyclitol at -1, ring 0 = 4-amino-4,6-dideoxyglucose at +1) and the protein
    residues mapped to the listed reference positions in contact_types.tsv.
 2. Hydrogen bonds: side-chain N/O of those residues within 3.5 A of a ligand O (same criterion as a40).
 3. Viewing direction: search over 20 000 directions; keep directions where the -1 ring is at least 40 % and the
    +1 ring at least 35 % face-on, then maximise the smallest on-page distance between residue side-chain centroids and ring centroids, with a
    penalty for residue atoms drawn in front of (or over) ring atoms.
 4. Rotate the coordinates so that direction is +z, write a PDB of residues within 12 A (protein cartoon omitted:
    the trimmed chain drew as broken segments easily confused with H-bond dashes), ray-trace with PyMOL (orthoscopic, identity view) and place labels with matplotlib using
    the exact orthoscopic projection.
Run:  ~/miniconda3/bin/python code/figures/fig_contact_classes.py
Writes figures/fig2_contact_classes.{png,pdf}; intermediates in data/figures/fig2_work/.
"""
import csv
import math
import sqlite3
import subprocess
from pathlib import Path

import gemmi
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

P = Path(__file__).resolve().parents[2]
CIF = P / 'data/census/cif'
WORK = P / 'data/figures/fig2_work'
WORK.mkdir(parents=True, exist_ok=True)
OUTDIR = P / 'figures'
CT = P / 'data/analysis/40_contact_types/contact_types.tsv'
DB = P / 'data/db/enzimas_experimental.sqlite'

CLASS = {  # reference position -> class (short draft)
    'GH13': {'D179': 'H', 'D289': 'H', 'H288': 'H', 'R177': 'H', 'H92': 'H', 'Y51': 'Ar',
             'E204': 'H', 'G183': 'H', 'F180': 'Al'},
    'GH31': {'D327': 'H', 'H600': 'H', 'D443': 'H', 'D542': 'H', 'W441': 'Ar', 'W539': 'Ar', 'F575': 'Ar',
             'I328': 'Al', 'I364': 'Al'},
}
SUBSITE_OF = {'GH13': {'E204': '+1', 'G183': '+1', 'F180': '+1'}, 'GH31': {}}
COL = {'H': '#1f5fbf', 'Ar': '#e07b00', 'Al': '#7f7f7f'}
PYCOL = {'H': [0.122, 0.373, 0.749], 'Ar': [0.878, 0.482, 0.0], 'Al': [0.55, 0.55, 0.55]}
SUGAR_C = [0.18, 0.55, 0.34]
ONE3 = {'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys', 'Q': 'Gln', 'E': 'Glu', 'G': 'Gly', 'H': 'His',
        'I': 'Ile', 'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro', 'S': 'Ser', 'T': 'Thr', 'W': 'Trp',
        'Y': 'Tyr', 'V': 'Val'}
RING_M1 = ['C1B', 'C2B', 'C3B', 'C4A', 'C5B', 'C7B']   # valienamine cyclitol, subsite -1
RING_P1 = ['C1', 'C2', 'C3', 'C4', 'C5', 'O5']         # 4-amino-4,6-dideoxyglucose, subsite +1
BACKBONE = {'N', 'C', 'O', 'OXT'}
HB = 3.5

PANELS = [
    dict(tag='a', family='GH13', pdb='6SAO', lig=('B', 3, 'AC1'), ref='1AMY'),
    dict(tag='b', family='GH31', pdb='3WEO', lig=('D', 7, 'AC1'), ref='2QMJ'),
]
PX = 2000          # ray-trace size (square), px
FOV = 20.0


def select_entries():
    """Rank entries per family: number of listed -1/+1 positions contacted (desc), resolution (asc)."""
    con = sqlite3.connect(f'file:{DB}?mode=ro', uri=True)
    res = {r[0]: (r[1], r[2]) for r in con.execute('select pdb, resolution_A, method from entry')}
    pos = {}
    for r in csv.DictReader(open(CT), delimiter='\t'):
        if r['position'] in CLASS.get(r['family'], {}):
            pos.setdefault((r['family'], r['pdb']), set()).add(r['position'])
    for fam in CLASS:
        L = sorted((-len(v), float(res[p][0]), p) for (f, p), v in pos.items()
                   if f == fam and 'X-RAY' in res[p][1])
        print(fam, 'top entries (positions contacted, resolution A):',
              ', '.join(f'{p} ({-n}, {r})' for n, r, p in L[:5]))
    return res


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
    return f'{present}{ref[1:]}' if present == refname else f'{present}({refname}{ref[1:]})'


def fib_dirs(n=20000):
    i = np.arange(n) + 0.5
    ph = np.arccos(1 - 2 * i / n)
    th = np.pi * (1 + 5 ** 0.5) * i
    return np.column_stack([np.cos(th) * np.sin(ph), np.sin(th) * np.sin(ph), np.cos(ph)])


def basis(ez):
    a = np.array([0, 1, 0]) if abs(ez[1]) < 0.9 else np.array([1, 0, 0])
    ex = np.cross(a, ez); ex /= np.linalg.norm(ex)
    ey = np.cross(ez, ex)
    return np.array([ex, ey, ez])


def prepare(panel):
    fam, pdb = panel['family'], panel['pdb']
    st = gemmi.read_structure(str(CIF / f'{pdb}.cif'))
    st.remove_hydrogens()
    m = st[0]
    lch, lseq, lname = panel['lig']
    lig = find(m, lch, lseq, lname)
    latoms = first_alt_heavy(lig)
    rows = [r for r in csv.DictReader(open(CT), delimiter='\t')
            if r['pdb'] == pdb and r['position'] in CLASS[fam] and r['comp'] == lname]
    residues = []
    for r in rows:
        pr = find(m, r['chain'], r['seq'], r['resname'])
        atoms = first_alt_heavy(pr)
        residues.append(dict(pos=r['position'], cls=CLASS[fam][r['position']], chain=r['chain'],
                             seq=int(r['seq']), resname=r['resname'], kind=r['kind'], atoms=atoms,
                             label=label_text(r['resname'], r['position'])))
    # hydrogen bonds
    lo = [a for a in latoms if a.element.name == 'O']
    hbonds = []
    for rd in residues:
        for a in rd['atoms']:
            if a.name in BACKBONE or a.element.name not in ('N', 'O'):
                continue
            for b in lo:
                d = a.pos.dist(b.pos)
                if d <= HB:
                    hbonds.append((rd['chain'], rd['seq'], a.name, b.name, d))
    # geometry for the view search
    L = {a.name: np.array(a.pos.tolist()) for a in latoms}
    c_m1 = np.mean([L[n] for n in RING_M1], 0)
    c_p1 = np.mean([L[n] for n in RING_P1], 0)
    ring_pts = np.array([L[n] for n in RING_M1 + RING_P1])
    r1 = np.array([L[n] for n in RING_M1]); n_m1 = np.linalg.svd(r1 - r1.mean(0))[2][2]
    r0 = np.array([L[n] for n in RING_P1]); n_p1 = np.linalg.svd(r0 - r0.mean(0))[2][2]
    anchors = [c_m1, c_p1]
    resatoms = []
    for rd in residues:
        sc = [np.array(a.pos.tolist()) for a in rd['atoms'] if a.name not in BACKBONE]
        rd['centroid'] = np.mean(sc, 0)
        anchors.append(rd['centroid'])
        resatoms += sc
    anchors = np.array(anchors)
    resatoms = np.array(resatoms)
    center = np.vstack([ring_pts, anchors]).mean(0)
    best, bestscore = None, -1e9
    for ez in fib_dirs():
        if abs(n_m1 @ ez) < 0.4 or abs(n_p1 @ ez) < 0.35:
            continue
        B = basis(ez)
        A = (anchors - center) @ B.T
        D = np.linalg.norm(A[:, None, :2] - A[None, :, :2], axis=2) + 1e3 * np.eye(len(A))
        R = (ring_pts - center) @ B.T
        Q = (resatoms - center) @ B.T
        dd = np.linalg.norm(Q[:, None, :2] - R[None, :, :2], axis=2)
        front = (Q[:, None, 2] > R[None, :, 2] + 0.5) & (dd < 1.2)
        occl = int(front.any(axis=0).sum())
        overl = int((dd < 1.0).any(axis=0).sum())   # ring atoms covered by a residue atom, any depth
        score = D.min() - 0.6 * occl - 0.3 * overl + 0.5 * abs(n_m1 @ ez) + 0.3 * abs(n_p1 @ ez)
        if score > bestscore:
            bestscore, best = score, (ez, D.min(), occl, abs(n_m1 @ ez), overl, abs(n_p1 @ ez))
    ez = best[0]
    print(f'{pdb}: view min anchor distance {best[1]:.2f} A, ring atoms occluded {best[2]}, '
          f'-1 ring face-on {best[3]:.2f}, +1 ring face-on {best[5]:.2f}, ring atoms overlapped {best[4]}')
    B = basis(ez)
    # rotate the whole model and write residues within 12 A of the ligand
    keep = set()
    ns = gemmi.NeighborSearch(m, st.cell, 6).populate()
    for a in latoms:
        for mk in ns.find_atoms(a.pos, '\0', radius=12):
            cra = mk.to_cra(m)
            keep.add((cra.chain.name, cra.residue.seqid.num))
    out = gemmi.Structure()
    out.cell = gemmi.UnitCell()
    mm = gemmi.Model('1')
    for ch in m:
        nc = gemmi.Chain(ch.name)
        for r in ch:
            if (ch.name, r.seqid.num) not in keep and not (ch.name == lch and r.seqid.num == lseq):
                continue
            nr = gemmi.Residue(); nr.name = r.name; nr.seqid = r.seqid; nr.het_flag = r.het_flag
            for a in first_alt_heavy(r):
                v = (np.array(a.pos.tolist()) - center) @ B.T
                na = gemmi.Atom(); na.name = a.name; na.element = a.element; na.b_iso = a.b_iso
                na.occ = 1.0; na.pos = gemmi.Position(*v)
                nr.add_atom(na)
            nc.add_residue(nr)
        if len(nc):
            mm.add_chain(nc)
    out.add_model(mm)
    out.setup_entities()
    pdbfile = WORK / f'{pdb}_rot.pdb'
    out.write_pdb(str(pdbfile))
    rot = lambda x: (np.asarray(x) - center) @ B.T
    shown = [rot(L[n]) for n in L]
    for rd in residues:
        rd['c2'] = rot(rd['centroid'])
        rd['tip2'] = rd['c2']
        shown += [rot(a.pos.tolist()) for a in rd['atoms'] if a.name not in ('N', 'C', 'O')]
    shown = np.array(shown)
    half = np.abs(shown[:, :2]).max() + 3.2
    dist = half / math.tan(math.radians(FOV / 2))
    return dict(panel=panel, residues=residues, hbonds=hbonds, pdbfile=pdbfile, dist=dist, half=half,
                c_m1=rot(c_m1), c_p1=rot(c_p1), shown=shown, lig=panel['lig'])


def render(info, png):
    lch, lseq, lname = info['lig']
    lsel = f'chain {lch} and resi {lseq} and resn {lname}'
    lines = ['reinitialize', f'load {info["pdbfile"]}, m', 'hide everything',
             'bg_color white', 'set ray_opaque_background, 1', 'set orthoscopic, 1',
             f'set field_of_view, {FOV}', 'set ray_shadows, 0', 'set antialias, 2',
             'set ambient, 0.45', 'set specular, 0.15', 'set direct, 0.6', 'set depth_cue, 0',
             'set stick_radius, 0.17', 'set dash_gap, 0.25', 'set dash_length, 0.25', 'set dash_radius, 0.06',
             'set dash_color, black', 'set cartoon_transparency, 0.82', 'set cartoon_color, grey70',
             'set ray_trace_mode, 1', 'set ray_trace_gain, 0.6', 'set ray_trace_color, grey25',
             f'show sticks, {lsel}', f'set stick_radius, 0.22, {lsel}',
             f'set_color sugc, {SUGAR_C}', f'color sugc, {lsel} and elem C',
             f'color red, {lsel} and elem O', f'color blue, {lsel} and elem N']
    for rd in info['residues']:
        s = f'chain {rd["chain"]} and resi {rd["seq"]} and not resn {lname}'
        cname = f'cls_{rd["cls"]}'
        lines += [f'set_color {cname}, {PYCOL[rd["cls"]]}',
                  f'show sticks, {s} and not name N+C+O', f'color {cname}, {s} and elem C',
                  f'color red, {s} and elem O', f'color blue, {s} and elem N']
    for k, (ch, seq, an, bn, d) in enumerate(info['hbonds']):
        lines.append(f'distance hb{k}, chain {ch} and resi {seq} and name {an} and not resn {lname}, '
                     f'{lsel} and name {bn}')
    lines += ['hide labels', 'set dash_round_ends, 1']
    D = info['dist']
    lines += [f'set_view (1,0,0, 0,1,0, 0,0,1, 0,0,{-D:.3f}, 0,0,0, {D-25:.3f},{D+25:.3f}, 1)',
              f'ray {PX},{PX}', f'png {png}, dpi=300']
    pml = WORK / (png.stem + '.pml')
    pml.write_text('\n'.join(lines) + '\n')
    subprocess.run(['pymol', '-cq', str(pml)], check=True)


def to_px(xy, info):
    s = PX / (2 * info['half'])
    return PX / 2 + xy[0] * s, PX / 2 - xy[1] * s


def place_labels(info, ax, fs):
    """Start each label beyond its side-chain tip, away from the -1 ring, then relax overlaps."""
    s = PX / (2 * info['half'])
    res = info['residues']
    atoms_px = np.array([to_px(p[:2], info) for p in info['shown']])
    c0 = np.array(to_px(info['c_m1'][:2], info)) * 0.5 + np.array(to_px(info['c_p1'][:2], info)) * 0.5
    anc = np.array([to_px(rd['tip2'][:2], info) for rd in res])
    pos = []
    for a in anc:
        v = a - c0; v /= np.linalg.norm(v) + 1e-9
        pos.append(a + v * 1.6 * s)
    pos = np.array(pos, float)
    w = np.array([len(rd['label']) * fs * 0.55 * PX / 300 * 1.25 for rd in res])  # approx width px at render scale
    h = fs * 1.6 * PX / 300
    for _ in range(600):
        moved = False
        for i in range(len(pos)):
            f = np.zeros(2)
            for j in range(len(pos)):
                if i == j:
                    continue
                dx, dy = pos[i] - pos[j]
                ox = (w[i] + w[j]) / 2 + 12 - abs(dx); oy = h + 8 - abs(dy)
                if ox > 0 and oy > 0:
                    if ox < oy * 1.5:
                        f[0] += np.sign(dx or 1) * ox * 0.5
                    else:
                        f[1] += np.sign(dy or 1) * oy * 0.5
            # keep off atoms
            d = atoms_px - pos[i]
            near = (np.abs(d[:, 0]) < w[i] / 2 + 10) & (np.abs(d[:, 1]) < h / 2 + 10)
            if near.any():
                f -= d[near].mean(0) * 0.15
            pos[i] += f
            pos[i][0] = np.clip(pos[i][0], w[i] / 2 + 10, PX - w[i] / 2 - 10)
            pos[i][1] = np.clip(pos[i][1], h + 10, PX - h - 10)
            moved |= np.abs(f).max() > 0.5
        if not moved:
            break
    for rd, a, p in zip(res, anc, pos):
        if np.linalg.norm(p - a) > 0.9 * s:
            ax.plot([a[0], p[0]], [a[1], p[1]], color=COL[rd['cls']], lw=0.5, alpha=0.8, zorder=2)
        ax.text(p[0], p[1], rd['label'], color=COL[rd['cls']], fontsize=fs, ha='center', va='center',
                fontweight='bold', zorder=3,
                bbox=dict(boxstyle='round,pad=0.12', fc='white', ec='none', alpha=0.85))
    for c, t in ((info['c_m1'], '−1'), (info['c_p1'], '+1')):
        x, y = to_px(c[:2], info)
        ax.text(x, y, t, fontsize=fs, ha='center', va='center', color='black', fontweight='bold', zorder=3,
                bbox=dict(boxstyle='round,pad=0.1', fc='white', ec='none', alpha=0.75))


def main():
    res = select_entries()
    infos = []
    for p in PANELS:
        info = prepare(p)
        png = WORK / f'{p["pdb"]}_ray.png'
        render(info, png)
        info['png'] = png
        info['resol'] = float(res[p['pdb']][0])
        infos.append(info)
        print(p['pdb'], 'hydrogen bonds:', ', '.join(f'{ch}{seq}:{a}-{b} {d:.2f}' for ch, seq, a, b, d in info['hbonds']))
        for rd in info['residues']:
            print(f"   {rd['label']:<14} {rd['resname']}{rd['seq']} class {rd['cls']} contact {rd['kind']}")
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.9))
    fs = 6.5
    for ax, info in zip(axes, infos):
        img = plt.imread(info['png'])
        ax.imshow(img, extent=(0, PX, PX, 0), interpolation='lanczos')
        place_labels(info, ax, fs)
        p = info['panel']
        ax.set_title(f"({p['tag']}) {p['family']}, PDB {p['pdb']} ({info['resol']:.2f} Å)", fontsize=8, loc='left')
        ax.set_xlim(0, PX); ax.set_ylim(PX, 0); ax.axis('off')
    handles = [Line2D([], [], color=COL['H'], lw=3, label='H: side chain H-bonds the sugar'),
               Line2D([], [], color=COL['Ar'], lw=3, label='Ar: aromatic ring stacks'),
               Line2D([], [], color=COL['Al'], lw=3, label='Al: aliphatic packing'),
               Line2D([], [], color='black', lw=1, ls=(0, (2, 2)), label='H-bond ≤ 3.5 Å')]
    fig.legend(handles=handles, loc='lower center', ncol=4, fontsize=6.5, frameon=False, handlelength=1.8)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.93, bottom=0.07, wspace=0.03)
    OUTDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTDIR / 'fig2_contact_classes.png', dpi=300, facecolor='white')
    fig.savefig(OUTDIR / 'fig2_contact_classes.pdf', dpi=300, facecolor='white')
    print('wrote', OUTDIR / 'fig2_contact_classes.png')


if __name__ == '__main__':
    main()
