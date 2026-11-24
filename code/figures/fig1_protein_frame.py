"""Data panels of the Figure 1 candidates: sugar-unit ring centroids of every density-supported unit at
subsites -3..+3 of GH13 and GH31, in the family frame.

1. Units: data/analysis/44_rule_test/units_used.tsv (the units the paper's tables use: subsites -3..+3,
   RSCC >= 0.80; GH13 308 units, GH31 82 units at the time of writing; counted at run time).
2. Ring centroid: the ring of the unit (ring_units() of code/analysis/a33_subsites.py, matched by the
   reference atom C1/C1<letter>) in the deposited mmCIF; moved into the family frame with the superposition of
   the unit's frame chain (data/calc/03_frames/superposition.tsv, x_frame = R x + t).
   Check: the reference atom moved the same way must reproduce x_frame,y_frame,z_frame of units_used.tsv
   (largest deviation printed and written to work/frame_check.txt).
3. Page: the frame coordinates are shown in the orientation of the whole-protein panel of the same family
   (6JQB for GH13, 3WEO for GH31): x_entry = R^T (x_frame - t) with the entry's own superposition, then the
   panel's basis and centre (work/whole_<pdb>.json).
4. Outline: the reference structure of the frame (1AMY chain A for GH13, 2QMJ chain A for GH31; family
   references of data/calc/03_frames/families.json) cut by a slab of +-SLAB A around the mean depth of the
   centroids; heavy atoms within the slab drawn as discs of 1.7 A on a 0.25 A grid, smoothed, contoured at 0.5.
   Reference catalytic carboxylates (Asp/Glu OD*/OE* centroid of 1AMY D179, E204; 2QMJ D443, D542) are marked.
Output: work/centroids_<family>.tsv
Run: ~/miniconda3/bin/python code/figures/fig1_protein_frame.py
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
import gemmi

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fig1_protein_common as C

FRAME = {'GH13': dict(ref='1AMY', chain='A', panel='6JQB', cat=[(179, 'ASP'), (204, 'GLU')]),
         'GH31': dict(ref='2QMJ', chain='A', panel='3WEO', cat=[(443, 'ASP'), (542, 'ASP')])}
SLAB = 5.0


def centroids():
    sup = C.superpositions()
    rows = list(csv.DictReader(open(C.UNITS44), delimiter='\t'))
    cache = {}
    out, dev = [], []
    for r in rows:
        if r['pdb'] not in cache:
            cache.clear()
            cache[r['pdb']] = C.load(r['pdb'])[0]
        m = cache[r['pdb']]
        atoms, ring = C.unit_atoms(m, r)
        R, t = sup[(r['family'], r['pdb'], r['frame_chain'])]
        cen = R @ np.mean(list(ring.values()), 0) + t
        ref = R @ ring[r['ref_atom']] + t
        dev.append(np.linalg.norm(ref - np.array([float(r['x_frame']), float(r['y_frame']), float(r['z_frame'])])))
        out.append(dict(family=r['family'], pdb=r['pdb'], chain=r['chain'], seq=r['seq'], comp=r['comp'],
                        ref_atom=r['ref_atom'], subsite=r['subsite'], rscc=r['rscc_f'],
                        xc=round(cen[0], 3), yc=round(cen[1], 3), zc=round(cen[2], 3)))
    msg = f'reference atom vs units_used.tsv frame coordinates: max deviation {max(dev):.4f} A over {len(dev)} units'
    print(msg)
    (C.WORK / 'frame_check.txt').write_text(msg + '\n')
    for fam in FRAME:
        with open(C.WORK / f'centroids_{fam}.tsv', 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter='\t')
            w.writeheader()
            w.writerows([o for o in out if o['family'] == fam])
    return out


def page_transform(fam):
    """frame -> page of the whole-protein panel"""
    f = FRAME[fam]
    meta = json.load(open(C.WORK / f'whole_{f["panel"]}.json'))
    B = np.array(meta['basis']); c = np.array(meta['center'])
    R, t = C.superpositions()[(fam, f['panel'], 'A')]
    return lambda X: ((np.asarray(X) - t) @ R - c) @ B.T      # (X - t) @ R == R^T (X - t) row-wise


def load_page(fam):
    rows = list(csv.DictReader(open(C.WORK / f'centroids_{fam}.tsv'), delimiter='\t'))
    T = page_transform(fam)
    X = T(np.array([[float(r['xc']), float(r['yc']), float(r['zc'])] for r in rows]))
    return rows, X


def reference_outline(fam, zmid, lim, grid=0.25, rad=1.7):
    """section of the reference protein at depth zmid +- SLAB, on the page; returns (xs, ys, field, cat_xy)"""
    f = FRAME[fam]
    T = page_transform(fam)
    st = C.load(f['ref'])
    pts, cat = [], []
    for ch in st[0]:
        if ch.name != f['chain']:
            continue
        for r in ch:
            if r.het_flag != 'A':
                continue
            for a in r:
                pts.append(a.pos.tolist())
            for seq, nm in f['cat']:
                if r.seqid.num == seq and r.name == nm:
                    cat.append(np.mean([a.pos.tolist() for a in r if a.name[:2] in ('OD', 'OE')], 0))
    Pp = T(np.array(pts))
    sel = np.abs(Pp[:, 2] - zmid) < SLAB
    Pp = Pp[sel]
    (x0, x1), (y0, y1) = lim
    xs = np.arange(x0, x1, grid); ys = np.arange(y0, y1, grid)
    F = np.zeros((len(ys), len(xs)))
    k = int(np.ceil(rad / grid))
    for x, y, _ in Pp:
        i, j = int((y - y0) / grid), int((x - x0) / grid)
        for di in range(-k, k + 1):
            for dj in range(-k, k + 1):
                ii, jj = i + di, j + dj
                if 0 <= ii < len(ys) and 0 <= jj < len(xs) and (di * di + dj * dj) * grid * grid <= rad * rad:
                    F[ii, jj] = 1.0
    from scipy.ndimage import gaussian_filter
    F = gaussian_filter(F, 1.5)
    return xs, ys, F, T(np.array(cat))


if __name__ == '__main__':
    centroids()
