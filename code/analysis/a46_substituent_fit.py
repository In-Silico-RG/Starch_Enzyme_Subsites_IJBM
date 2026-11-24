"""46 - Steric fit of common starch substituents at each glucose oxygen (O2, O3, O6) and subsite of GH13, GH31,
GH14 and GH15, on the deposited coordinates of the density-supported units of 45.

Question: for a hydroxyl that 45 calls held or free, would a substituent (methyl, hydroxyethyl, 2-hydroxypropyl,
acetyl, carboxymethyl, phosphate, succinyl, octenylsuccinyl) attached to that oxygen find room in the crystal
structure without overlapping a protein heavy atom or another sugar atom?

Inputs (read only):
  data/analysis/45_hydroxyl_map/units_labelled.tsv  439 units with graph-based labels C1..C6, O1..O6
  data/analysis/45_hydroxyl_map/atoms.tsv           held flag and polar partners per unit x atom; cl50
  data/analysis/45_hydroxyl_map/cell_summary.tsv    45 state per family x subsite x atom (INV/VAR/MC/EXP/BUR)
  data/census/cif/*.cif                             deposited coordinates
Method:
  Substituent fragments are built as Cprev-Cring-O-R (ethyl ether/ester of R) with RDKit: ETKDGv3 conformers,
  MMFF94 optimised, hydrogens dropped. Each conformer is placed on the deposited Cprev, Cring, O by matching the local
  frame at O (origin O, e1 along O->Cring, e2 toward Cprev in the Cprev-Cring-O plane), so that the deposited O is the
  substituent's attachment oxygen and the fragment's own C-O-X angle and internal geometry are kept; the pose is then
  rotated about the Cring-O bond in 30 deg steps (12 values), which samples the Cring-O torsion systematically.
  Cprev: C1 for O2, C2 for O3, C5 for O6.
  O6 rotamers: O6 is also modelled on C6 at omega(O5-C5-C6-O6) = -60 (gg), +60 (gt), 180 (tg), C6-O6 1.43 A,
  C5-C6-O6 111 deg; the deposited O6 is removed from the sugar set for those poses.
  Per pose: minimum distance of any substituent heavy atom to protein heavy atoms (all chains of the asymmetric unit,
  waters removed, symmetry mates not generated) and to sugar heavy atoms (the covalently connected non-protein,
  non-water component containing the unit; atoms within three bonds of the first substituent atom - the O itself,
  Cring and the ring neighbours of Cring - excluded; for modelled O6 also the deposited O6).
  A pose fits at cut c if d_prot >= c (c = 2.8, 3.0 main, 3.2 A) and d_sugar >= 2.8 A.
  Control: the hydroxyl hydrogen itself (O-H 0.97 A, Cring-O-H 108.5 deg, 12 torsions), fits if H...protein heavy
  >= 1.8 A and H...sugar heavy >= 1.8 A (1.6 and 2.0 also reported).
Aggregation: unit fits = any pose fits; entry fits if >= half of its units at that subsite fit; 50 % group fits if
  >= half of its entries fit (as 45).
Run: ~/miniconda3/bin/python code/analysis/a46_substituent_fit.py
"""
import sys, time, subprocess
from collections import defaultdict
from multiprocessing import Pool
from pathlib import Path
import numpy as np
import pandas as pd
import gemmi
import networkx as nx
from scipy.spatial import cKDTree
from rdkit import Chem
from rdkit.Chem import AllChem

P = Path(__file__).resolve().parents[2]
D45 = P / 'data/analysis/45_hydroxyl_map'
CIF = P / 'data/census/cif'
OUT = P / 'data/analysis/46_substituent_fit'
OUT.mkdir(parents=True, exist_ok=True)

FAMS = ['GH13', 'GH31', 'GH14', 'GH15']
SUBS = ['-3', '-2', '-1', '+1', '+2', '+3']
CUTS = [2.8, 3.0, 3.2]
MAIN = 3.0
CUT_SUGAR = 2.8
H_CUTS = [1.6, 1.8, 2.0]
UPPER = 8.0
ROT = np.deg2rad(np.arange(0, 360, 30))
OMEGA = {'O6gg': -60.0, 'O6gt': 60.0, 'O6tg': 180.0}
SEED = 4646

# name -> (list of SMILES with Cprev-Cring-O prefix, conformers per SMILES)
SUBST = {
    'methyl':        (['CCOC'], 30),
    'hydroxyethyl':  (['CCOCCO'], 300),
    'hydroxypropyl': (['CCOC[C@H](O)C', 'CCOC[C@@H](O)C'], 300),          # both C2 configurations (racemic)
    'acetyl':        (['CCOC(C)=O'], 100),
    'carboxymethyl': (['CCOCC(=O)[O-]'], 300),
    'phosphate':     (['CCOP(=O)([O-])[O-]'], 100),
    'succinyl':      (['CCOC(=O)CCC(=O)[O-]'], 400),
    'octenylsuccinyl': (['CCOC(=O)C[C@H](C(=O)[O-])C/C=C/CCCCC',
                         'CCOC(=O)C[C@@H](C(=O)[O-])C/C=C/CCCCC'], 400),  # both configurations at the branch C
}
ORDER = ['methyl', 'hydroxyethyl', 'hydroxypropyl', 'acetyl', 'carboxymethyl', 'phosphate', 'succinyl',
         'octenylsuccinyl']
ABBR = dict(methyl='Me', hydroxyethyl='HE', hydroxypropyl='HP', acetyl='Ac', carboxymethyl='CM', phosphate='P',
            succinyl='Suc', octenylsuccinyl='OSA')


def is_aa(name):
    t = gemmi.find_tabulated_residue(name)
    return bool(t and t.is_amino_acid())


def read(pdb):
    st = gemmi.read_structure(str(CIF / f'{pdb}.cif'))
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    st.setup_entities()
    return st


def frame(cprev, cring, o):
    e1 = cring - o
    e1 = e1 / np.linalg.norm(e1)
    v = cprev - cring
    e2 = v - v.dot(e1) * e1
    e2 = e2 / np.linalg.norm(e2)
    return np.stack([e1, e2, np.cross(e1, e2)])


def rotx(t):
    c, s = np.cos(t), np.sin(t)
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])


def place(a, b, c, bond, angle, tor):
    """NeRF: atom d with |cd| = bond, angle(b,c,d) = angle, dihedral(a,b,c,d) = tor (degrees)."""
    angle, tor = np.deg2rad(angle), np.deg2rad(tor)
    bc = (c - b) / np.linalg.norm(c - b)
    n = np.cross(b - a, bc)
    n /= np.linalg.norm(n)
    m = np.cross(n, bc)
    d2 = np.array([-bond * np.cos(angle), bond * np.sin(angle) * np.cos(tor), bond * np.sin(angle) * np.sin(tor)])
    return c + d2[0] * bc + d2[1] * m + d2[2] * n


def dihedral(a, b, c, d):
    b0, b1, b2 = a - b, c - b, d - c
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - b0.dot(b1n) * b1n
    w = b2 - b2.dot(b1n) * b1n
    return np.degrees(np.arctan2(np.cross(b1n, v).dot(w), v.dot(w)))


# ------------------------------------------------------------------ fragment library (local frame coordinates)
def build_library():
    lib, info = {}, []
    for name in ORDER:
        smis, n = SUBST[name]
        loc, els = [], None
        for smi in smis:
            mol = Chem.AddHs(Chem.MolFromSmiles(smi))
            ps = AllChem.ETKDGv3()
            ps.randomSeed = SEED
            ps.pruneRmsThresh = 0.05
            cids = list(AllChem.EmbedMultipleConfs(mol, numConfs=n, params=ps))
            res = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=2000)
            heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
            assert heavy[:3] == [0, 1, 2]
            e = [mol.GetAtomWithIdx(i).GetSymbol() for i in heavy[3:]]
            els = els or e
            nconv = sum(1 for r in res if r[0] == 0)
            for cid in cids:
                X = mol.GetConformer(cid).GetPositions()[heavy]
                M = frame(X[0], X[1], X[2])
                loc.append((X[3:] - X[2]) @ M.T)
            ang = np.degrees(np.arccos(np.dot((X[1] - X[2]) / np.linalg.norm(X[1] - X[2]),
                                              (X[3] - X[2]) / np.linalg.norm(X[3] - X[2]))))
            info.append(dict(substituent=name, smiles=smi, conformers=len(cids), mmff_converged=nconv,
                             n_heavy_substituent=len(e), C_O_X_angle_last=round(ang, 1),
                             O_X_bond_last=round(float(np.linalg.norm(X[3] - X[2])), 3)))
        L = np.array(loc)                                  # (nconf, nsub, 3)
        poses = np.concatenate([L @ rotx(t).T for t in ROT])   # (nconf*12, nsub, 3)
        lib[name] = (poses, els)
    # hydroxyl hydrogen control
    a = np.deg2rad(108.5)
    lib['OH_H'] = (np.array([[[0.97 * np.cos(a), 0.97 * np.sin(a) * np.cos(t), 0.97 * np.sin(a) * np.sin(t)]]
                             for t in ROT]), ['H'])
    return lib, pd.DataFrame(info)


LIB = None


def parse_labels(s):
    out = {}
    for tok in s.split():
        k, v = tok.split('=', 1)
        an, el = v[:-3], v[-2]
        out[k] = (an, el)
    return out


def neighbours_of_cring(lab):
    return {'O2': ['C1', 'C3'], 'O3': ['C2', 'C4'], 'O6': ['C5']}[lab]


def eval_poses(poses_local, M, O, ptree, stree, return_best=False):
    W = poses_local @ M + O                         # (np, ns, 3)
    flat = W.reshape(-1, 3)
    dp = ptree.query(flat, distance_upper_bound=UPPER)[0].reshape(W.shape[:2]).min(1)
    dp = np.minimum(dp, UPPER)
    if stree is not None:
        ds = stree.query(flat, distance_upper_bound=UPPER)[0].reshape(W.shape[:2]).min(1)
        ds = np.minimum(ds, UPPER)
    else:
        ds = np.full(len(W), UPPER)
    if return_best:
        return dp, ds, W
    return dp, ds


HB_POLAR = 2.5


def qmin(tree, X):
    if tree is None:
        return np.full(len(X), UPPER)
    return np.minimum(tree.query(X, distance_upper_bound=UPPER)[0], UPPER)


def eval_hb(poses_local, els, M, O, pp, pa, sp, sa, cut=MAIN):
    """Polar-aware criterion: substituent O/N to protein/sugar N/O >= 2.5 A (hydrogen-bond or salt-bridge distance);
    any other pair >= cut (protein) or >= 2.8 A (sugar)."""
    W = poses_local @ M + O
    ns = W.shape[1]
    flat = W.reshape(-1, 3)
    pol = np.tile(np.array([e in ('O', 'N') for e in els]), len(W))
    d_pp, d_pa = qmin(pp, flat), qmin(pa, flat)
    d_sp, d_sa = qmin(sp, flat), qmin(sa, flat)
    ok = (d_pa >= cut) & (d_sa >= CUT_SUGAR) & np.where(pol, (d_pp >= HB_POLAR) & (d_sp >= HB_POLAR),
                                                         (d_pp >= cut) & (d_sp >= CUT_SUGAR))
    return ok.reshape(len(W), ns).all(1)


def process_entry(args):
    pdb, rows = args
    st = read(pdb)
    m = st[0]
    PX, HX, Hkey, Hel, Hres, Pel = [], [], [], [], [], []
    for ch in m:
        for res in ch:
            if res.name in ('HOH', 'DOD', 'WAT'):
                continue
            prot = is_aa(res.name)
            for a in res:
                if prot:
                    PX.append(a.pos.tolist())
                    Pel.append(a.element.name.upper())
                else:
                    HX.append(a.pos.tolist())
                    Hkey.append((ch.name, str(res.seqid.num), res.name, a.name))
                    Hel.append(a.element.name.upper())
                    Hres.append((ch.name, str(res.seqid.num), res.name))
    PX, HX = np.array(PX), np.array(HX)
    ptree = cKDTree(PX)
    Ppol = np.isin(Pel, ['N', 'O'])
    pp, pa = cKDTree(PX[Ppol]), cKDTree(PX[~Ppol])
    Hpol = np.isin(Hel, ['N', 'O'])
    hidx = {k: i for i, k in enumerate(Hkey)}
    # covalent components of the non-protein atoms
    G = nx.Graph()
    G.add_nodes_from(set(Hres))
    ht = cKDTree(HX)
    for i, j in ht.query_pairs(1.95):
        lim = 1.95 if 'S' in (Hel[i], Hel[j]) else 1.75
        if np.linalg.norm(HX[i] - HX[j]) < lim and Hres[i] != Hres[j]:
            G.add_edge(Hres[i], Hres[j])
    comp_of = {}
    for c in nx.connected_components(G):
        for r in c:
            comp_of[r] = c
    Hres_arr = np.array(['|'.join(r) for r in Hres])
    out = []
    for r in rows:
        rk = (r['chain'], str(r['seq']), r['comp'])
        if rk not in comp_of:
            out.append(dict(**r['key'], label='', status='residue not found'))
            continue
        comp = {'|'.join(x) for x in comp_of[rk]}
        in_comp = np.isin(Hres_arr, list(comp))
        L = parse_labels(r['labels'])
        pos = {}
        for k, (an, el) in L.items():
            i = hidx.get((r['chain'], str(r['seq']), r['comp'], an))
            if i is not None:
                pos[k] = (HX[i], i, el)
        # overlapping copies (alpha/beta anomers or alternative ligands deposited as separate residues at the same
        # place): any other residue with >= 3 atoms within 1.2 A of this unit's atoms is removed from the sugar set
        own = np.array([v[0] for v in pos.values()])
        rkey = '|'.join(rk)
        dup = set()
        for i in set().union(*ht.query_ball_point(own, 1.2)):
            if Hres_arr[i] != rkey:
                dup.add(Hres_arr[i])
        dup = {d for d in dup if sum(1 for i in set().union(*ht.query_ball_point(own, 1.2)) if Hres_arr[i] == d) >= 3}
        in_comp = in_comp & ~np.isin(Hres_arr, list(dup))
        sites = []
        for lab, prev in (('O2', 'C1'), ('O3', 'C2'), ('O6', 'C5')):
            cr = {'O2': 'C2', 'O3': 'C3', 'O6': 'C6'}[lab]
            if lab in pos and pos[lab][2] == 'O' and cr in pos and prev in pos:
                excl = [pos[lab][1], pos[cr][1]] + [pos[x][1] for x in neighbours_of_cring(lab) if x in pos]
                sites.append((lab, 'deposited', pos[lab][0], pos[cr][0], pos[prev][0], excl))
        if all(k in pos for k in ('O5', 'C5', 'C6')):
            for lab, w in OMEGA.items():
                o6 = place(pos['O5'][0], pos['C5'][0], pos['C6'][0], 1.43, 111.0, w)
                excl = [pos['C6'][1], pos['C5'][1]] + ([pos['O6'][1]] if 'O6' in pos else [])
                sites.append((lab, 'modelled', o6, pos['C6'][0], pos['C5'][0], excl))
        for lab, src, O, Cr, Cp, excl in sites:
            mask = in_comp.copy()
            mask[excl] = False
            stree = cKDTree(HX[mask]) if mask.any() else None
            sp = cKDTree(HX[mask & Hpol]) if (mask & Hpol).any() else None
            sa = cKDTree(HX[mask & ~Hpol]) if (mask & ~Hpol).any() else None
            M = frame(Cp, Cr, O)
            o_prot = float(ptree.query(O)[0])
            o_sug = float(stree.query(O)[0]) if stree is not None else np.nan
            base = dict(**r['key'], label=lab, o_source=src, status='ok', o_min_prot=round(o_prot, 2),
                        o_min_sugar=round(o_sug, 2), overlapping_copies_removed=' '.join(sorted(dup)))
            if src == 'modelled':
                base['omega_check'] = round(dihedral(pos['O5'][0], pos['C5'][0], pos['C6'][0], O), 1)
            for name in ['OH_H'] + ORDER:
                poses, els = LIB[name]
                dp, ds = eval_poses(poses, M, O, ptree, stree)
                row = dict(base, substituent=name, n_poses=len(poses))
                if name == 'OH_H':
                    for c in H_CUTS:
                        ok = (dp >= c) & (ds >= c)
                        row[f'fits_{c}'] = bool(ok.any())
                        row[f'frac_{c}'] = round(float(ok.mean()), 4)
                    sok = ds >= 1.8
                else:
                    sok = ds >= CUT_SUGAR
                    row['frac_sugar_ok'] = round(float(sok.mean()), 4)
                    for c in CUTS:
                        ok = (dp >= c) & sok
                        row[f'fits_{c}'] = bool(ok.any())
                        row[f'frac_{c}'] = round(float(ok.mean()), 4)
                    ok = eval_hb(poses, els, M, O, pp, pa, sp, sa)
                    row['fits_hb'] = bool(ok.any())
                    row['frac_hb'] = round(float(ok.mean()), 4)
                    row['best_pose_hb'] = int(np.argmax(ok)) if ok.any() else -1
                row['best_clear_prot'] = round(float(dp[sok].max()), 2) if sok.any() else np.nan
                row['best_clear_prot_nosugar'] = round(float(dp.max()), 2)
                row['best_pose'] = int(np.argmax(np.where(sok, dp, -1.0)))
                out.append(row)
    return out


# ------------------------------------------------------------------ aggregation (as 45)
def agg(D, cut, tag):
    rows = []
    col = f'fits_{cut}'
    for (f, s, lab, sub), g in D.groupby(['family', 'subsite', 'label', 'substituent']):
        e = g.groupby('pdb').agg(fit=(col, 'mean'), frac=(f'frac_{cut}', 'mean'), cl50=('cl50', 'first'),
                                 clear=('best_clear_prot', 'median'), held=('held', 'mean'))
        e['fit'] = e.fit >= 0.5
        c = e.groupby('cl50').agg(fit=('fit', 'mean'))
        c['fit'] = c.fit >= 0.5
        fc = c.fit.mean()
        anyu = bool(g[col].any())
        cat = 'F' if fc >= 0.5 else ('s' if anyu else 'N')
        rows.append(dict(set=tag, cut=cut, family=f, subsite=s, label=lab, substituent=sub, n_units=len(g),
                         n_entries=len(e), n_cl50=len(c), frac_units_fit=round(g[col].mean(), 3),
                         frac_entries_fit=round(e.fit.mean(), 3), frac_cl_fit=round(fc, 3),
                         n_cl_fit=int(c.fit.sum()), mean_pose_frac=round(g[f'frac_{cut}'].mean(), 4),
                         median_best_clear=round(float(g.best_clear_prot.median()), 2),
                         frac_units_held=round(float(g.held.mean()), 3) if g.held.notna().any() else np.nan,
                         category=cat))
    return pd.DataFrame(rows)


def main():
    global LIB
    t0 = time.time()
    LIB, linfo = build_library()
    linfo.to_csv(OUT / 'fragment_library.tsv', sep='\t', index=False)
    print(linfo.to_string())
    print({k: v[0].shape for k, v in LIB.items()}, f'{time.time() - t0:.0f} s')

    U = pd.read_csv(D45 / 'units_labelled.tsv', sep='\t', dtype=str, keep_default_na=False)
    A = pd.read_csv(D45 / 'atoms.tsv', sep='\t', dtype=str, keep_default_na=False)
    CL = A.groupby(['family', 'pdb']).cl50.first().to_dict()
    HELD = {(r.pdb, r.instance_id, r.ring, r.label): (r.held == 'True', r.polar_partners) for r in A.itertuples()}
    jobs = defaultdict(list)
    for r in U.itertuples():
        key = dict(family=r.family, pdb=r.pdb, cl50=CL.get((r.family, r.pdb), '') or f'{r.family}_entry_{r.pdb}',
                   instance_id=r.instance_id, ring=r.ring, comp=r.comp, chain=r.chain, seq=r.seq, subsite=r.subsite,
                   kind=r.kind)
        jobs[r.pdb].append(dict(key=key, chain=r.chain, seq=r.seq, comp=r.comp, labels=r.labels))
    items = sorted(jobs.items())
    print('units', len(U), 'entries', len(items))
    res = []
    with Pool(7) as pool:
        for k, out in enumerate(pool.imap_unordered(process_entry, items, chunksize=1)):
            res.extend(out)
            if k % 20 == 0:
                print(k, len(items), f'{time.time() - t0:.0f} s', flush=True)
    D = pd.DataFrame(res)
    bad = D[D.status != 'ok']
    print('not ok', len(bad))
    D = D[D.status == 'ok'].copy()
    hl = []
    for r in D.itertuples():
        lab = 'O6' if r.label.startswith('O6') else r.label
        hl.append(HELD.get((r.pdb, r.instance_id, r.ring, lab), (np.nan, ''))[0])
    D['held'] = hl
    D['held_label_ref'] = np.where(D.label.str.startswith('O6'), 'O6 (deposited)', D.label)
    D.to_csv(OUT / 'site_results.tsv', sep='\t', index=False)
    # O6* = any of deposited + gg/gt/tg
    O6 = D[D.label.str.startswith('O6')]
    keycols = ['family', 'pdb', 'cl50', 'instance_id', 'ring', 'comp', 'chain', 'seq', 'subsite', 'kind', 'substituent']
    fcols = [c for c in D.columns if c.startswith('fits_')]
    pcols = [c for c in D.columns if c.startswith('frac_') and c != 'frac_sugar_ok']
    agg_spec = {c: 'any' for c in fcols}
    agg_spec.update({c: 'max' for c in pcols})
    agg_spec.update(best_clear_prot='max', held='first', o_min_prot='min', n_poses='sum')
    O6s = O6.groupby(keycols, dropna=False).agg(agg_spec).reset_index()
    O6s['label'] = 'O6*'
    for c in fcols:
        O6s[c] = O6s[c].astype(bool)
    DD = pd.concat([D, O6s], ignore_index=True)
    DD.to_pickle(OUT / '_dd.pkl')
    post(DD, linfo)
    print(f'done {time.time() - t0:.0f} s')


# ------------------------------------------------------------------ post-processing, tables, report
LIT = [
    # family, subsite, label, substituent, literature text, verdict (T tolerated / N not / P partial / M mixed)
    ('GH13', '-1', 'O3', 'phosphate', 'not productive (PPA, B. subtilis, A. oryzae; Takeda 1983/1986/1987)', 'N'),
    ('GH13', '-1', 'O6*', 'phosphate', 'not productive (PPA, B. subtilis, A. oryzae; Takeda 1983/1986/1987)', 'N'),
    ('GH13', '-1', 'O2', 'hydroxyethyl', 'not tolerated at -1 (PPA; Yoshida 1984 II; position mixed O2/O3/O6)', 'N'),
    ('GH13', '-1', 'O3', 'hydroxyethyl', 'not tolerated at -1 (PPA; Yoshida 1984 II; position mixed)', 'N'),
    ('GH13', '-1', 'O6*', 'hydroxyethyl', 'not tolerated at -1 (PPA; Yoshida 1984 II; position mixed)', 'N'),
    ('GH13', '+1', 'O6*', 'phosphate', 'tolerated (PPA, A. oryzae); not productive (B. subtilis)', 'M'),
    ('GH13', '+1', 'O3', 'phosphate', 'not productive (PPA, B. subtilis, A. oryzae)', 'N'),
    ('GH13', '-2', 'O3', 'phosphate', 'tolerated (PPA, B. subtilis); not productive (A. oryzae)', 'M'),
    ('GH13', '-2', 'O6*', 'phosphate', 'tolerated (B. subtilis); not productive (PPA, A. oryzae)', 'M'),
    ('GH13', '+2', 'O6*', 'phosphate', 'tolerated (B. subtilis)', 'T'),
    ('GH13', '+2', 'O3', 'phosphate', 'tolerated (A. oryzae); not productive (PPA, B. subtilis)', 'M'),
    ('GH13', '-3', 'O6*', 'phosphate', 'tolerated (PPA, A. oryzae)', 'T'),
    ('GH13', '-3', 'O3', 'phosphate', 'tolerated (A. oryzae; liquefying amylase, Kamasaka 1997)', 'T'),
    ('GH13', '+3', 'O3', 'phosphate', 'tolerated (PPA, B. subtilis)', 'T'),
    ('GH13', '+1', 'O2', 'methyl', 'terminal residues of B. subtilis limit dextrins >90 % unsubstituted (van der Burgt 2000b; position mixed)', 'N'),
    ('GH13', '+1', 'O3', 'methyl', 'idem (position mixed)', 'N'),
    ('GH13', '+1', 'O6*', 'methyl', 'idem (position mixed)', 'N'),
    ('GH14', '-1', 'O6*', 'phosphate', 'not productive (sweet-potato beta-amylase; Takeda 1987)', 'N'),
    ('GH14', '-2', 'O6*', 'phosphate', 'not productive (sweet-potato beta-amylase; Takeda 1987)', 'N'),
    ('GH14', '+1', 'O6*', 'phosphate', 'tolerated (sweet-potato beta-amylase; Takeda 1987)', 'T'),
    ('GH14', '+1', 'O3', 'phosphate', 'not productive (Takeda 1987)', 'N'),
    ('GH14', '-1', 'O3', 'phosphate', 'not productive (Takeda 1987)', 'N'),
    ('GH14', '-2', 'O3', 'phosphate', 'not productive (Takeda 1987)', 'N'),
    ('GH14', '-2', 'O2', 'methyl', 'not tolerated at -2/-1 (B. cereus beta-amylase; van der Burgt 2000b; position mixed)', 'N'),
    ('GH14', '-1', 'O2', 'methyl', 'idem (position mixed)', 'N'),
    ('GH14', '+1', 'O2', 'methyl', 'not tolerated at +1 (van der Burgt 2000b; position mixed)', 'N'),
    ('GH15', '-1', 'O6*', 'phosphate', 'not productive (A. niger GA; Abe 1982, Takeda 1983/1987)', 'N'),
    ('GH15', '-1', 'O3', 'phosphate', 'not productive (A. niger GA)', 'N'),
    ('GH15', '+1', 'O6*', 'phosphate', 'tolerated (A. niger GA; Takeda 1983/1987; Kamasaka 1997)', 'T'),
    ('GH15', '+1', 'O3', 'phosphate', 'not productive (A. niger GA)', 'N'),
    ('GH15', '+2', 'O3', 'phosphate', 'tolerated (implied: one Glc left on non-reducing side)', 'T'),
    ('GH15', '-1', 'O2', 'methyl', 'O-methyl blocks (Lemieux 1996, methyl isomaltoside)', 'N'),
    ('GH15', '-1', 'O3', 'methyl', 'O-methyl blocks (Lemieux 1996)', 'N'),
    ('GH15', '-1', 'O6*', 'methyl', 'O-methyl blocks (Lemieux 1996)', 'N'),
    ('GH15', '+1', 'O2', 'methyl', 'O-methyl tolerated, kinetics unaffected (Lemieux 1996; isomaltoside +1)', 'T'),
    ('GH15', '+1', 'O3', 'methyl', 'O-methyl costs 1.9 kcal/mol (~25-fold kcat/KM) (Lemieux 1996)', 'P'),
]


def cell_txt(r):
    if r is None:
        return '—'
    return f'{r.category} {r.frac_cl_fit:.2f} ({r.n_cl_fit}/{r.n_cl50}; ent {r.frac_entries_fit:.2f}/{r.n_entries})'


def post(DD, linfo):
    CS = pd.read_csv(D45 / 'cell_summary.tsv', sep='\t', dtype=str, keep_default_na=False)
    STATE = {(r.set, r.family, r.subsite, r.label): r.state for r in CS.itertuples()}
    SETS = {'all_glucoselike': DD, 'glc_units': DD[DD.kind == 'Glc']}
    S = []
    for tag, d in SETS.items():
        for c in CUTS + ['hb']:
            S.append(agg(d[d.substituent != 'OH_H'], c, tag))
        for c in H_CUTS:
            S.append(agg(d[d.substituent == 'OH_H'], c, tag))
    S = pd.concat(S, ignore_index=True)
    S['cut'] = S.cut.astype(str)
    lab_state = lambda tag, f, s, lab: STATE.get((tag, f, s, 'O6' if lab.startswith('O6') else lab), '—')
    S['state45'] = [lab_state(r.set, r.family, r.subsite, r.label) for r in S.itertuples()]
    so = {s: i for i, s in enumerate(SUBS)}
    lo = {l: i for i, l in enumerate(['O2', 'O3', 'O6', 'O6gg', 'O6gt', 'O6tg', 'O6*'])}
    S = S.sort_values(['set', 'cut', 'family', 'substituent', 'subsite', 'label'],
                      key=lambda x: x.map(so) if x.name == 'subsite' else (x.map(lo) if x.name == 'label' else x))
    S.to_csv(OUT / 'cell_fit_summary.tsv', sep='\t', index=False)
    IDX = {(r.set, r.cut, r.family, r.subsite, r.label, r.substituent): r for r in S.itertuples()}

    md = []
    w = md.append
    w('# 46 - Steric fit of starch substituents at each glucose oxygen and subsite\n')
    w(f'Generated by `code/analysis/a46_substituent_fit.py` on {time.strftime("%Y-%m-%d %H:%M")}.\n')
    w('METHODS_PLACEHOLDER\n')
    w('## Fragment library (RDKit ETKDGv3 + MMFF94, heavy atoms only; each conformer x 12 rotations about Cring-O)\n')
    w(linfo.to_markdown(index=False) + '\n')
    w('Pose counts per site: ' + ', '.join(f'{k} {LIB[k][0].shape[0]}' for k in ['OH_H'] + ORDER) + '.\n')

    # per-substituent tables, main cut, both sets
    labs = ['O2', 'O3', 'O6', 'O6*']
    CK = str(MAIN)
    for tag in ['all_glucoselike', 'glc_units']:
        w(f'\n## Fit tables at {MAIN} A, set {tag}\n')
        w('Cell: category (F = fits in >= 50 % of 50 % groups; s = fits in some unit but < 50 % of groups; '
          'N = fits in no unit), fraction of groups (n fit / n groups; fraction of entries / n entries). '
          'O6 = deposited O6 position; O6* = any of deposited, gg, gt, tg. 45 state in the first column of each '
          'oxygen.\n')
        for f in FAMS:
            for sub in ORDER:
                rows_ = []
                for s in SUBS:
                    cells = [IDX.get((tag, CK, f, s, l, sub)) for l in labs]
                    if all(c is None for c in cells):
                        continue
                    rows_.append(f'| {s} | ' + ' | '.join(
                        f'{lab_state(tag, f, s, l)} | {cell_txt(c)}' for l, c in zip(labs, cells)) + ' |')
                if rows_:
                    w(f'\n### {f} - {sub} ({tag})\n')
                    w('| subsite | 45 O2 | O2 | 45 O3 | O3 | 45 O6 | O6 (dep.) | 45 O6 | O6* |')
                    w('|---|---|---|---|---|---|---|---|---|')
                    md.extend(rows_)

    # combined compact table
    for tag, CK in [('all_glucoselike', '3.0'), ('glc_units', '3.0'), ('all_glucoselike', 'hb'), ('glc_units', 'hb')]:
        crit = f'{MAIN} A flat' if CK != 'hb' else 'polar-aware (O/N...N/O >= 2.5 A, other pairs >= 3.0 A)'
        w(f'\n## Combined table, {crit} ({tag}): 45 state, then category per substituent in the order '
          + ' '.join(ABBR[x] for x in ORDER) + '\n')
        for f in FAMS:
            w(f'\n### {f}\n')
            w('| subsite | O2 | O3 | O6 (dep.) | O6* |')
            w('|---|---|---|---|---|')
            for s in SUBS:
                cells = []
                ok = False
                for l in labs:
                    cs = [IDX.get((tag, CK, f, s, l, sub)) for sub in ORDER]
                    if all(c is None for c in cs):
                        cells.append('—')
                        continue
                    ok = True
                    cells.append(f'**{lab_state(tag, f, s, l)}** ' + ' '.join(c.category if c is not None else '·'
                                                                         for c in cs))
                if ok:
                    w(f'| {s} | ' + ' | '.join(cells) + ' |')

    # check (a): hydroxyl hydrogen and methyl vs exposure, GH13 glc units
    w('\n## Check (a): hydroxyl hydrogen control and methyl, GH13, glucose units (kind Glc)\n')
    w('OH-H fits at 1.8 A (H...protein heavy and H...sugar heavy >= 1.8 A); O itself: median minimum distance of the '
      'deposited O to a protein heavy atom and share of units with that distance < 3.0 A; methyl at 3.0 A.\n')
    w('| subsite | label | 45 state | n units | OH-H fits (units) | O-prot median A | O < 3.0 A (units) | methyl fits (units) | methyl F/s/N |')
    w('|---|---|---|---|---|---|---|---|---|')
    g = DD[(DD.family == 'GH13') & (DD.kind == 'Glc')]
    for s in SUBS:
        for l in ['O2', 'O3', 'O6', 'O6*']:
            h = g[(g.subsite == s) & (g.label == l) & (g.substituent == 'OH_H')]
            me = g[(g.subsite == s) & (g.label == l) & (g.substituent == 'methyl')]
            if h.empty:
                continue
            cat = IDX.get(('glc_units', '3.0', 'GH13', s, l, 'methyl'))
            w(f'| {s} | {l} | {lab_state("glc_units", "GH13", s, l)} | {len(h)} | {h["fits_1.8"].mean():.2f} | '
              f'{h.o_min_prot.median():.2f} | {(h.o_min_prot < 3.0).mean():.2f} | {me["fits_3.0"].mean():.2f} | '
              f'{cat.category if cat is not None else "—"} |')
    hh = DD[DD.substituent == 'OH_H']
    w('\nOH-H over all units, all families (deposited O only): fits at 1.6/1.8/2.0 A in '
      + ', '.join(f'{hh[hh.o_source == "deposited"][f"fits_{c}"].mean():.3f}' for c in H_CUTS)
      + f' of {int((hh.o_source == "deposited").sum())} unit-oxygen sites.\n')
    me = DD[(DD.substituent == 'methyl') & DD.label.isin(['O2', 'O3', 'O6'])].copy()
    me['state'] = [lab_state('all_glucoselike', r.family, r.subsite, r.label) for r in me.itertuples()]
    w('Methyl fit (3.0 A) by 45 state of the cell (all families, all units, deposited O2/O3/O6):\n')
    w('| 45 state | n unit-sites | methyl fits | held units among them |')
    w('|---|---|---|---|')
    for stt, gg in me.groupby('state'):
        w(f'| {stt} | {len(gg)} | {gg["fits_3.0"].mean():.2f} | {gg.held.mean():.2f} |')
    w('\nMethyl fit by the unit\'s own held flag (45, deposited O2/O3/O6): ' + '; '.join(
        f'held={k}: {v["fits_3.0"].mean():.2f} of {len(v)}' for k, v in me.groupby('held')) + '.\n')

    # check (b): cut-off sensitivity
    w('\n## Check (b): sensitivity to the protein clash cut-off (2.8 / 3.0 / 3.2 A)\n')
    w('Number of family x subsite x oxygen cells (O2, O3, O6, O6*; all_glucoselike) in each category, and number of '
      'cells whose category differs from the 3.0 A result.\n')
    w('| substituent | F/s/N at 2.8 | F/s/N at 3.0 | F/s/N at 3.2 | F/s/N polar-aware | changed 2.8 vs 3.0 | changed 3.2 vs 3.0 | changed polar-aware vs 3.0 | mean unit fit fraction 2.8/3.0/3.2/polar-aware |')
    w('|---|---|---|---|---|---|---|---|---|')
    SS = S[(S.set == 'all_glucoselike') & S.label.isin(labs) & (S.substituent != 'OH_H')]
    for sub in ORDER:
        x = SS[SS.substituent == sub]
        cc = ['2.8', '3.0', '3.2', 'hb']
        cnt = {c: x[x.cut == c].category.value_counts() for c in cc}
        piv = x.pivot_table(index=['family', 'subsite', 'label'], columns='cut', values='category', aggfunc='first')
        fr = [x[x.cut == c].frac_units_fit.mean() for c in cc]
        w(f'| {sub} | ' + ' | '.join(f'{cnt[c].get("F", 0)}/{cnt[c].get("s", 0)}/{cnt[c].get("N", 0)}' for c in cc)
          + f' | {(piv["2.8"] != piv["3.0"]).sum()} | {(piv["3.2"] != piv["3.0"]).sum()} | '
          + f'{(piv["hb"] != piv["3.0"]).sum()} | ' + '/'.join(f'{v:.2f}' for v in fr) + ' |')

    # literature comparison
    w('\n## Literature comparison (substituted_starch_hydrolysis_literature.md, §7, §8, §11)\n')
    w('Our columns: all_glucoselike set, 3.0 A (category, groups); 2.8 and 3.2 A categories; 45 state; share of '
      'units in which the oxygen is held. Agreement rule (stated, mechanical): literature "tolerated" agrees with F; '
      '"not tolerated" agrees with N; s or a mixed literature verdict is "partial"; otherwise "no".\n')
    w('| cell | substituent | literature | our fit 3.0 A | 2.8 / 3.2 / polar-aware | glc units only (3.0) | 45 state | held | agree? (3.0) | agree? (polar-aware) |')
    w('|---|---|---|---|---|---|---|---|---|---|')
    lrows = []
    for f, s, l, sub, txt, v in LIT:
        r = IDX.get(('all_glucoselike', '3.0', f, s, l, sub))
        rg = IDX.get(('glc_units', '3.0', f, s, l, sub))
        r28 = IDX.get(('all_glucoselike', '2.8', f, s, l, sub))
        r32 = IDX.get(('all_glucoselike', '3.2', f, s, l, sub))
        rhb = IDX.get(('all_glucoselike', 'hb', f, s, l, sub))

        def agree(r):
            if r is None:
                return 'no structural data'
            if v in ('M', 'P') or r.category == 's':
                return 'partial'
            if (v == 'T' and r.category == 'F') or (v == 'N' and r.category == 'N'):
                return 'yes'
            return 'no'
        ag, aghb = agree(r), agree(rhb)
        lrows.append(dict(family=f, subsite=s, label=l, substituent=sub, literature=txt, verdict=v,
                          fit_cat=r.category if r is not None else '', frac_cl=r.frac_cl_fit if r is not None else np.nan,
                          agree=ag, fit_cat_hb=rhb.category if rhb is not None else '', agree_hb=aghb))
        w(f'| {f} {s} {l} | {sub} | {txt} | {cell_txt(r)} | {r28.category if r28 is not None else "—"} / '
          f'{r32.category if r32 is not None else "—"} / {rhb.category if rhb is not None else "—"} | {cell_txt(rg)} | {lab_state("all_glucoselike", f, s, l)} | '
          f'{(f"{r.frac_units_held:.2f}" if r is not None and r.frac_units_held == r.frac_units_held else "—")} | {ag} | {aghb} |')
    pd.DataFrame(lrows).to_csv(OUT / 'literature_comparison.tsv', sep='\t', index=False)
    (OUT / '_report_body.md').write_text('\n'.join(md) + '\n')


# ------------------------------------------------------------------ spot checks: best pose written as PDB + PyMOL image
SPOTS = [  # family, pdb, subsite, label, substituent, tag, criterion (chosen after the first run; see REPORT.md)
    ('GH13', '5TD4', '-1', 'O3', 'phosphate', 'GH13_m1_O3_phosphate_5TD4', '3.0'),
    ('GH13', '1MFU', '+2', 'O6', 'phosphate', 'GH13_p2_O6_phosphate_1MFU', '3.0'),
    ('GH15', '1GAH', '+1', 'O6gt', 'phosphate', 'GH15_p1_O6gt_phosphate_1GAH', '3.0'),
]


def spot(fam, pdb, subsite, label, sub, tag, crit='3.0'):
    U = pd.read_csv(D45 / 'units_labelled.tsv', sep='\t', dtype=str, keep_default_na=False)
    r = U[(U.pdb == pdb) & (U.subsite == subsite) & (U.family == fam)].iloc[0]
    st = read(pdb)
    m = st[0]
    ch = m.find_chain(r.chain)
    res = [x for x in ch if str(x.seqid.num) == r.seq and x.name == r.comp][0]
    L = parse_labels(r.labels)
    pos = {k: np.array(res.find_atom(an, '*').pos.tolist()) for k, (an, el) in L.items() if res.find_atom(an, '*')}
    cr, prev = {'O2': ('C2', 'C1'), 'O3': ('C3', 'C2')}.get(label, ('C6', 'C5'))
    if label in OMEGA:
        O = place(pos['O5'], pos['C5'], pos['C6'], 1.43, 111.0, OMEGA[label])
    else:
        O = pos[label]
    D = pd.read_csv(OUT / 'site_results.tsv', sep='\t', dtype=str, keep_default_na=False)
    row = D[(D.pdb == pdb) & (D.instance_id == r.instance_id) & (D.ring == r.ring) & (D.label == label) &
            (D.substituent == sub)].iloc[0]
    poses, els = LIB[sub]
    bp = int(row.best_pose_hb) if crit == 'hb' and int(row.best_pose_hb) >= 0 else int(row.best_pose)
    W = poses[bp] @ frame(pos[prev], pos[cr], O) + O
    # neighbourhood
    lines, k = [], 1
    PX = []
    for c in m:
        for rr in c:
            if rr.name in ('HOH', 'DOD'):
                continue
            near = any(a.pos.dist(gemmi.Position(*O)) < 10.0 for a in rr)
            if not near:
                continue
            het = 'ATOM  ' if is_aa(rr.name) else 'HETATM'
            for a in rr:
                p = a.pos
                lines.append(f'{het}{k:5d} {a.name:<4s} {rr.name:>3s} {c.name[:1]}{rr.seqid.num:4d}    '
                             f'{p.x:8.3f}{p.y:8.3f}{p.z:8.3f}  1.00 20.00          {a.element.name:>2s}')
                k += 1
                if is_aa(rr.name):
                    PX.append((p.x, p.y, p.z, f'{rr.name}{rr.seqid.num}:{c.name}.{a.name}'))
    if label in OMEGA:
        lines.append(f'HETATM{k:5d} O6M  MOD Y   1    {O[0]:8.3f}{O[1]:8.3f}{O[2]:8.3f}  1.00 20.00           O')
        k += 1
    for i, (x, e) in enumerate(zip(W, els)):
        lines.append(f'HETATM{k:5d} {e + str(i + 1):<4s} SUB Z   1    {x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}  1.00 20.00'
                     f'          {e:>2s}')
        k += 1
    fn = OUT / 'spot' / f'{tag}.pdb'
    fn.parent.mkdir(exist_ok=True)
    fn.write_text('\n'.join(lines) + '\nEND\n')
    P_ = np.array([p[:3] for p in PX])
    dd = np.linalg.norm(W[:, None, :] - P_[None], axis=2)
    i, j = np.unravel_index(dd.argmin(), dd.shape)
    near_txt = f'{els[i]}{i + 1} - {PX[j][3]} {dd[i, j]:.2f} A'
    sel_o = f'resn MOD' if label in OMEGA else f'chain {r.chain[:1]} and resi {r.seq} and name {L[label][0]}'
    pml = f"""load {fn}, cx
hide everything
set ray_opaque_background, 1
bg_color white
show sticks, resn SUB or (byres (resn SUB around 6)) or (not polymer and not resn SUB)
hide sticks, polymer and not (byres (resn SUB around 6))
show cartoon, polymer
set cartoon_transparency, 0.75
color grey70, polymer and elem C
color yellow, (not polymer) and elem C and not resn SUB
color magenta, resn SUB and elem C
color magenta, resn SUB
show spheres, resn SUB
set sphere_scale, 0.3
distance d1, resn SUB and name {els[i]}{i + 1}, /cx//{PX[j][3].split(':')[1].split('.')[0][:1]}/{PX[j][3].split(':')[0][3:]}/{PX[j][3].split('.')[-1]}
color black, d1
orient resn SUB or ({sel_o})
zoom resn SUB, 6
set label_size, 20
ray 1200, 900
png {OUT / 'spot' / (tag + '.png')}, dpi=150
"""
    (OUT / 'spot' / f'{tag}.pml').write_text(pml)
    subprocess.run(['pymol', '-cq', str(OUT / 'spot' / f'{tag}.pml')], check=False)
    return dict(tag=tag, family=fam, pdb=pdb, subsite=subsite, label=label, substituent=sub, comp=r.comp,
                crit=crit, pose=bp, fits_3=row['fits_3.0'], frac_3=row['frac_3.0'], fits_hb=row['fits_hb'],
                frac_hb=row['frac_hb'], best_clear_prot=row.best_clear_prot,
                closest=near_txt)


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'spot':
        LIB, _ = build_library()
        out = [spot(*s) for s in SPOTS]
        pd.DataFrame(out).to_csv(OUT / 'spot' / 'spot_checks.tsv', sep='\t', index=False)
        print(pd.DataFrame(out).to_string())
    elif len(sys.argv) > 1 and sys.argv[1] == 'post':
        LIB, linfo = build_library()
        post(pd.read_pickle(OUT / '_dd.pkl'), linfo)
    else:
        main()
