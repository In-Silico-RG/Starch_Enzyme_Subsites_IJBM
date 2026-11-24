"""47 - Which glucose oxygens are held at each subsite in solution: MD of 12 enzyme-glucan complexes (calc 8).

Problem: the crystal hydroxyl map (45) at subsites -1/+1 rests mostly on inhibitors (acarbose cyclitols; no O6 at +1).
The 12 MD systems of calculation 8 carry wild-type, glucose-only glucans; their partial trajectories are read here.

Inputs (read only):
  data/calc/08_md/selection/selection.tsv           system -> PDB, label_asym of the glucan, family
  data/calc/08_md/systems/<S>/system.gro, map.tsv   start structure (hydrogens present), MD resid -> PDB residue
  data/analysis/47_md_hydroxyl/traj/<S>_pocket.xtc  reduced trajectory: glucan + whole residues with any atom
        within 10 A of the glucan in system.gro (index <S>_pocket.ndx, 1-based system.gro atom numbers); every 100 ps;
        written on GUANE with  gmx_mpi trjconv -f md.xtc -n <S>_pocket.ndx -dt 100  (gromacs/titan, login node)
  data/analysis/33_subsites/units.tsv               subsite labels of the deposited units (catalytic-pair rule)
  data/analysis/44_rule_test/posmap.tsv             GH13 residue -> reference position (1AMY)
  data/analysis/45_hydroxyl_map/posmap_GH14_GH15.tsv GH14 residue -> reference position (1BYA)
  data/analysis/45_hydroxyl_map/position_identity.tsv identity (m50, n_cl50) of the crystal partner positions
  data/analysis/45_hydroxyl_map/cell_summary.tsv    crystal state per family x subsite x oxygen
  data/calc/03_frames/catalytic_residues_mapped.tsv catalytic pair (nucleophile/base, acid) per entry

Method, per frame (and once for the start structure system.gro):
  - every atom placed at the minimum image of the first glucan atom (truncated octahedron, MDAnalysis
    minimize_vectors), then the pocket C-alpha atoms superposed on system.gro (glucan not fitted);
  - polar partner of a glucan oxygen = protein N or O within 3.5 A (as 45); 'held' = at least one;
    H-bond with angle = a polar pair whose donor hydrogen gives D-H...A >= 150 deg (30 deg from linear): protein
    N-H/O-H donating to the glucan oxygen, or a glucan hydroxyl H (H2O, H3O, H4O, H6O, ROH HO1) donating to a protein
    O or to a protein N that carries no hydrogen; ring O5 and glycosidic O4/O6 act only as acceptors;
  - SASA: Shrake-Rupley, 480 points, probe 1.4 A, Bondi radii, hydrogens removed (as 45): sasa_cplx with the glucan +
    pocket protein + ions, sasa_lig with the glucan alone; rel = sasa_cplx / sasa_lig;
  - partner residue -> PDB residue (map.tsv) -> reference position (GH13: 44 posmap; GH14: 45 posmap); GH57, GH49,
    AA13, GT35 have no position map (residue named by its PDB number).
Cell state (the 45 rule with frames in place of 50 % clusters): held in >= 0.5 of frames ->
  INV if a side chain at a lineage-invariant position (45 position_identity: m50 = 1, n_cl50 >= min(3, clusters of
  the family)) holds the oxygen in >= 0.5 of frames; else VAR if frames held by a variable side chain >= frames held by
  main chain only; else MC. Not held -> EXP if median rel >= 0.25 and median sasa_cplx >= 5 A2, else BUR.
  Families without a position map: held cells are SC (side chain, identity unknown) or MC.
Binding: RMSD of the glucan heavy atoms from system.gro in the pocket frame; per-unit ring-centroid displacement;
  distance of the -1 unit C1 to the midpoint of the catalytic carboxylate centroids (33 D_MAX = 4.5 A).
Run: ~/miniconda3/bin/python code/analysis/a47_md_hydroxyl_contacts.py
"""
import re, warnings
from pathlib import Path
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib.distances import minimize_vectors
from scipy.spatial import cKDTree

warnings.filterwarnings('ignore')
P = Path(__file__).resolve().parents[2]
MD = P / 'data/calc/08_md'
OUT = P / 'data/analysis/47_md_hydroxyl'
TR = OUT / 'traj'
A33 = P / 'data/analysis/33_subsites'
A44 = P / 'data/analysis/44_rule_test'
A45 = P / 'data/analysis/45_hydroxyl_map'
FR = P / 'data/calc/03_frames'

HB, ANG = 3.5, 150.0
REL_EXPOSED, SASA_EXPOSED, D_MAX = 0.25, 5.0, 4.5
PROBE = 1.4
RADII = {'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80, 'P': 1.80, 'CA': 2.31, 'NA': 2.27, 'CL': 1.75, 'ZN': 1.39}
SUGAR_RES = {'0GA', '4GA', '6GA', '0GB', '4GB', '6GB', 'ROH'}
ION_RES = {'CA', 'NA', 'CL', 'ZN', 'MG', 'K'}
COFAC = {'HP4'}
BACKBONE = {'N', 'O', 'OXT', 'C', 'CA', 'H', 'HA'}
OXY = ['O1', 'O2', 'O3', 'O4', 'O5', 'O6']
SUBS = ['-3', '-2', '-1', '+1', '+2', '+3']


def sphere(n=480):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    th = np.pi * (1 + 5 ** 0.5) * i
    return np.c_[np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)]


SPH = sphere()


def sasa(idx, X, R, mask):
    """SASA of atoms idx against atoms where mask is True (R includes the probe)."""
    sub = np.where(mask)[0]
    tree = cKDTree(X[sub])
    out = []
    for i in idx:
        pts = X[i] + SPH * R[i]
        nb = [sub[j] for j in tree.query_ball_point(X[i], R[i] + R.max()) if sub[j] != i]
        if not nb:
            out.append(4 * np.pi * R[i] ** 2)
            continue
        Y, RY = X[nb], R[nb]
        d2 = ((pts[:, None, :] - Y[None, :, :]) ** 2).sum(-1)
        out.append(4 * np.pi * R[i] ** 2 * (d2 > RY[None, :] ** 2).all(1).mean())
    return np.array(out)


def element(name, resname):
    if resname in ION_RES:
        return resname
    return name[0]


def kabsch(Pm, Qm):
    """rotation R and translation so that (Pm @ R.T + t) ~ Qm"""
    pc, qc = Pm.mean(0), Qm.mean(0)
    H = (Pm - pc).T @ (Qm - qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    return R, qc - pc @ R.T


# ------------------------------------------------------------------ reference tables
sel = pd.read_csv(MD / 'selection/selection.tsv', sep='\t')
units33 = pd.read_csv(A33 / 'units.tsv', sep='\t', dtype=str)
pm13 = pd.read_csv(A44 / 'posmap.tsv', sep='\t', dtype=str)
pm14 = pd.read_csv(A45 / 'posmap_GH14_GH15.tsv', sep='\t', dtype=str)
PM = pd.concat([pm13, pm14])
PMAP = {(r.family, r.pdb, r.chain, str(r.seq)): r.position for r in PM.itertuples()}
PI = pd.read_csv(A45 / 'position_identity.tsv', sep='\t')
NCLF = PI.groupby('family').n_cl50.max().to_dict()
INVPOS = {(r.family, r.position) for r in PI.itertuples()
          if r.m50 >= 1.0 and r.n_cl50 >= min(3, NCLF.get(r.family, 3))}
CRYPOS = {(r.family, r.position) for r in PI.itertuples()}
CS = pd.read_csv(A45 / 'cell_summary.tsv', sep='\t', dtype={'subsite': str})
CAT = pd.read_csv(FR / 'catalytic_residues_mapped.tsv', sep='\t', dtype=str)

frame_rows, start_rows, bind_rows, unit_rows, sys_rows = [], [], [], [], []

for srow in sel.itertuples():
    S, pdb, fam, lasym = srow.system, srow.pdb, srow.family, srow.label_asym
    sd = MD / 'systems' / S
    mp = pd.read_csv(sd / 'map.tsv', sep='\t', dtype={'orig_resnum': str})
    ndx = np.array(re.search(r'\[ pocket \]\n([^\[]*)', (TR / f'{S}_pocket.ndx').read_text()).group(1).split(), int) - 1
    u0 = mda.Universe(str(sd / 'system.gro'))
    ag0 = u0.atoms[ndx]
    u = mda.Merge(ag0)
    u.load_new(str(TR / f'{S}_pocket.xtc'))
    A = u.atoms
    names = A.names.astype(str)
    resn = A.resnames.astype(str)
    resid = A.resids
    el = np.array([element(n, r) for n, r in zip(names, resn)])
    heavy = el != 'H'
    is_sug = np.isin(resn, list(SUGAR_RES))
    is_ion = np.isin(resn, list(ION_RES))
    is_cof = np.isin(resn, list(COFAC))
    is_prot = ~is_sug & ~is_ion & ~is_cof
    rmap = {int(r.new_resid): r for r in mp.itertuples()}
    # sugar units: residues of kind sugar, orig_resnum 1 = reducing end; ROH carries O1 of unit 1
    sug_res = mp[mp.kind == 'sugar']
    roh = mp[mp.kind == 'sugar_cap']
    lab33 = {r.seq: r.subsite if isinstance(r.subsite, str) else '' for r in
             units33[(units33.pdb == pdb) & (units33.chain == lasym)].itertuples()}
    X0 = ag0.positions.astype(float)  # start = system.gro coordinates of the pocket atoms, made whole around the glucan
    _s0 = np.where(np.isin(resn, list(SUGAR_RES)))[0][0]
    X0 = X0[_s0] + minimize_vectors((X0 - X0[_s0]).astype(np.float32), u0.dimensions).astype(float)
    units = []
    for r in sug_res.itertuples():
        m = resid == r.new_resid
        ff = r.ff_name
        lab = lab33.get(str(r.orig_resnum), '')
        oxy = {}
        for o in ['O2', 'O3', 'O4', 'O5', 'O6']:
            w = np.where(m & (names == o))[0]
            if len(w):
                oxy[o] = int(w[0])
        if str(r.orig_resnum) == '1' and len(roh):
            w = np.where((resid == int(roh.new_resid.iloc[0])) & (names == 'O1'))[0]
            if len(w):
                oxy['O1'] = int(w[0])
        # role of each oxygen: glycosidic if an O4/O6 bonded to a C1 of another unit (start structure)
        c1s = np.where(is_sug & (names == 'C1') & ~m)[0]
        role = {}
        for o, i in oxy.items():
            if o == 'O5':
                role[o] = 'ring'
            elif o in ('O4', 'O6') and len(c1s) and np.min(np.linalg.norm(X0[c1s] - X0[i], axis=1)) < 1.7:
                role[o] = 'glycosidic'
            else:
                role[o] = 'OH'
        ring = [int(np.where(m & (names == a))[0][0]) for a in ['C1', 'C2', 'C3', 'C4', 'C5', 'O5']]
        c1 = int(np.where(m & (names == 'C1'))[0][0])
        units.append(dict(seq=str(r.orig_resnum), resid=int(r.new_resid), ff=ff, subsite=lab, oxy=oxy, role=role,
                          ring=ring, c1=c1, anomer='beta' if ff.endswith('B') else 'alpha'))
    sug_heavy = np.where(is_sug & heavy)[0]
    # hydrogens bonded to each polar atom (start structure, 1.25 A)
    Hidx = np.where(el == 'H')[0]
    tH = cKDTree(X0[Hidx])
    polar = np.where(np.isin(el, ['N', 'O']) & ~is_ion)[0]
    hmap = {int(i): [int(Hidx[j]) for j in tH.query_ball_point(X0[i], 1.25)] for i in polar}
    prot_polar = np.array([i for i in polar if is_prot[i]])
    # partner annotation
    pinfo = {}
    for i in prot_polar:
        r = rmap[int(resid[i])]
        key = (fam, pdb, r.orig_chain, str(r.orig_resnum))
        pos = PMAP.get(key, '')
        part = 'mc' if names[i] in BACKBONE else 'sc'
        if pos:
            idc = 'inv' if (fam, pos) in INVPOS else 'var'
        elif fam in ('GH13', 'GH14'):
            idc = 'var'   # frame chain residue without a structural equivalent in the reference
        else:
            idc = 'unk'
        pinfo[int(i)] = dict(res=f'{r.orig_name}{r.orig_resnum}', pos=pos, part=part, idc=idc, atom=names[i])
    # catalytic midpoint
    cat = CAT[(CAT.pdb == pdb) & (CAT.chain == 'A') & CAT.role.isin(['nuc_or_base', 'acid'])]
    cat_atoms = []
    for c in cat.itertuples():
        nr = mp[(mp.orig_chain == 'A') & (mp.orig_resnum == str(c.auth_seq)) & (mp.kind == 'protein')]
        if not len(nr):
            continue
        rr = int(nr.new_resid.iloc[0])
        w = np.where((resid == rr) & np.isin(names, ['OD1', 'OD2', 'OE1', 'OE2', 'NE2' if c.resname == 'GLN' else 'XX']))[0]
        if len(w):
            cat_atoms.append(w)
    has_cat = len(cat_atoms) == 2
    # glucans left unlabelled by 33 although the catalytic pair is present (7DCG: acid position is Gln, so 33 had no
    # carboxylate midpoint): same rule on the start structure, Gln OE1/NE2 centroid in place of a carboxylate;
    # linear alpha(1-4) chains only (seq 1 = reducing end)
    if has_cat and not any(uu['subsite'] for uu in units) and all(uu['ff'][0] in '04' for uu in units):
        mid0_ = np.mean([X0[w].mean(0) for w in cat_atoms], axis=0)
        dd = [np.linalg.norm(X0[uu['c1']] - mid0_) for uu in units]
        k = int(np.argmin(dd))
        if dd[k] <= D_MAX:
            for j, uu in enumerate(units):
                q = j - k
                uu['subsite'] = '-1' if q == 0 else (f'+{-q}' if q < 0 else f'-{q + 1}')
                uu['subsite_source'] = 'a47 (catalytic-pair rule on system.gro)'
    for uu in units:
        uu.setdefault('subsite_source', '33' if uu['subsite'] else '')
    # radii for SASA
    Rr = np.array([RADII.get(e, 1.8) for e in el]) + PROBE
    cplx_mask = heavy.copy()
    lig_mask = is_sug & heavy
    ca_idx = np.where(is_prot & (names == 'CA'))[0]
    all_ox = [(uu, o, i) for uu in units for o, i in uu['oxy'].items()]

    def analyse(X, tag, fi, t):
        rows = []
        for uu, o, i in all_ox:
            d = np.linalg.norm(X[prot_polar] - X[i], axis=1)
            close = prot_polar[d <= HB]
            parts, hb_any = [], False
            for j, dj in zip(close, d[d <= HB]):
                ok = False
                # protein donor -> glucan O
                for h in hmap.get(int(j), []):
                    v1, v2 = X[j] - X[h], X[i] - X[h]
                    if np.degrees(np.arccos(np.clip(v1 @ v2 / np.linalg.norm(v1) / np.linalg.norm(v2), -1, 1))) >= ANG:
                        ok = True
                # glucan OH donor -> protein O or H-less N
                if not ok and uu['role'][o] == 'OH' and (el[j] == 'O' or not hmap.get(int(j))):
                    for h in hmap.get(int(i), []):
                        v1, v2 = X[i] - X[h], X[j] - X[h]
                        if np.degrees(np.arccos(np.clip(v1 @ v2 / np.linalg.norm(v1) / np.linalg.norm(v2), -1, 1))) >= ANG:
                            ok = True
                hb_any |= ok
                pi = pinfo[int(j)]
                parts.append((pi['pos'] or pi['res'], pi['res'], pi['atom'], pi['part'], pi['idc'], round(float(dj), 2), ok))
            rows.append(dict(system=S, family=fam, pdb=pdb, frame=fi, time_ps=t, unit_seq=uu['seq'],
                             subsite=uu['subsite'], label=o, role=uu['role'][o], held=len(parts) > 0, hb_angle=hb_any,
                             n_polar=len(parts),
                             sc_inv=any(p[3] == 'sc' and p[4] == 'inv' for p in parts),
                             sc_var=any(p[3] == 'sc' and p[4] in ('var', 'unk') for p in parts),
                             mc=any(p[3] == 'mc' for p in parts),
                             partners=' '.join(sorted({f'{p[0]}[{p[3]}]' for p in parts})),
                             partners_hb=' '.join(sorted({f'{p[0]}[{p[3]}]' for p in parts if p[6]})),
                             partner_atoms=' '.join(f'{p[1]}.{p[2]}({p[5]})' for p in parts)))
        oi = [i for _, _, i in all_ox]
        sc_ = sasa(oi, X, Rr, cplx_mask)
        sl_ = sasa(oi, X, Rr, lig_mask)
        for r, a, b in zip(rows, sc_, sl_):
            r['sasa_cplx'] = round(float(a), 2)
            r['sasa_lig'] = round(float(b), 2)
            r['rel'] = round(float(a / b), 3) if b > 0.5 else np.nan
        return rows

    # start structure (deposited + hydrogens)
    start_rows += analyse(X0, 'start', -1, -1.0)
    ref_ca = X0[ca_idx]
    ref_sug = X0[sug_heavy]
    ref_cent = {uu['seq']: X0[uu['ring']].mean(0) for uu in units}
    nfr = 0
    for ts in u.trajectory:
        X = A.positions.astype(float)
        box = ts.dimensions
        ref = X[sug_heavy[0]].copy()
        X = ref + minimize_vectors((X - ref).astype(np.float32), box).astype(float)
        R, tt = kabsch(X[ca_idx], ref_ca)
        X = X @ R.T + tt
        frame_rows += analyse(X, 'md', ts.frame, float(ts.time))
        rmsd = float(np.sqrt(((X[sug_heavy] - ref_sug) ** 2).sum(1).mean()))
        b = dict(system=S, frame=ts.frame, time_ps=float(ts.time), glucan_rmsd_A=round(rmsd, 2),
                 pocket_ca_rmsd_A=round(float(np.sqrt(((X[ca_idx] - ref_ca) ** 2).sum(1).mean())), 2))
        mid = np.mean([X[w].mean(0) for w in cat_atoms], axis=0) if has_cat else None
        for uu in units:
            b[f'disp_{uu["seq"]}'] = round(float(np.linalg.norm(X[uu['ring']].mean(0) - ref_cent[uu['seq']])), 2)
            if mid is not None:
                b[f'dmid_{uu["seq"]}'] = round(float(np.linalg.norm(X[uu['c1']] - mid)), 2)
        bind_rows.append(b)
        nfr += 1
    mid0 = np.mean([X0[w].mean(0) for w in cat_atoms], axis=0) if has_cat else None
    for uu in units:
        unit_rows.append(dict(system=S, family=fam, pdb=pdb, unit_seq=uu['seq'], md_resid=uu['resid'], ff=uu['ff'],
                              anomer=uu['anomer'], subsite=uu['subsite'], subsite_source=uu['subsite_source'],
                              d_c1_catmid_start_A=round(float(np.linalg.norm(X0[uu['c1']] - mid0)), 2) if has_cat else np.nan,
                              oxygens=' '.join(f'{o}:{uu["role"][o]}' for o in sorted(uu['oxy']))))
    sys_rows.append(dict(system=S, family=fam, pdb=pdb, glucan=srow.glucan, pocket_atoms=len(ndx), frames=nfr,
                         t_first_ns=0.0, t_last_ns=round(float(u.trajectory[-1].time) / 1000, 1),
                         catalytic_pair_found=has_cat))
    print(S, nfr, 'frames', u.trajectory[-1].time / 1000, 'ns', flush=True)

F = pd.DataFrame(frame_rows)
S0 = pd.DataFrame(start_rows)
B = pd.DataFrame(bind_rows)
UN = pd.DataFrame(unit_rows)
SY = pd.DataFrame(sys_rows)
F.to_csv(OUT / 'frames_oxygen.tsv.gz', sep='\t', index=False)
S0.to_csv(OUT / 'start_oxygen.tsv', sep='\t', index=False)
B.to_csv(OUT / 'binding_timeseries.tsv', sep='\t', index=False)
UN.to_csv(OUT / 'units.tsv', sep='\t', index=False)


# ------------------------------------------------------------------ cell summaries
def state_of(held, inv, var, mco, relm, sasam, mapped):
    if held >= 0.5:
        if not mapped:
            return 'SC' if var >= mco else 'MC'
        if inv >= 0.5:
            return 'INV'
        return 'VAR' if var >= mco else 'MC'
    return 'EXP' if (relm >= REL_EXPOSED and sasam >= SASA_EXPOSED) else 'BUR'


rows = []
for (s, useq, lab), g in F.groupby(['system', 'unit_seq', 'label'], sort=False):
    fam = g.family.iloc[0]
    mapped = fam in ('GH13', 'GH14')
    n = len(g)
    held = g.held.mean()
    inv = g.sc_inv.mean()
    var = (g.sc_var & ~g.sc_inv).mean()
    mco = (g.mc & ~g.sc_inv & ~g.sc_var).mean()
    relm, sasam = float(np.nanmedian(g.rel)), float(g.sasa_cplx.median())
    pc = {}
    for p in g.partners:
        for q in set(p.split()):
            pc[q] = pc.get(q, 0) + 1
    top = sorted(pc.items(), key=lambda x: -x[1])[:4]
    h = n // 2
    st0 = S0[(S0.system == s) & (S0.unit_seq == useq) & (S0.label == lab)].iloc[0]
    rows.append(dict(system=s, family=fam, pdb=g.pdb.iloc[0], unit_seq=useq, subsite=g.subsite.iloc[0], label=lab,
                     role=g.role.iloc[0], n_frames=n, occ_held=round(held, 3), occ_hb_angle=round(g.hb_angle.mean(), 3),
                     occ_held_1st_half=round(g.held.iloc[:h].mean(), 3), occ_held_2nd_half=round(g.held.iloc[h:].mean(), 3),
                     occ_sc_inv=round(inv, 3), occ_sc_var=round(var, 3), occ_mc_only=round(mco, 3),
                     median_sasa_cplx=round(sasam, 1), median_rel=round(relm, 3),
                     md_state=state_of(held, inv, var, mco, relm, sasam, mapped),
                     main_partners='; '.join(f'{q} {c / n:.2f}' for q, c in top),
                     start_held=bool(st0.held), start_partners=st0.partners, start_sasa=st0.sasa_cplx,
                     start_rel=st0.rel,
                     start_state=state_of(float(st0.held), float(st0.sc_inv), float(st0.sc_var and not st0.sc_inv),
                                          float(st0.mc and not st0.sc_inv and not st0.sc_var),
                                          st0.rel if st0.rel == st0.rel else 0.0, st0.sasa_cplx, mapped)))
C = pd.DataFrame(rows)
# crystal comparison
for setn in ['all_glucoselike', 'glc_units']:
    cs = CS[CS.set == setn].set_index(['family', 'subsite', 'label'])
    C[f'crystal_{setn}'] = [cs.state.get((f, s, l), '') if s in SUBS else '' for f, s, l in
                            zip(C.family, C.subsite, C.label)]
    C[f'crystal_{setn}_held'] = [cs.frac_cl_held.get((f, s, l), np.nan) if s in SUBS else np.nan for f, s, l in
                                 zip(C.family, C.subsite, C.label)]
C['same_state_all'] = [(a == b) if b else np.nan for a, b in zip(C.md_state, C.crystal_all_glucoselike)]
C['same_held_all'] = [((a >= 0.5) == (b >= 0.5)) if b == b else np.nan for a, b in
                      zip(C.occ_held, C.crystal_all_glucoselike_held)]
C.to_csv(OUT / 'md_cells.tsv', sep='\t', index=False)

# binding summary per system
bs = []
for s, g in B.groupby('system', sort=False):
    d = dict(system=s, frames=len(g), t_last_ns=g.time_ps.max() / 1000,
             glucan_rmsd_median=g.glucan_rmsd_A.median(), glucan_rmsd_last2ns=g.glucan_rmsd_A.iloc[-20:].median(),
             glucan_rmsd_max=g.glucan_rmsd_A.max(), pocket_ca_rmsd_median=g.pocket_ca_rmsd_A.median())
    for c in g.columns:
        if c.startswith('disp_'):
            d[f'{c}_median'] = g[c].median()
            d[f'{c}_last2ns'] = g[c].iloc[-20:].median()
        if c.startswith('dmid_'):
            d[f'{c}_median'] = g[c].median()
            d[f'{c}_frac_le_4.5'] = (g[c] <= D_MAX).mean()
    bs.append(d)
BS = pd.DataFrame(bs)
BS.to_csv(OUT / 'binding_summary.tsv', sep='\t', index=False)
SY.to_csv(OUT / 'systems.tsv', sep='\t', index=False)

# GH13 units without a subsite label: nearest crystal subsite centroid in the family frame (33 frame coordinates)
u33 = units33[(units33.family == 'GH13') & (units33.in_tables == 'True')].copy()
for c in ['x_frame', 'y_frame', 'z_frame']:
    u33[c] = u33[c].astype(float)
cent = u33[u33.subsite.isin(SUBS)].groupby('subsite')[['x_frame', 'y_frame', 'z_frame']].mean()
nr = []
for srow in sel[sel.family.isin(['GH13'])].itertuples():
    q = units33[(units33.pdb == srow.pdb) & (units33.chain == srow.label_asym)]
    for r in q.itertuples():
        if r.x_frame != r.x_frame or r.x_frame is None:
            continue
        x = np.array([float(r.x_frame), float(r.y_frame), float(r.z_frame)])
        dd = np.linalg.norm(cent.values - x, axis=1)
        k = int(np.argmin(dd))
        nr.append(dict(system=srow.system, pdb=srow.pdb, seq=r.seq, subsite_33=r.subsite if isinstance(r.subsite, str) else '',
                       d_catmid_33=r.d_ref_to_midpoint_A, nearest_crystal_subsite=cent.index[k],
                       d_to_centroid_A=round(float(dd[k]), 2)))
pd.DataFrame(nr).to_csv(OUT / 'unlabelled_units_nearest_subsite.tsv', sep='\t', index=False)
print('done')
