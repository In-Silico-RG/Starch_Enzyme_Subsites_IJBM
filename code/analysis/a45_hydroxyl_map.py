"""45 - Atom-level map of the glucose oxygens held by the protein at subsites -3..+3 (GH13, GH31; GH14, GH15 apart).

Question (thesis on enzymatic processing of chemically substituted starch; substituents sit on O2, O3, O6): at which
subsite and which oxygen is there a protein hydrogen-bond partner, is that partner a lineage-invariant residue, and
where is an oxygen not held and open to solvent.

Inputs (read only):
  data/analysis/33_subsites/units.tsv          sugar ring units with subsite, RSCC, frame chain
  data/analysis/44_rule_test/posmap.tsv         residue -> reference position (Calpha within 2.0 A), GH13/GH31
  data/analysis/44_rule_test/clusters_bound.tsv 50 % / 30 % clusters of ligand-bound GH13/GH31 entries
  data/analysis/44_rule_test/mutations.tsv      engineered mutations (dropped from identity votes)
  data/analysis/44_rule_test/C_identity_T20.tsv modal fractions of a44 (cross-check)
  data/calc/03_frames/superposition.tsv         frame matrices (used here to map GH14/GH15 residues)
  data/calc/03_frames/catalytic_residues_mapped.tsv
  data/census/cif/*.cif

Unit selection: as a44 (in_tables, RSCC >= 0.80, subsite -3..+3), families GH13, GH31, GH14, GH15.
Atom labels: assigned from the bonded graph of the ring, not from atom names (multi-ring residues such as AC1, DAF,
AAO carry suffixed names): C1 = the unit's reference atom from 33 (anomeric or pseudo-anomeric carbon); the ring
neighbour of C1 that is O or N is position 5 (O5); in an all-carbon ring it is the neighbour without an exocyclic
heteroatom; the ring is then walked C1-C2-C3-C4-C5; C6 = exocyclic carbon on C5; O2/O3/O4 = exocyclic heteroatom
(O, N, F or S) on C2/C3/C4; O6 = heteroatom on C6; O1 = exocyclic heteroatom on C1 that has no other heavy-atom
neighbour (a free anomeric OH); a C1 heteroatom bonded to another ring is the glycosidic atom and is counted as
O4 of the acceptor unit, as for glycosidic bonds between residues.
Configuration check: sign of the volume (ring-previous, ring-next, exocyclic heavy substituent) at C2, C3, C4, C5
compared with alpha-D-glucose (GLC); a unit with a sign opposite to GLC at any of C2-C4 is an epimer and is excluded.
Contacts: protein heavy atoms of the deposited asymmetric unit (symmetry copies not generated, waters not partners);
polar partner = protein N or O within 3.5 A of the sugar atom; contact = any protein heavy atom within 4.0 A.
SASA: Shrake-Rupley, 480 points per atom, probe 1.4 A, Bondi radii, hydrogens absent, waters removed;
  sasa_cplx  = with all non-water heavy atoms of the asymmetric unit (protein, all ligands, other sugar units)
  sasa_lig   = with all non-protein, non-water heavy atoms only (the ligand without the protein)
  rel        = sasa_cplx / sasa_lig  (share of the atom's ligand-state surface the protein leaves open)
Identity: residue at the reference position in the entry's frame chain; one vote per 50 % cluster (modal residue
of the cluster's entries, ties share the vote, as a44b); over all ligand-bound entries of the family; mutants dropped.
Run: ~/miniconda3/bin/python code/analysis/a45_hydroxyl_map.py
"""
import itertools, re
from collections import defaultdict, Counter
from pathlib import Path
import numpy as np
import pandas as pd
import gemmi
import networkx as nx
from scipy.spatial import cKDTree
from Bio import Align
from Bio.Align import substitution_matrices

P = Path(__file__).resolve().parents[2]
S33 = P / 'data/analysis/33_subsites'
D44 = P / 'data/analysis/44_rule_test'
FR = P / 'data/calc/03_frames'
CIF = P / 'data/census/cif'
OUT = P / 'data/analysis/45_hydroxyl_map'
OUT.mkdir(parents=True, exist_ok=True)

FAMS = ['GH13', 'GH31', 'GH14', 'GH15']
SUBS = ['-3', '-2', '-1', '+1', '+2', '+3']
OXY = ['O1', 'O2', 'O3', 'O4', 'O5', 'O6']
CARB = ['C1', 'C2', 'C3', 'C4', 'C5', 'C6']
REF = {'GH13': ('1AMY', 'A'), 'GH31': ('2QMJ', 'A'), 'GH14': ('1BYA', 'A'), 'GH15': ('1AGM', 'A')}
CAT = {'GH13': {'D179', 'E204', 'D289'}, 'GH31': {'D443', 'D542'}, 'GH14': {'E380', 'E186'},
       'GH15': {'E400', 'E179'}}
HB, CONT, BUR = 3.5, 4.0, 4.5
REL_EXPOSED = 0.25
SASA_EXPOSED = 5.0   # A^2, median SASA in the complex; both criteria needed for 'exposed'
RADII = {'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80, 'F': 1.47, 'P': 1.80, 'CL': 1.75, 'BR': 1.85}
BACKBONE = {'N', 'O', 'OXT', 'CA', 'C'}
AROM_SC = {'PHE': {'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'}, 'TYR': {'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'},
           'TRP': {'CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2', 'CG', 'CD1', 'NE1'},
           'HIS': {'CG', 'ND1', 'CD2', 'CE1', 'NE2'}}
NONGLC_NAMED = {'GAL': 'galactose (named non-gluco; excluded by name)', 'MAN': 'mannose (named non-gluco)',
                '5DI': 'ketose (2,5-diulose) fluoride, C1 exocyclic', 'AFR': 'ketose (2,5-diulose) fluoride',
                '93Z': 'sulfated cyclitol, no unambiguous O5 equivalent'}
ONE = dict(ALA='A', ARG='R', ASN='N', ASP='D', CYS='C', GLN='Q', GLU='E', GLY='G', HIS='H', ILE='I', LEU='L',
           LYS='K', MET='M', PHE='F', PRO='P', SER='S', THR='T', TRP='W', TYR='Y', VAL='V')


def is_aa(name):
    t = gemmi.find_tabulated_residue(name)
    return bool(t and t.is_amino_acid())


def one(name):
    t = gemmi.find_tabulated_residue(name)
    return t.one_letter_code.upper() if t else 'X'


def read(pdb):
    st = gemmi.read_structure(str(CIF / f'{pdb}.cif'))
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    st.setup_entities()
    return st


# ------------------------------------------------------------------ sphere points for SASA
def sphere(n=480):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    th = np.pi * (1 + 5 ** 0.5) * i
    return np.c_[np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)]


SPH = sphere()


def sasa(idx, X, R, tree, mask):
    """SASA of atoms idx against atoms where mask is True (X coords, R radii incl. probe)."""
    out = []
    for i in idx:
        pts = X[i] + SPH * R[i]
        nb = [j for j in tree.query_ball_point(X[i], R[i] + R.max()) if j != i and mask[j]]
        if not nb:
            out.append(4 * np.pi * R[i] ** 2)
            continue
        Y, RY = X[nb], R[nb]
        d2 = ((pts[:, None, :] - Y[None, :, :]) ** 2).sum(-1)
        free = (d2 > RY[None, :] ** 2).all(1)
        out.append(4 * np.pi * R[i] ** 2 * free.mean())
    return np.array(out)


# ------------------------------------------------------------------ ring labelling
def res_graph(res):
    G = nx.Graph()
    at = {a.name: a for a in res}
    names = list(at)
    for n in names:
        G.add_node(n, el=at[n].element.name.upper())
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            lim = 1.95 if 'S' in (G.nodes[a]['el'], G.nodes[b]['el']) else 1.75
            if at[a].pos.dist(at[b].pos) < lim:
                G.add_edge(a, b)
    return G, at


HET = {'O', 'N', 'F', 'S'}
EXT_EL = {}


def label_unit(res, ref, ext=None):
    ext = ext or {}
    """-> (labels dict pos->atom name, info dict) or (None, reason)."""
    G, at = res_graph(res)
    if ref not in G:
        return None, 'reference atom absent'
    rings = [list(c) for c in nx.minimum_cycle_basis(G) if 5 <= len(c) <= 7]
    ring = [r for r in rings if ref in r]
    if not ring:
        return None, 'reference atom not in a ring'
    ring = set(ring[0])
    if len(ring) != 6:
        return None, f'{len(ring)}-membered ring'
    el = lambda n: G.nodes[n]['el']
    exo = lambda n: [m for m in G[n] if m not in ring]
    exo_het = lambda n: [m for m in exo(n) if el(m) in HET]
    nb = [m for m in G[ref] if m in ring]
    if len(nb) != 2:
        return None, 'C1 ring neighbours != 2'
    a, b = nb
    if el(a) in ('O', 'N') and el(b) == 'C':
        p5, c2 = a, b
    elif el(b) in ('O', 'N') and el(a) == 'C':
        p5, c2 = b, a
    elif el(a) == 'C' and el(b) == 'C':
        ha, hb = bool(exo_het(a)), bool(exo_het(b))
        if ha == hb:
            return None, 'carba ring: O5 equivalent ambiguous'
        p5, c2 = (b, a) if ha else (a, b)
        if any(el(m) == 'C' for m in exo(p5)):
            return None, 'carba ring: O5 equivalent carries a carbon substituent'
    else:
        return None, 'C1 neighbours both heteroatoms'
    order = [ref, c2]
    while len(order) < 5:
        nxt = [m for m in G[order[-1]] if m in ring and m not in order and m != p5]
        if len(nxt) != 1:
            return None, 'ring walk failed'
        order.append(nxt[0])
    if p5 not in G[order[-1]]:
        return None, 'ring walk did not close at position 5'
    L = {f'C{k + 1}': order[k] for k in range(5)}
    L['O5'] = p5  # may be N (imino) or C (carba); element kept in info
    c6 = [m for m in exo(L['C5']) if el(m) == 'C']
    if c6:
        L['C6'] = c6[0]
        o6 = [m for m in G[c6[0]] if m != L['C5'] and el(m) in HET]
        if o6:
            L['O6'] = o6[0]
        elif [x for x in ext.get(c6[0], []) if x[4] in HET]:
            x = [x for x in ext[c6[0]] if x[4] in HET][0]
            L['O6'] = 'EXT|' + '|'.join(x[:4])
            EXT_EL[L['O6']] = x[4]
    for k in (2, 3, 4):
        h = exo_het(L[f'C{k}'])
        if h:
            L[f'O{k}'] = h[0]
        elif [x for x in ext.get(L[f'C{k}'], []) if x[4] in HET]:
            x = [x for x in ext[L[f'C{k}']] if x[4] in HET][0]
            L[f'O{k}'] = 'EXT|' + '|'.join(x[:4])
            EXT_EL[L[f'O{k}']] = x[4]
    o1 = [m for m in exo_het(ref) if len([x for x in G[m] if x != ref]) == 0 and not ext.get(m)]
    if o1:
        L['O1'] = o1[0]
    # extra exocyclic heteroatoms (F, S, a second substituent) -> modified
    expected = {L.get(k) for k in ('O1', 'O2', 'O3', 'O4', 'O6')}
    extra = []
    for n in order + [p5]:
        for m in exo(n):
            if el(m) in HET and m not in expected and not (n == ref):
                extra.append(m)
    # chirality signs at C2..C5
    signs = {}
    for k in (2, 3, 4, 5):
        c = L[f'C{k}']
        prev = L[f'C{k - 1}']
        nxt = L[f'C{k + 1}'] if k < 5 else p5
        sub = [L['C6']] if k == 5 and 'C6' in L else ([L[f'O{k}']] if k < 5 and f'O{k}' in L and not L[f'O{k}'].startswith('EXT') else [])
        if not sub:
            continue
        v = [np.array(at[x].pos.tolist()) - np.array(at[c].pos.tolist()) for x in (prev, nxt, sub[0])]
        vol = float(np.dot(v[0], np.cross(v[1], v[2])))
        if abs(vol) > 1.0:
            signs[k] = int(np.sign(vol))
    info = dict(ring_el=el(p5), extra=','.join(extra), signs=signs,
                elements={k: (EXT_EL[v] if v.startswith('EXT|') else el(v)) for k, v in L.items()})
    return L, info


# ------------------------------------------------------------------ units
u = pd.read_csv(S33 / 'units.tsv', sep='\t', dtype=str, keep_default_na=False)
u['rscc_f'] = pd.to_numeric(u.rscc, errors='coerce')
base = u[u.family.isin(FAMS) & (u.in_tables == 'True') & u.subsite.isin(SUBS)].copy()
units = base[base.rscc_f >= 0.80].copy()
excl_rscc = base[~(base.rscc_f >= 0.80)]
print('units', len(units), 'entries', units.pdb.nunique(), '| without RSCC >= 0.80:', len(excl_rscc))

# ------------------------------------------------------------------ position maps (GH13/GH31 from 44; GH14/GH15 here)
pm = pd.read_csv(D44 / 'posmap.tsv', sep='\t', dtype=str, keep_default_na=False)
mut = pd.read_csv(D44 / 'mutations.tsv', sep='\t', dtype=str, keep_default_na=False)
cb = pd.read_csv(D44 / 'clusters_bound.tsv', sep='\t', dtype=str, keep_default_na=False)
sup = pd.read_csv(FR / 'superposition.tsv', sep='\t')
sup = sup[sup.status.str.startswith('accepted')]
mats = {}
for r in sup.itertuples():
    v = [float(x) for x in str(r[-1]).split()]
    mats[(r.family, r.pdb, r.chain)] = (np.array(v[:9]).reshape(3, 3), np.array(v[9:12]))

frame_chain = units.groupby(['family', 'pdb']).frame_chain.agg(lambda s: Counter(s).most_common(1)[0][0]).to_dict()

refs = {}
for f in ['GH14', 'GH15']:
    rp, rc = REF[f]
    st = read(rp)
    lab, xyz = [], []
    for res in st[0][rc]:
        ca = res.find_atom('CA', '*')
        if ca is not None and is_aa(res.name):
            lab.append(f'{one(res.name)}{res.seqid.num}{res.seqid.icode.strip()}')
            xyz.append(ca.pos.tolist())
    refs[f] = (lab, np.array(xyz))

pm_rows, mut_rows, seqs = [], [], {}
for f in ['GH14', 'GH15']:
    for pdb in sorted(set(units[units.family == f].pdb)):
      st = read(pdb)
      for ch in sorted(c for (ff, pp, c) in mats if ff == f and pp == pdb):
        R, t = mats[(f, pdb, ch)]
        lab, RX = refs[f]
        best = {}
        for res in st[0][ch]:
            ca = res.find_atom('CA', '*')
            if ca is None or not is_aa(res.name):
                continue
            x = R @ np.array(ca.pos.tolist()) + t
            d = np.linalg.norm(RX - x, axis=1)
            k = int(d.argmin())
            if d[k] <= 2.0 and (lab[k] not in best or d[k] < best[lab[k]][0]):
                best[lab[k]] = (float(d[k]), res.name, res.seqid.num, res.seqid.icode.strip())
        for q, (d, rn, sq, ic) in best.items():
            pm_rows.append(dict(family=f, pdb=pdb, chain=ch, position=q, resname=rn, seq=str(sq), icode=ic,
                                d_ca=round(d, 2)))
        if ch == frame_chain[(f, pdb)]:
            seqs[(f, pdb)] = gemmi.one_letter_code([r.name for r in st[0][ch].get_polymer()]).upper()
      doc = gemmi.cif.read(str(CIF / f'{pdb}.cif'))
      tab = doc.sole_block().find('_struct_ref_seq_dif.', ['pdbx_pdb_strand_id', 'pdbx_auth_seq_num', 'mon_id',
                                                           'db_mon_id', 'details'])
      for row in tab:
          det = gemmi.cif.as_string(row[4])
          if 'mutation' in det.lower() and 'engineered' in det.lower():
              mut_rows.append(dict(family=f, pdb=pdb, chain=gemmi.cif.as_string(row[0]),
                                   seq=gemmi.cif.as_string(row[1]), mutant=gemmi.cif.as_string(row[2]),
                                   wild_type=gemmi.cif.as_string(row[3]), details=det))
pm2 = pd.DataFrame(pm_rows)
pm2.to_csv(OUT / 'posmap_GH14_GH15.tsv', sep='\t', index=False)
mut2 = pd.DataFrame(mut_rows, columns=['family', 'pdb', 'chain', 'seq', 'mutant', 'wild_type', 'details'])
mut2.to_csv(OUT / 'mutations_GH14_GH15.tsv', sep='\t', index=False)

# 50 % clusters for GH14/GH15 bound entries (same aligner settings and single linkage as a44a)
al = Align.PairwiseAligner()
al.mode = 'global'
al.substitution_matrix = substitution_matrices.load('BLOSUM62')
al.open_gap_score, al.extend_gap_score, al.end_gap_score = -10, -0.5, 0


def ident(a, b):
    c = al.align(a, b)[0].counts()
    n = c.identities + c.mismatches
    return (c.identities / n if n else 0.0, c.identities / min(len(a), len(b)))


cl_rows = []
for f in ['GH14', 'GH15']:
    names = sorted(p for (ff, p) in seqs if ff == f)
    par50 = {n: n for n in names}
    par30 = {n: n for n in names}

    def root(par, x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    for a_, b_ in itertools.combinations(names, 2):
        i50, i30 = ident(seqs[(f, a_)], seqs[(f, b_)])
        if i50 >= 0.50:
            par50[root(par50, a_)] = root(par50, b_)
        if i30 >= 0.30:
            par30[root(par30, a_)] = root(par30, b_)
    r50, r30 = {}, {}
    for n in names:
        cl_rows.append(dict(family=f, pdb=n, frame_chain=frame_chain[(f, n)], protein_key='',
                            cl50=f'{f}_b50_{r50.setdefault(root(par50, n), len(r50) + 1)}',
                            cl30=f'{f}_b30_{r30.setdefault(root(par30, n), len(r30) + 1)}'))
cb2 = pd.DataFrame(cl_rows)
cb2.to_csv(OUT / 'clusters_GH14_GH15.tsv', sep='\t', index=False)
CB = pd.concat([cb, cb2], ignore_index=True)
PM = pd.concat([pm, pm2.astype(str)], ignore_index=True)
MUTS = pd.concat([mut, mut2.astype(str)], ignore_index=True)
CL50 = {(r.family, r.pdb): r.cl50 for r in CB.itertuples()}
CL30 = {(r.family, r.pdb): r.cl30 for r in CB.itertuples()}
FCH = {(r.family, r.pdb): r.frame_chain for r in CB.itertuples()}
MUTSET = {(r.pdb, r.chain, str(r.seq)) for r in MUTS.itertuples()}

# residue -> position and position identity over all bound entries
POS = {}
RESAT = defaultdict(dict)    # (fam, pdb) -> position -> resname (frame chain, mutants removed)
MUTPOS = set()
for r in PM.itertuples():
    POS[(r.pdb, r.chain, str(r.seq), r.icode)] = r.position
    if FCH.get((r.family, r.pdb)) != r.chain:
        continue
    if (r.pdb, r.chain, str(r.seq)) in MUTSET:
        MUTPOS.add((r.family, r.pdb, r.position))
        continue
    RESAT[(r.family, r.pdb)][r.position] = r.resname


def mode_vote(values):
    c = Counter(values)
    top = max(c.values())
    w = [k for k, v in c.items() if v == top]
    return {k: 1 / len(w) for k in w}


ID_CACHE = {}


def identity(fam, q):
    if (fam, q) in ID_CACHE:
        return ID_CACHE[(fam, q)]
    byc = defaultdict(list)
    for (f, p), cl in CL50.items():
        if f != fam or (f, p, q) in MUTPOS:
            continue
        byc[cl].append(RESAT.get((f, p), {}).get(q, 'gap'))
    acc, gaps = Counter(), 0.0
    for c, v in byc.items():
        mv = mode_vote(v)
        if 'gap' in mv:
            gaps += mv['gap']
            mv = {k: x for k, x in mv.items() if k != 'gap'}
            if not mv:
                continue
            s = sum(mv.values())
            mv = {k: x / s for k, x in mv.items()}
        for k, x in mv.items():
            acc[k] += x
    n = sum(acc.values())
    if n == 0:
        res = dict(n_cl50=0, modal='', m50=np.nan, gaps=gaps, types='')
    else:
        mod = max(acc, key=acc.get)
        res = dict(n_cl50=round(n, 2), modal=ONE.get(mod, mod), m50=round(acc[mod] / n, 3), gaps=round(gaps, 2),
                   types=', '.join(f'{ONE.get(k, k)} {acc[k]:.3g}' for k in sorted(acc, key=lambda z: -acc[z])))
    ID_CACHE[(fam, q)] = res
    return res


# ------------------------------------------------------------------ per unit: labels, contacts, SASA
atom_rows, partner_rows, unit_rows, skip_rows = [], [], [], []
REF_SIGNS = None
cache = {}
for uu in units.sort_values(['family', 'pdb']).itertuples():
    if uu.pdb not in cache:
        cache.clear()
        cache[uu.pdb] = read(uu.pdb)[0]
        NS = {uu.pdb: gemmi.NeighborSearch(cache[uu.pdb], gemmi.UnitCell(), 3).populate()}
    m = cache[uu.pdb]
    ch = m.find_chain(uu.chain)
    res = None
    if ch is not None:
        for r in ch:
            if str(r.seqid.num) == str(uu.seq) and r.name == uu.comp:
                res = r
                break
    key = dict(family=uu.family, pdb=uu.pdb, instance_id=uu.instance_id, ring=uu.ring, comp=uu.comp, chain=uu.chain,
               seq=uu.seq, subsite=uu.subsite, rscc=uu.rscc, unit_frame_chain=uu.frame_chain)
    if res is None:
        skip_rows.append({**key, 'reason': 'residue not found'})
        continue
    if uu.comp in NONGLC_NAMED:
        skip_rows.append({**key, 'reason': NONGLC_NAMED[uu.comp]})
        continue
    refat = uu.ref_atom or ('C1' if res.find_atom('C1', '*') else '')
    ext = defaultdict(list)
    for a_ in res:
        for mk in NS[uu.pdb].find_atoms(a_.pos, '\0', radius=1.75):
            cra = mk.to_cra(m)
            if mk.image_idx != 0 or cra.residue.name in ('HOH', 'DOD') or is_aa(cra.residue.name):
                continue
            if cra.chain.name == ch.name and cra.residue.seqid == res.seqid:
                continue
            if cra.atom.pos.dist(a_.pos) < 1.75:
                ext[a_.name].append((cra.chain.name, str(cra.residue.seqid.num), cra.residue.name, cra.atom.name,
                                     cra.atom.element.name.upper()))
    L, info = label_unit(res, refat, ext)
    if L is None:
        skip_rows.append({**key, 'reason': info})
        continue
    unit_rows.append(dict(**key, ref_atom=refat, ring_el=info['ring_el'], extra=info['extra'],
                          signs=';'.join(f'C{k}:{v:+d}' for k, v in sorted(info['signs'].items())),
                          labels=' '.join(f'{k}={v}({info["elements"][k]})' for k, v in sorted(L.items())),
                          Lmap=L, uinfo=info))

U = pd.DataFrame(unit_rows)
# reference signs from GLC units (modal)
glc_signs = Counter()
for r in U[U.comp == 'GLC'].itertuples():
    glc_signs[tuple(sorted(r.uinfo['signs'].items()))] += 1
ref_sign = dict(glc_signs.most_common(1)[0][0])
print('GLC reference chirality signs', ref_sign, 'from', glc_signs.most_common(3))


def kind_of(r):
    s, inf, L = r.uinfo['signs'], r.uinfo, r.Lmap
    epi = [k for k in (2, 3, 4) if k in s and k in ref_sign and s[k] != ref_sign[k]]
    if epi:
        return 'epimer', 'opposite configuration at ' + ','.join(f'C{k}' for k in epi)
    if 5 in s and 5 in ref_sign and s[5] != ref_sign[5]:
        return 'epimer', 'opposite configuration at C5 (L-ido)'
    if inf['ring_el'] == 'C':
        return 'carba', 'all-carbon ring (cyclitol, e.g. valienamine)'
    if inf['ring_el'] == 'N':
        return 'imino', 'ring N (iminosugar)'
    el = inf['elements']
    full = all(el.get(k) == 'O' for k in ('O2', 'O3', 'O4', 'O6')) and not inf['extra']
    if full:
        return 'Glc', ''
    miss = [k for k in ('O2', 'O3', 'O4', 'O6') if k not in el]
    sub = [f'{k}={el[k]}' for k in ('O2', 'O3', 'O4', 'O6') if k in el and el[k] != 'O']
    return 'Glc-modified', ' '.join(['missing ' + ','.join(miss)] * bool(miss) + sub +
                                     (['extra ' + inf['extra']] if inf['extra'] else []))


kk = U.apply(kind_of, axis=1)
U['kind'] = [k[0] for k in kk]
U['kind_note'] = [k[1] for k in kk]
for r in U[U.kind == 'epimer'].itertuples():
    skip_rows.append(dict(family=r.family, pdb=r.pdb, instance_id=r.instance_id, ring=r.ring, comp=r.comp,
                          chain=r.chain, seq=r.seq, subsite=r.subsite, rscc=r.rscc, reason='epimer: ' + r.kind_note))
U = U[U.kind != 'epimer'].copy()

# entry chemistry at -3..+3 (all in_tables units, any RSCC): glucose-only?
kind_by_comp = {}
for r in U.itertuples():
    kind_by_comp.setdefault(r.comp, Counter())[r.kind] += 1
ent_glc_only = {}
for (f, p), g in base.groupby(['family', 'pdb']):
    ok = True
    for r in g.itertuples():
        hit = U[(U.pdb == p) & (U.instance_id == r.instance_id) & (U.ring == r.ring)]
        if len(hit):
            ok &= hit.kind.iloc[0] == 'Glc'
        else:
            ok &= r.comp in ('GLC', 'BGC')   # unit not analysed (RSCC < 0.80): judged by component
    ent_glc_only[(f, p)] = ok

# contacts and SASA
cache = {}
for r in U.itertuples():
    if r.pdb not in cache:
        st = read(r.pdb)
        cache.clear()
        mm = st[0]
        Xl, Rl, metal = [], [], []
        for ch in mm:
            for res in ch:
                if res.name in ('HOH', 'DOD', 'WAT'):
                    continue
                prot = is_aa(res.name)
                for a in res:
                    e = a.element.name.upper()
                    Xl.append(a.pos.tolist())
                    Rl.append(RADII.get(e, 1.80) + 1.4)
                    metal.append((ch.name, str(res.seqid.num), res.seqid.icode.strip(), res.name, a.name, e, prot))
        X, R, meta = np.array(Xl), np.array(Rl), metal
        PROT = np.array([x[6] for x in meta])
        tree = cKDTree(X)
        idx_of = {(x[0], x[1], x[2], x[3], x[4]): i for i, x in enumerate(meta)}
        cache[r.pdb] = (mm, X, R, meta, PROT, tree, idx_of)
    m, X, R, meta, PROT, tree, idx_of = cache[r.pdb]
    L = r.Lmap
    ids = {}
    for lab, an in L.items():
        if an.startswith('EXT|'):
            _, c_, sq_, cp_, nm_ = an.split('|')
            k = (c_, sq_, '', cp_, nm_)
        else:
            k = (r.chain, str(r.seq), '', r.comp, an)
        if k in idx_of:
            ids[lab] = idx_of[k]
    labs = list(ids)
    ii = [ids[k] for k in labs]
    s_c = sasa(ii, X, R, tree, np.ones(len(X), bool))
    s_l = sasa(ii, X, R, tree, ~PROT)
    fch = r.unit_frame_chain   # chain whose catalytic pair defines this unit's subsite (33)
    for lab, i, sc, sl in zip(labs, ii, s_c, s_l):
        nb = tree.query_ball_point(X[i], BUR)
        pol, anyc, nbur, posset_pol, posset_any = [], [], 0, set(), set()
        for j in nb:
            if not PROT[j]:
                continue
            d = float(np.linalg.norm(X[j] - X[i]))
            c_, sq, ic, rn, an, e, _ = meta[j]
            nbur += 1
            if d > CONT:
                continue
            q = POS.get((r.pdb, c_, sq, ic), '') if c_ == fch else ''
            part = 'mc' if an in BACKBONE else 'sc'
            arom = an in AROM_SC.get(rn, ())
            ispol = e in ('N', 'O') and d <= HB
            partner_rows.append(dict(family=r.family, pdb=r.pdb, instance_id=r.instance_id, ring=r.ring, comp=r.comp,
                                     kind=r.kind, subsite=r.subsite, label=lab, ligand_atom=L[lab],
                                     ligand_el=r.uinfo['elements'][lab], p_chain=c_, p_seq=sq, p_icode=ic, p_resname=rn,
                                     p_atom=an, p_el=e, dist=round(d, 2), part=part, polar_hb=ispol, aromatic_atom=arom,
                                     frame_chain=c_ == fch, position=q))
        atom_rows.append(dict(family=r.family, pdb=r.pdb, cl50=CL50.get((r.family, r.pdb), ''),
                              cl30=CL30.get((r.family, r.pdb), ''), instance_id=r.instance_id, ring=r.ring,
                              comp=r.comp, kind=r.kind, subsite=r.subsite, rscc=r.rscc,
                              glc_only_entry=ent_glc_only.get((r.family, r.pdb), False), label=lab, atom=L[lab],
                              element=r.uinfo['elements'][lab], sasa_cplx=round(float(sc), 2),
                              sasa_lig=round(float(sl), 2),
                              rel=round(float(sc / sl), 3) if sl > 0.5 else np.nan, n_prot_4p5=nbur))

A = pd.DataFrame(atom_rows)
PT = pd.DataFrame(partner_rows)
U.drop(columns=['Lmap', 'uinfo']).to_csv(OUT / 'units_labelled.tsv', sep='\t', index=False)
pd.DataFrame(skip_rows).to_csv(OUT / 'units_excluded.tsv', sep='\t', index=False)
PT.to_csv(OUT / 'partners_4A.tsv', sep='\t', index=False)

# per atom: polar partners, invariant side-chain partners
INVM = 1.0
NCL_FAM = CB.groupby('family').cl50.nunique().to_dict()
NMIN = {f: min(3, NCL_FAM.get(f, 3)) for f in FAMS}
print('50% clusters per family (bound entries)', NCL_FAM, 'NMIN', NMIN)


def pos_info(fam, q):
    if not q:
        return dict(n_cl50=0, modal='', m50=np.nan, gaps=0, types='')
    return identity(fam, q)


pol = PT[PT.polar_hb]
g = pol.groupby(['pdb', 'instance_id', 'ring', 'label'])
held = g.size().rename('n_polar')
A = A.merge(held, left_on=['pdb', 'instance_id', 'ring', 'label'], right_index=True, how='left')
A['n_polar'] = A.n_polar.fillna(0).astype(int)
A['held'] = A.n_polar > 0


def partner_summary(df):
    out = []
    for r in df.itertuples():
        out.append(f'{r.position or (r.p_resname + r.p_seq + ":" + r.p_chain)}.{r.p_atom}({r.dist})')
    return ' '.join(sorted(set(out)))


desc = pol.groupby(['pdb', 'instance_id', 'ring', 'label']).apply(partner_summary).rename('polar_partners')
A = A.merge(desc, left_on=['pdb', 'instance_id', 'ring', 'label'], right_index=True, how='left')


def inv_flags(df):
    sc_inv = sc_var = mc = other = False
    for r in df.itertuples():
        if r.part == 'mc':
            mc = True
            continue
        if not r.frame_chain:
            other = True
            continue
        if not r.position:          # frame chain, no structural equivalent in the reference -> variable
            sc_var = True
            continue
        info = identity(r.family, r.position)
        if info['n_cl50'] >= NMIN[r.family] and info['m50'] >= INVM:
            sc_inv = True
        else:
            sc_var = True
    return pd.Series(dict(sc_inv=sc_inv, sc_var=sc_var, mc=mc, other_chain=other))


fl = pol.groupby(['pdb', 'instance_id', 'ring', 'label']).apply(inv_flags)
A = A.merge(fl, left_on=['pdb', 'instance_id', 'ring', 'label'], right_index=True, how='left')
for c in ['sc_inv', 'sc_var', 'mc', 'other_chain']:
    A[c] = A[c].fillna(False).astype(bool)
anyc = PT.groupby(['pdb', 'instance_id', 'ring', 'label']).size().rename('n_any_4A')
A = A.merge(anyc, left_on=['pdb', 'instance_id', 'ring', 'label'], right_index=True, how='left')
A['n_any_4A'] = A.n_any_4A.fillna(0).astype(int)
A.to_csv(OUT / 'atoms.tsv', sep='\t', index=False)

# ------------------------------------------------------------------ position identity table
prow = []
for (f, q), gg in pol[pol.position != ''].groupby(['family', 'position']):
    info = identity(f, q)
    prow.append(dict(family=f, position=q, catalytic=q in CAT[f], **info,
                     labels=' '.join(sorted(set(gg.subsite + ':' + gg.label))),
                     parts=' '.join(sorted(set(gg.part)))))
PI = pd.DataFrame(prow)
c20 = pd.read_csv(D44 / 'C_identity_T20.tsv', sep='\t')
PI = PI.merge(c20[['family', 'position', 'm50', 'n_cl50']].rename(columns={'m50': 'm50_a44_T20', 'n_cl50': 'n_cl50_a44'}),
              on=['family', 'position'], how='left')
PI.to_csv(OUT / 'position_identity.tsv', sep='\t', index=False)


# ------------------------------------------------------------------ aggregation per subsite x atom
def agg(Asub, Psub, tag):
    rows = []
    for (f, s, lab), gg in Asub.groupby(['family', 'subsite', 'label']):
        # entry level: share of the entry's units with the property; entry counts if >= 0.5
        e = gg.groupby('pdb').agg(held=('held', 'mean'), sc_inv=('sc_inv', 'mean'), sc_var=('sc_var', 'mean'),
                                  mc=('mc', 'mean'), anyc=('n_any_4A', lambda x: (x > 0).mean()),
                                  rel=('rel', 'median'), sasa=('sasa_cplx', 'median'), nbur=('n_prot_4p5', 'median'),
                                  cl50=('cl50', 'first'))
        for c in ['held', 'sc_inv', 'sc_var', 'mc', 'anyc']:
            e[c] = e[c] >= 0.5
        e['mc_only'] = e.mc & ~e.sc_inv & ~e.sc_var
        e['cl50'] = [c if c else f'{f}_entry_{p}' for p, c in zip(e.index, e.cl50)]
        c = e.groupby('cl50').agg(held=('held', 'mean'), sc_inv=('sc_inv', 'mean'), sc_var=('sc_var', 'mean'),
                                  mc_only=('mc_only', 'mean'), anyc=('anyc', 'mean'), rel=('rel', 'median'))
        for k in ['held', 'sc_inv', 'sc_var', 'mc_only', 'anyc']:
            c[k] = c[k] >= 0.5
        nE, nC = len(e), len(c)
        # positions providing a polar partner: share of clusters in which they hold this atom in >= half the
        # cluster's entries
        pp = Psub[(Psub.family == f) & (Psub.subsite == s) & (Psub.label == lab) & Psub.polar_hb]
        posstat = []
        gkey = [q if q else (('unmapped:' if fc else 'other-chain:') + rn) for q, fc, rn in
                zip(pp.position, pp.frame_chain, pp.p_resname)]
        for q, gq in pp.groupby(pd.Series(gkey, index=pp.index)):
            ent = set(gq.pdb)
            ecl = e.assign(h=[p in ent for p in e.index]).groupby('cl50').h.mean() >= 0.5
            parts = ''.join(sorted(set(gq.part)))
            info = identity(f, q) if ':' not in q else dict(m50=np.nan, modal='', n_cl50=0)
            posstat.append((ecl.sum() / nC, q, parts, info['modal'], info['m50'], info['n_cl50'],
                            q in CAT[f], len(ent) / nE))
        posstat.sort(key=lambda x: -x[0])
        ptxt = '; '.join(f'{q}[{pt}]{"*" if cat else ""} cl {fc:.2f} ent {fe:.2f} m50 {m50:.2f}/{ncl:g}'
                         if ':' not in q else f'{q}[{pt}] cl {fc:.2f} ent {fe:.2f}'
                         for fc, q, pt, mod, m50, ncl, cat, fe in posstat)
        relm = float(np.nanmedian(e.rel)) if e.rel.notna().any() else np.nan
        held_c, inv_c, var_c, mco_c = c.held.mean(), c.sc_inv.mean(), c.sc_var.mean(), c.mc_only.mean()
        exp_c = (c.rel >= REL_EXPOSED).mean() if c.rel.notna().any() else np.nan
        if held_c >= 0.5:
            if inv_c >= 0.5:
                state = 'INV'
            elif var_c >= mco_c:
                state = 'VAR'
            else:
                state = 'MC'
        else:
            state = 'EXP' if (relm >= REL_EXPOSED and float(e.sasa.median()) >= SASA_EXPOSED) else 'BUR'
        rows.append(dict(set=tag, family=f, subsite=s, label=lab, n_units=len(gg), n_entries=nE, n_cl50=nC,
                         frac_entries_held=round(e.held.mean(), 3), frac_cl_held=round(held_c, 3),
                         frac_cl_inv_sc=round(inv_c, 3), frac_cl_var_sc=round(var_c, 3),
                         frac_cl_mc_only=round(mco_c, 3), frac_cl_any4A=round(c.anyc.mean(), 3),
                         median_sasa_cplx=round(float(e.sasa.median()), 1), median_rel=round(relm, 3),
                         frac_cl_exposed=round(float(exp_c), 3) if exp_c == exp_c else np.nan,
                         median_nprot4p5=float(e.nbur.median()), state=state, positions=ptxt))
    return pd.DataFrame(rows)


SETS = {
    'all_glucoselike': A,
    'glc_units': A[A.kind == 'Glc'],
    'glc_only_entries': A[A.glc_only_entry & (A.kind == 'Glc')],
}
summ = []
for tag, Asub in SETS.items():
    keys = set(zip(Asub.pdb, Asub.instance_id, Asub.ring))
    Psub = PT[[k in keys for k in zip(PT.pdb, PT.instance_id, PT.ring)]]
    summ.append(agg(Asub, Psub, tag))
SUM = pd.concat(summ, ignore_index=True)
order_lab = {k: i for i, k in enumerate(OXY + CARB)}
SUM['_o'] = SUM.label.map(order_lab)
SUM['_s'] = SUM.subsite.map({s: i for i, s in enumerate(SUBS)})
SUM = SUM.sort_values(['set', 'family', '_s', '_o']).drop(columns=['_o', '_s'])
SUM.to_csv(OUT / 'cell_summary.tsv', sep='\t', index=False)

# state tables (oxygens only)
with open(OUT / 'state_tables.md', 'w') as fh:
    for tag in SETS:
        for f in FAMS:
            sub = SUM[(SUM.set == tag) & (SUM.family == f) & SUM.label.isin(OXY)]
            if sub.empty:
                continue
            ncl = sub.groupby('subsite').n_cl50.max()
            nen = sub.groupby('subsite').n_entries.max()
            fh.write(f'\n### {f} - {tag}\n\n| subsite | n entries / 50% clusters | ' + ' | '.join(OXY) + ' |\n')
            fh.write('|---|---|' + '---|' * len(OXY) + '\n')
            for s in SUBS:
                if s not in ncl.index:
                    continue
                cells = []
                for o in OXY:
                    x = sub[(sub.subsite == s) & (sub.label == o)]
                    if x.empty:
                        cells.append('—')
                        continue
                    x = x.iloc[0]
                    cells.append(f'**{x.state}** {x.frac_cl_held:.2f}/{x.frac_cl_inv_sc:.2f}; {x.median_sasa_cplx:.0f} Å², '
                                 f'{x.median_rel:.2f}')
                fh.write(f'| {s} | {nen[s]} / {ncl[s]} | ' + ' | '.join(cells) + ' |\n')

# apolar face per subsite
face = []
for tag, Asub in SETS.items():
    keys = set(zip(Asub.pdb, Asub.instance_id, Asub.ring))
    Psub = PT[[k in keys for k in zip(PT.pdb, PT.instance_id, PT.ring)] & PT.label.isin(CARB)]
    for (f, s), gg in Asub[Asub.label.isin(CARB)].groupby(['family', 'subsite']):
        unitk = gg.groupby(['pdb', 'instance_id', 'ring']).cl50.first()
        stack, anyc = {}, {}
        for k in unitk.index:
            pk = Psub[(Psub.pdb == k[0]) & (Psub.instance_id == k[1]) & (Psub.ring == k[2])]
            anyc[k] = len(pk) > 0
            ar = pk[pk.aromatic_atom]
            best = ar.groupby(['p_chain', 'p_seq']).label.nunique().max() if len(ar) else 0
            stack[k] = best >= 3
        e = pd.DataFrame(dict(pdb=[k[0] for k in unitk.index], cl50=[c or 'e' + k[0] for k, c in unitk.items()],
                              anyc=[anyc[k] for k in unitk.index], stack=[stack[k] for k in unitk.index]))
        ee = e.groupby('pdb').agg(anyc=('anyc', 'mean'), stack=('stack', 'mean'), cl50=('cl50', 'first'))
        ee = ee.assign(anyc=ee['anyc'] >= 0.5, stack=ee['stack'] >= 0.5)
        cc = ee.groupby('cl50').agg(anyc=('anyc', 'mean'), stack=('stack', 'mean'))
        ar = Psub[(Psub.family == f) & (Psub.subsite == s) & Psub.aromatic_atom]
        apos = ar.groupby('position').pdb.nunique().sort_values(ascending=False)
        face.append(dict(set=tag, family=f, subsite=s, n_entries=len(ee), n_cl50=len(cc),
                         frac_cl_ring_carbon_contact=round((cc.anyc >= 0.5).mean(), 3),
                         frac_cl_aromatic_face=round((cc['stack'] >= 0.5).mean(), 3),
                         aromatic_positions_entries=' '.join(f'{q or "unmapped"}:{n}' for q, n in apos.items())))
pd.DataFrame(face).to_csv(OUT / 'apolar_face.tsv', sep='\t', index=False)

# kinds per component
kc = U.groupby(['family', 'comp', 'kind', 'kind_note']).size().rename('n_units').reset_index()
kc.to_csv(OUT / 'unit_kinds.tsv', sep='\t', index=False)
print(kc.to_string())
print('analysed units', len(U), 'entries', U.pdb.nunique(), '| excluded', len(skip_rows))
print('done')

# ------------------------------------------------------------------ substituted glucose units present in the census
# Components with an O-substituent other than a glycosidic bond to another sugar, found in GH13/GH14/GH15/GH31
# entries (local census scan of _chem_comp names + RCSB search 2026-11-24, see REPORT.md). None of them received a
# subsite in 33 (non-standard atom names or reference atom > 4.5 A from the catalytic midpoint). Here each pyranose
# ring is placed in the family frame and assigned to the nearest subsite (mean frame position of the C1 atoms of
# the density-supported units of that subsite), and its substituent atoms' protein contacts are listed.
SUBST = [('GH13', '1U33', 'A', 'LM2'), ('GH13', '5E0F', 'A', '5J7'), ('GH13', '3IJ7', 'A', 'B8D'),
         ('GH13', '1JFH', 'A', 'MA3'), ('GH13', '1JFH', 'A', 'MA2'), ('GH13', '1P6W', 'A', 'GTM'),
         ('GH14', '1J11', 'A', 'EPG'), ('GH14', '1J12', 'A', 'EBG'), ('GH14', '2XGB', 'A', 'EPG')]
cen = units.groupby(['family', 'subsite'])[['x_frame', 'y_frame', 'z_frame']].agg(
    lambda s: pd.to_numeric(s, errors='coerce').mean())
srows = []
for f, pdb, pch, comp in SUBST:
    st = read(pdb)
    m = st[0]
    chs = sorted(c for (ff, pp, c) in mats if ff == f and pp == pdb)
    prot = [(a, r, c) for c in m for r in c if is_aa(r.name) for a in r]
    PX = np.array([a.pos.tolist() for a, _, _ in prot])
    ptree = cKDTree(PX)
    cc = cen.loc[f]
    for c in m:
        for res in c:
            if res.name != comp:
                continue
            G, at = res_graph(res)
            rings6 = [set(x) for x in nx.minimum_cycle_basis(G) if len(x) == 6]
            sugar = []
            for rg in rings6:
                ro = [n for n in rg if G.nodes[n]['el'] == 'O']
                if len(ro) != 1:
                    continue
                # C1 = ring-O neighbour without an exocyclic carbon (C5 carries C6)
                c1 = [n for n in G[ro[0]] if n in rg and not any(G.nodes[x]['el'] == 'C' for x in G[n] if x not in rg)]
                if len(c1) == 1:
                    sugar.append((rg, c1[0], None))
            sugar_atoms = set().union(*[x[0] for x in sugar]) if sugar else set()
            for k, (rg, c1, c5) in enumerate(sugar):
                # place: chain matrix giving the nearest subsite
                best = None
                for pc in chs:
                    R, t = mats[(f, pdb, pc)]
                    x = R @ np.array(at[c1].pos.tolist()) + t
                    d = np.linalg.norm(cc.values.astype(float) - x, axis=1)
                    jj = int(np.nanargmin(d))
                    if best is None or d[jj] < best[0]:
                        best = (d[jj], cc.index[jj], pc)
                # substituents: heteroatoms on ring carbons or on their exocyclic carbons, with atoms beyond them
                anchors = []
                for n in rg:
                    for o in G[n]:
                        if o in rg:
                            continue
                        if G.nodes[o]['el'] in HET:
                            anchors.append((n, o))
                        elif G.nodes[o]['el'] == 'C':
                            anchors += [(o, h) for h in G[o] if h != n and G.nodes[h]['el'] in HET]
                subst = []
                for n, o in anchors:
                    beyond = [y for y in G[o] if y != n]
                    if not beyond:
                        continue
                    H = G.subgraph(set(G) - rg - {o, n})
                    comp_atoms = nx.node_connected_component(H, beyond[0]) if beyond[0] in H else set()
                    if comp_atoms & sugar_atoms:
                        kind = 'glycosidic (to another sugar ring)'
                    else:
                        kind = 'substituent'
                    near = set()
                    if kind == 'substituent':
                        for an in comp_atoms:
                            for q_ in ptree.query_ball_point(at[an].pos.tolist(), 4.0):
                                a_, r_, c_ = prot[q_]
                                near.add(f'{r_.name}{r_.seqid.num}')
                    subst.append((n, o, kind, sorted(comp_atoms), near))
                ss = [x for x in subst if x[2] == 'substituent'] or [('', '', 'none', [], set())]
                for n, o, kind, sa, near in ss:
                    srows.append(dict(family=f, pdb=pdb, comp=comp, lig_chain=c.name, lig_seq=res.seqid.num,
                                      ring_c1=c1, frame_chain_used=best[2], nearest_subsite=best[1],
                                      d_to_subsite_C1_mean_A=round(float(best[0]), 2), carbon=n, substituted_atom=o,
                                      n_substituent_atoms=len(sa), substituent_atoms=' '.join(sa[:12]),
                                      substituent_protein_residues_4A=' '.join(sorted(near))))
pd.DataFrame(srows).to_csv(OUT / 'substituted_glucose_in_census.tsv', sep='\t', index=False)
print(pd.DataFrame(srows).to_string())
