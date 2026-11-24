"""Rule test, step 1: build the tables that the tests of a44b_tests.py read.

Read-only on every input. Inputs: data/analysis/33_subsites/{units,contacts_mapped}.tsv,
data/calc/03_frames/superposition.tsv (accepted chain -> reference-frame matrices, GH13 ref 1AMY:A,
GH31 ref 2QMJ:A), data/census/cif/*.cif, data/analysis/grip_independence_check/{entries_minus1.tsv,
ident.json}, data/db/enzimas_experimental.sqlite (opened read-only).

Outputs (data/analysis/44_rule_test/):
  units_used.tsv        density-supported units (RSCC >= 0.80, in_tables) at subsites -3..+3 of GH13/GH31
  posmap.tsv            residue at every reference position for every accepted GH13/GH31 chain
                        (residue Calpha within 2.0 A of the reference Calpha after the frame superposition;
                        nearest residue per reference position)
  contacts_geom.tsv     per unit x contacted reference position: minimum side-chain and main-chain N/O to
                        ligand N/O distance, with and without a heavy-atom angle criterion (both angles
                        neighbour-X...Y and neighbour-Y...X >= 90 deg; neighbours = heavy atoms within 1.75 A
                        in the same residue)
  mutations.tsv         engineered mutations (_struct_ref_seq_dif.details) of every entry used
  clusters_bound.tsv    ligand-bound entries -> protein key, 50 % cluster (I_al, single linkage), 30 % cluster
                        (I_sh); protein-level identities from ident.json, extra entries aligned here
  clusters_heldout.tsv  superposed entries whose chain is < 50 % (I_al) identical to every ligand-bound protein,
                        clustered among themselves (50 % I_al, 30 % I_sh)
Run: ~/miniconda3/bin/python code/analysis/a44a_build.py
"""
import csv, json, itertools, sqlite3
from collections import defaultdict, Counter
from multiprocessing import Pool
from pathlib import Path
import numpy as np
import pandas as pd
import gemmi
from Bio import Align
from Bio.Align import substitution_matrices

P = Path(__file__).resolve().parents[2]
S = P / 'data/analysis/33_subsites'
FR = P / 'data/calc/03_frames'
CIF = P / 'data/census/cif'
GI = P / 'data/analysis/grip_independence_check'
OUT = P / 'data/analysis/44_rule_test'
OUT.mkdir(parents=True, exist_ok=True)
FAMS = ['GH13', 'GH31']
SUBS = ['-3', '-2', '-1', '+1', '+2', '+3']
REF = {'GH13': ('1AMY', 'A'), 'GH31': ('2QMJ', 'A')}
MAP_MAX = 2.0
BACKBONE = {'N', 'O', 'OXT'}


def read(pdb):
    st = gemmi.read_structure(str(CIF / f'{pdb}.cif'))
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    st.setup_entities()
    return st


def one(name):
    t = gemmi.find_tabulated_residue(name)
    return t.one_letter_code.upper() if t else 'X'


def is_aa(name):
    t = gemmi.find_tabulated_residue(name)
    return bool(t and t.is_amino_acid())


# ---------------------------------------------------------------- units
u = pd.read_csv(S / 'units.tsv', sep='\t', dtype=str, keep_default_na=False)
u['rscc_f'] = pd.to_numeric(u.rscc, errors='coerce')
units = u[u.family.isin(FAMS) & (u.in_tables == 'True') & (u.rscc_f >= 0.80) & u.subsite.isin(SUBS)].copy()
units.to_csv(OUT / 'units_used.tsv', sep='\t', index=False)
print('units', len(units), 'entries', units.pdb.nunique())

# ---------------------------------------------------------------- superposition matrices
sup = pd.read_csv(FR / 'superposition.tsv', sep='\t')
sup = sup[sup.status.str.startswith('accepted') & sup.family.isin(FAMS)]
mats = {}
for r in sup.itertuples():
    v = [float(x) for x in str(r[-1]).split()]
    mats[(r.family, r.pdb, r.chain)] = (np.array(v[:9]).reshape(3, 3), np.array(v[9:12]))

refs = {}
for f, (rp, rc) in REF.items():
    st = read(rp)
    lab, xyz = [], []
    for res in st[0][rc]:
        ca = res.find_atom('CA', '*')
        if ca is not None and is_aa(res.name):
            lab.append(f'{one(res.name)}{res.seqid.num}{res.seqid.icode.strip()}')
            xyz.append(ca.pos.tolist())
    refs[f] = (lab, np.array(xyz))

# ---------------------------------------------------------------- position map, sequences, mutations
chains_by_pdb = defaultdict(list)
for (f, p, c) in mats:
    chains_by_pdb[(f, p)].append(c)

pm_rows, seq_rows, mut_rows = [], [], []
for (f, pdb), chs in sorted(chains_by_pdb.items()):
    st = read(pdb)
    m = st[0]
    lab, RX = refs[f]
    for ch in chs:
        if m.find_chain(ch) is None:
            continue
        R, t = mats[(f, pdb, ch)]
        best = {}
        for res in m[ch]:
            ca = res.find_atom('CA', '*')
            if ca is None or not is_aa(res.name):
                continue
            x = R @ np.array(ca.pos.tolist()) + t
            d = np.linalg.norm(RX - x, axis=1)
            k = int(d.argmin())
            if d[k] <= MAP_MAX and (lab[k] not in best or d[k] < best[lab[k]][0]):
                best[lab[k]] = (float(d[k]), res.name, res.seqid.num, res.seqid.icode.strip())
        for q, (d, rn, sq, ic) in best.items():
            pm_rows.append(dict(family=f, pdb=pdb, chain=ch, position=q, resname=rn, seq=sq, icode=ic, d_ca=round(d, 2)))
        ent = st.get_entity_of(m[ch].get_polymer()) if hasattr(st, 'get_entity_of') else None
        seq = ''
        if ent is not None and ent.full_sequence:
            seq = ''.join(one(gemmi.Entity.first_mon(x)) for x in ent.full_sequence)
        if not seq:
            seq = gemmi.one_letter_code([r.name for r in m[ch].get_polymer()]).upper()
        seq_rows.append(dict(family=f, pdb=pdb, chain=ch, seq=seq))
    doc = gemmi.cif.read(str(CIF / f'{pdb}.cif'))
    tab = doc.sole_block().find('_struct_ref_seq_dif.', ['pdbx_pdb_strand_id', 'pdbx_auth_seq_num', 'mon_id',
                                                         'db_mon_id', 'details'])
    for row in tab:
        det = gemmi.cif.as_string(row[4])
        if 'mutation' in det.lower() and 'engineered' in det.lower():
            mut_rows.append(dict(family=f, pdb=pdb, chain=gemmi.cif.as_string(row[0]), seq=gemmi.cif.as_string(row[1]),
                                 mutant=gemmi.cif.as_string(row[2]), wild_type=gemmi.cif.as_string(row[3]), details=det))
pm = pd.DataFrame(pm_rows)
pm.to_csv(OUT / 'posmap.tsv', sep='\t', index=False)
seqs = pd.DataFrame(seq_rows)
seqs.to_csv(OUT / 'chain_sequences.tsv', sep='\t', index=False)
pd.DataFrame(mut_rows).to_csv(OUT / 'mutations.tsv', sep='\t', index=False)
print('posmap rows', len(pm), 'chains', len(seqs), 'mutation rows', len(mut_rows))

# ---------------------------------------------------------------- contact geometry per unit and position
lig_atoms, res_at = defaultdict(set), defaultdict(set)
for r in csv.DictReader(open(S / 'contacts_mapped.tsv'), delimiter='\t'):
    if r['family'] not in FAMS or r['asu'] != 'True' or r['same_chain_as_frame'] != 'True':
        continue
    for ring in str(r['rings']).split(','):
        key = (r['pdb'], r['instance_id'], ring.strip())
        lig_atoms[key].add(r['ligand_atom'])
        if r['ref_position']:
            res_at[key].add((r['ref_position'], r['chain'], str(r['seq']), r['ins_code'], r['resname']))


def find_res(m, chain, seq, icode, name):
    ch = m.find_chain(chain)
    if ch is None:
        return None
    for res in ch:
        if str(res.seqid.num) == str(seq) and res.seqid.icode.strip() == (icode or '').strip() and res.name == name:
            return res
    return None


def neighbours(res, atom):
    return [b for b in res if b.name != atom.name and b.element.name != 'H' and b.pos.dist(atom.pos) <= 1.75]


def ang_ok(pres, X, lres, Y):
    for nb in neighbours(pres, X):
        if np.degrees(gemmi.calculate_angle(nb.pos, X.pos, Y.pos)) < 90:
            return False
    for nb in neighbours(lres, Y):
        if np.degrees(gemmi.calculate_angle(nb.pos, Y.pos, X.pos)) < 90:
            return False
    return True


def geometry(units_df, outname):
    geo = []
    cache = {}
    for uu in units_df.itertuples():
        if uu.pdb not in cache:
            cache.clear()
            cache[uu.pdb] = read(uu.pdb)[0]
        m = cache[uu.pdb]
        key = (uu.pdb, uu.instance_id, uu.ring)
        lres = find_res(m, uu.chain, uu.seq, '', uu.comp)
        if lres is None:
            continue
        names = lig_atoms.get(key, set())
        lpol = [a for a in lres if a.name in names and a.element.name in ('O', 'N')]
        for (q, ch, sq, ic, rn) in sorted(res_at.get(key, ())):
            pres = find_res(m, ch, sq, ic, rn)
            if pres is None:
                continue
            rec = dict(family=uu.family, subsite=uu.subsite, position=q, pdb=uu.pdb, instance_id=uu.instance_id,
                       ring=uu.ring, comp=uu.comp, chain=ch, seq=sq, icode=ic, resname=rn)
            for part, sel in (('sc', lambda a: a.name not in BACKBONE), ('bb', lambda a: a.name in BACKBONE)):
                d_all, d_ang = np.inf, np.inf
                for X in pres:
                    if X.element.name not in ('O', 'N') or not sel(X):
                        continue
                    for Y in lpol:
                        d = X.pos.dist(Y.pos)
                        if d < d_all:
                            d_all = d
                        if d < d_ang and d <= 4.0 and ang_ok(pres, X, lres, Y):
                            d_ang = d
                rec[f'd_{part}'] = round(d_all, 3) if np.isfinite(d_all) else ''
                rec[f'd_{part}_ang'] = round(d_ang, 3) if np.isfinite(d_ang) else ''
            geo.append(rec)
    pd.DataFrame(geo).to_csv(OUT / outname, sep='\t', index=False)
    print(outname, 'rows', len(geo))


geometry(units, 'contacts_geom.tsv')
# reference-protein entries for the class read from the reference alone (GH13: 1BG9, barley AMY2 = protein of
# the reference 1AMY, which holds no sugar; GH31: 2QMJ itself, whose units have no RSCC and are kept here only)
ref_units = u[u.pdb.isin(['1BG9', '2QMJ']) & (u.in_tables == 'True') & u.subsite.isin(SUBS)]
geometry(ref_units, 'contacts_geom_ref.tsv')

# ---------------------------------------------------------------- clusters
al = Align.PairwiseAligner()
al.mode = 'global'
al.substitution_matrix = substitution_matrices.load('BLOSUM62')
al.open_gap_score, al.extend_gap_score, al.end_gap_score = -10, -0.5, 0


def ident(pair):
    a, b = pair
    c = al.align(a, b)[0].counts()
    n = c.identities + c.mismatches
    return (c.identities / n if n else 0.0, c.identities / min(len(a), len(b)))


def single_link(names, sim, thr):
    par = {n: n for n in names}

    def f(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    for (i, j), s in sim.items():
        if s >= thr:
            par[f(i)] = f(j)
    roots = {}
    return {n: roots.setdefault(f(n), len(roots) + 1) for n in names}


J = json.load(open(GI / 'ident.json'))
e1 = pd.read_csv(GI / 'entries_minus1.tsv', sep='\t', keep_default_na=False)
frame_chain = units.groupby('pdb').frame_chain.agg(lambda s: Counter(s).most_common(1)[0][0]).to_dict()
seqd = {(r.pdb, r.chain): r.seq for r in seqs.itertuples()}

bound_rows = []
heldout_rows = []
with Pool(8) as pool:
    for fam in FAMS:
        o = J['out'][fam]
        names = list(o['names'])
        M, Ms = np.array(o['M']), np.array(o['Ms'])
        sim50, sim30 = {}, {}
        for i, j in itertools.combinations(range(len(names)), 2):
            sim50[(names[i], names[j])] = M[i, j]
            sim30[(names[i], names[j])] = Ms[i, j]
        ee = e1[e1.family == fam]
        key_of = dict(zip(ee.pdb, ee.protein_key))
        rep_seq = {}
        for r in ee.itertuples():
            rep_seq.setdefault(r.protein_key, []).append(r.chain_seq)
        # extra ligand-bound entries (at -3..+3 but not at -1)
        extra = sorted(set(units[units.family == fam].pdb) - set(key_of))
        for p in extra:
            s = seqd.get((p, frame_chain[p]), '')
            res = pool.map(ident, [(s, rep_seq[k][0]) for k in names])
            best = max(range(len(names)), key=lambda i: res[i][0])
            if res[best][0] >= 0.995:
                key_of[p] = names[best]
                print(f'{fam} extra {p}: same protein as {names[best]} (I_al {res[best][0]:.3f})')
                continue
            k = 'PDB:' + p
            for i, n in enumerate(names):
                sim50[(k, n)] = res[i][0]
                sim30[(k, n)] = res[i][1]
            names.append(k)
            rep_seq[k] = [s]
            key_of[p] = k
            print(f'{fam} extra {p}: new protein key, max I_al {res[best][0]:.3f} to {names[best]}')
        c50 = single_link(names, sim50, 0.50)
        c30 = single_link(names, sim30, 0.30)
        for p, k in sorted(key_of.items()):
            if p in set(units.pdb):
                bound_rows.append(dict(family=fam, pdb=p, frame_chain=frame_chain.get(p, ''), protein_key=k,
                                       cl50=f'{fam}_b50_{c50[k]}', cl30=f'{fam}_b30_{c30[k]}'))
        print(fam, 'bound protein keys', len(set(key_of.values())), '50% clusters', len(set(c50[key_of[p]] for p in key_of)),
              '30% clusters', len(set(c30[key_of[p]] for p in key_of)))
        # held-out: superposed entries not ligand-bound at -3..+3, one chain per entry (first accepted, sorted)
        bound_seqs = [rep_seq[k][0] for k in names]
        cand = sorted({p for (f, p, c) in mats if f == fam} - set(key_of))
        cseq = {}
        for p in cand:
            ch = sorted(c for (f, pp, c) in mats if f == fam and pp == p)[0]
            cseq[p] = (ch, seqd.get((p, ch), ''))
        useq = sorted({v[1] for v in cseq.values() if v[1]})
        res = pool.map(ident, [(s, b) for s in useq for b in bound_seqs], chunksize=50)
        nb = len(bound_seqs)
        maxb = {s: max(res[i * nb + j][0] for j in range(nb)) for i, s in enumerate(useq)}
        held = [s for s in useq if maxb[s] < 0.50]
        pairs = list(itertools.combinations(held, 2))
        res2 = pool.map(ident, pairs, chunksize=50)
        s50 = {pr: r[0] for pr, r in zip(pairs, res2)}
        s30 = {pr: r[1] for pr, r in zip(pairs, res2)}
        h50, h30 = single_link(held, s50, 0.50), single_link(held, s30, 0.30)
        for p in cand:
            ch, s = cseq[p]
            heldout_rows.append(dict(family=fam, pdb=p, chain=ch, max_ident_to_bound=round(maxb.get(s, np.nan), 3),
                                     heldout=s in h50, cl50=f'{fam}_h50_{h50[s]}' if s in h50 else '',
                                     cl30=f'{fam}_h30_{h30[s]}' if s in h30 else ''))
        print(fam, 'non-bound superposed entries', len(cand), 'unique seqs', len(useq), 'held-out seqs', len(held),
              'held-out 50% clusters', len(set(h50.values())), '30%', len(set(h30.values())))

pd.DataFrame(bound_rows).to_csv(OUT / 'clusters_bound.tsv', sep='\t', index=False)
pd.DataFrame(heldout_rows).to_csv(OUT / 'clusters_heldout.tsv', sep='\t', index=False)
print('done')
