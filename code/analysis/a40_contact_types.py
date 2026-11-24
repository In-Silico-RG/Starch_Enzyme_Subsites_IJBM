"""Contact type (side-chain polar / backbone polar / apolar only) at conserved subsite positions.

For each density-supported sugar unit at a subsite (the units of a33_subsites.py), and each
reference position of interest, the protein residue mapped to that position in the same entry is
read from the deposited mmCIF (model 1, first altloc), together with the ligand atoms assigned to
that unit. A contact is 'side-chain polar' if a side-chain N/O lies within 3.5 A of a ligand N/O,
'backbone polar' if only a main-chain N/O does, and 'apolar' otherwise.
Usage: python a40_contact_types.py            GH13/GH31 positions of the short draft
       python a40_contact_types.py GH14 GH15  every position contacted in >= half the entries
                                              at subsites -2, -1, +1 of those families
Output in data/analysis/40_contact_types/ (prefix 'all_' for the second form).
"""
import sqlite3
import sys
import csv
from collections import defaultdict, Counter
from pathlib import Path
import gemmi
import numpy as np

P = Path(__file__).resolve().parents[2]
S = P / 'data/analysis/33_subsites'
CIF = P / 'data/census/cif'
OUT = P / 'data/analysis/40_contact_types'
OUT.mkdir(parents=True, exist_ok=True)
HB = 3.5
BACKBONE = {'N', 'O', 'OXT'}

TARGETS = {
    'GH13': {'-1': ['D289', 'H92', 'D179', 'R177', 'Y51', 'H288'],
             '+1': ['F180', 'G183', 'E204'], '+2': ['W206']},
    'GH31': {'-1': ['D327', 'H600', 'D443', 'D542', 'W539', 'W441', 'I328', 'F575', 'I364']},
}


PREFIX = ''
if len(sys.argv) > 1:
    PREFIX = 'all_'
    TARGETS = defaultdict(lambda: defaultdict(list))
    for r in csv.DictReader(open(S / 'subsite_residues_reference_numbering.tsv'), delimiter='\t'):
        if r['family'] in sys.argv[1:] and r['subsite'] in ('-2', '-1', '+1') and float(r['fraction']) >= 0.5:
            TARGETS[r['family']][r['subsite']].append(r['ref_position'])
    TARGETS = {f: dict(v) for f, v in TARGETS.items()}

con = sqlite3.connect(f"file:{P / 'data/db/enzimas_experimental.sqlite'}?mode=ro", uri=True)
UNIPROT = {pdb: (acc or pdb) for pdb, acc in con.execute('select pdb, uniprot from entry')}


def rscc_ok(v):
    try:
        return float(v) >= 0.80
    except ValueError:
        return False


units = [u for u in csv.DictReader(open(S / 'units.tsv'), delimiter='\t')
         if u['family'] in TARGETS and u['in_tables'] == 'True' and rscc_ok(u['rscc'])
         and u['subsite'] in TARGETS[u['family']]]

# ligand atoms per unit and the protein residue per (entry, position), from the asymmetric-unit contacts
lig_atoms = defaultdict(set)
res_at = defaultdict(set)
for r in csv.DictReader(open(S / 'contacts_mapped.tsv'), delimiter='\t'):
    if r['family'] not in TARGETS or r['asu'] != 'True':
        continue
    for ring in str(r['rings']).split(','):
        key = (r['pdb'], r['instance_id'], ring.strip())
        lig_atoms[key].add(r['ligand_atom'])
        if r['ref_position']:
            res_at[key + (r['ref_position'],)].add((r['chain'], str(r['seq']), r['ins_code'], r['resname']))

cache = {}


def model(pdb):
    if pdb not in cache:
        st = gemmi.read_structure(str(CIF / f'{pdb}.cif'))
        st.setup_entities()
        cache[pdb] = st[0]
    return cache[pdb]


def residue(m, chain, seq, icode, name=None):
    for ch in m:
        if ch.name != chain:
            continue
        for res in ch:
            if str(res.seqid.num) == str(seq) and (res.seqid.icode.strip() == (icode or '').strip()) \
                    and (name is None or res.name == name):
                return res
    return None


def first_alt(res):
    out, seen = [], set()
    for a in res:
        if a.name in seen:
            continue
        seen.add(a.name)
        out.append(a)
    return out


rows = []
for u in units:
    fam, sub = u['family'], u['subsite']
    key = (u['pdb'], u['instance_id'], u['ring'])
    m = model(u['pdb'])
    lres = residue(m, u['chain'], u['seq'], '', u['comp'])
    if lres is None:
        continue
    names = lig_atoms.get(key, set())
    latoms = [a for a in first_alt(lres) if a.name in names]
    lpol = np.array([a.pos.tolist() for a in latoms if a.element.name in ('O', 'N')]) if latoms else np.empty((0, 3))
    for pos in TARGETS[fam][sub]:
        for (ch, seq, ic, rn) in res_at.get(key + (pos,), ()):
            pres = residue(m, ch, seq, ic, rn)
            if pres is None:
                continue
            pa = [a for a in first_alt(pres) if a.element.name in ('O', 'N')]
            sc = [a for a in pa if a.name not in BACKBONE]
            bb = [a for a in pa if a.name in BACKBONE]

            def close(atoms):
                if not atoms or len(lpol) == 0:
                    return False
                X = np.array([a.pos.tolist() for a in atoms])
                return bool((np.linalg.norm(X[:, None, :] - lpol[None, :, :], axis=2) <= HB).any())
            kind = 'sidechain_polar' if close(sc) else ('backbone_polar' if close(bb) else 'apolar')
            rows.append(dict(family=fam, subsite=sub, position=pos, pdb=u['pdb'],
                             protein=UNIPROT.get(u['pdb'], u['pdb']), comp=u['comp'],
                             chain=ch, seq=seq, resname=rn, kind=kind))

with open(OUT / f'{PREFIX}contact_types.tsv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
    w.writeheader()
    w.writerows(rows)

# one vote per entry and position: the strongest kind seen in that entry
rank = {'sidechain_polar': 2, 'backbone_polar': 1, 'apolar': 0}
best = {}
types = defaultdict(set)
prots = defaultdict(lambda: defaultdict(set))
for r in rows:
    prots[(r['family'], r['subsite'], r['position'])][r['resname']].add(r['protein'])
    k = (r['family'], r['subsite'], r['position'], r['pdb'])
    if k not in best or rank[r['kind']] > rank[best[k]]:
        best[k] = r['kind']
    types[k].add(r['resname'])
summ = defaultdict(Counter)
rtypes = defaultdict(Counter)
for k, v in best.items():
    summ[k[:3]][v] += 1
    for t in types[k]:
        rtypes[k[:3]][t] += 1
with open(OUT / f'{PREFIX}summary.tsv', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['family', 'subsite', 'position', 'entries', 'sidechain_polar', 'backbone_polar', 'apolar', 'residue_types', 'proteins_per_type'])
    for k in sorted(summ):
        c = summ[k]
        n = sum(c.values())
        w.writerow([*k, n, c['sidechain_polar'], c['backbone_polar'], c['apolar'],
                    ', '.join(f'{t} {x}' for t, x in rtypes[k].most_common()),
                    ', '.join(f'{t} {len(v)}' for t, v in sorted(prots[k].items(), key=lambda z: -len(z[1])))])
print(open(OUT / f'{PREFIX}summary.tsv').read())
