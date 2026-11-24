"""Control for the contact-type rule: sugars bound away from the catalytic site.

Density-supported sugar units with no subsite and more than 12 A from the midpoint of the catalytic
carboxylates (units.tsv of a33_subsites.py). For each reference position they contact (asymmetric
unit, contacts_mapped.tsv), the residue and the kind of contact are read from the mmCIF exactly as
in a40_contact_types.py. Positions contacted in >= MIN_ENTRIES entries from >= MIN_PROT proteins
are kept. For each, the fraction of entries holding its most common residue type is reported,
next to the same quantity for the catalytic-subsite positions of a40.
Usage: python a41_control_noncatalytic.py [family ...]   (default GH13)
Output: data/analysis/41_control_noncatalytic/
"""
import csv
import sqlite3
import sys
from collections import defaultdict, Counter
from pathlib import Path
import gemmi
import numpy as np

P = Path(__file__).resolve().parents[2]
S = P / 'data/analysis/33_subsites'
CIF = P / 'data/census/cif'
OUT = P / 'data/analysis/41_control_noncatalytic'
OUT.mkdir(parents=True, exist_ok=True)
FAMS = sys.argv[1:] or ['GH13']
import os
DMIN, HB = 12.0, 3.5
MIN_ENTRIES = int(os.environ.get("MIN_ENTRIES", 10))
MIN_PROT = int(os.environ.get("MIN_PROT", 3))
BACKBONE = {'N', 'O', 'OXT'}
AROM = {'PHE', 'TRP', 'TYR', 'HIS'}

con = sqlite3.connect(f"file:{P / 'data/db/enzimas_experimental.sqlite'}?mode=ro", uri=True)
UNIPROT = {pdb: ((acc or '').split(';')[0].strip() or pdb) for pdb, acc in con.execute('select pdb, uniprot from entry')}


def fnum(v):
    try:
        return float(v)
    except ValueError:
        return float('nan')


units = [u for u in csv.DictReader(open(S / 'units.tsv'), delimiter='\t')
         if u['family'] in FAMS and u['subsite'] == '' and u['ligand_class'] == 'sugar'
         and fnum(u['rscc']) >= 0.80 and fnum(u['d_ref_to_midpoint_A']) > DMIN]

lig_atoms, res_at = defaultdict(set), defaultdict(set)
for r in csv.DictReader(open(S / 'contacts_mapped.tsv'), delimiter='\t'):
    if r['family'] not in FAMS or r['asu'] != 'True':
        continue
    for ring in str(r['rings']).split(','):
        key = (r['pdb'], r['instance_id'], ring.strip())
        lig_atoms[key].add(r['ligand_atom'])
        if r['ref_position']:
            res_at[key].add((r['ref_position'], r['chain'], str(r['seq']), r['ins_code'], r['resname']))

# a unit is non-catalytic only if it touches no position that occurs anywhere in the subsite tables
cleft = defaultdict(set)
for r in csv.DictReader(open(S / 'subsite_residues_reference_numbering.tsv'), delimiter='\t'):
    cleft[r['family']].add(r['ref_position'])
n_dist = len(units)
units = [u for u in units
         if not any(p in cleft[u['family']] for (p, *_r) in res_at.get((u['pdb'], u['instance_id'], u['ring']), ()))]
print(f'units > {DMIN} A from the catalytic midpoint: {n_dist}; of these touching no subsite-table position: {len(units)}')

cache = {}


def model(pdb):
    if pdb not in cache:
        cache[pdb] = gemmi.read_structure(str(CIF / f'{pdb}.cif'))[0]
    return cache[pdb]


def residue(m, chain, seq, icode, name):
    for ch in m:
        if ch.name == chain:
            for res in ch:
                if str(res.seqid.num) == str(seq) and res.seqid.icode.strip() == (icode or '').strip() and res.name == name:
                    return res
    return None


def first_alt(res):
    seen, out = set(), []
    for a in res:
        if a.name not in seen:
            seen.add(a.name)
            out.append(a)
    return out


def close(atoms, lpol):
    if not atoms or len(lpol) == 0:
        return False
    X = np.array([a.pos.tolist() for a in atoms])
    return bool((np.linalg.norm(X[:, None, :] - lpol[None, :, :], axis=2) <= HB).any())


rank = {'sidechain_polar': 2, 'backbone_polar': 1, 'apolar': 0}
best, rtype = {}, {}
for u in units:
    key = (u['pdb'], u['instance_id'], u['ring'])
    m = model(u['pdb'])
    lres = residue(m, u['chain'], u['seq'], '', u['comp'])
    if lres is None:
        continue
    names = lig_atoms.get(key, set())
    lpol = np.array([a.pos.tolist() for a in first_alt(lres) if a.name in names and a.element.name in ('O', 'N')])
    for (pos, ch, seq, ic, rn) in res_at.get(key, ()):
        pres = residue(m, ch, seq, ic, rn)
        if pres is None:
            continue
        pa = [a for a in first_alt(pres) if a.element.name in ('O', 'N')]
        kind = ('sidechain_polar' if close([a for a in pa if a.name not in BACKBONE], lpol)
                else 'backbone_polar' if close([a for a in pa if a.name in BACKBONE], lpol) else 'apolar')
        k = (u['family'], pos, u['pdb'])
        if k not in best or rank[kind] > rank[best[k]]:
            best[k] = kind
        rtype[k] = rn

with open(OUT / 'control_entries.tsv', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['family', 'position', 'pdb', 'protein', 'resname', 'kind'])
    for (fam, pos, pdb), kind in sorted(best.items()):
        w.writerow([fam, pos, pdb, UNIPROT.get(pdb, pdb), rtype[(fam, pos, pdb)], kind])

agg = defaultdict(lambda: {'kinds': Counter(), 'types': Counter(), 'prots': set()})
for (fam, pos, pdb), kind in best.items():
    a = agg[(fam, pos)]
    a['kinds'][kind] += 1
    a['types'][rtype[(fam, pos, pdb)]] += 1
    a['prots'].add(UNIPROT.get(pdb, pdb))

rows = []
for (fam, pos), a in agg.items():
    n = sum(a['kinds'].values())
    if n < MIN_ENTRIES or len(a['prots']) < MIN_PROT:
        continue
    kind = a['kinds'].most_common(1)[0][0]
    top, ntop = a['types'].most_common(1)[0]
    cls = 'H' if kind == 'sidechain_polar' else ('Ar' if top in AROM else ('BB' if kind == 'backbone_polar' else 'Al'))
    rows.append(dict(family=fam, position=pos, entries=n, proteins=len(a['prots']), contact=kind, cls=cls,
                     dominant=top, dominant_fraction=round(ntop / n, 3),
                     types=', '.join(f'{t} {c}' for t, c in a['types'].most_common())))
rows.sort(key=lambda r: (r['family'], r['cls'], -r['entries']))
with open(OUT / 'control_positions.tsv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
    w.writeheader()
    w.writerows(rows)
print(f'units used: {len(units)}; positions kept: {len(rows)}')
for r in rows:
    print('\t'.join(str(r[k]) for k in rows[0]))
