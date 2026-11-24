"""Data for Figure 2 (v2) of the short draft: residue types per conserved position, one vote per protein.

Problem: Table 1 counts residue types once per PDB entry; proteins deposited many times weigh more. The figure
data panel counts one vote per protein (UniProt accession of the chain that carries the residue, from SIFTS).

Data: units at subsites (data/analysis/33_subsites/units.tsv, RSCC >= 0.80, in_tables), their contacts
(33_subsites/contacts_mapped.tsv), the deposited mmCIF files (data/census/cif), and the SIFTS residue map
(data/db/enzimas_experimental.sqlite, table sifts_residue; fallback entry.uniprot).

Method: identical to code/analysis/a40_contact_types.py (same units, same 3.5 A N/O-N/O criterion, model 1,
first altloc), re-run here for the position set of the figure:
  Table 1 positions minus Gly183 (+1; contacted in 34/44 = 0.773 of +1 entries, below the 0.80 threshold of
  33_subsites/subsite_residues_reference_numbering.tsv) plus Lys182 (+2; contacted in 37/41 = 0.902).
The rows of positions already in 40_contact_types/contact_types.tsv are checked against that file.
Engineered mutants are excluded from the votes: 1QI5, 5TD4, 9UNS, 9XG6 (D->N at Asp289), 7WJC, 7WJF (D542A).
Per protein and position, the residue type is the one seen in that protein's entries (if a protein shows two
types, each type gets a fractional vote 1/k); the contact kind of the protein is its majority kind over entries.

Writes data/figures/fig2_v2_work/{contact_rows.tsv, per_protein.tsv, per_protein_counts.tsv}.
Run: ~/miniconda3/bin/python code/figures/fig2_v2_data.py
"""
import csv
import sqlite3
from collections import defaultdict, Counter
from pathlib import Path

import gemmi
import numpy as np

P = Path(__file__).resolve().parents[2]
S = P / 'data/analysis/33_subsites'
CIF = P / 'data/census/cif'
CT = P / 'data/analysis/40_contact_types/contact_types.tsv'
DB = P / 'data/db/enzimas_experimental.sqlite'
WORK = P / 'data/figures/fig2_v2_work'
WORK.mkdir(parents=True, exist_ok=True)
HB = 3.5
BACKBONE = {'N', 'O', 'OXT'}
MUTANTS = {'1QI5', '5TD4', '9UNS', '9XG6', '7WJC', '7WJF'}

TARGETS = {
    'GH13': {'-1': ['D179', 'D289', 'H288', 'R177', 'H92', 'Y51'],
             '+1': ['E204', 'F180'], '+2': ['K182', 'W206']},
    'GH31': {'-1': ['D327', 'H600', 'D443', 'D542', 'W441', 'W539', 'F575', 'I328', 'I364']},
}
# contact class of each position, as in Table 1 (majority kind); K182 assigned from its own data below
CLASS = {'D179': 'H', 'D289': 'H', 'H288': 'H', 'R177': 'H', 'H92': 'H', 'E204': 'H', 'Y51': 'Ar',
         'F180': 'Al', 'W206': 'Al', 'D327': 'H', 'H600': 'H', 'D443': 'H', 'D542': 'H',
         'W441': 'Ar', 'W539': 'Ar', 'F575': 'Ar', 'I328': 'Al', 'I364': 'Al'}


def rscc_ok(v):
    try:
        return float(v) >= 0.80
    except ValueError:
        return False


def first_alt(res):
    out, seen = [], set()
    for a in res:
        if a.name in seen:
            continue
        seen.add(a.name)
        out.append(a)
    return out


def residue(m, chain, seq, icode, name=None):
    for ch in m:
        if ch.name != chain:
            continue
        for res in ch:
            if str(res.seqid.num) == str(seq) and res.seqid.icode.strip() == (icode or '').strip() \
                    and (name is None or res.name == name):
                return res
    return None


def compute_rows():
    units = [u for u in csv.DictReader(open(S / 'units.tsv'), delimiter='\t')
             if u['family'] in TARGETS and u['in_tables'] == 'True' and rscc_ok(u['rscc'])
             and u['subsite'] in TARGETS[u['family']]]
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
    rows = []
    for u in units:
        fam, sub = u['family'], u['subsite']
        key = (u['pdb'], u['instance_id'], u['ring'])
        if u['pdb'] not in cache:
            cache[u['pdb']] = gemmi.read_structure(str(CIF / f"{u['pdb']}.cif"))[0]
        m = cache[u['pdb']]
        lres = residue(m, u['chain'], u['seq'], '', u['comp'])
        if lres is None:
            continue
        names = lig_atoms.get(key, set())
        latoms = [a for a in first_alt(lres) if a.name in names]
        lpol = np.array([a.pos.tolist() for a in latoms if a.element.name in ('O', 'N')]) \
            if latoms else np.empty((0, 3))
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
                rows.append(dict(family=fam, subsite=sub, position=pos, pdb=u['pdb'], comp=u['comp'],
                                 chain=ch, seq=seq, resname=rn, kind=kind))
    return rows


def main():
    rows = compute_rows()
    # check against the a40 output for the positions it holds
    old = Counter((r['family'], r['subsite'], r['position'], r['pdb'], r['chain'], r['seq'], r['resname'], r['kind'])
                  for r in csv.DictReader(open(CT), delimiter='\t') if r['position'] != 'G183')
    new = Counter((r['family'], r['subsite'], r['position'], r['pdb'], r['chain'], r['seq'], r['resname'], r['kind'])
                  for r in rows if r['position'] != 'K182')
    print(f'check vs contact_types.tsv: {sum(old.values())} old rows, {sum(new.values())} new rows, '
          f'identical = {old == new}')
    if old != new:
        print('  only old:', list((old - new).items())[:5]); print('  only new:', list((new - old).items())[:5])

    con = sqlite3.connect(f'file:{DB}?mode=ro', uri=True)
    ent = {p: u for p, u in con.execute('select pdb, uniprot from entry')}
    sifts = {}
    for p, ch, s, u in con.execute('select pdb, chain, auth_seq, uniprot from sifts_residue'):
        sifts[(p, ch, str(s))] = u
    for r in rows:
        r['uniprot'] = sifts.get((r['pdb'], r['chain'], str(r['seq']))) or ent.get(r['pdb']) or r['pdb']
        r['mutant'] = r['pdb'] in MUTANTS
    with open(WORK / 'contact_rows.tsv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader(); w.writerows(rows)

    # one vote per protein
    rank = {'sidechain_polar': 2, 'backbone_polar': 1, 'apolar': 0}
    per = defaultdict(lambda: dict(types=Counter(), kinds=Counter(), entries=set()))
    entries_all = defaultdict(set)
    for r in rows:
        k = (r['family'], r['subsite'], r['position'])
        entries_all[k].add(r['pdb'])
        if r['mutant']:
            continue
        d = per[k + (r['uniprot'],)]
        d['types'][r['resname']] += 1
        d['entries'].add(r['pdb'])
    # contact kind per entry = strongest kind in that entry (as a40), then majority per protein
    best = {}
    for r in rows:
        if r['mutant']:
            continue
        ke = (r['family'], r['subsite'], r['position'], r['uniprot'], r['pdb'])
        if ke not in best or rank[r['kind']] > rank[best[ke]]:
            best[ke] = r['kind']
    for ke, kd in best.items():
        per[ke[:4]]['kinds'][kd] += 1
    with open(WORK / 'per_protein.tsv', 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['family', 'subsite', 'position', 'uniprot', 'n_entries', 'residue_types', 'kinds_per_entry'])
        for k in sorted(per):
            d = per[k]
            w.writerow([*k, len(d['entries']), ','.join(f'{t}:{n}' for t, n in d['types'].most_common()),
                        ','.join(f'{t}:{n}' for t, n in d['kinds'].most_common())])
    agg = defaultdict(lambda: dict(votes=Counter(), kinds=Counter(), n=0))
    for k, d in per.items():
        pk = k[:3]
        types = list(d['types'])
        for t in types:
            agg[pk]['votes'][t] += 1.0 / len(types)
        agg[pk]['kinds'][d['kinds'].most_common(1)[0][0]] += 1
        agg[pk]['n'] += 1
    with open(WORK / 'per_protein_counts.tsv', 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['family', 'subsite', 'position', 'class', 'n_entries_all', 'n_entries_native', 'n_proteins',
                    'votes', 'protein_kind_majority'])
        for pk in sorted(agg):
            a = agg[pk]
            nat = len({e for kk, d in per.items() if kk[:3] == pk for e in d['entries']})
            cls = CLASS.get(pk[2])
            if cls is None:   # K182: class from its own majority contact kind
                top = a['kinds'].most_common(1)[0][0]
                cls = 'H' if top == 'sidechain_polar' else 'Al'
            w.writerow([*pk, cls, len(entries_all[pk]), nat, a['n'],
                        ', '.join(f'{t} {v:g}' for t, v in a['votes'].most_common()),
                        ', '.join(f'{t} {v}' for t, v in a['kinds'].most_common())])
    print(open(WORK / 'per_protein_counts.tsv').read())


if __name__ == '__main__':
    main()
