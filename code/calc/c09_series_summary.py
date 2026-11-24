#!/usr/bin/env python3
"""Calculation 9, series: MM/GBSA Delta G of the GAA iminosugar complexes against their measured potency.

pChEMBL per inhibitor (database read-only): activities in chembl_activity on target CHEMBL2608 (lysosomal
alpha-glucosidase, human) with a pchembl_value, keeping only assays whose description names the human enzyme
("human" together with "lysosomal" or "GAA") and dropping cell assays ("Caco-2"); the median of what remains.
MM/GBSA: DELTA TOTAL of data/calc/09_mmgbsa/series/<PDB>_<COMP>/FINAL_RESULTS_MMPBSA.dat (c09_series.py).
Writes data/calc/09_mmgbsa/series/series_summary.tsv and prints Pearson r and Spearman rho with n.
"""
import os, re, sys, statistics
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from common import db, CALC

SER = [('5NN5', 'NOJ'), ('5NN6', 'MIG'), ('9GTC', 'A1IO7'), ('9GTW', 'A1IPC')]
c = db()
rows = []
for pdb, comp in SER:
    mid = c.execute('select molecule_chembl_id from chembl_pdb_ligand where comp_id=?', (comp,)).fetchone()[0]
    acts = c.execute("select activity_id, standard_type, pchembl_value, document_chembl_id, assay_description from chembl_activity "
                     "where molecule_chembl_id=? and target_chembl_id='CHEMBL2608' and pchembl_value is not null", (mid,)).fetchall()
    keep = [a for a in acts if 'human' in (a[4] or '').lower() and ('lysosomal' in (a[4] or '').lower() or 'gaa' in (a[4] or '').lower())
            and 'caco-2' not in (a[4] or '').lower()]
    p = statistics.median(float(a[2]) for a in keep) if keep else None
    f = os.path.join(CALC, '09_mmgbsa', 'series', f'{pdb}_{comp}', 'FINAL_RESULTS_MMPBSA.dat')
    dg = None
    if os.path.exists(f):
        m = re.search(r'DELTA TOTAL\s+(-?\d+\.\d+)', open(f).read())
        dg = float(m.group(1)) if m else None
    rows.append((pdb, comp, mid, p, ';'.join(f'{a[1]} {a[2]} {a[3]}' for a in keep), dg))
out = os.path.join(CALC, '09_mmgbsa', 'series', 'series_summary.tsv')
with open(out, 'w') as f:
    f.write('pdb\tcomp\tmolecule_chembl_id\tpchembl_median_GAA\tactivities_used(type value document)\tdeltaG_MMGBSA_kcal_mol\n')
    for r in rows:
        f.write('\t'.join('' if x is None else str(x) for x in r) + '\n')
ok = [(r[3], r[5]) for r in rows if r[3] is not None and r[5] is not None]
print(open(out).read())
if len(ok) >= 3:
    x, y = np.array(ok).T
    pr = np.corrcoef(x, y)[0, 1]
    rk = lambda v: np.argsort(np.argsort(v))
    sr = np.corrcoef(rk(x), rk(y))[0, 1]
    print(f'n = {len(ok)}; Pearson r(pChEMBL, DeltaG) = {pr:.2f}; Spearman rho = {sr:.2f} (a better binder has larger pChEMBL and more negative DeltaG, so agreement is a negative r)')
