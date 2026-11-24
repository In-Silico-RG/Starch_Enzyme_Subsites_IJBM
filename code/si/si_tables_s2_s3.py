"""Supplementary Tables S2 and S3 and the water-site numbers of Section S2, from the package data.

Written for this package (2026-11-24): in the project these numbers were computed interactively; this script
reproduces them from the same tables.
Inputs:  data/calc/03_frames/conserved_clusters.tsv, conserved_cluster_members.tsv, superposition.tsv (calc. 3);
         data/calc/01_propka/propka_all.tsv (calc. 1); data/db/enzimas_experimental.sqlite (entry.resolution_A).
Output:  printed tables.
S2 sites: water clusters occupied in >= half of the superposed entries, nearest catalytic oxygen <= 3.6 A.
S2 resolution table: GH13 W3/W4, every GH13 entry of superposition.tsv binned by resolution. As in the SI (412
entries), this includes 4OKD, whose chains were rejected by the superposition (it falls in the 2.2-2.6 A band); with
the 411 accepted entries that band reads 82 entries, W3 56.1 %, W4 47.6 %. Set ACCEPTED_ONLY = True for that count.
S3: catalytic Asp/Glu (catalytic_role nuc_or_base / acid, no '?'), 'with ligand' = a sugar or 'other' molecule
listed within 4.0 A of the residue (column sugar_other_ligands_within_4A non-empty); median pKa, number of chains.
"""
import sqlite3
from pathlib import Path
import pandas as pd

P = Path(__file__).resolve().parents[2]
FR = P / 'data/calc/03_frames'

print('== Section S2: conserved water sites at the catalytic residues')
cl = pd.read_csv(FR / 'conserved_clusters.tsv', sep='\t', keep_default_na=False)
w = cl[(cl.kind == 'water') & (cl.fraction_entries.astype(float) >= 0.5)]
for _, r in w.iterrows():
    cat = [x.split(':') for x in r.distance_to_catalytic_A.split(';') if x.startswith(('nuc_or_base_', 'acid_'))]
    near = [(n, float(dd)) for n, dd in cat if float(dd) <= 3.6]
    for n, dd in near:
        print(f'{r.family}\t{r.cluster_id}\t{n}\t{dd:.2f} A\t{r.n_entries} of {r.n_entries_superposed}')

print('\n== Table S2: GH13 W3/W4 occupancy by resolution band')
con = sqlite3.connect(f"file:{P / 'data/db/enzimas_experimental.sqlite'}?mode=ro", uri=True)
res = dict(con.execute('select pdb, resolution_A from entry'))
sup = pd.read_csv(FR / 'superposition.tsv', sep='\t', keep_default_na=False)
ACCEPTED_ONLY = False
g13 = sup[sup.family == 'GH13']
ents = sorted(set(g13[g13.status.str.startswith('accepted')].pdb if ACCEPTED_ONLY else g13.pdb))
mem = pd.read_csv(FR / 'conserved_cluster_members.tsv', sep='\t', keep_default_na=False)
occ = {c: set(mem[mem.cluster_id == c].pdb) for c in ('GH13_W3', 'GH13_W4')}
bands = [('<1.6', 0, 1.6), ('1.6-1.9', 1.6, 1.9), ('1.9-2.2', 1.9, 2.2), ('2.2-2.6', 2.2, 2.6), ('>=2.6', 2.6, 99)]
print(f'{len(ents)} GH13 entries in superposition.tsv' + (' (accepted only)' if ACCEPTED_ONLY else ''))
for lab, lo, hi in bands:
    e = [p for p in ents if res.get(p) is not None and lo <= float(res[p]) < hi]
    print(f'{lab}\t{len(e)}\tW3 {100 * len(occ["GH13_W3"] & set(e)) / len(e):.1f}%\tW4 {100 * len(occ["GH13_W4"] & set(e)) / len(e):.1f}%')

print('\n== Table S3: median computed pKa (PROPKA) of the catalytic residues')
d = pd.read_csv(P / 'data/calc/01_propka/propka_all.tsv', sep='\t', dtype=str, keep_default_na=False)
d = d[d.catalytic_role.isin(['nuc_or_base', 'acid']) & d.family.isin(['GH13', 'GH31', 'GH57'])].copy()
d['pKa'] = d.pKa.astype(float)
d['with_ligand'] = d.sugar_other_ligands_within_4A != ''
for (fam, role), g in d.groupby(['family', 'catalytic_role']):
    wo, wi = g[~g.with_ligand].pKa, g[g.with_ligand].pKa
    print(f'{fam}\t{role}\t{g.resname.mode()[0]}\twithout {wo.median():.2f} ({len(wo)})\twith {wi.median():.2f} ({len(wi)})')
