"""Ring puckering of density-supported sugar units by subsite (GH13, GH31, GH14, GH15).

Each census linkage (data/census/bound_linkages_checked.tsv) gives the Cremer-Pople theta of its
non-reducing (C1 side) and reducing (O4 side) ring. A ring is matched to a subsite unit of
a33_subsites.py by (pdb, chain, component, residue number); only units with RSCC >= 0.80 that enter
the subsite tables are used. theta is folded to min(theta, 180 - theta), the distance from the nearest
chair; a ring is called distorted above 20 degrees. Counts are given per ring, per entry and per protein.
Usage: python a42_ring_distortion_subsite.py   -> data/analysis/42_ring_distortion/
"""
import csv
import sqlite3
from collections import defaultdict
from pathlib import Path

P = Path(__file__).resolve().parents[2]
OUT = P / 'data/analysis/42_ring_distortion'
OUT.mkdir(parents=True, exist_ok=True)
CUT = 20.0

con = sqlite3.connect(f"file:{P / 'data/db/enzimas_experimental.sqlite'}?mode=ro", uri=True)
UNIPROT = {pdb: ((acc or '').split(';')[0].strip() or pdb) for pdb, acc in con.execute('select pdb, uniprot from entry')}


def fnum(v):
    try:
        return float(v)
    except ValueError:
        return None


units = {}
for u in csv.DictReader(open(P / 'data/analysis/33_subsites/units.tsv'), delimiter='\t'):
    r = fnum(u['rscc'])
    if u['subsite'] and u['in_tables'] == 'True' and r is not None and r >= 0.80:
        units[(u['pdb'], u['chain'], u['comp'], str(u['seq']))] = (u['family'], u['subsite'])


def parse(label):          # 'D:GLC2' -> ('D', 'GLC', '2')
    ch, rest = label.split(':')
    i = next(k for k, c in enumerate(rest) if c.isdigit() or c == '-')
    return ch, rest[:i], rest[i:]


rows = []
for L in csv.DictReader(open(P / 'data/census/bound_linkages_checked.tsv'), delimiter='\t'):
    for side, lab, th in (('nonreducing', 'nonreducing_ring(C1)', 'theta_nonred_deg'),
                          ('reducing', 'reducing_ring(O4)', 'theta_red_deg')):
        try:
            ch, comp, seq = parse(L[lab])
        except (ValueError, StopIteration):
            continue
        key = (L['pdb'], ch, comp, seq)
        t = fnum(L[th])
        if key in units and t is not None:
            fam, sub = units[key]
            tf = min(t, 180.0 - t)
            rows.append(dict(family=fam, subsite=sub, side=side, pdb=L['pdb'], protein=UNIPROT.get(L['pdb'], L['pdb']),
                             ring=L[lab], theta=t, theta_folded=round(tf, 1), distorted=tf > CUT))

# one ring is counted once even if it appears in two linkages (as reducing ring of one, non-reducing of the next)
seen, uniq = set(), []
for r in rows:
    k = (r['pdb'], r['ring'])
    if k not in seen:
        seen.add(k)
        uniq.append(r)
with open(OUT / 'rings.tsv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(uniq[0]), delimiter='\t')
    w.writeheader()
    w.writerows(uniq)

g = defaultdict(list)
for r in uniq:
    g[(r['family'], '-1' if r['subsite'] == '-1' else 'other')].append(r)
with open(OUT / 'summary.tsv', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['family', 'subsite', 'rings', 'distorted', 'entries', 'entries_with_distorted', 'proteins',
                'proteins_with_distorted', 'distorted_entries'])
    for k in sorted(g):
        v = g[k]
        d = [r for r in v if r['distorted']]
        w.writerow([*k, len(v), len(d), len({r['pdb'] for r in v}), len({r['pdb'] for r in d}),
                    len({r['protein'] for r in v}), len({r['protein'] for r in d}),
                    ' '.join(sorted({r['pdb'] for r in d}))])
print(open(OUT / 'summary.tsv').read())
