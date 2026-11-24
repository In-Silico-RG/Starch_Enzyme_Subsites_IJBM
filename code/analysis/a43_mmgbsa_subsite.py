#!/usr/bin/env python3
"""Per-residue MM/GBSA energy of the four GAA iminosugar complexes at the conserved GH31 subsite -1 positions.

Input:  data/calc/09_mmgbsa/series/<PDB>_<COMP>/{complex_in.pdb, FINAL_DECOMP_MMPBSA.dat, FINAL_RESULTS_MMPBSA.dat}
        (c09_series.py, job 33445 on GUANE <cluster-account>, 2026-11-24), and
        data/analysis/40_contact_types/contact_types.tsv (a40_contact_types.py: 2QMJ reference position ->
        5NN5 author residue number and contact class).
Numbering: complex_in.pdb numbers the protein 1..N in file order (tleap). For 5NN5 the author number is the
        sequential number + 93 over the subsite region (checked: the residue name at seq+93 equals the author
        residue name at all nine positions of contact_types.tsv). The other three complexes are mapped to 5NN5
        by aligning their residue-name sequences (difflib), then to author numbers through 5NN5.
Output: data/analysis/43_mmgbsa_subsite/{per_position.tsv, per_class.tsv, series_vs_potency.tsv}
"""
import os, csv, difflib, statistics
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
SER = os.path.join(ROOT, 'data/calc/09_mmgbsa/series')
OUT = os.path.join(ROOT, 'data/analysis/43_mmgbsa_subsite'); os.makedirs(OUT, exist_ok=True)
CX = [('5NN5', 'NOJ', 6.82), ('5NN6', 'MIG', 6.46), ('9GTC', 'A1IO7', 7.72), ('9GTW', 'A1IPC', 5.68)]
CLASS = {'D327': 'H', 'H600': 'H', 'D443': 'H', 'D542': 'H', 'W441': 'Ar', 'W539': 'Ar', 'F575': 'Ar',
         'I328': 'Al', 'I364': 'Al'}  # as in Table tab:rule of grip_finding.tex
THREE = {'HID': 'HIS', 'HIE': 'HIS', 'HIP': 'HIS', 'GLH': 'GLU', 'ASH': 'ASP', 'CYX': 'CYS', 'LYN': 'LYS'}

def seq(pdb):
    out, last = [], None
    for l in open(pdb):
        if l.startswith('ATOM'):
            k = (l[21], int(l[22:26]))
            if k != last: out.append((int(l[22:26]), THREE.get(l[17:20], l[17:20]))); last = k
    return out

def decomp(f):
    e, on = {}, False
    for l in open(f):
        if l.startswith('Total Energy'): on = True; continue
        if l.startswith('Sidechain') or l.startswith('Backbone'): on = False
        p = l.split(',')
        if on and len(p) > 17 and p[0][:3].isalpha() and p[0][3:4] == ' ':
            e[int(p[0].split()[1])] = (p[0].split()[0], float(p[17]))
    return e

ref = {}
for r in csv.DictReader(open(os.path.join(ROOT, 'data/analysis/40_contact_types/contact_types.tsv')), delimiter='\t'):
    if r.get('pdb', r.get(list(r)[3])) == '5NN5':
        pass
rows = [l.rstrip('\n').split('\t') for l in open(os.path.join(ROOT, 'data/analysis/40_contact_types/contact_types.tsv'))]
for p in rows:
    if len(p) > 7 and p[0] == 'GH31' and p[3] == '5NN5':
        ref[p[2]] = (int(p[6]), p[7])
assert set(ref) == set(CLASS), ref

s5 = seq(os.path.join(SER, '5NN5_NOJ/complex_in.pdb'))
OFF = 93
for pos, (auth, name) in ref.items():
    i = auth - OFF
    assert dict(s5)[i] == name, (pos, auth, dict(s5)[i], name)

per_pos, per_cls, tot = [], [], []
for pdb, comp, pch in CX:
    d = os.path.join(SER, f'{pdb}_{comp}')
    s = seq(os.path.join(d, 'complex_in.pdb'))
    sm = difflib.SequenceMatcher(None, [x[1] for x in s5], [x[1] for x in s], autojunk=False)
    m = {}
    for a, b, n in sm.get_matching_blocks():
        for k in range(n): m[s5[a + k][0]] = s[b + k][0]
    e = decomp(os.path.join(d, 'FINAL_DECOMP_MMPBSA.dat'))
    dg = float([l.split()[2] for l in open(os.path.join(d, 'FINAL_RESULTS_MMPBSA.dat')) if l.startswith('DELTA TOTAL')][0])
    lig = [v[1] for v in e.values() if v[0] == 'LIG'][0]
    prot = sum(v[1] for v in e.values() if v[0] != 'LIG')
    cls = {'H': 0.0, 'Ar': 0.0, 'Al': 0.0}
    for pos in CLASS:
        auth, name = ref[pos]
        j = m.get(auth - OFF)
        nm, en = e.get(j, ('-', 0.0)) if j else ('-', 0.0)
        nm = THREE.get(nm, nm)
        per_pos.append([pdb, comp, pos, CLASS[pos], auth, j, nm, f'{en:.2f}'])
        cls[CLASS[pos]] += en
    per_cls.append([pdb, comp, f'{dg:.2f}', f'{lig:.2f}', f'{prot:.2f}'] + [f'{cls[c]:.2f}' for c in ('H', 'Ar', 'Al')])
    tot.append((pdb, comp, pch, dg))

w = lambda n, h, r: open(os.path.join(OUT, n), 'w').write('\t'.join(h) + '\n' + ''.join('\t'.join(map(str, x)) + '\n' for x in r))
w('per_position.tsv', ['pdb', 'comp', 'pos_2QMJ', 'class', 'author_resnum_5NN5', 'seq_resnum', 'resname', 'dG_kcal_mol'], per_pos)
w('per_class.tsv', ['pdb', 'comp', 'dG_total', 'dG_ligand_half', 'dG_protein_half', 'H', 'Ar', 'Al'], per_cls)
w('series_vs_potency.tsv', ['pdb', 'comp', 'pchembl_median', 'dG_MMGBSA'], tot)

def rank(v): o = sorted(range(len(v)), key=lambda i: v[i]); r = [0] * len(v); [r.__setitem__(i, k) for k, i in enumerate(o)]; return r
x = [t[2] for t in tot]; y = [-t[3] for t in tot]
def pear(a, b):
    ma, mb = statistics.mean(a), statistics.mean(b)
    return sum((i - ma) * (j - mb) for i, j in zip(a, b)) / (sum((i - ma) ** 2 for i in a) * sum((j - mb) ** 2 for j in b)) ** .5
print('Pearson r(pChEMBL, -dG) =', round(pear(x, y), 3), '| Spearman rho =', round(pear(rank(x), rank(y)), 3), '| n = 4')
for r in per_cls: print('\t'.join(r))
for r in per_pos: print('\t'.join(map(str, r)))
