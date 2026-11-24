#!/usr/bin/env python3
"""LaTeX tables for the short paper from the rule test (a44).

Input:  data/analysis/44_rule_test/C_identity_T20.tsv (per position, one vote per 50% sequence group),
        data/analysis/44_rule_test/D_tests.tsv, B2_split_half_summary.tsv, E_inner_vs_outer_H.tsv.
Output: data/analysis/44_rule_test/tab_rule.tex, tab_tests.tex
"""
import csv, os
D = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'analysis', '44_rule_test')
rd = lambda f: list(csv.DictReader(open(os.path.join(D, f)), delimiter='\t'))
THREE = dict(A='Ala', R='Arg', N='Asn', D='Asp', C='Cys', Q='Gln', E='Glu', G='Gly', H='His', I='Ile', L='Leu',
             K='Lys', M='Met', F='Phe', P='Pro', S='Ser', T='Thr', W='Trp', Y='Tyr', V='Val')
fmt = lambda x: ('%g' % float(x))

def pos(p):
    return THREE[p[0]] + p[1:]

rows = rd('C_identity_T20.tsv')
order = {'H': 0, 'Ar': 1, 'Al': 2}
rows.sort(key=lambda r: (order.get(r['class_cluster50'], 3), r['family'], -float(r['m50']), r['position']))
L = [r'\begin{table}[h]', r'\centering\small',
     r'\caption{Every reference position contacted in at least 20\% of the entries with a sugar unit at subsites $-3$ to $+3$ '
     r'(GH13 reference 1AMY/1BG9, GH31 reference 2QMJ). Class: contact kind made in most sequence groups (H, side-chain hydrogen '
     r'bond; Ar, aromatic apolar; Al, aliphatic apolar). Groups: proteins sharing $\ge$50\% sequence identity; residues: votes '
     r'per group (a mixed group splits its vote); $m$: modal fraction (eq~\ref{eq:modal}); Arom.: aromatic fraction. '
     r'$^\ast$catalytic residue. Engineered mutants removed. Source: \texttt{44\_rule\_test/C\_identity\_T20.tsv}.}',
     r'\label{tab:rule}', r'\begin{tabular}{lllrp{4.6cm}rr}', r'\toprule',
     r'Class & Family & Position (subsites) & Groups & Residues (group votes) & $m$ & Arom.\\', r'\midrule']
last = None
for r in rows:
    c = r['class_cluster50']
    c = 'none' if c.startswith('tie') else c
    if last is not None and c != last:
        L.append(r'\midrule')
    last = c
    res = ', '.join(THREE[t.split()[0]] + ' ' + fmt(t.split()[1]) for t in r['types50'].split(', '))
    star = r'$^\ast$' if r['catalytic'] == 'True' else ''
    sub = r['subsites'].replace('-', '$-$').replace('+', '$+$')
    L.append(f"{c} & {r['family']} & {pos(r['position'])}{star} ({sub}) & {int(float(r['n_cl50']))} & {res} & "
             f"{float(r['m50']):.2f} & {float(r['arom50']):.2f}\\\\")
L += [r'\bottomrule', r'\end{tabular}', r'\end{table}']
open(os.path.join(D, 'tab_rule.tex'), 'w').write('\n'.join(L) + '\n')

T = {(r['setting'], r['statistic']): r for r in rd('D_tests.tsv')}
def line(label, setting, stat, ci=True):
    r = T[(setting, stat)]
    c = f"[{float(r['ci_lo']):.2f}, {float(r['ci_hi']):.2f}]" if ci and r['ci_lo'] else '--'
    p = float(r['p_one_sided'])
    ps = '$<$0.001' if p < 0.001 else f'{p:.3f}'
    return f"{label} & {stat.replace('-', '$-$')} & {r['n_H']}/{r['n_Ar']}/{r['n_Al']} & {float(r['observed']):+.2f} & {c} & {ps}\\\\"
S = [r'\begin{table}[h]', r'\centering\small',
     r'\caption{Difference in mean modal fraction $m$ between classes (H$-$Al, H$-$Ar) and between the aromatic fraction of Ar '
     r'positions and $m$ of Al positions (Ar(arom)$-$Al); one-sided $p$ from 10\,000 permutations of the class labels over positions; '
     r'95\% interval from 2\,000 bootstrap resamples of sequence groups (conservative, see text). Positions contacted in $\ge$20\% of '
     r'entries; one vote per 50\% sequence group. Source: \texttt{44\_rule\_test/D\_tests.tsv}.}',
     r'\label{tab:tests}', r'\begin{tabular}{llrrlr}', r'\toprule',
     r'Set & Statistic & $n$ H/Ar/Al & Difference & 95\% interval & $p$\\', r'\midrule']
for lab, st in [('All positions', 'T20 cl50 class=cluster all'), ('Without catalytic', 'T20 cl50 class=cluster no-catalytic'),
                ('Held-out groups', 'HELDOUT T20 cl50 class=cluster(bound) all'),
                ('Outer subsites', 'OUTER T20 cl50 class=cluster(outer contacts) all')]:
    for stat in ('H-Al', 'H-Ar', 'Ar(arom)-Al'):
        S.append(line(lab, st, stat))
    S.append(r'\addlinespace')
sp = {(r['T'], r['subset'], r['statistic']): r for r in rd('B2_split_half_summary.tsv')}
r = sp[('0.2000', 'all', 'H-Al')]
S.append(f"Split halves (1\\,000) & H$-$Al & -- & {float(r['median']):+.2f} (median) & [{float(r['p2_5']):.2f}, {float(r['p97_5']):.2f}] & "
         f"$>$0 in {100*float(r['frac_positive']):.0f}\\%\\\\")
e = rd('E_inner_vs_outer_H.tsv')[0]
S.append(f"H at $-1$/$+1$ vs outer & H$-$H & {e['n_inner']}/{e['n_outer']} & {float(e['diff']):+.2f} & -- & {float(e['p_one_sided']):.4f}\\\\")
S += [r'\bottomrule', r'\end{tabular}', r'\end{table}']
open(os.path.join(D, 'tab_tests.tex'), 'w').write('\n'.join(S) + '\n')
print('\n'.join(L[-12:])); print('\n'.join(S[8:]))
