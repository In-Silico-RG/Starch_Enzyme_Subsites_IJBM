"""Figure 3 of the IJBM manuscript: the map subsite x glucose oxygen per family, with literature results.

Data (read only):
  data/analysis/45_hydroxyl_map/cell_summary.tsv, set 'all_glucoselike' (code/analysis/a45_hydroxyl_map.py):
    state  INV  held (>= half of the 50 % sequence groups) by a side chain invariant over the family's groups
                (m50 = 1 over >= 3 groups; GH14/GH15: over their 2 groups)
           VAR  held, mostly by side chains that are not invariant
           MC   held only by main-chain N/O
           EXP  not held; median open surface in the complex >= 5 A^2 and >= 25 % of the ligand-alone surface
           BUR  not held and not exposed
    number in the cell = frac_cl_held (share of 50 % groups in which the oxygen has a protein N/O within 3.5 A)
  Literature marks: data/analysis/substituted_starch_hydrolysis_literature_v2.md, Sections 6-8 and 11
    (phosphate: Takeda 1983/1986/1987, Abe 1982, Kamasaka 1997; hydroxyethyl: Yoshida 1984; O-methyl: Lemieux 1996,
    van der Burgt 1998, Mischnick 2001; 6-O-Gal / 3-O-GlcNAc: Okada 2000; 6-pyridylamino: Omichi 1988;
    6-iodo/6-amino: Omichi 1992; 4,6-O-benzylidene: Kandra 1997/2000; C6 replacements: Bock & Pedersen 1988;
    deoxy (OH removed, not a substituent; GH31 -1 only): Nishio 2002, Ogawa 2004). Same cells as Table 3.
    mark '+' tolerated, 'x' not tolerated, '+-' differs between enzymes of the family, '-' no data;
    deoxy marks are drawn in italics and are not counted in the agreement (v2 update 2026-11-24).
Rows: O2, O3, O6 (substitution sites) and O4 of the +1 unit (the scissile glycosidic oxygen).
Output: manuscript_ijbm/manuscript/figures/fig3_map.pdf/.png
Run: ~/miniconda3/bin/python code/figures/fig3_map.py
"""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
import pandas as pd

P = Path(__file__).resolve().parents[2]
CS = P / 'data/analysis/45_hydroxyl_map/cell_summary.tsv'
OUT = P / 'manuscript_ijbm/manuscript/figures'
OUT.mkdir(parents=True, exist_ok=True)

SUBS = ['-3', '-2', '-1', '+1', '+2', '+3']
SUBLAB = ['−3', '−2', '−1', '+1', '+2', '+3']
ROWS = ['O2', 'O3', 'O6', 'O4']
# distinct hues (IJBM editorial review 2026-11-24); INV and VAR as in the map-on-molecule figure (orange fixed, blue variable)
FILL = {'INV': '#B5540B', 'VAR': '#5AB4E5', 'MC': '#E9DCA6', 'EXP': '#FFFFFF', 'BUR': '#8C8C8C'}
TXT = {'INV': 'white', 'VAR': 'black', 'MC': 'black', 'EXP': 'black', 'BUR': 'white'}
NAME = {'INV': 'fixed (invariant side chain)', 'VAR': 'variable side chain', 'MC': 'main chain only',
        'EXP': 'free, exposed', 'BUR': 'free, buried'}
TICK, CROSS, MIXED = '✓', '✗', '±'
DISAGREE = '#CC0066'

# (family, subsite, oxygen): (symbol, substituents); Table 3 of the manuscript, same sources
LIT = {
    ('GH13', '-3', 'O3'): (TICK, 'P,Hx'), ('GH13', '-3', 'O6'): (TICK, 'P,Hx,I…'),
    ('GH13', '-2', 'O3'): (MIXED, 'P,Hx'), ('GH13', '-2', 'O6'): (MIXED, 'P,Hx'),
    ('GH13', '-1', 'O2'): (CROSS, 'HE,Me'), ('GH13', '-1', 'O3'): (CROSS, 'P,HE,Me'),
    ('GH13', '-1', 'O6'): (CROSS, 'P,HE,Me'),
    ('GH13', '+1', 'O3'): (CROSS, 'P'), ('GH13', '+1', 'O6'): (MIXED, 'P'),
    ('GH13', '+2', 'O3'): (MIXED, 'P'), ('GH13', '+2', 'O6'): (TICK, 'P,NH₂'), ('GH13', '+3', 'O3'): (TICK, 'P'),
    ('GH14', '-2', 'O3'): (CROSS, 'P'), ('GH14', '-2', 'O6'): (CROSS, 'P'), ('GH14', '-1', 'O3'): (CROSS, 'P'),
    ('GH14', '-1', 'O6'): (CROSS, 'P'), ('GH14', '+1', 'O3'): (CROSS, 'P'), ('GH14', '+1', 'O6'): (TICK, 'P'),
    ('GH15', '-1', 'O2'): (CROSS, 'Me'), ('GH15', '-1', 'O3'): (CROSS, 'P,Me'), ('GH15', '-1', 'O6'): (CROSS, 'P,Me,C6'),
    ('GH15', '+1', 'O2'): (TICK, 'Me'), ('GH15', '+1', 'O3'): (CROSS, 'P,Me'), ('GH15', '+1', 'O6'): (TICK, 'P,C6'),
    ('GH15', '+2', 'O3'): (TICK, 'P'),
    # deoxy (hydroxyl removed): tick = OH dispensable, cross = OH needed, +- = differs between enzymes
    ('GH31', '-1', 'O2'): (TICK, 'd'), ('GH31', '-1', 'O3'): (MIXED, 'd'), ('GH31', '-1', 'O6'): (CROSS, 'd'),
}
DEOXY = {('GH31', '-1', 'O2'), ('GH31', '-1', 'O3'), ('GH31', '-1', 'O6')}
# cells where a definite map prediction (INV -> not tolerated, EXP -> tolerated) contradicts the literature
DIS = {('GH14', '+1', 'O6'), ('GH15', '+1', 'O2')}

d = pd.read_csv(CS, sep='\t', dtype={'subsite': str})
d = d[d.set == 'all_glucoselike']
CELL = {(r.family, r.subsite, r.label): r for r in d.itertuples()}

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7.5, 'pdf.fonttype': 42, 'ps.fonttype': 42})
fig = plt.figure(figsize=(7.4, 5.6))
# axes: [left, bottom, width, height] in figure fraction
AX = {"GH13": [0.06, 0.53, 0.42, 0.39], "GH31": [0.555, 0.53, 0.42, 0.39],
      "GH14": [0.06, 0.20, 0.42, 0.26], "GH15": [0.555, 0.20, 0.42, 0.26]}
TITLE = {'GH13': 'GH13 (29 groups)', 'GH31': 'GH31 (6 groups)',
         'GH14': 'GH14 (2 groups)', 'GH15': 'GH15 (2 groups)'}

for fam, box in AX.items():
    ax = fig.add_axes(box)
    small = fam in ('GH14', 'GH15')
    fs_num = 6.5 if small else 7.5
    fs_lit = 6.2 if small else 7.0
    for j, s in enumerate(SUBS):
        for i, o in enumerate(ROWS):
            y = len(ROWS) - 1 - i
            if o == 'O4' and s != '+1':
                continue
            r = CELL.get((fam, s, o))
            lit = LIT.get((fam, s, o))
            if r is None:
                ax.add_patch(Rectangle((j, y), 1, 1, fill=False, ec='#BBBBBB', lw=0.5, ls=(0, (2, 2))))
                if lit:
                    ax.text(j + 0.5, y + 0.5, f'{lit[0]}{lit[1]}', ha='center', va='center', fontsize=fs_lit)
                continue
            st = r.state
            ax.add_patch(Rectangle((j, y), 1, 1, fc=FILL[st], ec='black', lw=0.6,
                                   hatch='////' if st == 'EXP' else None))
            if st == 'EXP':   # white box behind the text so the hatching does not cross it
                ax.add_patch(Rectangle((j + 0.22, y + 0.47), 0.56, 0.36, fc='white', ec='none'))
            ax.text(j + 0.5, y + 0.66, f'{r.frac_cl_held:.2f}', ha='center', va='center', fontsize=fs_num,
                    color=TXT[st])
            if st == 'EXP' and o != 'O4':
                ax.add_patch(Rectangle((j + 0.08, y + 0.1), 0.84, 0.33, fc='white', ec='none'))
            if o == 'O4':
                ax.text(j + 0.5, y + 0.27, '✂', ha='center', va='center', fontsize=fs_num + 1, color=TXT[st])
            elif lit:
                ax.text(j + 0.5, y + 0.27, f'{lit[0]}{lit[1]}', ha='center', va='center', fontsize=fs_lit,
                        color=TXT[st], fontweight='normal' if (fam, s, o) in DEOXY else 'bold',
                        style='italic' if (fam, s, o) in DEOXY else 'normal')
            else:
                ax.text(j + 0.5, y + 0.27, '–', ha='center', va='center', fontsize=fs_lit,
                        color='#BBBBBB' if st in ('INV', 'VAR', 'BUR') else '#999999')
            if (fam, s, o) in DIS:
                ax.add_patch(Rectangle((j + 0.04, y + 0.04), 0.92, 0.92, fill=False, ec=DISAGREE, lw=2.0))
    ax.set_xlim(0, len(SUBS))
    ax.set_ylim(0, len(ROWS))
    ax.set_xticks([j + 0.5 for j in range(len(SUBS))])
    ax.set_xticklabels(SUBLAB)
    ax.set_yticks([len(ROWS) - 1 - i + 0.5 for i in range(len(ROWS))])
    ax.set_yticklabels(ROWS)
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.axvline(3, color='black', lw=1.6, ymin=0, ymax=1)   # the scissile bond between -1 and +1
    ax.set_title(TITLE[fam], fontsize=8, loc='left', pad=3, fontweight='bold')
    ax.set_xlabel('subsite', fontsize=7, labelpad=1)

handles = [Patch(fc=FILL[k], ec='black', lw=0.6, hatch='////' if k == 'EXP' else None, label=NAME[k])
           for k in ['INV', 'VAR', 'MC', 'EXP', 'BUR']]
handles += [Patch(fc='white', ec=DISAGREE, lw=2.0, label='map and experiment disagree'),
            Patch(fc='white', ec='#BBBBBB', ls=(0, (2, 2)), lw=0.6, label='no density-supported unit or atom')]
fig.legend(handles=handles, loc='lower left', bbox_to_anchor=(0.04, -0.005), ncol=4, frameon=False, fontsize=6.6,
           handlelength=1.6, columnspacing=1.2)
fig.text(0.04, 0.068,
         f'{TICK} tolerated   {CROSS} not tolerated   {MIXED} differs between enzymes   \u2013 no data   \u2702 scissile O4   '
         'd deoxy (OH removed; not counted)\n'
         'P phosphate   HE hydroxyethyl   Me methyl   Hx 6-O-Gal or 3-O-GlcNAc   I\u2026 6-iodo, 6-pyridylamino, 4,6-benzylidene   '
         'NH\u2082 6-amino   C6 C6 replaced',
         fontsize=6.6, linespacing=1.6)
for ext in ('pdf', 'png'):
    fig.savefig(OUT / f'fig3_map.{ext}', dpi=400)
print('written', OUT / 'fig3_map.pdf')
