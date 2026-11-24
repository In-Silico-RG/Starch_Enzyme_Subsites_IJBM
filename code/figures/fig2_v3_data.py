"""Figure 2 v3, panel c: how fixed each contact position is across proteins.

Data: data/analysis/44_rule_test/C_identity_T20.tsv, 36 positions (GH13 17, GH31 19), reference numbering
(GH13 1AMY, GH31 2QMJ). Per position: m50 = fraction of the 50 %-identity sequence groups (one vote per group)
whose residue at the position is the modal type (modal50); n_cl50 = number of groups with a residue there;
class_cluster50 = contact class (H, Ar, Al) by the same group vote; subsites = subsites contacted.

Drawing: one dot per position at height m50 (axis 0-1), with a stem from 0. Positions grouped by family, then by
subsite group (-1/+1 = contacts subsite -1 or +1; outer = only -2, -3, +2, +3) and by class (H, Ar, Al; one
tie at GH31 Thr544 drawn grey); within a group by m50 then by number. Catalytic residues drawn as open rings
around the dot; positions with fewer than 5 sequence groups (n_cl50 < 5) drawn as open dots. Below each dot: reference residue and number; above: modal residue type (one-letter) when it
differs from the reference residue. Okabe-Ito colours as in panels a and b.

Writes figures/candidates/fig2_v3/fig2_v3_c.{png,pdf} (7 in wide, 300 dpi).
Run: ~/miniconda3/bin/python code/figures/fig2_v3_data.py
"""
import csv
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fig2_v3_common as C

FONT = 'Liberation Sans'
plt.rcParams.update({'font.family': FONT, 'font.size': 7, 'pdf.fonttype': 42, 'axes.linewidth': 0.6,
                     'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
                     'mathtext.fontset': 'custom', 'mathtext.it': 'Liberation Sans:italic',
                     'mathtext.rm': 'Liberation Sans'})
ONE = {v.upper(): k for k, v in C.ONE3.items()}
CLS_ORDER = ['H', 'Ar', 'Al', 'tie']
NMIN = 5


def load():
    rows = list(csv.DictReader(open(C.D44 / 'C_identity_T20.tsv'), delimiter='\t'))
    for r in rows:
        subs = r['subsites'].split(',')
        r['inner'] = any(s in ('-1', '+1') for s in subs)
        r['cls'] = r['class_cluster50'] if r['class_cluster50'] in ('H', 'Ar', 'Al') else 'tie'
        r['m'] = float(r['m50'])
        r['n'] = int(float(r['n_cl50']))
        r['num'] = int(r['position'][1:])
        r['cat'] = r['catalytic'] == 'True'
    return rows


def layout(rows):
    """x positions with gaps between groups; returns list of (x, row) and group spans."""
    xs, spans = [], []
    x = 0.0
    for fam in ('GH13', 'GH31'):
        fam_start = x
        for inner in (True, False):
            g_start = x
            for cls in CLS_ORDER:
                sel = [r for r in rows if r['family'] == fam and r['inner'] == inner and r['cls'] == cls]
                sel.sort(key=lambda r: (-r['m'], r['num']))
                if not sel:
                    continue
                for r in sel:
                    xs.append((x, r)); x += 1
                x += 0.35
            spans.append((fam, inner, g_start, x - 1.35))
            x += 0.6
        x += 0.9
    return xs, spans


H_C = 2.05      # inches


def draw(fig):
    """Draw panel c into a figure or subfigure 7 in wide and H_C in high."""
    rows = load()
    assert len(rows) == 36, len(rows)
    xs, spans = layout(rows)
    ax = fig.add_axes((0.065, 0.30, 0.925, 0.52))
    for y in (0.25, 0.5, 0.75):
        ax.axhline(y, color='0.90', lw=0.5, zorder=0)
    ax.axhline(1.0, color='0.75', lw=0.6, ls=(0, (2, 2)), zorder=0)
    for x, r in xs:
        col = C.COL.get(r['cls'], '0.55')
        ax.plot([x, x], [0, r['m']], color=col, lw=1.4, solid_capstyle='butt', zorder=1)
        if r['n'] >= NMIN:
            ax.scatter([x], [r['m']], s=22, color=col, zorder=3, edgecolor='white', linewidth=0.4)
        else:                                       # fewer than NMIN sequence groups: open dot
            ax.scatter([x], [r['m']], s=20, facecolor='white', edgecolor=col, linewidth=1.0, zorder=3)
        if r['cat']:
            ax.scatter([x], [r['m']], s=58, facecolor='none', edgecolor=C.CAT_USED, linewidth=0.9, zorder=4)
        mod = r['modal50']
        ref1 = r['position'][0]
        if mod != ref1:
            if r['m'] > 0.9:
                ax.text(x + 0.32, r['m'], mod, ha='left', va='center', fontsize=6.5, color='0.25')
            else:
                ax.text(x, r['m'] + 0.07, mod, ha='center', va='bottom', fontsize=6.5, color='0.25')
        ax.text(x, -0.06, f'{C.ONE3[ref1]}{r["num"]}', rotation=90, ha='center', va='top', fontsize=6.5,
                color='0.15')
    ax.set_ylim(0, 1.2)
    ax.set_xlim(xs[0][0] - 0.8, xs[-1][0] + 0.8)
    ax.set_yticks([0, 0.5, 1]); ax.set_yticklabels(['0', '0.5', '1'])
    ax.set_xticks([])
    ax.set_ylabel('modal fraction $m$', fontsize=7, labelpad=2)
    for s in ('top', 'right', 'bottom'):
        ax.spines[s].set_visible(False)
    ax.axhline(0, color='0.3', lw=0.6, zorder=2)
    # group headers
    for fam, inner, a, b in spans:
        yy = 1.17
        ax.plot([a - 0.3, b + 0.3], [yy, yy], color='0.35', lw=0.6, clip_on=False)
        ax.text((a + b) / 2, yy + 0.03, 'subsites −1/+1' if inner else 'outer subsites', ha='center',
                va='bottom', fontsize=6.5, color='0.25', clip_on=False)
    for fam in ('GH13', 'GH31'):
        fx = [x for x, r in xs if r['family'] == fam]
        ax.text((min(fx) + max(fx)) / 2, 1.45, fam, ha='center', va='bottom', fontsize=7.5, weight='bold',
                transform=ax.transData, clip_on=False)
    # legend
    from matplotlib.lines import Line2D
    h = [Line2D([], [], marker='o', ls='', color=C.COL['H'], ms=4.5, label='H (hydrogen bond)'),
         Line2D([], [], marker='o', ls='', color=C.COL['Ar'], ms=4.5, label='Ar (aromatic residue)'),
         Line2D([], [], marker='o', ls='', color=C.COL['Al'], ms=4.5, label='Al (aliphatic)'),
         Line2D([], [], marker='o', ls='', color='0.55', ms=4.5, label='tie'),
         Line2D([], [], marker='o', ls='', mfc='white', mec='0.4', ms=4.5, mew=1.0, label='< 5 sequence groups'),
         Line2D([], [], marker='o', ls='', mfc='none', mec=C.CAT_USED, ms=6.5, mew=0.9, label='catalytic')]
    fig.legend(handles=h, loc='lower center', bbox_to_anchor=(0.53, -0.005), ncol=6, frameon=False, fontsize=6.5,
               handletextpad=0.3, columnspacing=1.2)
    fig.text(0.005, 0.97, 'c', fontsize=9, weight='bold', va='top')
    return xs


def main():
    fig = plt.figure(figsize=(7.0, H_C))
    xs = draw(fig)
    fig.savefig(C.OUT / 'fig2_v3_c.png', dpi=300)
    fig.savefig(C.OUT / 'fig2_v3_c.pdf')
    plt.close(fig)
    for x, r in xs:
        print(r['family'], r['position'], r['subsites'], r['cls'], r['modal50'], r['m'], r['n'])


if __name__ == '__main__':
    main()
