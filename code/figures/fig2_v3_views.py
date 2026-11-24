"""Figure 2 v3: viewing-direction search for the two active-site panels, with a montage for the choice by eye.

Method (fig2_v3_common.view_search): 30 000 directions on a Fibonacci sphere. Kept: -1 ring normal at
|n.z| 0.45-0.97 (neither edge-on nor flat), +1 and +2 rings |n.z| >= 0.35, aromatic rings |n.z| >= 0.35.
Score: -3 x on-page crossings of bonds / H-bond dashes of different groups, -1 x atom pairs of different groups
closer than 1.3 A on the page, -2.5 x residue atoms or bond midpoints inside a sugar ring outline (+0.6 A),
+1 x smallest distance between label anchors (capped 4 A), small terms for ring tilts and for a page aspect
near 1.15. The -1 -> +2 vector is put along the page x axis. The 12 best directions at least 20 degrees apart
are drawn as 2D projections in work/views_<fam>.png.
Run: ~/miniconda3/bin/python code/figures/fig2_v3_views.py
"""
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parent))
import fig2_v3_common as C

for fam in ('GH13', 'GH31'):
    info = C.prepare(fam)
    print(fam, [(rd['pos'], rd['resname'], rd['seq'], rd['cls']) for rd in info['residues']])
    print('  hbonds', [(h['res']['pos'], h['a'], h['unit']['sub'], h['b'], round(h['d'], 2)) for h in info['hbonds']])
    picked, center, right = C.view_search(info)
    fig, axs = plt.subplots(3, 4, figsize=(16, 12))
    for ax, (k, (s, ez, rep)) in zip(axs.ravel(), enumerate(picked)):
        B = C.basis(ez, right)
        for rd in info['residues']:
            col = C.CAT_USED if rd['role'] else C.COL.get(rd['cls'], 'grey')
            for a, b in C.bonds(rd['X']):
                p = (np.array([rd['X'][a], rd['X'][b]]) - center) @ B.T
                ax.plot(p[:, 0], p[:, 1], color=col, lw=2.5)
            t = (C.tip(rd) - center) @ B.T
            ax.text(t[0], t[1], rd['pos'], fontsize=8)
        SX = info['sugarX']; ks = list(SX)
        for i in range(len(ks)):
            for j in range(i + 1, len(ks)):
                if np.linalg.norm(SX[ks[i]] - SX[ks[j]]) < 1.9:
                    p = (np.array([SX[ks[i]], SX[ks[j]]]) - center) @ B.T
                    ax.plot(p[:, 0], p[:, 1], color=C.SUGAR_C, lw=4)
        for sub, (c, nn) in info['rings'].items():
            p = (c - center) @ B.T
            ax.text(p[0], p[1], sub, fontsize=10, weight='bold', ha='center', va='center')
        for h in info['hbonds']:
            p = (np.array([h['res']['X'][h['a']], h['unit']['X'][h['b']]]) - center) @ B.T
            ax.plot(p[:, 0], p[:, 1], 'k:', lw=1)
        ax.set_aspect('equal'); ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f'{k}: {np.round(ez, 3).tolist()}\n{rep}', fontsize=6)
        print(f'  {k}: ez={np.round(ez, 4).tolist()} {rep}')
    fig.tight_layout()
    fig.savefig(C.WORK / f'views_{fam}.png', dpi=70)
    plt.close(fig)
