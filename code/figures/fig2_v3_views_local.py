"""Figure 2 v3: local refinement of a viewing direction (directions within 40 degrees of a start direction),
with an extra penalty on residue-residue atom pairs closer than 2 A on the page (w_rr), to separate side
chains that the global search leaves stacked. Montage in work/views_local_<fam>.png.
Run: ~/miniconda3/bin/python code/figures/fig2_v3_views_local.py GH31 -0.7836 0.5235 0.3346
"""
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parent))
import fig2_v3_common as C

fam = sys.argv[1]
start = np.array([float(v) for v in sys.argv[2:5]]); start /= np.linalg.norm(start)
info = C.prepare(fam)
D = C.fib_dirs(60000)
D = D[D @ start > np.cos(np.radians(40))]
picked, center, right = C.view_search(info, dirs=D, top=9, sep=8, w_rr=0.6)
fig, axs = plt.subplots(3, 3, figsize=(13, 12))
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
    ax.set_aspect('equal'); ax.set_xticks([]); ax.set_yticks([])
    ang = np.degrees(np.arccos(np.clip(ez @ start, -1, 1)))
    ax.set_title(f'{k}: {np.round(ez, 4).tolist()} ({ang:.0f} deg)\n{rep}', fontsize=6)
    print(k, np.round(ez, 4).tolist(), round(ang), rep)
fig.tight_layout(); fig.savefig(C.WORK / f'views_local_{fam}.png', dpi=70)
