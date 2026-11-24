"""Candidate A: local refinement of the viewing direction (in the family frame) for the hydroxyl-map panels,
with fig2_v3_common.view_search (directions within 40 deg of the fig2_v3 direction carried into the frame, extra
penalty on residue-residue overlaps w_rr = 0.6), then the -2 ring tilt |n.z| added to the score (GH13 panel
shows -2, which the fig2_v3 search did not include). Montage of the 9 best in descriptive/work/A_views_<fam>.png.
Run: ~/miniconda3/bin/python code/figures/figA_views.py GH13
"""
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parent))
import fig2_v3_common as C
import figdesc_common as F
import figA_map_on_molecule as A


def shim(info):
    res = []
    for r in info['residues']:
        res.append(dict(X=r['X'], role=None, cls='H', resname=r['resname'] if r['sc'] else 'MC', pos=r['label']))
    units, rings = [], {}
    for u in info['units']:
        X = {n: x for n, x in zip(u['ring'], u['ringX'])}
        units.append(dict(sub=u['sub'], X=X, ring=u['ring']))
        rings[u['sub']] = C.plane(u['ringX'])
    hb = [dict(res=dict(X={'p': h['Y']}), a='p', unit=dict(X={'o': h['o']['X']}), b='o') for h in info['hbonds']]
    sx = {k: v for k, v in info['sugarX'].items()}
    return dict(residues=res, units=units, rings=rings, sugarX=sx, hbonds=hb)


fam = sys.argv[1]
info = A.prepare(fam)
S = shim(info)
start = F.view_dir(fam)
D = C.fib_dirs(60000)
ANG = float(sys.argv[2]) if len(sys.argv) > 2 else 40.0
D = D[D @ start > np.cos(np.radians(ANG))]
picked, center, right = C.view_search(S, dirs=D, top=40, sep=6, w_rr=0.6)
rows = []
for s, ez, rep in picked:
    t2 = abs(S['rings']['-2'][1] @ ez) if '-2' in S['rings'] else 1.0
    rows.append((s + 3.0 * min(t2, 0.6), ez, rep, t2))
rows.sort(key=lambda t: -t[0])
fig, axs = plt.subplots(3, 3, figsize=(13, 12))
for ax, (k, (s, ez, rep, t2)) in zip(axs.ravel(), enumerate(rows[:9])):
    B = C.basis(ez, right)
    for r in S['residues']:
        for a, b in C.bonds(r['X']):
            p = (np.array([r['X'][a], r['X'][b]]) - center) @ B.T
            ax.plot(p[:, 0], p[:, 1], color='tab:blue', lw=2)
        t = (C.tip(r) - center) @ B.T
        ax.text(t[0], t[1], r['pos'], fontsize=7)
    ks = list(S['sugarX'])
    for i in range(len(ks)):
        for j in range(i + 1, len(ks)):
            if np.linalg.norm(S['sugarX'][ks[i]] - S['sugarX'][ks[j]]) < 1.9:
                p = (np.array([S['sugarX'][ks[i]], S['sugarX'][ks[j]]]) - center) @ B.T
                ax.plot(p[:, 0], p[:, 1], color=C.SUGAR_C, lw=4)
    for o in info['oxy']:
        p = (o['X'] - center) @ B.T
        ax.scatter([p[0]], [p[1]], s=60, color=A.STATE_COL[o['state']], edgecolor='k', zorder=5)
    for sub, (cc, nn) in S['rings'].items():
        p = (cc - center) @ B.T
        ax.text(p[0], p[1], sub, fontsize=10, weight='bold', ha='center')
    ax.set_aspect('equal'); ax.set_xticks([]); ax.set_yticks([])
    ang = np.degrees(np.arccos(np.clip(ez @ start, -1, 1)))
    ax.set_title(f'{k}: {np.round(ez, 4).tolist()} ({ang:.0f} deg) t-2={t2:.2f}\n{rep}', fontsize=6)
    print(k, np.round(ez, 4).tolist(), round(ang), round(t2, 2), rep)
fig.tight_layout(); fig.savefig(F.WORK / f'A_views_{fam}.png', dpi=70)
