"""Figure 2 v3, assembled: a (GH13, PDB 1G94) and b (GH31, PDB 3WEO) active sites side by side, a shared
legend, and c (modal fraction m of the 36 contact positions). 7 in wide, 300 dpi PNG and PDF.

Panels a and b come from fig2_v3_struct.py (renders cached in fig2_v3/work/ by a hash of the settings),
panel c from fig2_v3_data.py. Also writes the single panels (fig2_v3_a/b/c).
Run: ~/miniconda3/bin/python code/figures/fig2_v3_figure.py [--rerender]
"""
import sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parent))
import fig2_v3_common as C
import fig2_v3_struct as S
import fig2_v3_data as D

W = 7.0
GAP = W - 2 * S.PANEL_W
HP = S.PANEL_W * S.ASPECT + S.TITLE_H


def main():
    panels = S.main()                    # also writes fig2_v3_a, fig2_v3_b
    D.main()                             # writes fig2_v3_c
    H = HP + S.LEG_H + 0.08 + D.H_C
    fig = plt.figure(figsize=(W, H))
    top, leg, bot = fig.subfigures(3, 1, height_ratios=[HP, S.LEG_H + 0.08, D.H_C])
    fa, fb = top.subfigures(1, 2, width_ratios=[1, 1], wspace=GAP / W)
    for sf, fam, letter in ((fa, 'GH13', 'a'), (fb, 'GH31', 'b')):
        info, zp, tp = panels[fam]
        S.panel_into(sf, info, zp, tp, letter, HP)
    S.struct_legend(leg, 0.5, 0.55)
    D.draw(bot)
    fig.savefig(C.OUT / 'fig2_v3.png', dpi=300)
    fig.savefig(C.OUT / 'fig2_v3.pdf', dpi=300)
    plt.close(fig)
    print('figure', W, 'x', round(H, 2), 'in')


if __name__ == '__main__':
    main()
