"""Graphical abstract of the IJBM manuscript, composed from the approved figure material (no new structure renders).

Material reused
  Left, functionalised starch: an alpha-1,4 glucan of 4 units built and depicted by RDKit with the functions of
      fig1_protein_chem.py (build_chain, depict; CoordGen 2D coordinates, nothing placed by hand). Rings filled gold
      as in Fig. 1e; substituents marked on chosen O2/O3/O6 as coloured tags (HP hydroxypropyl, Ac acetyl,
      OSA octenyl succinyl, CM carboxymethyl), the four groups named in the abstract.
  Centre, enzyme cleft: the cached ChimeraX layers of descriptive candidate A (figA_map_on_molecule.py), GH13
      maltogenic alpha-amylase PDB 1QHO in the 1AMY frame, render hash 90e1d703 (cartoon layer faded by 0.70 as in
      figA; stick layer: gold sugar at -2..+2, O2/O3/O6 spheres coloured by their state in the hydroxyl map of
      analysis 45 (vermillion held by a fixed residue, sky blue held by a variable residue, white free), holding
      residues vermillion (fixed) / sky blue (variable)). Subsite numbers and the scissile bond (-1/+1) placed
      at pixel positions read off that render.
  Right, the rule, from the Conclusions of the manuscript: -1/+1 fixed residues, no substituent; -2/+2
      variable residues, depends on the enzyme; outer subsites and O6 mostly free, substituent fits; O6
      tolerated better than O2 or O3 (design rule (i)).
Output: 13.3 x 5.32 cm (5:2), 600 dpi (3142 x 1257 px, above the Elsevier minimum 1328 x 531 px); fonts >= 8 pt.
  candidates -> figures/candidates/graphical_abstract/ga_{A,B,C}.{png,pdf}
  --pick X   -> manuscript_ijbm/submission/graphical_abstract.{png,pdf,tiff}
Run: ~/miniconda3/bin/python code/figures/graphical_abstract.py [--pick B]
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, FancyBboxPatch, Circle, FancyArrowPatch
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fig1_protein_chem as K          # noqa: E402

P = Path(__file__).resolve().parents[2]
WORKA = P / 'figures/candidates/descriptive/work'
OUT = P / 'figures/candidates/graphical_abstract'
SUB = P / 'manuscript_ijbm/submission'
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({'font.family': 'Liberation Sans', 'font.size': 8, 'pdf.fonttype': 42})
CM = 1 / 2.54
WID, HEI = 13.3 * CM, 5.32 * CM
DPI = 600

SUGAR = '#F4B942'
SUGAR_FILL = '#F9DFA6'
O_COL = '#E03C31'
FIXED = '#D55E00'
VAR = '#56B4E9'
FREE = '#FFFFFF'
TXT_FIXED, TXT_VAR = '#B84A00', '#2A86BF'
SUBST = {'HP': '#7B5EA7', 'Ac': '#2C7F7A', 'OSA': '#8C5A2B', 'CM': '#3C3C3C'}
# substituent tags: (unit, oxygen, tag)
TAGS = [(0, 'O6', 'HP'), (1, 'O2', 'Ac'), (1, 'O6', 'OSA'), (2, 'O3', 'CM')]

# ---------------------------------------------------------------- cleft (cached figA layers, hash 90e1d703)
HASH = '90e1d703'
CROP = (470, 600, 1780, 1560)       # x0, y0, x1, y1 in the 2000 x 2000 render
SUBPOS = {'−2': (715, 790), '−1': (995, 1088), '+1': (1395, 1200), '+2': (1615, 1085)}
CUT = (1095, 985)                   # scissors; cut bar across the -1 C1 - N4 (+1) bond
CUTBAR = ((1160, 1092), (1196, 1160))   # across the bond to the bridging N4 of acarbose (scissile position)


def cleft_image():
    cart = plt.imread(WORKA / f'A_cart_GH13_{HASH}.png')[..., :3]
    stk = plt.imread(WORKA / f'A_stk_GH13_{HASH}.png')
    base = 1 - 0.70 * (1 - cart)
    a = stk[..., 3:4]
    img = base * (1 - a) + stk[..., :3] * a
    x0, y0, x1, y1 = CROP
    return img[y0:y1, x0:x1]


def draw_cleft(ax, labels=True, scissors=True, fs_sub=9):
    img = cleft_image()
    x0, y0, x1, y1 = CROP
    ax.imshow(img, extent=(x0, x1, y1, y0), interpolation='lanczos')
    ax.set_xlim(x0, x1); ax.set_ylim(y1, y0); ax.axis('off')
    if labels:
        for t, (x, y) in SUBPOS.items():
            ax.text(x, y, t, fontsize=fs_sub, weight='bold', ha='center', va='center', zorder=5,
                    bbox=dict(boxstyle='round,pad=0.12', fc='white', ec='none', alpha=0.85))
    if scissors:
        (a, b), (c, d) = CUTBAR
        ax.plot([a, c], [b, d], color='black', lw=1.8, zorder=6, solid_capstyle='round')
        ax.text(CUT[0], CUT[1], '✂', fontsize=12, family='DejaVu Sans', ha='center', va='center',
                rotation=-60, color='black', zorder=6)


# ---------------------------------------------------------------- functionalised chain
def chain_geometry():
    m = K.depict(K.build_chain(4), straighten_units=4)
    conf = m.GetConformer()
    X = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y] for i in range(m.GetNumAtoms())])
    nm = {a.GetIdx(): (a.GetIntProp('unit'), a.GetProp('nm')) for a in m.GetAtoms()}
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in m.GetBonds()]
    el = {a.GetIdx(): a.GetSymbol() for a in m.GetAtoms()}
    return m, X, nm, bonds, el


def draw_chain(ax, tag_fs=8, title=None):
    m, X, nm, bonds, el = chain_geometry()
    idx = {v: k for k, v in nm.items()}
    for k in range(4):
        ring = [idx[(k, n)] for n in ('C1', 'C2', 'C3', 'C4', 'C5', 'O5')]
        ax.add_patch(Polygon(X[ring], closed=True, fc=SUGAR_FILL, ec='none', zorder=1))
    for i, j in bonds:
        ax.plot(X[[i, j], 0], X[[i, j], 1], color='black', lw=1.0, solid_capstyle='round', zorder=2)
    for i, e in el.items():
        if e == 'O':
            ax.add_patch(Circle(X[i], 0.17, fc=O_COL, ec='none', zorder=3))
    blen = np.median([np.linalg.norm(X[i] - X[j]) for i, j in bonds])
    tp = []
    for k, o, tag in TAGS:
        io = idx[(k, o)]
        nb = [j for a, b in bonds for j in ((b,) if a == io else (a,) if b == io else ())]
        v = X[io] - X[nb[0]]
        v /= np.linalg.norm(v)
        p = X[io] + v * blen * 1.25
        tp.append(p)
        ax.plot([X[io, 0], p[0]], [X[io, 1], p[1]], color=SUBST[tag], lw=1.6, zorder=2.5)
        ax.text(p[0], p[1], tag, fontsize=tag_fs, weight='bold', color='white', ha='center', va='center', zorder=4,
                bbox=dict(boxstyle='round,pad=0.22', fc=SUBST[tag], ec='none'))
    A = np.vstack([X, np.array(tp)])
    lo, hi = A.min(0), A.max(0)
    ax.set_xlim(lo[0] - 1.1, hi[0] + 0.6); ax.set_ylim(lo[1] - 0.8, hi[1] + 0.8)
    ax.set_aspect('equal'); ax.axis('off')
    if title:
        ax.set_title(title, fontsize=8.5, weight='bold', pad=1)


# ---------------------------------------------------------------- rule
def ball(ax, x, y, col, ms=10):
    ax.plot([x], [y], 'o', ms=ms, mfc=col, mec='black', mew=0.6, transform=ax.transAxes, clip_on=False, zorder=3)


def rule_rows(ax, y0, dy, fs=8, with_marks=True):
    rows = [(FIXED, '−1 / +1: fixed', 'substituent blocked', '✗', '#C0392B'),
            (VAR, '−2 / +2: variable', 'enzyme-dependent', '?', TXT_VAR),
            (FREE, 'outer sites, O6: free', 'substituent fits', '✓', '#1E8449')]
    for k, (col, l1, l2, mk, mc) in enumerate(rows):
        y = y0 - k * dy
        ball(ax, 0.07, y, col)
        ax.text(0.17, y + 0.035, l1, fontsize=fs, weight='bold', va='center', transform=ax.transAxes)
        ax.text(0.17, y - 0.06, l2, fontsize=fs, va='center', transform=ax.transAxes)
        if with_marks:
            ax.text(0.97, y, mk, fontsize=fs + 4, family='DejaVu Sans', weight='bold', color=mc, ha='right',
                    va='center', transform=ax.transAxes)


# ---------------------------------------------------------------- layouts
def arrow(fig, x0, x1, y):
    fig.add_artist(FancyArrowPatch((x0, y), (x1, y), transform=fig.transFigure, arrowstyle='-|>',
                                   mutation_scale=10, lw=1.2, color='#555555'))


def layout_A():
    """chain | cleft | rule rows with tick/cross marks"""
    fig = plt.figure(figsize=(WID, HEI))
    axc = fig.add_axes((0.005, 0.12, 0.27, 0.76))
    draw_chain(axc, title='functionalised starch')
    axm = fig.add_axes((0.30, 0.02, 0.36, 0.96))
    draw_cleft(axm)
    arrow(fig, 0.265, 0.305, 0.5)
    axr = fig.add_axes((0.675, 0.04, 0.32, 0.92)); axr.axis('off')
    axr.text(0.0, 0.95, 'Enzyme subsites', fontsize=8.5, weight='bold', transform=axr.transAxes, va='center')
    rule_rows(axr, 0.74, 0.27)
    axr.text(0.0, 0.06, 'Tolerance: O6 > O2, O3', fontsize=8, weight='bold', transform=axr.transAxes,
             va='center')
    return fig


def layout_B():
    """chain top-left over a small arrow; larger cleft; rule rows as the colour key (no marks)"""
    fig = plt.figure(figsize=(WID, HEI))
    axc = fig.add_axes((0.0, 0.30, 0.25, 0.62))
    draw_chain(axc, title='functionalised starch')
    axm = fig.add_axes((0.255, 0.0, 0.42, 1.0))
    draw_cleft(axm)
    arrow(fig, 0.215, 0.275, 0.62)
    fig.text(0.125, 0.17, 'Where can a substituted\nglucose sit?', fontsize=8, ha='center', va='center',
             style='italic')
    axr = fig.add_axes((0.69, 0.04, 0.31, 0.92)); axr.axis('off')
    rule_rows(axr, 0.86, 0.29, with_marks=True)
    axr.add_patch(FancyBboxPatch((0.0, 0.0), 0.98, 0.14, boxstyle='round,pad=0.0,rounding_size=0.04',
                                 fc='#F2F2F2', ec='none', transform=axr.transAxes))
    axr.text(0.49, 0.07, 'Best site: O6 > O2, O3', fontsize=8, weight='bold', ha='center', va='center',
             transform=axr.transAxes)
    return fig


def layout_C():
    """three framed boxes: problem -> structure -> rule"""
    fig = plt.figure(figsize=(WID, HEI))
    boxes = [(0.005, 0.30), (0.315, 0.36), (0.685, 0.31)]
    for x, w in boxes:
        fig.add_artist(FancyBboxPatch((x, 0.02), w - 0.01, 0.96, boxstyle='round,pad=0,rounding_size=0.015',
                                      transform=fig.transFigure, fc='none', ec='#BBBBBB', lw=0.6))
    fig.text(0.155, 0.91, 'Functionalised starch', fontsize=8.5, weight='bold', ha='center', va='center')
    axc = fig.add_axes((0.015, 0.10, 0.28, 0.74))
    draw_chain(axc)
    fig.text(0.49, 0.91, '637 PDB structures', fontsize=8.5, weight='bold', ha='center',
             va='center')
    axm = fig.add_axes((0.32, 0.03, 0.34, 0.82))
    draw_cleft(axm, fs_sub=8.5)
    for x in (0.297, 0.667):
        arrow(fig, x, x + 0.03, 0.5)
    axr = fig.add_axes((0.695, 0.05, 0.29, 0.90)); axr.axis('off')
    axr.text(0.5, 0.955, 'Design rule', fontsize=8.5, weight='bold', ha='center', va='center',
             transform=axr.transAxes)
    rule_rows(axr, 0.76, 0.25)
    axr.text(0.5, 0.06, 'Best site: O6 > O2, O3', fontsize=8.5, weight='bold', ha='center', va='center',
             transform=axr.transAxes)
    return fig


def save(fig, stem, folder, tiff=False):
    fig.savefig(folder / f'{stem}.png', dpi=DPI)
    fig.savefig(folder / f'{stem}.pdf', dpi=DPI)
    if tiff:
        im = Image.open(folder / f'{stem}.png').convert('RGB')
        im.save(folder / f'{stem}.tiff', compression='tiff_lzw', dpi=(DPI, DPI))
    plt.close(fig)


def main():
    lay = {'A': layout_A, 'B': layout_B, 'C': layout_C}
    if '--pick' in sys.argv:
        k = sys.argv[sys.argv.index('--pick') + 1]
        save(lay[k](), 'graphical_abstract', SUB, tiff=True)
        print('picked', k)
        return
    for k, f in lay.items():
        save(f(), f'ga_{k}', OUT)
    print('candidates written to', OUT)


if __name__ == '__main__':
    main()
