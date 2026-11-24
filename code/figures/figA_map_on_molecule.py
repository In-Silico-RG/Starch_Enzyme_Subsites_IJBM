"""Descriptive candidate A, "The map on the molecule": the hydroxyl map of 45 drawn on one crystal structure.

Data
  data/analysis/45_hydroxyl_map/cell_summary.tsv  state of each family x subsite x oxygen, set
      all_glucoselike: INV = held (polar partner <= 3.5 A in >= half the 50 % groups) by a lineage-invariant side
      chain; VAR = held by a side chain at a variable position; MC = held by main chain only; EXP = not held and
      open to solvent (median rel SASA >= 0.25 and median SASA >= 5 A^2); BUR = not held and not open.
  data/analysis/45_hydroxyl_map/units_labelled.tsv  ring units with graph-based labels (O2, O3, O6 ->
      atom names).
  data/analysis/45_hydroxyl_map/partners_4A.tsv     polar partners (protein N/O <= 3.5 A) of each oxygen
      in the entry drawn; these residues are the "holding residues" shown.
  data/analysis/45_hydroxyl_map/position_identity.tsv  fixed residue = 45's criterion (>= 3 groups,
      m50 = 1.0); others, and residues without a reference position, are variable.
  data/census/cif/{1QHO,3WEO}.cif, model 1, first altloc.
Entries
  GH13: 1QHO (Bacillus stearothermophilus maltogenic alpha-amylase Novamyl, 1.70 A, wild type, no engineered
        mutation in 44_rule_test/mutations.tsv): density-supported units at -3..+3 with a real glucose (GLC) at
        -1, -2 and +2; +1 is the 4-amino-4,6-dideoxy-glucose of acarbose (AGL), which has no O6. The only
        GH13 entry with glucose units at -2..+2 that is native is this one or entries without +1/+2
        (5TD4 has glucose at -3..+2 but carries D300N at the reference Asp289 that holds O2/O3 at -1).
  GH31: 3WEO (sugar beet alpha-glucosidase, 1.45 A, acarviosyl-maltotriose): -1 valienamine ring, +1 AGL-like
        ring (no O6), +2 glucose; the entry of fig2_v3 b. GH31 has no -2 subsite.
Drawing: as fig2_v3 (fig2_v3_struct.render_zoom): faded cartoon layer and a stick layer, both in the 1AMY / 2QMJ
frame (03_frames superposition), viewing direction of fig2_v3 carried into the frame. Sugar units gold sticks
r = 0.22 A; each O2, O3, O6 a sphere r = 0.50 A coloured by its 45 state; holding residues thin sticks r = 0.14 A
(side chain from CA, or N-CA-C-O when the partner is main chain), carbons vermillion (fixed position) or sky blue
(variable position); thin grey dashes from each oxygen to its partner atoms.
Writes figures/candidates/descriptive/figA_map_on_molecule.{png,pdf}, 7 in wide, 300 dpi.
Run: ~/miniconda3/bin/python code/figures/figA_map_on_molecule.py [--rerender]
"""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fig2_v3_common as C            # noqa: E402
import figdesc_common as F            # noqa: E402

PX, TPX, MS, SS = 2000, 900, (64 if "--draft" in sys.argv else 256), (2 if "--draft" in sys.argv else 3)
PANEL_W, ASPECT, TITLE_H, LEG_H = 3.42, 0.92, 0.22, 0.62
MARGIN = 2.6
CART_ALPHA = 0.70
SPH_R = 0.50
STATE_COL = {'INV': '#D55E00', 'VAR': '#56B4E9', 'MC': '#009E73', 'EXP': '#FFFFFF', 'BUR': '#9A9A9A'}
STATE_TXT = {'INV': 'held by a fixed residue', 'VAR': 'held by a variable residue', 'MC': 'held by main chain',
             'EXP': 'free, open to solvent', 'BUR': 'free, buried'}
RES_COL = {'fixed': '#D55E00', 'variable': '#56B4E9', 'mc': '#8C8C8C'}
RES_TXT = {'#D55E00': '#B84A00', '#56B4E9': '#2A86BF', '#8C8C8C': '#555555'}
SUBS = ['-2', '-1', '+1', '+2']
ENTRY = {'GH13': dict(pdb='1QHO', chain='A', res=1.70, short='maltogenic α-amylase'),
         'GH31': dict(pdb='3WEO', chain='A', res=1.45, short='α-glucosidase')}
THUMB = 0.25
# GH13: direction refined in the frame by figA_views.py (candidate 0, 11 deg from the fig2_v3 direction; the -2
# ring is edge-on from every direction within 40 deg of it; the direction that shows -2 flat, 86 deg away, looks
# through the protein and puts +1/+2 edge-on, tried and dropped); GH31: the fig2_v3 direction unchanged.
VIEW_FRAME = {'GH13': [0.0215, 0.4419, 0.8968]}


def prepare(fam):
    e = ENTRY[fam]
    pdb, pch = e['pdb'], e['chain']
    st = F.structure(pdb)
    m = st[0]
    states = F.cell_states()
    ident = F.position_identity()
    pm = C.posmap()
    rows = [r for r in F.units_labelled() if r['pdb'] == pdb and r['subsite'] in SUBS]
    tf = lambda X: F.to_frame(fam, pdb, pch, X)
    units, sugarX, oxy = [], {}, []
    for r in rows:
        lab = F.unit_labels(r)
        X, res = F.residue_atoms(m, r['chain'], r['seq'], r['comp'])
        for n, x in X.items():
            sugarX[(r['chain'], int(r['seq']), n)] = tf(x)[0]
        ring = [lab[k] for k in ('C1', 'C2', 'C3', 'C4', 'C5', 'O5')]
        cen = tf(np.array([X[a] for a in ring])).mean(0)
        units.append(dict(sub=r['subsite'], chain=r['chain'], seq=int(r['seq']), comp=r['comp'], ring=ring,
                          cen=cen, inst=r['instance_id'], ringno=r['ring'],
                          ringX=tf(np.array([X[a] for a in ring]))))
        for o in ('O2', 'O3', 'O6'):
            if o not in lab:
                continue
            st_ = states[(fam, r['subsite'], o)]['state']
            oxy.append(dict(sub=r['subsite'], lab=o, atom=lab[o], chain=r['chain'], seq=int(r['seq']),
                            X=tf(X[lab[o]])[0], state=st_, inst=r['instance_id'], ringno=r['ring']))
    # holding residues: polar partners of the drawn O2/O3/O6 in this entry
    keys = {(o['inst'], o['ringno'], o['lab']): o for o in oxy}
    parts = [p for p in F.read_tsv(F.D45 / 'partners_4A.tsv') if p['pdb'] == pdb and p['polar_hb'] == 'True'
             and (p['instance_id'], p['ring'], p['label']) in keys]
    resd = {}
    hb = {}
    for p in parts:
        k = (p['p_chain'], int(p['p_seq']))
        d = resd.setdefault(k, dict(chain=p['p_chain'], seq=int(p['p_seq']), resname=p['p_resname'],
                                    pos=p['position'], parts=set()))
        d['parts'].add(p['part'])
        o = keys[(p['instance_id'], p['ring'], p['label'])]
        hk = (o['inst'], o['ringno'], o['lab'], k)
        if hk not in hb or float(p['dist']) < hb[hk][2]:
            hb[hk] = (o, p['p_atom'], float(p['dist']), k)
    residues = []
    for k, d in sorted(resd.items(), key=lambda kv: kv[0][1]):
        X, res = F.residue_atoms(m, d['chain'], d['seq'], d['resname'])
        sc = 'sc' in d['parts']
        keep = [n for n in X if n not in C.BACKBONE] if sc else [n for n in ('N', 'CA', 'C', 'O') if n in X]
        if sc and d['resname'] == 'GLY':
            keep = ['N', 'CA', 'C', 'O']
        fixed = bool(d['pos']) and F.is_fixed(fam, d['pos'], ident)
        kind = ('fixed' if fixed else 'variable') if sc else 'mc'
        if d['pos']:
            lab = C.label_text(d['resname'], d['pos'], None)
        else:
            lab = f'{d["resname"].capitalize()}{d["seq"]}*'
        residues.append(dict(d, X={n: F.to_frame(fam, pdb, pch, X[n])[0] for n in keep}, keep=keep, kind=kind,
                             col=RES_COL[kind], label=lab, sc=sc))
    hbonds = []
    for (inst, ringno, olab, k), (o, patom, dist, kk) in hb.items():
        rd = next(r for r in residues if (r['chain'], r['seq']) == k)
        X, _ = F.residue_atoms(m, rd['chain'], rd['seq'], rd['resname'])
        hbonds.append(dict(o=o, res=rd, a=patom, Y=F.to_frame(fam, pdb, pch, X[patom])[0], d=dist))
    # whole chain CA (thumbnail)
    ca = F.to_frame(fam, pdb, pch, np.array([r['CA'][0].pos.tolist() for r in m[pch] if r.find_atom('CA', '*')]))
    info = dict(fam=fam, pdb=pdb, chain=pch, e=e, units=units, sugarX=sugarX, oxy=oxy, residues=residues,
                hbonds=hbonds, ca=ca)
    cen = {u['sub']: u['cen'] for u in units}
    vd = np.array(VIEW_FRAME[fam]) if fam in VIEW_FRAME else F.view_dir(fam)
    B = C.basis(vd / np.linalg.norm(vd), cen['+2'] - cen['-1'])
    pts = np.array(list(sugarX.values()) + [x for r in residues for x in r['X'].values()])
    c, W, zmin, zmax = F.field(B, pts, MARGIN, ASPECT)
    info.update(B=B, c=c, W=W, zmin=zmin, zmax=zmax)
    print(fam, pdb, 'units', [(u['sub'], u['comp']) for u in units])
    print('  hydroxyls', [(o['sub'], o['lab'], o['state']) for o in oxy])
    print('  residues', [(r['label'], r['kind'], ''.join(sorted(r['parts']))) for r in residues])
    return info


def specs(info):
    sug = '/' + ','.join(sorted({u['chain'] for u in info['units']})) + ':' + \
          ','.join(sorted({str(u['seq']) for u in info['units']}))
    return sug


def render_zoom(info, png_cart, png_stk):
    fam, pdb, pch = info['fam'], info['pdb'], info['chain']
    mx = F.cx_matrix(fam, pdb, pch)
    cam = F.camera(info['B'], info['c'], info['W'], PX)
    head = [f'open {F.CIF / (pdb + ".cif")}', f'view matrix models #1,{mx}', 'delete solvent', 'delete H',
            'set bgColor white', 'graphics quality 6', 'hide atoms', 'hide cartoons']
    clip = (f'clip near {-(info["zmax"] + 1.5):.2f} far 10.0 position '
            f'{",".join(f"{v:.3f}" for v in info["c"])}')
    L1 = head + [f'cartoon /{pch}', f'cartoon style /{pch} width 1.6 thickness 0.35',
                 f'color /{pch} {C.CARTOON} target c', f'lighting soft multiShadow {MS}', 'lighting depthCue true',
                 'lighting depthCueStart 0.2', 'lighting depthCueEnd 1.0', 'material dull',
                 'graphics silhouettes true width 1.5 color #8C8C8C depthJump 0.03', *cam, clip,
                 f'save {png_cart} width {PX} height {PX} supersample {SS} transparentBackground false']
    F.run_chimerax(L1, Path(png_cart).stem)
    sug = specs(info)
    L2 = head + [f'name lig {sug}', 'show lig atoms', 'style lig stick', 'size lig stickRadius 0.22',
                 f'color lig & C {C.SUGAR_C} target a', f'color lig & O {C.O_COL} target a',
                 f'color lig & N {C.N_COL} target a']
    for r in info['residues']:
        spec = f'/{r["chain"]}:{r["seq"]}@' + ','.join(r['keep'])
        L2 += [f'show {spec}', f'style {spec} stick', f'size {spec} stickRadius 0.14',
               f'color {spec} & C {r["col"]} target a', f'color {spec} & O {C.O_COL} target a',
               f'color {spec} & N {C.N_COL} target a', f'color {spec} & S #E6C229 target a']
    for o in info['oxy']:
        spec = f'/{o["chain"]}:{o["seq"]}@{o["atom"]}'
        L2 += [f'style {spec} sphere', f'size {spec} atomRadius {SPH_R}', f'color {spec} {STATE_COL[o["state"]]} target a']
    L2 += ['hide pbonds', f'lighting soft multiShadow {MS}', 'lighting depthCue false', 'material dull',
           'graphics silhouettes true width 3.0 color black depthJump 0.01', *cam, 'clip off',
           f'save {png_stk} width {PX} height {PX} supersample {SS} transparentBackground true']
    F.run_chimerax(L2, Path(png_stk).stem)


def thumb_geometry(info):
    B = info['B']
    v = info['ca'] @ B.T
    lo, hi = v[:, :2].min(0), v[:, :2].max(0)
    mid = np.r_[(lo + hi) / 2, (info['c'] @ B.T)[2]]
    return mid @ B, max(hi - lo) + 9.0


def render_thumb(info, png):
    fam, pdb, pch = info['fam'], info['pdb'], info['chain']
    ct, Wt = thumb_geometry(info)
    sug = specs(info)
    L = [f'open {F.CIF / (pdb + ".cif")}', f'view matrix models #1,{F.cx_matrix(fam, pdb, pch)}', 'delete solvent',
         'delete H', 'set bgColor white', 'hide atoms', 'hide cartoons', f'surface /{pch} enclose /{pch}',
         f'color /{pch} #D5D9DE target s', f'show {sug}', f'style {sug} sphere', f'color {sug} {C.SUGAR_C} target a',
         f'color {sug} & O #E8833A target a', 'lighting soft multiShadow 128', 'lighting depthCue false',
         'material dull', 'graphics silhouettes true width 3 color black depthJump 0.01',
         *F.camera(info['B'], ct, Wt, TPX), 'clip off',
         f'save {png} width {TPX} height {TPX} supersample {SS} transparentBackground true']
    F.run_chimerax(L, Path(png).stem)


def renders(info, rerender):
    key = json.dumps(dict(pdb=info['pdb'], W=round(info['W'], 3), B=np.round(info['B'], 4).tolist(), ms=MS, ss=SS,
                          res=[r['label'] for r in info['residues']], sph=SPH_R, v=4), sort_keys=True)
    h = hashlib.md5(key.encode()).hexdigest()[:8]
    zc, zs = F.WORK / f'A_cart_{info["fam"]}_{h}.png', F.WORK / f'A_stk_{info["fam"]}_{h}.png'
    tp = F.WORK / f'A_thumb_{info["fam"]}_{h}.png'
    if rerender or not zs.exists():
        render_zoom(info, zc, zs)
    if rerender or not tp.exists():
        render_thumb(info, tp)
    return (zc, zs), tp


def obstacles(info):
    P_, W_, R_ = [], [], []
    for r in info['residues']:
        for x in F.segment_samples(r['X']):
            P_.append(x); W_.append(3.0); R_.append(0.4)
    for x in F.segment_samples(info['sugarX']):
        P_.append(x); W_.append(6.0); R_.append(0.55)
    for o in info['oxy']:
        P_.append(o['X']); W_.append(8.0); R_.append(SPH_R + 0.1)
    for h in info['hbonds']:
        for t in np.linspace(0.2, 0.8, 4):
            P_.append(h['o']['X'] * (1 - t) + h['Y'] * t); W_.append(1.5); R_.append(0.25)
    Q = C.to_px(np.array(P_), info['B'], info['c'], info['W'], PX)
    return Q[:, :2], np.array(W_), np.array(R_) * PX / info['W']


def thumb_place(info, crop, tp, OP, OW, OR):
    timg = plt.imread(tp)
    a = timg[..., 3] > 0.02
    rows, cols = np.where(a.any(1))[0], np.where(a.any(0))[0]
    v0, v1, u0, u1 = rows[0], rows[-1] + 1, cols[0], cols[-1] + 1
    sub = timg[v0:v1, u0:u1]
    x0, y0, x1, y1 = crop
    f = THUMB * (x1 - x0) / max(u1 - u0, v1 - v0)
    tw, th = (u1 - u0) * f, (v1 - v0) * f
    m = 0.01 * (x1 - x0)
    pos = {'tl': (x0 + m, y0 + m), 'tr': (x1 - tw - m, y0 + m), 'bl': (x0 + m, y1 - th - m),
           'br': (x1 - tw - m, y1 - th - m)}
    best = None
    for name, (bx, by) in pos.items():
        b = (bx, by, bx + tw, by + th)
        inside = ((OP[:, 0] + OR > b[0]) & (OP[:, 0] - OR < b[2]) & (OP[:, 1] + OR > b[1]) & (OP[:, 1] - OR < b[3]))
        s = OW[inside].sum()
        if best is None or s < best[0]:
            best = (s, name, b)
    s, name, b = best
    ct, Wt = thumb_geometry(info)

    def s2p(p):
        q = C.to_px(p, info['B'], ct, Wt, TPX)[:, :2]
        return np.column_stack([b[0] + (q[:, 0] - u0) * f, b[1] + (q[:, 1] - v0) * f])
    return sub, b, s2p


def draw_panel(info, zp, tp, ax):
    B, c, W = info['B'], info['c'], info['W']
    cart = 1 - CART_ALPHA * (1 - plt.imread(zp[0])[..., :3])
    stk = plt.imread(zp[1])
    hpx = PX * ASPECT
    crop = (0, (PX - hpx) / 2, PX, (PX + hpx) / 2)
    ax.imshow(cart, extent=(0, PX, PX, 0), interpolation='lanczos', zorder=0)
    for h in info['hbonds']:
        q = C.to_px(np.array([h['o']['X'], h['Y']]), B, c, W, PX)
        ax.plot(q[:, 0], q[:, 1], color=F.HB_COL, lw=0.7, ls=(0, (2.2, 1.6)), zorder=0.5, dash_capstyle='round')
    ax.imshow(stk, extent=(0, PX, PX, 0), interpolation='lanczos', zorder=1)
    ax.set_xlim(crop[0], crop[2]); ax.set_ylim(crop[3], crop[1]); ax.axis('off')
    OP, OW, OR = obstacles(info)
    sub, ib, s2p = thumb_place(info, crop, tp, OP, OW, OR)
    ax.imshow(sub, extent=(ib[0], ib[2], ib[3], ib[1]), interpolation='lanczos', zorder=3)
    cp = np.array([[crop[0], crop[1]], [crop[2], crop[1]], [crop[2], crop[3]], [crop[0], crop[3]]])
    sc = np.column_stack([(cp[:, 0] - PX / 2) * W / PX, -(cp[:, 1] - PX / 2) * W / PX, np.zeros(4)]) @ B + c
    ax.add_patch(plt.Polygon(s2p(sc), closed=True, fill=False, lw=0.9, ec='black', zorder=5))
    pt = (crop[2] - crop[0]) / PANEL_W / 72.0
    scale = PX / W
    items = []
    for r in info['residues']:
        names = [n for n in r['X'] if n != 'CA'] or list(r['X'])
        anc = C.to_px(np.array([r['X'][n] for n in names]), B, c, W, PX)[:, :2]
        col = RES_TXT[r['col']]
        w, h = F.text_size_px(r['label'], 7.0, pt)
        items.append(dict(kind='res', txt=r['label'], ancs=anc, col=col, fs=7.0, weight='normal', w=w, h=h))
    for u in info['units']:
        cen = C.to_px(u['cen'], B, c, W, PX)[0, :2]
        t = u['sub'].replace('-', '−')
        w, h = F.text_size_px(t, 8.0, pt, 'bold')
        items.append(dict(kind='sub', txt=t, ancs=cen[None], col='black', fs=8.0, weight='bold', w=w, h=h))
    for o in info['oxy']:
        anc = C.to_px(o['X'], B, c, W, PX)[0, :2]
        t = o['lab']
        w, h = F.text_size_px(t, 5.8, pt)
        items.append(dict(kind='oxy', txt=t, ancs=anc[None], col='#333333', fs=5.8, weight='normal', w=w, h=h,
                          rs=[0.62, 0.8, 1.0, 1.3]))
    labs = F.place_labels(items, OP, OW, OR, crop, scale, forbid=[ib])
    F.draw_labels(ax, labs)
    return crop


def title_of(info):
    e = info['e']
    return f'{info["fam"]} {e["short"]}, PDB {info["pdb"]} ({e["res"]:.2f} Å)'


def legend(sf, H, present, fs=6.5):
    h = [Line2D([], [], marker='o', ls='', mfc=STATE_COL[s], mec='black', mew=0.6, ms=6, label=STATE_TXT[s])
         for s in ('INV', 'VAR', 'MC', 'EXP', 'BUR') if s in present]
    h1 = [Line2D([], [], color=C.SUGAR_C, lw=3.0, label='sugar'),
          Line2D([], [], color=RES_COL['fixed'], lw=1.6, label='fixed residue'),
          Line2D([], [], color=RES_COL['variable'], lw=1.6, label='variable residue'),
          Line2D([], [], color=F.HB_COL, lw=0.75, ls=(0, (2.2, 1.6)), label='N/O–O ≤ 3.5 Å'),
          Line2D([], [], ls='', label='* no reference position (author number)')]
    l1 = sf.legend(handles=h, loc='center', bbox_to_anchor=(0.5, 0.40 / H), ncol=5, frameon=False, fontsize=fs,
                   handletextpad=0.2, columnspacing=0.9, title='O2, O3, O6, state of the cell in the family map:',
                   title_fontsize=fs)
    l1._legend_box.align = 'center'
    sf.add_artist(l1)
    sf.legend(handles=h1, loc='center', bbox_to_anchor=(0.5, 0.10 / H), ncol=5, frameon=False, fontsize=fs,
              handlelength=1.8, columnspacing=1.2, handletextpad=0.5)


def main():
    rer = '--rerender' in sys.argv
    panels = {}
    for fam in ('GH13', 'GH31'):
        info = prepare(fam)
        panels[fam] = (info, *renders(info, rer))
    W = 7.0
    HP = PANEL_W * ASPECT + TITLE_H
    H = 0.30 + HP + LEG_H
    fig = plt.figure(figsize=(W, H))
    # no in-figure title (IJBM editorial review 2026-11-24): the title is in the caption
    for k, (fam, letter) in enumerate((('GH13', 'a'), ('GH31', 'b'))):
        info, zp, tp = panels[fam]
        x0 = 0 if k == 0 else (W - PANEL_W) / W
        ax = fig.add_axes((x0, LEG_H / H, PANEL_W / W, PANEL_W * ASPECT / H))
        draw_panel(info, zp, tp, ax)
        ytit = (LEG_H + PANEL_W * ASPECT + 0.04) / H
        fig.text(x0 + 0.003, ytit, letter, fontsize=9, weight='bold', va='bottom')
        fig.text(x0 + 0.035, ytit, title_of(info), fontsize=7.5, va='bottom')
    legend(fig, H, {o['state'] for fam in panels for o in panels[fam][0]['oxy']})
    fig.savefig(F.OUT / 'figA_map_on_molecule.png', dpi=300)
    fig.savefig(F.OUT / 'figA_map_on_molecule.pdf', dpi=300)
    plt.close(fig)


if __name__ == '__main__':
    main()
