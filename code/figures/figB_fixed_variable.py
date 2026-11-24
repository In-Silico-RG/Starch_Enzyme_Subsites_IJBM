"""Descriptive candidate B, "Fixed versus variable": the residue at a fixed contact position and at a variable
one next to it in sequence, from one representative of every 50 % sequence group, superposed in the family frame.

Data
  data/analysis/44_rule_test/clusters_bound.tsv  50 % groups (cl50) of the ligand-bound GH13 / GH31 entries,
      with the frame chain of each entry.
  data/analysis/44_rule_test/posmap.tsv          residue of each entry chain at each reference position
      (C-alpha within 2.0 A after superposition), reference numbering 1AMY (GH13) / 2QMJ (GH31).
  data/analysis/44_rule_test/mutations.tsv       engineered mutations (an entry mutated at the position is
      not used as representative).
  data/calc/03_frames/superposition.tsv          R, t taking each chain into the family frame.
  data/census/cif/*.cif                          deposited coordinates (model 1, first altloc).
Representative of a group: among its entries with a residue at the position and no engineered mutation there,
  the one with the best resolution (database table entry.resolution_A). A group with no residue at the position
  (gap) is counted, not drawn.
Pairs (all sequence neighbours): GH13 Asp179 (fixed, catalytic nucleophile) / Phe180 (variable, +1);
  GH13 Glu204 (fixed, acid/base) / Trp206 (variable, +1/+2); GH31 Asp327 (fixed) / Ile328 (variable).
Drawing: ChimeraX, frame of 1AMY / 2QMJ. Faded reference cartoon; one thin stick (r = 0.10 A, side chain from
  C-alpha) per group, carbons by residue type: aromatic F/Y/W orange #E69F00, aliphatic A/V/L/I/M/P/C/G pink
  #CC79A7, polar S/T/N/Q/D/E/K/R/H blue #0072B2 (the colours of the contact classes; here assigned by residue type);
  the sugar of one entry (GH13 1QHO, GH31 3WEO, both in the frame) as gold sticks with subsite numbers.
Writes figures/candidates/descriptive/figB_fixed_variable.{png,pdf}, 7 in wide, 300 dpi.
Run: ~/miniconda3/bin/python code/figures/figB_fixed_variable.py [--rerender] [--draft]
"""
import hashlib
import json
import sqlite3
import sys
from collections import Counter, defaultdict
from pathlib import Path

import gemmi
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fig2_v3_common as C            # noqa: E402
import figdesc_common as F            # noqa: E402

DRAFT = '--draft' in sys.argv
PX, MS, SS = 1600, (64 if DRAFT else 256), (2 if DRAFT else 3)
PANEL_W, ASPECT, TITLE_H = 7.0 / 3, 1.0, 0.22
MARGIN = 2.4
CART_ALPHA = 0.60
TYPE_CLS = {**{k: 'Ar' for k in ('PHE', 'TYR', 'TRP')},
            **{k: 'Al' for k in ('ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PRO', 'CYS', 'GLY')},
            **{k: 'H' for k in ('SER', 'THR', 'ASN', 'GLN', 'ASP', 'GLU', 'LYS', 'ARG', 'HIS')}}
# aromatic carbons in the darker orange of the fig2_v3 Ar labels: #E69F00 is not told apart from the gold sugar
# where a Trp stacks under a ring
CLS_COL = {'H': C.COL['H'], 'Ar': '#B87800', 'Al': C.COL['Al']}
# weight of the aromatic face-on term per variable position: on for Trp206 (rings stacked under +2 read only
# face-on); off for Phe180, where it turned the view until Asp179 hid behind the -1/-2 units
AROM_W = {'F180': 0.0, 'W206': 0.4, 'I328': 0.0}
PAIRS = [('GH13', 'D179', 'F180', ['-1', '+1', '+2']),
         ('GH13', 'E204', 'W206', ['-1', '+1', '+2']),
         ('GH31', 'D327', 'I328', ['-1', '+1'])]
SUGAR = {'GH13': ('1QHO', 'A', ['-2', '-1', '+1', '+2']), 'GH31': ('3WEO', 'A', ['-1', '+1', '+2'])}
ONE = {v.upper(): k for k, v in C.ONE3.items()}


def resolutions():
    con = sqlite3.connect(F.P / 'data/db/enzimas_experimental.sqlite')
    out = {}
    for pdb, r in con.execute('select pdb, resolution_A from entry'):
        try:
            out[pdb] = float(r)
        except (TypeError, ValueError):
            out[pdb] = 9.9
    return out


def representatives(fam, pos):
    res = resolutions()
    pm = {(r['pdb'], r['chain'], r['position']): r for r in F.read_tsv(F.D44 / 'posmap.tsv')}
    mut = {(r['pdb'], r['chain'], r['seq']) for r in F.read_tsv(F.D44 / 'mutations.tsv')}
    groups = defaultdict(list)
    for r in F.read_tsv(F.D44 / 'clusters_bound.tsv'):
        if r['family'] == fam:
            groups[r['cl50']].append(r)
    reps, gaps = [], []
    for cl, ents in sorted(groups.items()):
        cand = []
        for e in ents:
            p = pm.get((e['pdb'], e['frame_chain'], pos))
            if not p:
                continue
            if (e['pdb'], p['chain'], p['seq']) in mut:
                continue
            if (fam, e['pdb'], p['chain']) not in F.mats():
                continue
            cand.append((res.get(e['pdb'], 9.9), e['pdb'], p))
        if not cand:
            gaps.append(cl)
            continue
        cand.sort(key=lambda t: (t[0], t[1]))
        reps.append(dict(cl=cl, pdb=cand[0][1], chain=cand[0][2]['chain'], seq=int(cand[0][2]['seq']),
                         resname=cand[0][2]['resname'], res=cand[0][0]))
    return reps, gaps


def side_chain(fam, rep):
    st = F.structure(rep['pdb'])
    X, _ = F.residue_atoms(st[0], rep['chain'], rep['seq'], rep['resname'])
    keep = [n for n in X if n not in C.BACKBONE]
    return {n: F.to_frame(fam, rep['pdb'], rep['chain'], X[n])[0] for n in keep}


def sugar(fam):
    pdb, ch, subs = SUGAR[fam]
    st = F.structure(pdb)
    rows = [r for r in F.units_labelled() if r['pdb'] == pdb and r['subsite'] in subs]
    units, X = [], {}
    for r in rows:
        lab = F.unit_labels(r)
        A, _ = F.residue_atoms(st[0], r['chain'], r['seq'], r['comp'])
        for n, x in A.items():
            X[(r['chain'], int(r['seq']), n)] = F.to_frame(fam, pdb, ch, x)[0]
        ring = [A[lab[k]] for k in ('C1', 'C2', 'C3', 'C4', 'C5', 'O5')]
        RX = F.to_frame(fam, pdb, ch, np.array(ring))
        units.append(dict(sub=r['subsite'], chain=r['chain'], seq=int(r['seq']), cen=RX.mean(0), ringX=RX))
    return units, X


def best_dir(fam, sets, units, SX, focus, right, max_angle=70.0, w_arom=0.4):
    """Viewing direction that separates the two bundles on the page: among Fibonacci directions within
    max_angle of the family (fig2_v3) direction with the -1 ring normal at |n.z| >= 0.4, minimise
    (fraction of atoms of one bundle within 1.2 A of the other on the page, both ways) + 0.5 x (fraction of
    bundle atoms within 1.0 A of a sugar atom on the page) + 0.3 x angle / max_angle + w_arom x (1 - mean |n.z| of the aromatic rings; AROM_W), so rings
    are seen face-on rather than edge-on."""
    from scipy.spatial import cKDTree
    start = F.view_dir(fam)
    A = np.array([x for r in sets[0]['reps'] for x in r['X'].values()])
    Bb = np.array([x for r in sets[1]['reps'] for x in r['X'].values()])
    keys = {(u['chain'], u['seq']) for u in units if u['sub'] in focus}
    S = np.array([x for k, x in SX.items() if k[:2] in keys])
    m1 = next(u for u in units if u['sub'] == '-1')
    ringn = C.plane(m1['ringX'])[1]
    arom = []
    for st_ in sets:
        for r in st_['reps']:
            if r['resname'] in ('PHE', 'TYR', 'TRP') and all(a in r['X'] for a in C.AROM[r['resname']]):
                arom.append(C.plane(np.array([r['X'][a] for a in C.AROM[r['resname']]]))[1])
    arom = np.array(arom) if arom else np.zeros((0, 3))
    D = C.fib_dirs(8000)
    D = D[D @ start > np.cos(np.radians(max_angle))]
    best = None
    for ez in D:
        if abs(ringn @ ez) < 0.4:
            continue
        Bm = C.basis(ez, right)
        a2, b2, s2 = (A @ Bm.T)[:, :2], (Bb @ Bm.T)[:, :2], (S @ Bm.T)[:, :2]
        ta, tb, ts = cKDTree(a2), cKDTree(b2), cKDTree(s2)
        fa = np.mean([len(x) > 0 for x in tb.query_ball_point(a2, 1.2)])
        fb = np.mean([len(x) > 0 for x in ta.query_ball_point(b2, 1.2)])
        fs = np.mean([len(x) > 0 for x in ts.query_ball_point(np.vstack([a2, b2]), 1.0)])
        ang = np.degrees(np.arccos(np.clip(ez @ start, -1, 1)))
        tar = float(np.mean(np.abs(arom @ ez))) if len(arom) else 1.0
        sc = fa + fb + 0.5 * fs + 0.3 * ang / max_angle + w_arom * (1 - tar)
        if best is None or sc < best[0]:
            best = (sc, ez, fa, fb, fs, ang)
    print(f'  view: score {best[0]:.3f} overlap {best[2]:.2f}/{best[3]:.2f} sugar {best[4]:.2f} '
          f'angle {best[5]:.0f} deg dir {np.round(best[1], 4).tolist()}')
    return best[1]


def write_overlay(path, groups):
    """groups: list of (chain id, list of (resname, {atom: xyz}))."""
    st = gemmi.Structure()
    md = gemmi.Model('1')
    for chid, items in groups:
        ch = gemmi.Chain(chid)
        for k, (rn, X) in enumerate(items):
            r = gemmi.Residue()
            r.name = rn
            r.seqid = gemmi.SeqId(k + 1, ' ')
            for n, x in X.items():
                a = gemmi.Atom()
                a.name = n
                a.element = gemmi.Element(n[0])
                a.pos = gemmi.Position(*x)
                a.occ = 1.0
                r.add_atom(a)
            ch.add_residue(r)
        md.add_chain(ch)
    st.add_model(md)
    st.setup_entities()
    st.write_pdb(str(path))


def prepare(fam, posA, posB, focus):
    out = dict(fam=fam, pos=[posA, posB])
    sets = []
    for pos in (posA, posB):
        reps, gaps = representatives(fam, pos)
        for r in reps:
            r['X'] = side_chain(fam, r)
            r['cls'] = TYPE_CLS.get(r['resname'], 'H')
        sets.append(dict(pos=pos, reps=reps, gaps=gaps))
        cnt = Counter(ONE[r['resname']] for r in reps)
        print(f'{fam} {pos}: {len(reps)} groups drawn, {len(gaps)} without residue; {dict(cnt.most_common())}')
        print('   ', [(r['cl'], r['pdb'], r['chain'], r['seq'], r['resname'], r['res']) for r in reps])
    units, SX = sugar(fam)
    out.update(sets=sets, units=units, sugarX=SX)
    cen = {u['sub']: u['cen'] for u in units}
    right = cen['+2'] - cen['-1'] if '+2' in cen else cen['+1'] - cen['-1']
    B = C.basis(best_dir(fam, sets, units, SX, focus, right, w_arom=AROM_W.get(posB, 0.0)), right)
    pts = np.array([x for s in sets for r in s['reps'] for x in r['X'].values()] +
                   [x for k, x in SX.items() if any(u['sub'] in focus and (u['chain'], u['seq']) == k[:2]
                                                    for u in units)])
    c, W, zmin, zmax = F.field(B, pts, MARGIN, ASPECT)
    out.update(B=B, c=c, W=W, zmin=zmin, zmax=zmax, focus=focus)
    return out


def render(info, png_cart, png_stk):
    fam = info['fam']
    rp, rc = F.REF[fam]
    spdb, sch, _ = SUGAR[fam]
    ov = F.WORK / f'B_overlay_{fam}_{info["pos"][0]}_{info["pos"][1]}.pdb'
    write_overlay(ov, [(chid, [(r['resname'], r['X']) for r in s['reps']])
                       for chid, s in zip(('P', 'Q'), info['sets'])])
    cam = F.camera(info['B'], info['c'], info['W'], PX)
    clip = (f'clip near {-(info["zmax"] + 1.5):.2f} far 10.0 position '
            f'{",".join(f"{v:.3f}" for v in info["c"])}')
    L1 = [f'open {F.CIF / (rp + ".cif")}', 'delete solvent', 'delete H', 'set bgColor white', 'graphics quality 6',
          'hide atoms', 'hide cartoons', f'cartoon /{rc}', f'cartoon style /{rc} width 1.6 thickness 0.35',
          f'color /{rc} {C.CARTOON} target c', f'lighting soft multiShadow {MS}', 'lighting depthCue true',
          'lighting depthCueStart 0.2', 'lighting depthCueEnd 1.0', 'material dull',
          'graphics silhouettes true width 1.5 color #8C8C8C depthJump 0.03', *cam, clip,
          f'save {png_cart} width {PX} height {PX} supersample {SS} transparentBackground false']
    F.run_chimerax(L1, Path(png_cart).stem)
    sug = '/' + ','.join(sorted({u['chain'] for u in info['units']})) + ':' + \
          ','.join(sorted({str(u['seq']) for u in info['units']}))
    L2 = [f'open {F.CIF / (spdb + ".cif")}', f'view matrix models #1,{F.cx_matrix(fam, spdb, sch)}',
          'delete solvent', 'delete H', 'hide #1 atoms', 'hide #1 cartoons', f'open {ov}', 'hide #2 cartoons',
          'set bgColor white', 'graphics quality 6',
          f'show #1{sug}', f'style #1{sug} stick', f'size #1{sug} stickRadius 0.20',
          f'color #1{sug} & C {C.SUGAR_C} target a', f'color #1{sug} & O {C.O_COL} target a',
          f'color #1{sug} & N {C.N_COL} target a',
          'show #2 atoms', 'style #2 stick', 'size #2 stickRadius 0.10']
    for chid, s in zip(('P', 'Q'), info['sets']):
        for k, r in enumerate(s['reps']):
            L2.append(f'color #2/{chid}:{k + 1} & C {CLS_COL[r["cls"]]} target a')
    L2 += [f'color #2 & O {C.O_COL} target a', f'color #2 & N {C.N_COL} target a', f'color #2 & S #E6C229 target a',
           'hide pbonds', f'lighting soft multiShadow {MS}', 'lighting depthCue false', 'material dull',
           'graphics silhouettes true width 2.2 color black depthJump 0.01', *cam, 'clip off',
           f'save {png_stk} width {PX} height {PX} supersample {SS} transparentBackground true']
    F.run_chimerax(L2, Path(png_stk).stem)


def renders(info, rer):
    key = json.dumps(dict(f=info['fam'], p=info['pos'], W=round(info['W'], 3), B=np.round(info['B'], 4).tolist(),
                          ms=MS, ss=SS, px=PX, v=4), sort_keys=True)
    h = hashlib.md5(key.encode()).hexdigest()[:8]
    zc = F.WORK / f'B_cart_{info["fam"]}_{info["pos"][0]}_{h}.png'
    zs = F.WORK / f'B_stk_{info["fam"]}_{info["pos"][0]}_{h}.png'
    if rer or not zs.exists():
        render(info, zc, zs)
    return zc, zs


def counts_text(s):
    cnt = Counter(ONE[r['resname']] for r in s['reps'])
    t = ' '.join(f'{k}{v}' for k, v in sorted(cnt.items(), key=lambda kv: (-kv[1], kv[0])))
    if s['gaps']:
        t += f' –{len(s["gaps"])}'
    return t


def draw_panel(info, zc, zs, ax, ptw):
    B, c, W = info['B'], info['c'], info['W']
    cart = 1 - CART_ALPHA * (1 - plt.imread(zc)[..., :3])
    ax.imshow(cart, extent=(0, PX, PX, 0), interpolation='lanczos', zorder=0)
    ax.imshow(plt.imread(zs), extent=(0, PX, PX, 0), interpolation='lanczos', zorder=1)
    crop = (0, 0, PX, PX)
    ax.set_xlim(0, PX); ax.set_ylim(PX, 0); ax.axis('off')
    P_, W_, R_ = [], [], []
    for s in info['sets']:
        for r in s['reps']:
            for x in F.segment_samples(r['X']):
                P_.append(x); W_.append(3.0); R_.append(0.35)
    for x in F.segment_samples(info['sugarX']):
        P_.append(x); W_.append(6.0); R_.append(0.5)
    Q = C.to_px(np.array(P_), B, c, W, PX)
    OP, OW, OR = Q[:, :2], np.array(W_), np.array(R_) * PX / W
    pt = PX / ptw / 72.0
    items = []
    for s in info['sets']:
        fam, pos = info['fam'], s['pos']
        name = f'{C.ONE3[pos[0]]}{pos[1:]}'
        txt = f'{name}  {counts_text(s)}'
        tips = np.array([np.mean([x for n, x in r['X'].items() if n not in ('CA',)], 0) for r in s['reps']])
        anc = C.to_px(np.array([np.median(tips, 0)]), B, c, W, PX)[:, :2]
        # several anchors: every tip, so the label may sit beside any member of the bundle
        ancs = np.vstack([anc, C.to_px(tips, B, c, W, PX)[:, :2]])
        w, h = F.text_size_px(txt, 7.0, pt)
        items.append(dict(kind='res', txt=txt, ancs=ancs, col='#222222', fs=7.0, weight='normal', w=w, h=h))
    for u in info['units']:
        cen = C.to_px(u['cen'], B, c, W, PX)[0, :2]
        if not (0 < cen[0] < PX and 0 < cen[1] < PX):
            continue
        t = u['sub'].replace('-', '−')
        w, h = F.text_size_px(t, 8.0, pt, 'bold')
        items.append(dict(kind='sub', txt=t, ancs=cen[None], col='black', fs=8.0, weight='bold', w=w, h=h))
    labs = F.place_labels(items, OP, OW, OR, crop, PX / W)
    F.draw_labels(ax, labs)
    ax.add_patch(plt.Rectangle((2, 2), PX - 4, PX - 4, fill=False, lw=0.6, ec='0.6', zorder=7))


def main():
    rer = '--rerender' in sys.argv
    W = 7.0
    HP = PANEL_W * ASPECT
    H = 0.30 + TITLE_H + HP + 0.40
    fig = plt.figure(figsize=(W, H))
    # no in-figure title (IJBM editorial review 2026-11-24): the title is in the caption
    for k, (fam, a, b, focus) in enumerate(PAIRS):
        info = prepare(fam, a, b, focus)
        zc, zs = renders(info, rer)
        x0 = k * PANEL_W / W
        ax = fig.add_axes((x0 + 0.004, 0.40 / H, PANEL_W / W - 0.008, HP / H))
        draw_panel(info, zc, zs, ax, PANEL_W)
        fig.text(x0 + 0.006, (0.40 + HP + 0.05) / H, 'abc'[k], fontsize=9, weight='bold', va='bottom')
        rp = F.REF[fam][0]
        fig.text(x0 + 0.04, (0.40 + HP + 0.05) / H,
                 f'{fam}: {C.ONE3[a[0]]}{a[1:]} and {C.ONE3[b[0]]}{b[1:]}', fontsize=7.5, va='bottom')
    h = [Line2D([], [], color=C.COL['H'], lw=2, label='polar (S T N Q D E K R H)'),
         Line2D([], [], color=CLS_COL['Ar'], lw=2, label='aromatic (F Y W)'),
         Line2D([], [], color=C.COL['Al'], lw=2, label='aliphatic (A V L I M P)'),
         Line2D([], [], color=C.SUGAR_C, lw=3, label='sugar (1QHO / 3WEO)'),
         ]
    fig.legend(handles=h, loc='center', bbox_to_anchor=(0.5, 0.26 / H), ncol=4, frameon=False, fontsize=6.5,
               handlelength=1.8, columnspacing=1.2, handletextpad=0.5)
    fig.text(0.5, 0.05 / H, 'A7 = Ala in 7 groups (one representative per 50 % sequence group); '
             '–n = n groups without a residue at the position', ha='center', va='bottom', fontsize=6.5, color='0.25')
    fig.savefig(F.OUT / 'figB_fixed_variable.png', dpi=300)
    fig.savefig(F.OUT / 'figB_fixed_variable.pdf', dpi=300)
    plt.close(fig)


if __name__ == '__main__':
    main()
