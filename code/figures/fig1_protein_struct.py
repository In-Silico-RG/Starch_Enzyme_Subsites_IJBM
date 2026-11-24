"""Structure panels of Figure 1 candidates: whole protein (cartoon by domain, transparent surface, bound
pseudo-maltooligosaccharide coloured by subsite, catalytic residues) and close-up of the catalytic centre.

Entries (reasons in the captions file):
  GH13  6JQB chain A  maltotetraose-forming amylase, Pseudomonas saccharophila STB07, 1.10 A; acarbose-derived
        heptasaccharide 7SA (one residue, 7 rings) at subsites -4..+3, RSCC 0.971 (whole residue).
  GH31  3WEO chain A  sugar beet alpha-glucosidase, 1.45 A; acarviosyl-maltohexaose (chain D, AC1 + 6 GLC) at
        subsites -1..+7, RSCC 0.805-0.991 per residue.
Domains:
  6JQB  CATH-B (PDBe API /mappings/cath_b/6jqb, read 2026-11-24): 6jqbA01 3.20.20.80 residues 1-360,
        6jqbA02 2.60.40.1180 residues 362-417. Domain B is not a CATH domain (it lies inside 6jqbA01); it is
        taken as the excursion between barrel strand beta3 and helix alpha3 (the usual definition of GH13
        domain B), located on the PyMOL dss secondary structure of 6JQB: beta3 = 108-113 (ends two residues
        before His117, the region I His mapped to 1AMY H92), alpha3 = 168-185 (the helix before beta4 189-192,
        which ends at the nucleophile Asp193) -> B = 114-167; A = 1-113 + 168-361; C = 362-418 (418 = last
        modelled residue).
  3WEO  CATH (database table sifts_domain, SIFTS): 3weoA01 2.60.40.1760 38-313 (N-terminal beta-sandwich),
        3weoA02 3.20.20.80 314-655 ((beta/alpha)8 barrel, domain A), 3weoA03 2.60.40.1180 656-756 (proximal C),
        3weoA04 2.60.40.1180 757-909 (distal C).
Catalytic residues: data/calc/03_frames/catalytic_residues_mapped.tsv (nucleophile, acid/base) and
data/analysis/44_rule_test/posmap.tsv (6JQB Asp294 at reference position D289, the transition-state
stabiliser of GH13).

View of the whole protein: a search over 4000 directions maximises the on-page length of the sugar chain
(subsites -3..+3 shown) and penalises (i) sugar atoms covered on the page by backbone atoms lying in front of
them and (ii) protein atoms in front of the sugar within 12 A on the page (the cleft is seen from its open side).
The chain is then turned to run from the lowest to the highest subsite left to right.
Output: work/<name>.png (ray-traced, orthoscopic) and work/<name>.json (page coordinates of the label anchors).
Run: ~/miniconda3/bin/python code/figures/fig1_protein_struct.py
"""
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import gemmi

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fig1_protein_common as C
WORK = C.WORK

ENTRIES = {
    '6JQB': dict(family='GH13', chain='A', resolution=1.101,
                 domains=[('A', 1, 113), ('B', 114, 167), ('A', 168, 361), ('C', 362, 418)],
                 cat=[('A', 193, 'ASP', 'nucleophile'), ('A', 219, 'GLU', 'acid/base'),
                      ('A', 294, 'ASP', 'TS stabiliser')],
                 view=None, close_view=None),
    '3WEO': dict(family='GH31', chain='A', resolution=1.45,
                 domains=[('N', 38, 313), ('A', 314, 655), ('C', 656, 756), ('C2', 757, 909)],
                 cat=[('A', 469, 'ASP', 'nucleophile'), ('A', 568, 'ASP', 'acid/base')],
                 view=None, close_view=None),
}
# chosen directions (deposited frame) after inspecting work/views_<pdb>.png; None -> best of the search
VIEWS = {'6JQB': [0.612, -0.389, -0.689], '3WEO': [-0.615, 0.585, 0.529]}
CLOSE_VIEWS = {'6JQB': [0.8358, -0.0189, 0.5487], '3WEO': [-0.0125, -0.153, 0.9882]}


def fib_dirs(n):
    i = np.arange(n) + 0.5
    ph = np.arccos(1 - 2 * i / n)
    th = np.pi * (1 + 5 ** 0.5) * i
    return np.column_stack([np.cos(th) * np.sin(ph), np.sin(th) * np.sin(ph), np.cos(ph)])


def basis(ez):
    ez = ez / np.linalg.norm(ez)
    a = np.array([0, 1, 0]) if abs(ez[1]) < 0.9 else np.array([1, 0, 0])
    ex = np.cross(a, ez); ex /= np.linalg.norm(ex)
    ey = np.cross(ez, ex)
    return np.array([ex, ey, ez])


def collect(pdb):
    spec = ENTRIES[pdb]
    st = C.load(pdb)
    m = st[0]
    prot = [(a, r) for ch in m if ch.name == spec['chain'] for r in ch if r.het_flag == 'A' for a in r]
    P_all = np.array([a.pos.tolist() for a, r in prot])
    P_bb = np.array([a.pos.tolist() for a, r in prot if a.name in ('N', 'CA', 'C', 'O')])
    units = []
    for row in C.units_of(pdb):
        atoms, ring = C.unit_atoms(m, row)
        units.append(dict(sub=row['subsite'], chain=row['chain'], seq=int(row['seq']), comp=row['comp'],
                          ref=row['ref_atom'], rscc=row['rscc'], atoms=atoms, ring=ring,
                          cen=np.mean(list(ring.values()), 0)))
    cats = []
    for ch, seq, name, role in spec['cat']:
        r = C.find_res(m, ch, seq, name)
        X = {a.name: np.array(a.pos.tolist()) for a in r}
        sc = [x for n, x in X.items() if n not in ('N', 'C', 'O', 'CA')]
        cats.append(dict(chain=ch, seq=seq, name=name, role=role, X=X, cen=np.mean(sc, 0)))
    return dict(pdb=pdb, spec=spec, st=st, P_all=P_all, P_bb=P_bb, units=units, cats=cats)


def orient_inplane(B, a, b):
    """rotate basis about z so that page vector a->b points right"""
    v = (b - a) @ B.T
    ang = math.atan2(v[1], v[0])
    Rz = np.array([[math.cos(-ang), -math.sin(-ang), 0], [math.sin(-ang), math.cos(-ang), 0], [0, 0, 1]])
    return Rz @ B


def search_whole(info, n=4000, keep=8):
    shown = [u for u in info['units'] if u['sub'] in C.SUBSITES]
    order = sorted(shown, key=lambda u: C.SUBSITES.index(u['sub']))
    S = np.vstack([np.array(list(u['atoms'].values())) for u in shown])
    cen = np.array([u['cen'] for u in order])
    true_len = sum(np.linalg.norm(cen[i + 1] - cen[i]) for i in range(len(cen) - 1))
    sc = S.mean(0)
    Pb = info['P_bb'] - sc
    Pa = info['P_all'][::3] - sc
    Sc = S - sc
    cands = []
    for ez in fib_dirs(n):
        B = basis(ez)
        c2 = (cen - sc) @ B.T
        span = sum(np.linalg.norm(c2[i + 1, :2] - c2[i, :2]) for i in range(len(c2) - 1)) / true_len
        s2 = Sc @ B.T
        b2 = Pb @ B.T
        a2 = Pa @ B.T
        d = np.linalg.norm(b2[:, None, :2] - s2[None, :, :2], axis=2)
        front = b2[:, None, 2] > s2[None, :, 2] + 1.0
        cover = float(((d < 2.0) & front).any(axis=0).mean())
        r = np.linalg.norm(a2[:, :2], axis=1)
        infront = float(((r < 12) & (a2[:, 2] > 2.0)).sum() / max(1, (r < 12).sum()))
        score = 1.5 * span - 2.0 * cover - 1.5 * infront
        cands.append(dict(ez=ez, span=span, cover=cover, infront=infront, score=score))
    alts = []
    for c in sorted(cands, key=lambda c: -c['score']):
        if all(np.degrees(np.arccos(np.clip(c['ez'] @ a['ez'], -1, 1))) > 25 for a in alts):
            alts.append(c)
        if len(alts) >= keep:
            break
    return alts, order


def search_close(info, subs=('-1', '+1'), n=20000):
    us = [u for u in info['units'] if u['sub'] in subs]
    cats = info['cats']
    anchors = np.array([u['cen'] for u in us] + [c['cen'] for c in cats])
    S = np.vstack([np.array(list(u['atoms'].values())) for u in us])
    R = np.vstack([np.array(list(c['X'].values())) for c in cats])
    rid = np.concatenate([[k] * len(c['X']) for k, c in enumerate(cats)])
    ring_m1 = np.array(list(next(u for u in us if u['sub'] == '-1')['ring'].values()))
    nrm = np.linalg.svd(ring_m1 - ring_m1.mean(0))[2][2]
    ring_p1 = np.array(list(next(u for u in us if u['sub'] == '+1')['ring'].values()))
    nrm1 = np.linalg.svd(ring_p1 - ring_p1.mean(0))[2][2]
    ctr = np.vstack([S, R]).mean(0)
    uid = np.concatenate([[k] * len(u['atoms']) for k, u in enumerate(us)])
    best = None
    for ez in fib_dirs(n):
        tilt = abs(nrm @ ez); tilt1 = abs(nrm1 @ ez)
        if not (0.6 <= tilt <= 0.97 and 0.6 <= tilt1 <= 0.97):    # neither ring edge-on
            continue
        B = basis(ez)
        A = (anchors - ctr) @ B.T
        if np.linalg.norm(A[0, :2] - A[1, :2]) < 0.85 * np.linalg.norm(anchors[0] - anchors[1]):
            continue                                  # -1 -> +1 step seen nearly in the page plane
        D = np.linalg.norm(A[:, None, :2] - A[None, :, :2], axis=2) + 1e3 * np.eye(len(A))
        s2 = (S - ctr) @ B.T; r2 = (R - ctr) @ B.T
        occl = int((np.linalg.norm(r2[:, None, :2] - s2[None, :, :2], axis=2) < 1.1).any(axis=0).sum())
        rr = np.linalg.norm(r2[:, None, :2] - r2[None, :, :2], axis=2)
        cross = int(((rr < 0.9) & (rid[:, None] != rid[None, :])).sum() // 2)
        ss = np.linalg.norm(s2[:, None, :2] - s2[None, :, :2], axis=2)
        sclash = int(((ss < 0.7) & (uid[:, None] != uid[None, :])).sum() // 2)
        score = min(D.min(), 4.0) - 1.0 * sclash - 1.2 * occl - 1.0 * cross + 0.3 * tilt + 0.3 * tilt1
        if best is None or score > best[0]:
            best = (score, ez, D.min(), occl, cross, tilt)
    print(f"{info['pdb']} close-up: min anchor distance {best[2]:.2f} A, sugar atoms covered {best[3]}, "
          f"residue overlaps {best[4]}, -1 ring |n.z| {best[5]:.2f}")
    return best[1]


def write_rotated(info, B, ctr, path):
    out = info['st'].clone()
    for ch in out[0]:
        for r in ch:
            for a in r:
                a.pos = gemmi.Position(*((np.array(a.pos.tolist()) - ctr) @ B.T))
    out.setup_entities()
    out.write_pdb(str(path))


def rot(x, B, ctr):
    return (np.asarray(x) - ctr) @ B.T


def pml_header(px_w, px_h, half_h, slab_near, slab_far, fov=20.0):
    dist = half_h / math.tan(math.radians(fov / 2))
    return dist, ['bg_color white', 'set ray_opaque_background, 1', 'set orthoscopic, 1',
                  f'set field_of_view, {fov}', 'set ray_shadows, 0', 'set antialias, 2',
                  'set ambient, 0.45', 'set direct, 0.45', 'set reflect, 0.35', 'set specular, 0.15',
                  'set spec_power, 80', 'set depth_cue, 0', 'set ray_trace_mode, 1', 'set ray_trace_gain, 0.08',
                  'set ray_trace_color, black', 'set stick_quality, 24', 'set valence, 0',
                  'set cartoon_fancy_helices, 1', 'set cartoon_loop_radius, 0.25',
                  'set cartoon_oval_length, 1.2', 'set cartoon_rect_length, 1.2',
                  'set surface_quality, 1', 'set transparency_mode, 1',
                  f'set_color ocol, {C.hexrgb(C.O_COL)}', f'set_color ncol, {C.hexrgb(C.N_COL)}',
                  f'set_color catc, {C.hexrgb(C.CAT_C)}', f'set_color outc, {C.hexrgb(C.OUTSIDE_C)}',
                  f'set_color surfc, {C.hexrgb(C.SURF_C)}'] + \
        [f'set_color sub{i}, {C.hexrgb(C.SUGAR_C)}' for i, s in enumerate(C.SUBSITES)] + \
        [f'set_color dom_{k}, {C.hexrgb(v)}' for k, v in C.DOM_COL.items()]


def unit_sel(u):
    return f'chain {u["chain"]} and resi {u["seq"]} and resn {u["comp"]} and name {"+".join(u["atoms"])}'


def render_whole(info, name, ez, px=3000, surface=True, transp=0.62, style='transparent'):
    spec = info['spec']
    order = sorted([u for u in info['units'] if u['sub'] in C.SUBSITES], key=lambda u: C.SUBSITES.index(u['sub']))
    B = basis(np.asarray(ez, float))
    B = orient_inplane(B, order[0]['cen'], order[-1]['cen'])
    ctr = info['P_all'].mean(0)
    Pp = rot(info['P_all'], B, ctr)
    lo, hi = Pp[:, :2].min(0) - 3, Pp[:, :2].max(0) + 3
    mid = np.array([(lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2, 0.0])
    ctr = ctr + mid @ B          # recentre on the projected extent
    Pp = rot(info['P_all'], B, ctr)
    half_w, half_h = (hi - lo) / 2
    zr = np.abs(Pp[:, 2]).max() + 6
    pdbfile = WORK / f'{name}_rot.pdb'
    write_rotated(info, B, ctr, pdbfile)
    # PyMOL field: vertical half-height half_h; width set by the aspect of the ray-traced image
    pw, ph = px, int(round(px * half_h / half_w))
    dist, L = pml_header(pw, ph, half_h, 0, 0)
    L = ['reinitialize', f'load {pdbfile}, m', 'remove hydro', 'remove resn HOH', 'hide everything'] + L
    L += ['set cartoon_transparency, 0.0', f'select prot, polymer.protein and chain {spec["chain"]}']
    for dname, a, b in spec['domains']:
        L += [f'color dom_{dname}, prot and resi {a}-{b}']
    if style == 'opaque':
        # opaque molecular surface coloured by domain (cartoon hidden): the cleft reads as a groove
        L += [f'set surface_color, dom_{dn}, prot and resi {a}-{b}' for dn, a, b in spec['domains']]
        L += ['set surface_quality, 1', 'show surface, prot', 'set transparency, 0.0', 'set surface_mode, 1',
              'set ambient, 0.55', 'set direct, 0.35']
    else:
        L += ['show cartoon, prot']
    if surface and style != 'opaque':
        L += ['create surf, prot', 'hide cartoon, surf', 'show surface, surf', 'color surfc, surf',
              f'set transparency, {transp}, surf', 'set surface_mode, 1']
    for u in info['units']:
        s = unit_sel(u)
        col = f'sub{C.SUBSITES.index(u["sub"])}' if u['sub'] in C.SUBSITES else 'outc'
        L += [f'show sticks, {s}', f'set stick_radius, 0.30, {s}', f'color {col}, {s} and elem C',
              f'color ocol, {s} and elem O', f'color ncol, {s} and elem N']
    for c in info['cats']:
        s = f'prot and resi {c["seq"]} and not name N+C+O'
        L += [f'show sticks, {s}', f'set stick_radius, 0.30, {s}', f'color catc, {s} and elem C',
              f'color ocol, {s} and elem O', f'color ncol, {s} and elem N']
    L += [f'set_view (1,0,0, 0,1,0, 0,0,1, 0,0,{-dist:.3f}, 0,0,0, {dist - zr:.3f},{dist + zr:.3f}, 1)',
          f'ray {pw},{ph}', f'png {WORK / (name + ".png")}, dpi=300']
    pml = WORK / f'{name}.pml'
    pml.write_text('\n'.join(L) + '\n')
    subprocess.run(['pymol', '-cq', str(pml)], check=True, stdout=subprocess.DEVNULL)
    meta = dict(pdb=info['pdb'], png=str(WORK / (name + '.png')), px=[pw, ph], half=[float(half_w), float(half_h)],
                basis=B.tolist(), center=ctr.tolist(),
                units=[dict(sub=u['sub'], comp=u['comp'], chain=u['chain'], seq=u['seq'], ref=u['ref'],
                            rscc=u['rscc'], xyz=rot(u['cen'], B, ctr).tolist(),
                            atoms=[rot(x, B, ctr).tolist() for x in u['atoms'].values()]) for u in info['units']],
                cats=[dict(seq=c['seq'], name=c['name'], role=c['role'], xyz=rot(c['cen'], B, ctr).tolist(),
                           atoms=[rot(x, B, ctr).tolist() for x in c['X'].values()]) for c in info['cats']],
                domains=[dict(name=d, start=a, end=b,
                              xyz=rot(np.mean([x for x, (at, r) in zip(info['P_all'], [(a_, r_) for ch in info['st'][0]
                                   if ch.name == spec['chain'] for r_ in ch if r_.het_flag == 'A' for a_ in r_])
                                   if a <= r.seqid.num <= b and at.name == 'CA'], 0), B, ctr).tolist())
                         for d, a, b in spec['domains']])
    json.dump(meta, open(WORK / f'{name}.json', 'w'), indent=1)
    return meta


def render_close(info, name, ez, subs=('-1', '+1'), px=2400, slab=9.0, context_subs=()):
    spec = info['spec']
    us = [u for u in info['units'] if u['sub'] in subs]
    ctx = [u for u in info['units'] if u['sub'] in context_subs]
    B = basis(np.asarray(ez, float))
    m1 = next(u for u in us if u['sub'] == '-1'); p1 = next(u for u in us if u['sub'] == '+1')
    B = orient_inplane(B, m1['cen'], p1['cen'])
    pts = np.vstack([np.array(list(u['atoms'].values())) for u in us + ctx] +
                    [np.array(list(c['X'].values())) for c in info['cats']])
    ctr = pts.mean(0)
    Pp = rot(pts, B, ctr)
    lo, hi = Pp[:, :2].min(0) - 2.0, Pp[:, :2].max(0) + 2.0
    ctr = ctr + np.array([(lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2, 0.0]) @ B
    half_w, half_h = (hi - lo) / 2
    pdbfile = WORK / f'{name}_rot.pdb'
    write_rotated(info, B, ctr, pdbfile)
    pw, ph = px, int(round(px * half_h / half_w))
    dist, L = pml_header(pw, ph, half_h, 0, 0)
    L = ['reinitialize', f'load {pdbfile}, m', 'remove hydro', 'remove resn HOH', 'hide everything'] + L
    L += [f'select prot, polymer.protein and chain {spec["chain"]}',
          'set cartoon_transparency, 0.75', 'set cartoon_color, grey80', 'show cartoon, prot',
          'set cartoon_loop_radius, 0.12', 'set cartoon_oval_width, 0.18', 'set cartoon_rect_width, 0.25']
    for u in ctx:
        s = unit_sel(u)
        L += [f'show sticks, {s}', f'set stick_radius, 0.10, {s}', f'color grey65, {s}']
    for u in us:
        s = unit_sel(u)
        col = f'sub{C.SUBSITES.index(u["sub"])}'
        L += [f'show sticks, {s}', f'set stick_radius, 0.15, {s}', f'color {col}, {s} and elem C',
              f'color ocol, {s} and elem O', f'color ncol, {s} and elem N']
    for c in info['cats']:
        s = f'prot and resi {c["seq"]} and not name N+C+O'
        L += [f'show sticks, {s}', f'set stick_radius, 0.13, {s}', f'color catc, {s} and elem C',
              f'color ocol, {s} and elem O', f'color ncol, {s} and elem N']
    L += [f'set_view (1,0,0, 0,1,0, 0,0,1, 0,0,{-dist:.3f}, 0,0,0, {dist - slab:.3f},{dist + slab:.3f}, 1)',
          f'ray {pw},{ph}', f'png {WORK / (name + ".png")}, dpi=300']
    pml = WORK / f'{name}.pml'
    pml.write_text('\n'.join(L) + '\n')
    subprocess.run(['pymol', '-cq', str(pml)], check=True, stdout=subprocess.DEVNULL)
    meta = dict(pdb=info['pdb'], png=str(WORK / (name + '.png')), px=[pw, ph], half=[float(half_w), float(half_h)],
                units=[dict(sub=u['sub'], xyz=rot(u['cen'], B, ctr).tolist(),
                            atoms=[rot(x, B, ctr).tolist() for x in u['atoms'].values()]) for u in us],
                cats=[dict(seq=c['seq'], name=c['name'], role=c['role'], xyz=rot(c['cen'], B, ctr).tolist(),
                           atoms=[rot(x, B, ctr).tolist() for x in c['X'].values()],
                           X={n: rot(x, B, ctr).tolist() for n, x in c['X'].items()}) for c in info['cats']])
    json.dump(meta, open(WORK / f'{name}.json', 'w'), indent=1)
    return meta


def montage(info, alts):
    """quick low-resolution renders of the candidate directions, for choosing the view"""
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    pngs = []
    for k, a in enumerate(alts):
        meta = render_whole(info, f'view_{info["pdb"]}_{k}', a['ez'], px=700)
        pngs.append((meta['png'], a))
    fig, axs = plt.subplots(2, 4, figsize=(16, 8))
    for ax, (p, a) in zip(axs.ravel(), pngs):
        ax.imshow(plt.imread(p)); ax.set_axis_off()
        ax.set_title(f"{np.round(a['ez'], 3).tolist()}\nspan {a['span']:.2f} cover {a['cover']:.2f} "
                     f"front {a['infront']:.2f}", fontsize=8)
    fig.savefig(WORK / f'views_{info["pdb"]}.png', dpi=80)
    plt.close(fig)


def main(args):
    for pdb in ENTRIES:
        info = collect(pdb)
        if 'search' in args:
            alts, _ = search_whole(info)
            for a in alts:
                print(pdb, np.round(a['ez'], 4).tolist(), f"span {a['span']:.2f} cover {a['cover']:.2f} front {a['infront']:.2f}")
            montage(info, alts)
            continue
        ez = VIEWS[pdb] if VIEWS[pdb] is not None else search_whole(info)[0][0]['ez']
        if 'close' not in args and 'only_close' not in args:
            render_whole(info, f'whole_{pdb}', ez)
            render_whole(info, f'whole_{pdb}_opaque', ez, style='opaque')
        if 'whole_only' in args:
            continue
        cz = CLOSE_VIEWS[pdb] if CLOSE_VIEWS[pdb] is not None else search_close(info)
        print(pdb, 'close view', np.round(cz, 4).tolist())
        render_close(info, f'close_{pdb}', cz)


if __name__ == '__main__':
    main(sys.argv[1:])
