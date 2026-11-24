"""Whole-protein panels of the Figure 1 candidates rendered with ChimeraX (Recipe A, cleft variant, of
figures/references/REFERENCES.md section C): soft lighting (ambient occlusion), black
silhouettes, cartoon in desaturated domain colours, light-grey molecular surface 65 % transparent (style
'transparent', candidate A) or opaque surface in the domain colours (style 'opaque', candidate B); the sugar
units ball-and-stick with carbons in the subsite colour (units outside -3..+3 grey), catalytic residues as sticks
with dark-grey carbons.

Geometry is the same as fig1_protein_struct.py: the coordinates are rotated into the chosen view (VIEWS there),
written to work/<name>_rot.pdb, and photographed by an orthographic camera looking down -z with field width
2*half_w, so a point (x, y) of the rotated frame lands at the same place as in the JSON written next to the image
(the label overlay uses those page coordinates).
Run: ~/miniconda3/bin/python code/figures/fig1_protein_cx.py      (ChimeraX 1.11.1, chimerax --nogui --offscreen)
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fig1_protein_common as C
import fig1_protein_struct as S

WORK = C.WORK


def hexs(h):
    return h if h.startswith('#') else '#' + h


def render(pdb, style, px=2400):
    info = S.collect(pdb)
    spec = info['spec']
    order = sorted([u for u in info['units'] if u['sub'] in C.SUBSITES], key=lambda u: C.SUBSITES.index(u['sub']))
    B = S.basis(np.asarray(S.VIEWS[pdb], float))
    B = S.orient_inplane(B, order[0]['cen'], order[-1]['cen'])
    ctr = info['P_all'].mean(0)
    Pp = S.rot(info['P_all'], B, ctr)
    lo, hi = Pp[:, :2].min(0) - 3, Pp[:, :2].max(0) + 3
    ctr = ctr + np.array([(lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2, 0.0]) @ B
    half_w, half_h = (hi - lo) / 2
    name = f'whole_{pdb}' + ('' if style == 'transparent' else '_' + style)
    pdbfile = WORK / f'{name}_rot.pdb'
    S.write_rotated(info, B, ctr, pdbfile)
    pw, ph = px, int(round(px * half_h / half_w))
    ch = spec['chain']
    prot = f'/{ch} & protein'
    L = [f'open {pdbfile}', 'delete solvent', 'delete H', 'set bgColor white', 'graphics quality 5',
         'lighting soft multiShadow 256', 'lighting depthCue false', 'material dull',
         'graphics silhouettes true width 3 color black depthJump 0.01',
         'hide atoms', 'hide cartoons', 'hide pbonds', 'label delete', f'cartoon {prot}', 'cartoon style width 1.4 thickness 0.35']
    for dn, a, b in spec['domains']:
        L.append(f'color /{ch}:{a}-{b} {C.DOM_COL[dn]} target c')
    L.append(f'surface {prot} enclose {prot}')
    if style == 'transparent':
        L += [f'color {prot} {C.SURF_C} target s', f'transparency {prot} 60 target s']
    else:
        for dn, a, b in spec['domains']:
            L.append(f'color /{ch}:{a}-{b} {C.DOM_COL[dn]} target s')
        L += [f'transparency {prot} 0 target s']
    for k, u in enumerate(info['units']):
        sel = f'/{u["chain"]}:{u["seq"]}@{",".join(u["atoms"])}'
        inside = u['sub'] in C.SUBSITES
        col = C.SUGAR_C if inside else C.CTX_C
        L += [f'name u{k} {sel}', f'show u{k} atoms', f'style u{k} stick',
              f'size u{k} stickRadius {0.34 if inside else 0.20}',
              f'color u{k} & C {col} target a', f'color u{k} & O {C.O_COL} target a',
              f'color u{k} & N {C.N_COL} target a']
    for c in info['cats']:
        sel = f'/{c["chain"]}:{c["seq"]} & sidechain'
        L += [f'show {sel} atoms', f'style {sel} stick', f'size {sel} stickRadius 0.28',
              f'color {sel} & C {C.CAT_C} target a',
              f'color {sel} & O {C.O_COL} target a']
    # 'camera ortho' fails offscreen in 1.11.1 (no stereo attribute): set the orthographic camera from Python
    py = ['from chimerax.core.commands import run', 'from chimerax.graphics import OrthographicCamera',
          'from chimerax.geometry import Place']
    py += [f'run(session, {c!r})' for c in L]
    py += ['cam = OrthographicCamera()', f'cam.field_width = {2 * half_w:.4f}',
           'cam.position = Place(origin=(0.0, 0.0, 400.0))', 'session.main_view.camera = cam',
           "run(session, 'clip off')",
           f"run(session, {('save ' + str(WORK / (name + '.png')) + f' width {pw} height {ph} supersample 2')!r})",
           "run(session, 'exit')"]
    cxc = WORK / f'{name}_cx.py'
    cxc.write_text('\n'.join(py) + '\n')
    subprocess.run(['chimerax', '--nogui', '--offscreen', '--exit', str(cxc)], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    meta = dict(pdb=pdb, style=style, renderer='ChimeraX 1.11.1', png=str(WORK / (name + '.png')), px=[pw, ph],
                half=[float(half_w), float(half_h)], basis=B.tolist(), center=ctr.tolist(),
                units=[dict(sub=u['sub'], comp=u['comp'], chain=u['chain'], seq=u['seq'], ref=u['ref'],
                            rscc=u['rscc'], xyz=S.rot(u['cen'], B, ctr).tolist(),
                            atoms=[S.rot(x, B, ctr).tolist() for x in u['atoms'].values()]) for u in info['units']],
                cats=[dict(seq=c['seq'], name=c['name'], role=c['role'], xyz=S.rot(c['cen'], B, ctr).tolist(),
                           atoms=[S.rot(x, B, ctr).tolist() for x in c['X'].values()]) for c in info['cats']])
    json.dump(meta, open(WORK / f'{name}.json', 'w'), indent=1)
    return meta


if __name__ == '__main__':
    pdbs = [a for a in sys.argv[1:] if a in S.ENTRIES] or list(S.ENTRIES)
    styles = [a for a in sys.argv[1:] if a in ('transparent', 'opaque')] or ['transparent', 'opaque']
    for pdb in pdbs:
        for st in styles:
            m = render(pdb, st)
            print('rendered', m['png'])
