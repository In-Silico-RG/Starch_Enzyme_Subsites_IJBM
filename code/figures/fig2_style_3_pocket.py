"""Style 3: pocket surface. ChimeraX solvent-excluded surface of the protein (chain A, probe 1.4 A) shown only
within 7 A of the sugar units, semi-transparent; surface patches of the figure residues coloured by their class,
the rest light grey; residues as thin sticks under the surface, sugar as sticks (ball-and-stick) inside.
Writes figures/candidates/styles/style_3*.png.
Run: ~/miniconda3/bin/python code/figures/fig2_style_3_pocket.py
"""
import fig2_style_common as C

info = C.setup()
png = C.WORK / 'style_3_render.png'
sug = ' | '.join(C.unit_sel_cx(u) for u in info['units'])
L = C.cx_atoms(info, stick=0.12, sugar_stick=0.12, sugar_ball=0.14)
L += ['surface /A & protein enclose /A & protein gridSpacing 0.3',
      f'surface zone /A nearAtoms {sug} distance 7.0',
      'color /A & protein #E6E6E6 target s']
for rd in info['residues']:
    L.append(f'color /{rd["chain"]}:{rd["seq"]} {C.COL[rd["cls"]]} target s')
L += ['transparency /A 62 target s',
      'lighting soft', 'lighting intensity 0.35 ambientIntensity 1.25', 'material dull',
      'graphics silhouettes true width 3 color black depthJump 0.02', 'graphics quality 20']
L += C.chimerax_camera(info, png)
C.run_chimerax(L, 'style_3')
C.compose(info, png, 'style_3')
