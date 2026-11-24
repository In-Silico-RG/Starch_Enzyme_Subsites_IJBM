"""Style 2: ChimeraX, soft ambient-occlusion lighting (lighting soft: 64 directional shadows), silhouette edges,
thin sticks, sugar ball-and-stick, no cartoon, slight perspective.
Writes figures/candidates/styles/style_2*.png.
Run: ~/miniconda3/bin/python code/figures/fig2_style_2_chimerax.py
"""
import fig2_style_common as C

info = C.setup()
png = C.WORK / 'style_2_render.png'
L = C.cx_atoms(info, stick=0.13, sugar_stick=0.095, sugar_ball=0.16)
L += ['lighting soft', 'lighting intensity 0.35 ambientIntensity 1.25', 'material dull',
      'graphics silhouettes true width 3.5 color black depthJump 0.01', 'graphics quality 20']
L += C.chimerax_camera(info, png)
C.run_chimerax(L, 'style_2')
C.compose(info, png, 'style_2')
