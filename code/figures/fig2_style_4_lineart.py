"""Style 4: minimalist line-art. Only the sugar units (ball-and-stick, neutral grey carbons) and the side chains
from CA (class-coloured carbons), hydrogen-bond dashes; nothing else. ChimeraX flat shading (no lighting
gradient) with thick silhouette outlines, orthographic camera: a journal line-drawing look.
Writes figures/candidates/styles/style_4*.png.
Run: ~/miniconda3/bin/python code/figures/fig2_style_4_lineart.py
"""
import fig2_style_common as C

info = C.setup()
info['ortho'] = True
png = C.WORK / 'style_4_render.png'
L = C.cx_atoms(info, stick=0.15, sugar_stick=0.12, sugar_ball=0.17)
L += ['lighting flat', 'graphics silhouettes true width 5 color black depthJump 0.01', 'graphics quality 20']
L += C.chimerax_camera(info, png)
C.run_chimerax(L, 'style_4')
C.compose(info, png, 'style_4')
