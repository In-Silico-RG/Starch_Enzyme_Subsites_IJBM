"""Style 1: PyMOL ray-traced, black outlines (ray_trace_mode 1), thin sticks, sugar ball-and-stick, faint thin
backbone trace, depth cue, slight perspective. Writes figures/candidates/styles/style_1*.png.
Run: ~/miniconda3/bin/python code/figures/fig2_style_1_pymol.py
"""
import fig2_style_common as C

info = C.setup()
png = C.WORK / 'style_1_render.png'
L = ['reinitialize', f'load {info["pdbfile"]}, m', 'remove hydro', 'remove solvent', 'hide everything',
     'bg_color white', 'set ray_opaque_background, 1', 'set antialias, 2', 'set hash_max, 300',
     # lighting: two lights, soft ambient, low specular, no hard shadows
     'set light_count, 3', 'set ambient, 0.42', 'set direct, 0.35', 'set reflect, 0.55', 'set specular, 0.25',
     'set spec_power, 80', 'set spec_reflect, 0.25', 'set ray_shadows, 0', 'set ray_light_scale, 0.6',
     'set depth_cue, 1', 'set fog_start, 0.55', 'set fog, 0.6',
     # outlines
     'set ray_trace_mode, 1', 'set ray_trace_gain, 0.22', 'set ray_trace_color, black',
     'set ray_trace_disco_factor, 1',
     # geometry
     'set stick_radius, 0.14', 'set stick_quality, 30', 'set sphere_quality, 4', 'set stick_ball, 0',
     'set valence, 0',
     'set dash_gap, 0.20', 'set dash_length, 0.16', 'set dash_radius, 0.045', 'set dash_color, grey20',
     'set dash_round_ends, 1',
     'set cartoon_transparency, 0.88', 'set cartoon_color, grey80', 'set cartoon_tube_radius, 0.12',
     f'set_color sugc, {C.hexrgb(C.SUGAR_C)}', f'set_color ocol, {C.hexrgb(C.O_COL)}',
     f'set_color ncol, {C.hexrgb(C.N_COL)}']
# faint thin backbone trace near the site only
near = ' or '.join(C.unit_sel_pymol(u) for u in info['units'])
L += [f'select ctx, polymer.protein and name CA and byres (polymer.protein within 13 of ({near}))',
      'set cartoon_trace_atoms, 0', 'show cartoon, byres ctx', 'cartoon tube, byres ctx']
for u in info['units']:
    s = C.unit_sel_pymol(u)
    L += [f'show sticks, {s}', f'show spheres, {s}', f'set sphere_scale, 0.19, {s}', f'set stick_radius, 0.095, {s}',
          f'color sugc, {s} and elem C', f'color ocol, {s} and elem O', f'color ncol, {s} and elem N']
# bonds between sugar units (inter-residue link inside DAF already covalent; GLC-DAF link is O4-C1)
for rd in info['residues']:
    s = C.res_sel_pymol(rd)
    L += [f'set_color cls_{rd["cls"]}, {C.hexrgb(C.COL[rd["cls"]])}', f'show sticks, {s}',
          f'color cls_{rd["cls"]}, {s} and elem C', f'color ocol, {s} and elem O', f'color ncol, {s} and elem N+S']
for k, ((rc, rs, ra), (uc, us, ucomp, ub), d) in enumerate(C.hbond_pairs(info)):
    L.append(f'distance hb{k}, polymer.protein and chain {rc} and resi {rs} and name {ra}, '
             f'chain {uc} and resi {us} and resn {ucomp} and name {ub}')
L += ['hide labels'] + C.pymol_camera(info) + [f'ray {C.PX},{C.PX}', f'png {png}, dpi=300']
C.run_pymol(L, 'style_1')
C.compose(info, png, 'style_1')
