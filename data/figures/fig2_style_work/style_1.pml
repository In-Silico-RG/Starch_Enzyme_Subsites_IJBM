reinitialize
load <PROJECT_ROOT>/resultados/figures/fig2_style_work/1G94_style_rot.pdb, m
remove hydro
remove solvent
hide everything
bg_color white
set ray_opaque_background, 1
set antialias, 2
set hash_max, 300
set light_count, 3
set ambient, 0.42
set direct, 0.35
set reflect, 0.55
set specular, 0.25
set spec_power, 80
set spec_reflect, 0.25
set ray_shadows, 0
set ray_light_scale, 0.6
set depth_cue, 1
set fog_start, 0.55
set fog, 0.6
set ray_trace_mode, 1
set ray_trace_gain, 0.22
set ray_trace_color, black
set ray_trace_disco_factor, 1
set stick_radius, 0.14
set stick_quality, 30
set sphere_quality, 4
set stick_ball, 0
set valence, 0
set dash_gap, 0.20
set dash_length, 0.16
set dash_radius, 0.045
set dash_color, grey20
set dash_round_ends, 1
set cartoon_transparency, 0.88
set cartoon_color, grey80
set cartoon_tube_radius, 0.12
set_color sugc, [0.8313725490196079, 0.8313725490196079, 0.8313725490196079]
set_color ocol, [0.8784313725490196, 0.1411764705882353, 0.16862745098039217]
set_color ncol, [0.1843137254901961, 0.30980392156862746, 0.8470588235294118]
select ctx, polymer.protein and name CA and byres (polymer.protein within 13 of ((chain C and resi 3 and resn DAF and name N4+C1H+C2H+O2H+C3H+O3H+C4H+O4H+C5H+C6H+C7H+O7H) or (chain C and resi 3 and resn DAF and name C1+C2+O2+C3+O3+C4+C5+O5+C6) or (chain C and resi 2 and resn GLC and name C1+C2+C3+C4+C5+C6+O2+O3+O4+O5+O6)))
set cartoon_trace_atoms, 0
show cartoon, byres ctx
cartoon tube, byres ctx
show sticks, (chain C and resi 3 and resn DAF and name N4+C1H+C2H+O2H+C3H+O3H+C4H+O4H+C5H+C6H+C7H+O7H)
show spheres, (chain C and resi 3 and resn DAF and name N4+C1H+C2H+O2H+C3H+O3H+C4H+O4H+C5H+C6H+C7H+O7H)
set sphere_scale, 0.19, (chain C and resi 3 and resn DAF and name N4+C1H+C2H+O2H+C3H+O3H+C4H+O4H+C5H+C6H+C7H+O7H)
set stick_radius, 0.095, (chain C and resi 3 and resn DAF and name N4+C1H+C2H+O2H+C3H+O3H+C4H+O4H+C5H+C6H+C7H+O7H)
color sugc, (chain C and resi 3 and resn DAF and name N4+C1H+C2H+O2H+C3H+O3H+C4H+O4H+C5H+C6H+C7H+O7H) and elem C
color ocol, (chain C and resi 3 and resn DAF and name N4+C1H+C2H+O2H+C3H+O3H+C4H+O4H+C5H+C6H+C7H+O7H) and elem O
color ncol, (chain C and resi 3 and resn DAF and name N4+C1H+C2H+O2H+C3H+O3H+C4H+O4H+C5H+C6H+C7H+O7H) and elem N
show sticks, (chain C and resi 3 and resn DAF and name C1+C2+O2+C3+O3+C4+C5+O5+C6)
show spheres, (chain C and resi 3 and resn DAF and name C1+C2+O2+C3+O3+C4+C5+O5+C6)
set sphere_scale, 0.19, (chain C and resi 3 and resn DAF and name C1+C2+O2+C3+O3+C4+C5+O5+C6)
set stick_radius, 0.095, (chain C and resi 3 and resn DAF and name C1+C2+O2+C3+O3+C4+C5+O5+C6)
color sugc, (chain C and resi 3 and resn DAF and name C1+C2+O2+C3+O3+C4+C5+O5+C6) and elem C
color ocol, (chain C and resi 3 and resn DAF and name C1+C2+O2+C3+O3+C4+C5+O5+C6) and elem O
color ncol, (chain C and resi 3 and resn DAF and name C1+C2+O2+C3+O3+C4+C5+O5+C6) and elem N
show sticks, (chain C and resi 2 and resn GLC and name C1+C2+C3+C4+C5+C6+O2+O3+O4+O5+O6)
show spheres, (chain C and resi 2 and resn GLC and name C1+C2+C3+C4+C5+C6+O2+O3+O4+O5+O6)
set sphere_scale, 0.19, (chain C and resi 2 and resn GLC and name C1+C2+C3+C4+C5+C6+O2+O3+O4+O5+O6)
set stick_radius, 0.095, (chain C and resi 2 and resn GLC and name C1+C2+C3+C4+C5+C6+O2+O3+O4+O5+O6)
color sugc, (chain C and resi 2 and resn GLC and name C1+C2+C3+C4+C5+C6+O2+O3+O4+O5+O6) and elem C
color ocol, (chain C and resi 2 and resn GLC and name C1+C2+C3+C4+C5+C6+O2+O3+O4+O5+O6) and elem O
color ncol, (chain C and resi 2 and resn GLC and name C1+C2+C3+C4+C5+C6+O2+O3+O4+O5+O6) and elem N
set_color cls_H, [0.0, 0.4470588235294118, 0.6980392156862745]
show sticks, (polymer.protein and chain A and resi 200 and not name N+C+O+OXT)
color cls_H, (polymer.protein and chain A and resi 200 and not name N+C+O+OXT) and elem C
color ocol, (polymer.protein and chain A and resi 200 and not name N+C+O+OXT) and elem O
color ncol, (polymer.protein and chain A and resi 200 and not name N+C+O+OXT) and elem N+S
set_color cls_Al, [0.8, 0.4745098039215686, 0.6549019607843137]
show sticks, (polymer.protein and chain A and resi 175 and not name N+C+O+OXT)
color cls_Al, (polymer.protein and chain A and resi 175 and not name N+C+O+OXT) and elem C
color ocol, (polymer.protein and chain A and resi 175 and not name N+C+O+OXT) and elem O
color ncol, (polymer.protein and chain A and resi 175 and not name N+C+O+OXT) and elem N+S
set_color cls_H, [0.0, 0.4470588235294118, 0.6980392156862745]
show sticks, (polymer.protein and chain A and resi 177 and not name N+C+O+OXT)
color cls_H, (polymer.protein and chain A and resi 177 and not name N+C+O+OXT) and elem C
color ocol, (polymer.protein and chain A and resi 177 and not name N+C+O+OXT) and elem O
color ncol, (polymer.protein and chain A and resi 177 and not name N+C+O+OXT) and elem N+S
set_color cls_Al, [0.8, 0.4745098039215686, 0.6549019607843137]
show sticks, (polymer.protein and chain A and resi 202 and not name N+C+O+OXT)
color cls_Al, (polymer.protein and chain A and resi 202 and not name N+C+O+OXT) and elem C
color ocol, (polymer.protein and chain A and resi 202 and not name N+C+O+OXT) and elem O
color ncol, (polymer.protein and chain A and resi 202 and not name N+C+O+OXT) and elem N+S
set_color cls_H, [0.0, 0.4470588235294118, 0.6980392156862745]
show sticks, (polymer.protein and chain A and resi 174 and not name N+C+O+OXT)
color cls_H, (polymer.protein and chain A and resi 174 and not name N+C+O+OXT) and elem C
color ocol, (polymer.protein and chain A and resi 174 and not name N+C+O+OXT) and elem O
color ncol, (polymer.protein and chain A and resi 174 and not name N+C+O+OXT) and elem N+S
set_color cls_H, [0.0, 0.4470588235294118, 0.6980392156862745]
show sticks, (polymer.protein and chain A and resi 264 and not name N+C+O+OXT)
color cls_H, (polymer.protein and chain A and resi 264 and not name N+C+O+OXT) and elem C
color ocol, (polymer.protein and chain A and resi 264 and not name N+C+O+OXT) and elem O
color ncol, (polymer.protein and chain A and resi 264 and not name N+C+O+OXT) and elem N+S
set_color cls_H, [0.0, 0.4470588235294118, 0.6980392156862745]
show sticks, (polymer.protein and chain A and resi 263 and not name N+C+O+OXT)
color cls_H, (polymer.protein and chain A and resi 263 and not name N+C+O+OXT) and elem C
color ocol, (polymer.protein and chain A and resi 263 and not name N+C+O+OXT) and elem O
color ncol, (polymer.protein and chain A and resi 263 and not name N+C+O+OXT) and elem N+S
set_color cls_H, [0.0, 0.4470588235294118, 0.6980392156862745]
show sticks, (polymer.protein and chain A and resi 172 and not name N+C+O+OXT)
color cls_H, (polymer.protein and chain A and resi 172 and not name N+C+O+OXT) and elem C
color ocol, (polymer.protein and chain A and resi 172 and not name N+C+O+OXT) and elem O
color ncol, (polymer.protein and chain A and resi 172 and not name N+C+O+OXT) and elem N+S
set_color cls_H, [0.0, 0.4470588235294118, 0.6980392156862745]
show sticks, (polymer.protein and chain A and resi 89 and not name N+C+O+OXT)
color cls_H, (polymer.protein and chain A and resi 89 and not name N+C+O+OXT) and elem C
color ocol, (polymer.protein and chain A and resi 89 and not name N+C+O+OXT) and elem O
color ncol, (polymer.protein and chain A and resi 89 and not name N+C+O+OXT) and elem N+S
set_color cls_Ar, [0.9019607843137255, 0.6235294117647059, 0.0]
show sticks, (polymer.protein and chain A and resi 50 and not name N+C+O+OXT)
color cls_Ar, (polymer.protein and chain A and resi 50 and not name N+C+O+OXT) and elem C
color ocol, (polymer.protein and chain A and resi 50 and not name N+C+O+OXT) and elem O
color ncol, (polymer.protein and chain A and resi 50 and not name N+C+O+OXT) and elem N+S
distance hb0, polymer.protein and chain A and resi 174 and name OD1, chain C and resi 3 and resn DAF and name O7H
distance hb1, polymer.protein and chain A and resi 264 and name OD1, chain C and resi 3 and resn DAF and name O3H
distance hb2, polymer.protein and chain A and resi 264 and name OD2, chain C and resi 3 and resn DAF and name O2H
distance hb3, polymer.protein and chain A and resi 200 and name OE1, chain C and resi 3 and resn DAF and name O3
distance hb4, polymer.protein and chain A and resi 200 and name OE2, chain C and resi 3 and resn DAF and name N4
distance hb5, polymer.protein and chain A and resi 263 and name NE2, chain C and resi 3 and resn DAF and name O3H
distance hb6, polymer.protein and chain A and resi 89 and name NE2, chain C and resi 3 and resn DAF and name O7H
distance hb7, polymer.protein and chain A and resi 177 and name NZ, chain C and resi 2 and resn GLC and name O3
distance hb8, polymer.protein and chain A and resi 172 and name NH1, chain C and resi 3 and resn DAF and name O2H
hide labels
set field_of_view, 16.0
set orthoscopic, 0
set_view (1,0,0, 0,1,0, 0,0,1, 0,0,-80.281, 0,0,0, 72.207,94.880, -16.0)
ray 2400,2400
png <PROJECT_ROOT>/resultados/figures/fig2_style_work/style_1_render.png, dpi=300
