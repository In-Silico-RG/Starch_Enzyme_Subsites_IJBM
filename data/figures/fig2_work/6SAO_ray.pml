reinitialize
load <PROJECT_ROOT>/resultados/figures/fig2_work/6SAO_rot.pdb, m
hide everything
bg_color white
set ray_opaque_background, 1
set orthoscopic, 1
set field_of_view, 20.0
set ray_shadows, 0
set antialias, 2
set ambient, 0.45
set specular, 0.15
set direct, 0.6
set depth_cue, 0
set stick_radius, 0.17
set dash_gap, 0.25
set dash_length, 0.25
set dash_radius, 0.06
set dash_color, black
set cartoon_transparency, 0.82
set cartoon_color, grey70
set ray_trace_mode, 1
set ray_trace_gain, 0.6
set ray_trace_color, grey25
show sticks, chain B and resi 3 and resn AC1
set stick_radius, 0.22, chain B and resi 3 and resn AC1
set_color sugc, [0.18, 0.55, 0.34]
color sugc, chain B and resi 3 and resn AC1 and elem C
color red, chain B and resi 3 and resn AC1 and elem O
color blue, chain B and resi 3 and resn AC1 and elem N
set_color cls_H, [0.122, 0.373, 0.749]
show sticks, chain A and resi 280 and not resn AC1 and not name N+C+O
color cls_H, chain A and resi 280 and not resn AC1 and elem C
color red, chain A and resi 280 and not resn AC1 and elem O
color blue, chain A and resi 280 and not resn AC1 and elem N
set_color cls_H, [0.122, 0.373, 0.749]
show sticks, chain A and resi 114 and not resn AC1 and not name N+C+O
color cls_H, chain A and resi 114 and not resn AC1 and elem C
color red, chain A and resi 114 and not resn AC1 and elem O
color blue, chain A and resi 114 and not resn AC1 and elem N
set_color cls_H, [0.122, 0.373, 0.749]
show sticks, chain A and resi 190 and not resn AC1 and not name N+C+O
color cls_H, chain A and resi 190 and not resn AC1 and elem C
color red, chain A and resi 190 and not resn AC1 and elem O
color blue, chain A and resi 190 and not resn AC1 and elem N
set_color cls_H, [0.122, 0.373, 0.749]
show sticks, chain A and resi 188 and not resn AC1 and not name N+C+O
color cls_H, chain A and resi 188 and not resn AC1 and elem C
color red, chain A and resi 188 and not resn AC1 and elem O
color blue, chain A and resi 188 and not resn AC1 and elem N
set_color cls_Ar, [0.878, 0.482, 0.0]
show sticks, chain A and resi 74 and not resn AC1 and not name N+C+O
color cls_Ar, chain A and resi 74 and not resn AC1 and elem C
color red, chain A and resi 74 and not resn AC1 and elem O
color blue, chain A and resi 74 and not resn AC1 and elem N
set_color cls_H, [0.122, 0.373, 0.749]
show sticks, chain A and resi 279 and not resn AC1 and not name N+C+O
color cls_H, chain A and resi 279 and not resn AC1 and elem C
color red, chain A and resi 279 and not resn AC1 and elem O
color blue, chain A and resi 279 and not resn AC1 and elem N
set_color cls_Al, [0.55, 0.55, 0.55]
show sticks, chain A and resi 191 and not resn AC1 and not name N+C+O
color cls_Al, chain A and resi 191 and not resn AC1 and elem C
color red, chain A and resi 191 and not resn AC1 and elem O
color blue, chain A and resi 191 and not resn AC1 and elem N
set_color cls_H, [0.122, 0.373, 0.749]
show sticks, chain A and resi 194 and not resn AC1 and not name N+C+O
color cls_H, chain A and resi 194 and not resn AC1 and elem C
color red, chain A and resi 194 and not resn AC1 and elem O
color blue, chain A and resi 194 and not resn AC1 and elem N
set_color cls_H, [0.122, 0.373, 0.749]
show sticks, chain A and resi 214 and not resn AC1 and not name N+C+O
color cls_H, chain A and resi 214 and not resn AC1 and elem C
color red, chain A and resi 214 and not resn AC1 and elem O
color blue, chain A and resi 214 and not resn AC1 and elem N
distance hb0, chain A and resi 280 and name OD1 and not resn AC1, chain B and resi 3 and resn AC1 and name O3B
distance hb1, chain A and resi 280 and name OD2 and not resn AC1, chain B and resi 3 and resn AC1 and name O2B
distance hb2, chain A and resi 280 and name OD2 and not resn AC1, chain B and resi 3 and resn AC1 and name O3B
distance hb3, chain A and resi 114 and name NE2 and not resn AC1, chain B and resi 3 and resn AC1 and name O6B
distance hb4, chain A and resi 190 and name OD2 and not resn AC1, chain B and resi 3 and resn AC1 and name O6B
distance hb5, chain A and resi 188 and name NH2 and not resn AC1, chain B and resi 3 and resn AC1 and name O2B
distance hb6, chain A and resi 279 and name NE2 and not resn AC1, chain B and resi 3 and resn AC1 and name O2B
distance hb7, chain A and resi 279 and name NE2 and not resn AC1, chain B and resi 3 and resn AC1 and name O3B
distance hb8, chain A and resi 194 and name NE2 and not resn AC1, chain B and resi 3 and resn AC1 and name O2
distance hb9, chain A and resi 214 and name OE1 and not resn AC1, chain B and resi 3 and resn AC1 and name O3
distance hb10, chain A and resi 214 and name OE2 and not resn AC1, chain B and resi 3 and resn AC1 and name O3
distance hb11, chain A and resi 214 and name OE2 and not resn AC1, chain B and resi 3 and resn AC1 and name O2B
hide labels
set dash_round_ends, 1
set_view (1,0,0, 0,1,0, 0,0,1, 0,0,-69.275, 0,0,0, 44.275,94.275, 1)
ray 2000,2000
png <PROJECT_ROOT>/resultados/figures/fig2_work/6SAO_ray.png, dpi=300
