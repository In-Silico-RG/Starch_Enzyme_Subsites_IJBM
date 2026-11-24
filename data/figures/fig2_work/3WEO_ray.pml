reinitialize
load <PROJECT_ROOT>/resultados/figures/fig2_work/3WEO_rot.pdb, m
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
show sticks, chain D and resi 7 and resn AC1
set stick_radius, 0.22, chain D and resi 7 and resn AC1
set_color sugc, [0.18, 0.55, 0.34]
color sugc, chain D and resi 7 and resn AC1 and elem C
color red, chain D and resi 7 and resn AC1 and elem O
color blue, chain D and resi 7 and resn AC1 and elem N
set_color cls_H, [0.122, 0.373, 0.749]
show sticks, chain A and resi 357 and not resn AC1 and not name N+C+O
color cls_H, chain A and resi 357 and not resn AC1 and elem C
color red, chain A and resi 357 and not resn AC1 and elem O
color blue, chain A and resi 357 and not resn AC1 and elem N
set_color cls_H, [0.122, 0.373, 0.749]
show sticks, chain A and resi 626 and not resn AC1 and not name N+C+O
color cls_H, chain A and resi 626 and not resn AC1 and elem C
color red, chain A and resi 626 and not resn AC1 and elem O
color blue, chain A and resi 626 and not resn AC1 and elem N
set_color cls_H, [0.122, 0.373, 0.749]
show sticks, chain A and resi 469 and not resn AC1 and not name N+C+O
color cls_H, chain A and resi 469 and not resn AC1 and elem C
color red, chain A and resi 469 and not resn AC1 and elem O
color blue, chain A and resi 469 and not resn AC1 and elem N
set_color cls_H, [0.122, 0.373, 0.749]
show sticks, chain A and resi 568 and not resn AC1 and not name N+C+O
color cls_H, chain A and resi 568 and not resn AC1 and elem C
color red, chain A and resi 568 and not resn AC1 and elem O
color blue, chain A and resi 568 and not resn AC1 and elem N
set_color cls_Ar, [0.878, 0.482, 0.0]
show sticks, chain A and resi 565 and not resn AC1 and not name N+C+O
color cls_Ar, chain A and resi 565 and not resn AC1 and elem C
color red, chain A and resi 565 and not resn AC1 and elem O
color blue, chain A and resi 565 and not resn AC1 and elem N
set_color cls_Ar, [0.878, 0.482, 0.0]
show sticks, chain A and resi 467 and not resn AC1 and not name N+C+O
color cls_Ar, chain A and resi 467 and not resn AC1 and elem C
color red, chain A and resi 467 and not resn AC1 and elem O
color blue, chain A and resi 467 and not resn AC1 and elem N
set_color cls_Al, [0.55, 0.55, 0.55]
show sticks, chain A and resi 358 and not resn AC1 and not name N+C+O
color cls_Al, chain A and resi 358 and not resn AC1 and elem C
color red, chain A and resi 358 and not resn AC1 and elem O
color blue, chain A and resi 358 and not resn AC1 and elem N
set_color cls_Ar, [0.878, 0.482, 0.0]
show sticks, chain A and resi 601 and not resn AC1 and not name N+C+O
color cls_Ar, chain A and resi 601 and not resn AC1 and elem C
color red, chain A and resi 601 and not resn AC1 and elem O
color blue, chain A and resi 601 and not resn AC1 and elem N
set_color cls_Al, [0.55, 0.55, 0.55]
show sticks, chain A and resi 396 and not resn AC1 and not name N+C+O
color cls_Al, chain A and resi 396 and not resn AC1 and elem C
color red, chain A and resi 396 and not resn AC1 and elem O
color blue, chain A and resi 396 and not resn AC1 and elem N
distance hb0, chain A and resi 357 and name OD1 and not resn AC1, chain D and resi 7 and resn AC1 and name O6B
distance hb1, chain A and resi 357 and name OD2 and not resn AC1, chain D and resi 7 and resn AC1 and name O4
distance hb2, chain A and resi 626 and name NE2 and not resn AC1, chain D and resi 7 and resn AC1 and name O3B
distance hb3, chain A and resi 626 and name NE2 and not resn AC1, chain D and resi 7 and resn AC1 and name O4
distance hb4, chain A and resi 568 and name OD1 and not resn AC1, chain D and resi 7 and resn AC1 and name O2B
distance hb5, chain A and resi 568 and name OD2 and not resn AC1, chain D and resi 7 and resn AC1 and name O2B
hide labels
set dash_round_ends, 1
set_view (1,0,0, 0,1,0, 0,0,1, 0,0,-73.597, 0,0,0, 48.597,98.597, 1)
ray 2000,2000
png <PROJECT_ROOT>/resultados/figures/fig2_work/3WEO_ray.png, dpi=300
