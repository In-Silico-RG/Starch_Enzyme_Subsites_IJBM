reinitialize
load <PROJECT_ROOT>/resultados/figures/fig2_v2_work/3BC9_rot.pdb, m
remove hydro
hide everything
bg_color white
set ray_opaque_background, 1
set orthoscopic, 1
set field_of_view, 20.0
set ray_shadows, 0
set antialias, 2
set ambient, 0.5
set specular, 0.2
set spec_power, 60
set direct, 0.55
set reflect, 0.3
set depth_cue, 0
set stick_radius, 0.2
set stick_quality, 24
set sphere_quality, 3
set ray_trace_mode, 1
set ray_trace_gain, 0.4
set ray_trace_color, grey30
set dash_gap, 0.22
set dash_length, 0.18
set dash_radius, 0.055
set dash_color, grey15
set dash_round_ends, 1
set cartoon_transparency, 0.72
set cartoon_color, grey75
set cartoon_loop_radius, 0.12
set cartoon_oval_length, 1.0
set cartoon_oval_width, 0.18
set cartoon_rect_length, 1.0
set cartoon_rect_width, 0.25
set cartoon_fancy_helices, 1
set_color sugc, [0.7843137254901961, 0.7843137254901961, 0.7843137254901961]
set_color ocol, [0.8392156862745098, 0.12549019607843137, 0.12549019607843137]
set_color ncol, [0.16862745098039217, 0.23137254901960785, 0.8392156862745098]
show sticks, chain A and resi 806 and resn ACI and name C1+C2+C3+C4+C5+C6+C7+N1+O2+O3+O4+O6
set stick_radius, 0.21, chain A and resi 806 and resn ACI and name C1+C2+C3+C4+C5+C6+C7+N1+O2+O3+O4+O6
color sugc, chain A and resi 806 and resn ACI and name C1+C2+C3+C4+C5+C6+C7+N1+O2+O3+O4+O6 and elem C
color ocol, chain A and resi 806 and resn ACI and name C1+C2+C3+C4+C5+C6+C7+N1+O2+O3+O4+O6 and elem O
color ncol, chain A and resi 806 and resn ACI and name C1+C2+C3+C4+C5+C6+C7+N1+O2+O3+O4+O6 and elem N
set_color cls_Ar, [0.9019607843137255, 0.6235294117647059, 0.0]
show sticks, polymer.protein and chain A and resi 184 and not name N+C+O+OXT
color cls_Ar, polymer.protein and chain A and resi 184 and not name N+C+O+OXT and elem C
color ocol, polymer.protein and chain A and resi 184 and not name N+C+O+OXT and elem O
color ncol, polymer.protein and chain A and resi 184 and not name N+C+O+OXT and elem N+S
hide labels
set label_size, 0
set_view (1,0,0, 0,1,0, 0,0,1, 0,0,-54.735, 0,0,0, 46.735,62.735, 1)
ray 1400,1400
png <PROJECT_ROOT>/resultados/figures/fig2_v2_work/3BC9_stack_A.png, dpi=300
