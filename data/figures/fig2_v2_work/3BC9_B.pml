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
set stick_radius, 0.17
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
show cartoon, polymer.protein
show sticks, chain A and resi 806 and resn ACI and name C1+C2+C3+C4+C5+C6+C7+N1+O2+O3+O4+O6
set stick_radius, 0.21, chain A and resi 806 and resn ACI and name C1+C2+C3+C4+C5+C6+C7+N1+O2+O3+O4+O6
color sugc, chain A and resi 806 and resn ACI and name C1+C2+C3+C4+C5+C6+C7+N1+O2+O3+O4+O6 and elem C
color ocol, chain A and resi 806 and resn ACI and name C1+C2+C3+C4+C5+C6+C7+N1+O2+O3+O4+O6 and elem O
color ncol, chain A and resi 806 and resn ACI and name C1+C2+C3+C4+C5+C6+C7+N1+O2+O3+O4+O6 and elem N
show sticks, chain D and resi 3 and resn GLD and name C1+C2+C3+C4+C5+C6+O2+O3+O5
set stick_radius, 0.21, chain D and resi 3 and resn GLD and name C1+C2+C3+C4+C5+C6+O2+O3+O5
color sugc, chain D and resi 3 and resn GLD and name C1+C2+C3+C4+C5+C6+O2+O3+O5 and elem C
color ocol, chain D and resi 3 and resn GLD and name C1+C2+C3+C4+C5+C6+O2+O3+O5 and elem O
color ncol, chain D and resi 3 and resn GLD and name C1+C2+C3+C4+C5+C6+O2+O3+O5 and elem N
show sticks, chain D and resi 2 and resn GLC and name C1+C2+C3+C4+C5+C6+O2+O3+O4+O5+O6
set stick_radius, 0.21, chain D and resi 2 and resn GLC and name C1+C2+C3+C4+C5+C6+O2+O3+O4+O5+O6
color sugc, chain D and resi 2 and resn GLC and name C1+C2+C3+C4+C5+C6+O2+O3+O4+O5+O6 and elem C
color ocol, chain D and resi 2 and resn GLC and name C1+C2+C3+C4+C5+C6+O2+O3+O4+O5+O6 and elem O
color ncol, chain D and resi 2 and resn GLC and name C1+C2+C3+C4+C5+C6+O2+O3+O4+O5+O6 and elem N
set_color cls_H, [0.0, 0.4470588235294118, 0.6980392156862745]
show sticks, polymer.protein and chain A and resi 353 and not name N+C+O+OXT
color cls_H, polymer.protein and chain A and resi 353 and not name N+C+O+OXT and elem C
color ocol, polymer.protein and chain A and resi 353 and not name N+C+O+OXT and elem O
color ncol, polymer.protein and chain A and resi 353 and not name N+C+O+OXT and elem N+S
set_color cls_Al, [0.8, 0.4745098039215686, 0.6549019607843137]
show sticks, polymer.protein and chain A and resi 382 and not name N+C+O+OXT
color cls_Al, polymer.protein and chain A and resi 382 and not name N+C+O+OXT and elem C
color ocol, polymer.protein and chain A and resi 382 and not name N+C+O+OXT and elem O
color ncol, polymer.protein and chain A and resi 382 and not name N+C+O+OXT and elem N+S
set_color cls_H, [0.0, 0.4470588235294118, 0.6980392156862745]
show sticks, polymer.protein and chain A and resi 380 and not name N+C+O+OXT
color cls_H, polymer.protein and chain A and resi 380 and not name N+C+O+OXT and elem C
color ocol, polymer.protein and chain A and resi 380 and not name N+C+O+OXT and elem O
color ncol, polymer.protein and chain A and resi 380 and not name N+C+O+OXT and elem N+S
set_color cls_Al, [0.8, 0.4745098039215686, 0.6549019607843137]
show sticks, polymer.protein and chain A and resi 351 and not name N+C+O+OXT
color cls_Al, polymer.protein and chain A and resi 351 and not name N+C+O+OXT and elem C
color ocol, polymer.protein and chain A and resi 351 and not name N+C+O+OXT and elem O
color ncol, polymer.protein and chain A and resi 351 and not name N+C+O+OXT and elem N+S
set_color cls_H, [0.0, 0.4470588235294118, 0.6980392156862745]
show sticks, polymer.protein and chain A and resi 350 and not name N+C+O+OXT
color cls_H, polymer.protein and chain A and resi 350 and not name N+C+O+OXT and elem C
color ocol, polymer.protein and chain A and resi 350 and not name N+C+O+OXT and elem O
color ncol, polymer.protein and chain A and resi 350 and not name N+C+O+OXT and elem N+S
set_color cls_H, [0.0, 0.4470588235294118, 0.6980392156862745]
show sticks, polymer.protein and chain A and resi 447 and not name N+C+O+OXT
color cls_H, polymer.protein and chain A and resi 447 and not name N+C+O+OXT and elem C
color ocol, polymer.protein and chain A and resi 447 and not name N+C+O+OXT and elem O
color ncol, polymer.protein and chain A and resi 447 and not name N+C+O+OXT and elem N+S
set_color cls_H, [0.0, 0.4470588235294118, 0.6980392156862745]
show sticks, polymer.protein and chain A and resi 446 and not name N+C+O+OXT
color cls_H, polymer.protein and chain A and resi 446 and not name N+C+O+OXT and elem C
color ocol, polymer.protein and chain A and resi 446 and not name N+C+O+OXT and elem O
color ncol, polymer.protein and chain A and resi 446 and not name N+C+O+OXT and elem N+S
set_color cls_H, [0.0, 0.4470588235294118, 0.6980392156862745]
show sticks, polymer.protein and chain A and resi 348 and not name N+C+O+OXT
color cls_H, polymer.protein and chain A and resi 348 and not name N+C+O+OXT and elem C
color ocol, polymer.protein and chain A and resi 348 and not name N+C+O+OXT and elem O
color ncol, polymer.protein and chain A and resi 348 and not name N+C+O+OXT and elem N+S
set_color cls_H, [0.0, 0.4470588235294118, 0.6980392156862745]
show sticks, polymer.protein and chain A and resi 233 and not name N+C+O+OXT
color cls_H, polymer.protein and chain A and resi 233 and not name N+C+O+OXT and elem C
color ocol, polymer.protein and chain A and resi 233 and not name N+C+O+OXT and elem O
color ncol, polymer.protein and chain A and resi 233 and not name N+C+O+OXT and elem N+S
set_color cls_Ar, [0.9019607843137255, 0.6235294117647059, 0.0]
show sticks, polymer.protein and chain A and resi 184 and not name N+C+O+OXT
color cls_Ar, polymer.protein and chain A and resi 184 and not name N+C+O+OXT and elem C
color ocol, polymer.protein and chain A and resi 184 and not name N+C+O+OXT and elem O
color ncol, polymer.protein and chain A and resi 184 and not name N+C+O+OXT and elem N+S
distance hb0, polymer.protein and chain A and resi 350 and name OD1, chain A and resi 806 and resn ACI and name O2
distance hb1, polymer.protein and chain A and resi 350 and name OD2, chain A and resi 806 and resn ACI and name O6
distance hb2, polymer.protein and chain A and resi 447 and name OD1, chain A and resi 806 and resn ACI and name O3
distance hb3, polymer.protein and chain A and resi 447 and name OD2, chain A and resi 806 and resn ACI and name O2
distance hb4, polymer.protein and chain A and resi 380 and name OE1, chain D and resi 3 and resn GLD and name O3
distance hb5, polymer.protein and chain A and resi 380 and name OE2, chain A and resi 806 and resn ACI and name N1
distance hb6, polymer.protein and chain A and resi 446 and name NE2, chain A and resi 806 and resn ACI and name O2
distance hb7, polymer.protein and chain A and resi 353 and name NZ, chain D and resi 2 and resn GLC and name O3
distance hb8, polymer.protein and chain A and resi 348 and name NH2, chain A and resi 806 and resn ACI and name O2
hide labels
set label_size, 0
set_view (1,0,0, 0,1,0, 0,0,1, 0,0,-67.746, 0,0,0, 56.746,78.746, 1)
ray 2400,2400
png <PROJECT_ROOT>/resultados/figures/fig2_v2_work/3BC9_B.png, dpi=300
