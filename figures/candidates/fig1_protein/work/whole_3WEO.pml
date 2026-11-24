reinitialize
load /media/aldo/Aldo/Documents/Investigacion/Amilosa/Proyecto_F_Enzimas_Amiloliticas/manuscript_ijbm/figures/candidates/fig1_protein/work/whole_3WEO_rot.pdb, m
remove hydro
remove resn HOH
hide everything
bg_color white
set ray_opaque_background, 1
set orthoscopic, 1
set field_of_view, 20.0
set ray_shadows, 0
set antialias, 2
set ambient, 0.45
set direct, 0.45
set reflect, 0.35
set specular, 0.15
set spec_power, 80
set depth_cue, 0
set ray_trace_mode, 1
set ray_trace_gain, 0.08
set ray_trace_color, black
set stick_quality, 24
set valence, 0
set cartoon_fancy_helices, 1
set cartoon_loop_radius, 0.25
set cartoon_oval_length, 1.2
set cartoon_rect_length, 1.2
set surface_quality, 1
set transparency_mode, 1
set_color ocol, [0.8392, 0.1255, 0.1255]
set_color ncol, [0.1686, 0.2314, 0.8392]
set_color catc, [0.1882, 0.1882, 0.1882]
set_color outc, [0.6039, 0.6039, 0.6039]
set_color surfc, [0.902, 0.902, 0.902]
set_color sub0, [0.3373, 0.7059, 0.9137]
set_color sub1, [0.0, 0.4471, 0.698]
set_color sub2, [0.0, 0.6196, 0.451]
set_color sub3, [0.902, 0.6235, 0.0]
set_color sub4, [0.8353, 0.3686, 0.0]
set_color sub5, [0.8, 0.4745, 0.6549]
set_color dom_A, [0.6863, 0.7647, 0.8549]
set_color dom_B, [0.8902, 0.8118, 0.651]
set_color dom_C, [0.7098, 0.8392, 0.6824]
set_color dom_N, [0.851, 0.7216, 0.8118]
set_color dom_C2, [0.8118, 0.8118, 0.8118]
set cartoon_transparency, 0.0
select prot, polymer.protein and chain A
color dom_N, prot and resi 38-313
color dom_A, prot and resi 314-655
color dom_C, prot and resi 656-756
color dom_C2, prot and resi 757-909
show cartoon, prot
create surf, prot
hide cartoon, surf
show surface, surf
color surfc, surf
set transparency, 0.62, surf
set surface_mode, 1
show sticks, chain D and resi 7 and resn AC1 and name C1+O2+C2+C3+O3+C4+N4A+C5+O5+C6
set stick_radius, 0.30, chain D and resi 7 and resn AC1 and name C1+O2+C2+C3+O3+C4+N4A+C5+O5+C6
color sub3, chain D and resi 7 and resn AC1 and name C1+O2+C2+C3+O3+C4+N4A+C5+O5+C6 and elem C
color ocol, chain D and resi 7 and resn AC1 and name C1+O2+C2+C3+O3+C4+N4A+C5+O5+C6 and elem O
color ncol, chain D and resi 7 and resn AC1 and name C1+O2+C2+C3+O3+C4+N4A+C5+O5+C6 and elem N
show sticks, chain D and resi 2 and resn GLC and name C1+C2+C3+C4+C5+C6+O2+O3+O4+O5+O6
set stick_radius, 0.30, chain D and resi 2 and resn GLC and name C1+C2+C3+C4+C5+C6+O2+O3+O4+O5+O6
color outc, chain D and resi 2 and resn GLC and name C1+C2+C3+C4+C5+C6+O2+O3+O4+O5+O6 and elem C
color ocol, chain D and resi 2 and resn GLC and name C1+C2+C3+C4+C5+C6+O2+O3+O4+O5+O6 and elem O
color ncol, chain D and resi 2 and resn GLC and name C1+C2+C3+C4+C5+C6+O2+O3+O4+O5+O6 and elem N
show sticks, chain D and resi 5 and resn GLC and name C1+C2+C3+C4+C5+C6+O2+O3+O4+O5+O6
set stick_radius, 0.30, chain D and resi 5 and resn GLC and name C1+C2+C3+C4+C5+C6+O2+O3+O4+O5+O6
color sub5, chain D and resi 5 and resn GLC and name C1+C2+C3+C4+C5+C6+O2+O3+O4+O5+O6 and elem C
color ocol, chain D and resi 5 and resn GLC and name C1+C2+C3+C4+C5+C6+O2+O3+O4+O5+O6 and elem O
color ncol, chain D and resi 5 and resn GLC and name C1+C2+C3+C4+C5+C6+O2+O3+O4+O5+O6 and elem N
show sticks, chain D and resi 3 and resn GLC and name C1+C2+C3+C4+C5+C6+O2+O3+O4+O5+O6
set stick_radius, 0.30, chain D and resi 3 and resn GLC and name C1+C2+C3+C4+C5+C6+O2+O3+O4+O5+O6
color outc, chain D and resi 3 and resn GLC and name C1+C2+C3+C4+C5+C6+O2+O3+O4+O5+O6 and elem C
color ocol, chain D and resi 3 and resn GLC and name C1+C2+C3+C4+C5+C6+O2+O3+O4+O5+O6 and elem O
color ncol, chain D and resi 3 and resn GLC and name C1+C2+C3+C4+C5+C6+O2+O3+O4+O5+O6 and elem N
show sticks, chain D and resi 1 and resn GLC and name C1+C2+C3+C4+C5+C6+O1+O2+O3+O4+O5+O6
set stick_radius, 0.30, chain D and resi 1 and resn GLC and name C1+C2+C3+C4+C5+C6+O1+O2+O3+O4+O5+O6
color outc, chain D and resi 1 and resn GLC and name C1+C2+C3+C4+C5+C6+O1+O2+O3+O4+O5+O6 and elem C
color ocol, chain D and resi 1 and resn GLC and name C1+C2+C3+C4+C5+C6+O1+O2+O3+O4+O5+O6 and elem O
color ncol, chain D and resi 1 and resn GLC and name C1+C2+C3+C4+C5+C6+O1+O2+O3+O4+O5+O6 and elem N
show sticks, chain D and resi 4 and resn GLC and name C1+C2+C3+C4+C5+C6+O2+O3+O4+O5+O6
set stick_radius, 0.30, chain D and resi 4 and resn GLC and name C1+C2+C3+C4+C5+C6+O2+O3+O4+O5+O6
color outc, chain D and resi 4 and resn GLC and name C1+C2+C3+C4+C5+C6+O2+O3+O4+O5+O6 and elem C
color ocol, chain D and resi 4 and resn GLC and name C1+C2+C3+C4+C5+C6+O2+O3+O4+O5+O6 and elem O
color ncol, chain D and resi 4 and resn GLC and name C1+C2+C3+C4+C5+C6+O2+O3+O4+O5+O6 and elem N
show sticks, chain D and resi 7 and resn AC1 and name C4A+N4A+C1B+C2B+O2B+C3B+O3B+O4+C5B+C7B+C6B+O6B
set stick_radius, 0.30, chain D and resi 7 and resn AC1 and name C4A+N4A+C1B+C2B+O2B+C3B+O3B+O4+C5B+C7B+C6B+O6B
color sub2, chain D and resi 7 and resn AC1 and name C4A+N4A+C1B+C2B+O2B+C3B+O3B+O4+C5B+C7B+C6B+O6B and elem C
color ocol, chain D and resi 7 and resn AC1 and name C4A+N4A+C1B+C2B+O2B+C3B+O3B+O4+C5B+C7B+C6B+O6B and elem O
color ncol, chain D and resi 7 and resn AC1 and name C4A+N4A+C1B+C2B+O2B+C3B+O3B+O4+C5B+C7B+C6B+O6B and elem N
show sticks, chain D and resi 6 and resn GLC and name C1+C2+C3+C4+C5+C6+O2+O3+O4+O5+O6
set stick_radius, 0.30, chain D and resi 6 and resn GLC and name C1+C2+C3+C4+C5+C6+O2+O3+O4+O5+O6
color sub4, chain D and resi 6 and resn GLC and name C1+C2+C3+C4+C5+C6+O2+O3+O4+O5+O6 and elem C
color ocol, chain D and resi 6 and resn GLC and name C1+C2+C3+C4+C5+C6+O2+O3+O4+O5+O6 and elem O
color ncol, chain D and resi 6 and resn GLC and name C1+C2+C3+C4+C5+C6+O2+O3+O4+O5+O6 and elem N
show sticks, prot and resi 469 and not name N+C+O
set stick_radius, 0.30, prot and resi 469 and not name N+C+O
color catc, prot and resi 469 and not name N+C+O and elem C
color ocol, prot and resi 469 and not name N+C+O and elem O
color ncol, prot and resi 469 and not name N+C+O and elem N
show sticks, prot and resi 568 and not name N+C+O
set stick_radius, 0.30, prot and resi 568 and not name N+C+O
color catc, prot and resi 568 and not name N+C+O and elem C
color ocol, prot and resi 568 and not name N+C+O and elem O
color ncol, prot and resi 568 and not name N+C+O and elem N
set_view (1,0,0, 0,1,0, 0,0,1, 0,0,-234.571, 0,0,0, 185.783,283.359, 1)
ray 3000,3437
png /media/aldo/Aldo/Documents/Investigacion/Amilosa/Proyecto_F_Enzimas_Amiloliticas/manuscript_ijbm/figures/candidates/fig1_protein/work/whole_3WEO.png, dpi=300
