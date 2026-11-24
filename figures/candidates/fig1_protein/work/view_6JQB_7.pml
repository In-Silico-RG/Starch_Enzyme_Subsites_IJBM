reinitialize
load /media/aldo/Aldo/Documents/Investigacion/Amilosa/Proyecto_F_Enzimas_Amiloliticas/manuscript_ijbm/figures/candidates/fig1_protein/work/view_6JQB_7_rot.pdb, m
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
color dom_A, prot and resi 1-113
color dom_B, prot and resi 114-167
color dom_A, prot and resi 168-361
color dom_C, prot and resi 362-418
show cartoon, prot
create surf, prot
hide cartoon, surf
show surface, surf
color surfc, surf
set transparency, 0.62, surf
set surface_mode, 1
show sticks, chain A and resi 605 and resn 7SA and name N4F+C4F+C3F+O3F+C2F+O2F+C5F+C6F+O5F+C1F+O4E
set stick_radius, 0.30, chain A and resi 605 and resn 7SA and name N4F+C4F+C3F+O3F+C2F+O2F+C5F+C6F+O5F+C1F+O4E
color sub0, chain A and resi 605 and resn 7SA and name N4F+C4F+C3F+O3F+C2F+O2F+C5F+C6F+O5F+C1F+O4E and elem C
color ocol, chain A and resi 605 and resn 7SA and name N4F+C4F+C3F+O3F+C2F+O2F+C5F+C6F+O5F+C1F+O4E and elem O
color ncol, chain A and resi 605 and resn 7SA and name N4F+C4F+C3F+O3F+C2F+O2F+C5F+C6F+O5F+C1F+O4E and elem N
show sticks, chain A and resi 605 and resn 7SA and name N4C+C1D+C7D+C5D+C6D+O6D+C4D+O4D+C3D+O3D+C2D+O2D
set stick_radius, 0.30, chain A and resi 605 and resn 7SA and name N4C+C1D+C7D+C5D+C6D+O6D+C4D+O4D+C3D+O3D+C2D+O2D
color sub2, chain A and resi 605 and resn 7SA and name N4C+C1D+C7D+C5D+C6D+O6D+C4D+O4D+C3D+O3D+C2D+O2D and elem C
color ocol, chain A and resi 605 and resn 7SA and name N4C+C1D+C7D+C5D+C6D+O6D+C4D+O4D+C3D+O3D+C2D+O2D and elem O
color ncol, chain A and resi 605 and resn 7SA and name N4C+C1D+C7D+C5D+C6D+O6D+C4D+O4D+C3D+O3D+C2D+O2D and elem N
show sticks, chain A and resi 605 and resn 7SA and name O4B+C1C+O5C+C5C+C6C+C4C+C3C+O3C+C2C+O2C+N4C
set stick_radius, 0.30, chain A and resi 605 and resn 7SA and name O4B+C1C+O5C+C5C+C6C+C4C+C3C+O3C+C2C+O2C+N4C
color sub3, chain A and resi 605 and resn 7SA and name O4B+C1C+O5C+C5C+C6C+C4C+C3C+O3C+C2C+O2C+N4C and elem C
color ocol, chain A and resi 605 and resn 7SA and name O4B+C1C+O5C+C5C+C6C+C4C+C3C+O3C+C2C+O2C+N4C and elem O
color ncol, chain A and resi 605 and resn 7SA and name O4B+C1C+O5C+C5C+C6C+C4C+C3C+O3C+C2C+O2C+N4C and elem N
show sticks, chain A and resi 605 and resn 7SA and name O6G+C6G+C5G+C4G+O4G+C3G+O3G+C2G+O2G+C7G+C1G+N4F
set stick_radius, 0.30, chain A and resi 605 and resn 7SA and name O6G+C6G+C5G+C4G+O4G+C3G+O3G+C2G+O2G+C7G+C1G+N4F
color outc, chain A and resi 605 and resn 7SA and name O6G+C6G+C5G+C4G+O4G+C3G+O3G+C2G+O2G+C7G+C1G+N4F and elem C
color ocol, chain A and resi 605 and resn 7SA and name O6G+C6G+C5G+C4G+O4G+C3G+O3G+C2G+O2G+C7G+C1G+N4F and elem O
color ncol, chain A and resi 605 and resn 7SA and name O6G+C6G+C5G+C4G+O4G+C3G+O3G+C2G+O2G+C7G+C1G+N4F and elem N
show sticks, chain A and resi 605 and resn 7SA and name O6A+C6A+C5A+O5A+C1A+O1A+C4A+C3A+O3A+C2A+O2A+O4A
set stick_radius, 0.30, chain A and resi 605 and resn 7SA and name O6A+C6A+C5A+O5A+C1A+O1A+C4A+C3A+O3A+C2A+O2A+O4A
color sub5, chain A and resi 605 and resn 7SA and name O6A+C6A+C5A+O5A+C1A+O1A+C4A+C3A+O3A+C2A+O2A+O4A and elem C
color ocol, chain A and resi 605 and resn 7SA and name O6A+C6A+C5A+O5A+C1A+O1A+C4A+C3A+O3A+C2A+O2A+O4A and elem O
color ncol, chain A and resi 605 and resn 7SA and name O6A+C6A+C5A+O5A+C1A+O1A+C4A+C3A+O3A+C2A+O2A+O4A and elem N
show sticks, chain A and resi 605 and resn 7SA and name O4A+C1B+O5B+C5B+C6B+O6B+C4B+C3B+O3B+C2B+O2B+O4B
set stick_radius, 0.30, chain A and resi 605 and resn 7SA and name O4A+C1B+O5B+C5B+C6B+O6B+C4B+C3B+O3B+C2B+O2B+O4B
color sub4, chain A and resi 605 and resn 7SA and name O4A+C1B+O5B+C5B+C6B+O6B+C4B+C3B+O3B+C2B+O2B+O4B and elem C
color ocol, chain A and resi 605 and resn 7SA and name O4A+C1B+O5B+C5B+C6B+O6B+C4B+C3B+O3B+C2B+O2B+O4B and elem O
color ncol, chain A and resi 605 and resn 7SA and name O4A+C1B+O5B+C5B+C6B+O6B+C4B+C3B+O3B+C2B+O2B+O4B and elem N
show sticks, chain A and resi 605 and resn 7SA and name O4D+O4E+C4E+C3E+O3E+C2E+O2E+C5E+C6E+O6E+O5E+C1E
set stick_radius, 0.30, chain A and resi 605 and resn 7SA and name O4D+O4E+C4E+C3E+O3E+C2E+O2E+C5E+C6E+O6E+O5E+C1E
color sub1, chain A and resi 605 and resn 7SA and name O4D+O4E+C4E+C3E+O3E+C2E+O2E+C5E+C6E+O6E+O5E+C1E and elem C
color ocol, chain A and resi 605 and resn 7SA and name O4D+O4E+C4E+C3E+O3E+C2E+O2E+C5E+C6E+O6E+O5E+C1E and elem O
color ncol, chain A and resi 605 and resn 7SA and name O4D+O4E+C4E+C3E+O3E+C2E+O2E+C5E+C6E+O6E+O5E+C1E and elem N
show sticks, prot and resi 193 and not name N+C+O
set stick_radius, 0.30, prot and resi 193 and not name N+C+O
color catc, prot and resi 193 and not name N+C+O and elem C
color ocol, prot and resi 193 and not name N+C+O and elem O
color ncol, prot and resi 193 and not name N+C+O and elem N
show sticks, prot and resi 219 and not name N+C+O
set stick_radius, 0.30, prot and resi 219 and not name N+C+O
color catc, prot and resi 219 and not name N+C+O and elem C
color ocol, prot and resi 219 and not name N+C+O and elem O
color ncol, prot and resi 219 and not name N+C+O and elem N
show sticks, prot and resi 294 and not name N+C+O
set stick_radius, 0.30, prot and resi 294 and not name N+C+O
color catc, prot and resi 294 and not name N+C+O and elem C
color ocol, prot and resi 294 and not name N+C+O and elem O
color ncol, prot and resi 294 and not name N+C+O and elem N
set_view (1,0,0, 0,1,0, 0,0,1, 0,0,-138.588, 0,0,0, 106.179,170.997, 1)
ray 700,449
png /media/aldo/Aldo/Documents/Investigacion/Amilosa/Proyecto_F_Enzimas_Amiloliticas/manuscript_ijbm/figures/candidates/fig1_protein/work/view_6JQB_7.png, dpi=300
