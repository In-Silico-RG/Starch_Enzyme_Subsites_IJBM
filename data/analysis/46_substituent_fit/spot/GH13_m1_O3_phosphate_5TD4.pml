load <PROJECT_ROOT>/resultados/analysis/46_substituent_fit/spot/GH13_m1_O3_phosphate_5TD4.pdb, cx
hide everything
set ray_opaque_background, 1
bg_color white
show sticks, resn SUB or (byres (resn SUB around 6)) or (not polymer and not resn SUB)
hide sticks, polymer and not (byres (resn SUB around 6))
show cartoon, polymer
set cartoon_transparency, 0.75
color grey70, polymer and elem C
color yellow, (not polymer) and elem C and not resn SUB
color magenta, resn SUB and elem C
color magenta, resn SUB
show spheres, resn SUB
set sphere_scale, 0.3
distance d1, resn SUB and name O4, /cx//A/62/CE2
color black, d1
orient resn SUB or (chain D and resi 3 and name O3)
zoom resn SUB, 6
set label_size, 20
ray 1200, 900
png <PROJECT_ROOT>/resultados/analysis/46_substituent_fit/spot/GH13_m1_O3_phosphate_5TD4.png, dpi=150
