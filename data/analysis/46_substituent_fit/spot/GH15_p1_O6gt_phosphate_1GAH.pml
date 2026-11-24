load <PROJECT_ROOT>/resultados/analysis/46_substituent_fit/spot/GH15_p1_O6gt_phosphate_1GAH.pdb, cx
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
distance d1, resn SUB and name O3, /cx//A/52/CH2
color black, d1
orient resn SUB or (resn MOD)
zoom resn SUB, 6
set label_size, 20
ray 1200, 900
png <PROJECT_ROOT>/resultados/analysis/46_substituent_fit/spot/GH15_p1_O6gt_phosphate_1GAH.png, dpi=150
