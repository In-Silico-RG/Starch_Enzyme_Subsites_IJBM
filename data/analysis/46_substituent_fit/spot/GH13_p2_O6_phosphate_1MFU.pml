load <PROJECT_ROOT>/resultados/analysis/46_substituent_fit/spot/GH13_p2_O6_phosphate_1MFU.pdb, cx
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
distance d1, resn SUB and name O4, /cx//A/163/OG
color black, d1
orient resn SUB or (chain C and resi 2 and name O6)
zoom resn SUB, 6
set label_size, 20
ray 1200, 900
png <PROJECT_ROOT>/resultados/analysis/46_substituent_fit/spot/GH13_p2_O6_phosphate_1MFU.png, dpi=150
