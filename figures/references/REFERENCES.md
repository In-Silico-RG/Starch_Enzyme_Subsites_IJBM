# Visual references for the active-site and cleft figures (GH13 / GH31)

Collected 2026-11-24. Purpose: set the visual style of two figures of the IJBM draft
(α-amylase GH13, α-glucosidase GH31) and of the graphical abstract. No manuscript was edited.

**How the images were obtained.** All article figures come from the PMC Open Access subset,
downloaded from the NCBI open-data bucket `https://pmc-oa-opendata.s3.amazonaws.com/<PMCID>.<ver>/<file>`.
(The PMC and Europe PMC web pages return 403/HTML to scripted requests, and ACS returns 410.)
The PDB-101 images come from `cdn.rcsb.org`. Each entry below gives the licence. Files under
CC BY-NC-ND or ACS AuthorChoice are kept **only as private style references**. Do not reuse
them in our figures; we imitate the style, not the image.

Tools installed locally (checked 2026-11-24):
- **UCSF ChimeraX 1.11.1**: `/usr/bin/chimerax`, deb `ucsf-chimerax`. Offscreen rendering works (`chimerax --nogui --offscreen`).
- **PyMOL open-source 2.5.0**: `/usr/bin/pymol`, deb. Also **PyMOL 3.1.0** as the flatpak `org.pymol.PyMOL`.
- **VMD**: `/usr/local/bin/vmd`.
- **gfortran**: `/usr/bin/gfortran`. With it, Goodsell's Illustrate (Fortran source) can be compiled if wanted.
- **Not installed:** Blender and Molecular Nodes, Illustrate, UCSF Chimera (legacy).

---

## A. Catalogue of examples

Column codes. **Cat.** is the category: 1 = zoomed active site, 2 = cleft overview with surface,
3 = graphical abstract, 4 = 2D interaction diagram, X = counter-example.
**Renderer** is the one stated in the paper's Methods; "PyMOL look" means my reading of the image
when the paper does not say.

### 01. Human lysosomal acid α-glucosidase (GAA, GH31) with DNJ, NHE-DNJ and acarbose
- File: `01_activesite_GAA_acarbose_NatCommun2017_Fig3.jpg`
- Source: Roig-Zamboni V. et al., *Nat. Commun.* 8, 1111 (2017), doi:10.1038/s41467-017-01263-3. Fig. 3. https://pmc.ncbi.nlm.nih.gov/articles/PMC5653652/
- Image URL: https://pmc-oa-opendata.s3.amazonaws.com/PMC5653652.1/41467_2017_1263_Fig3_HTML.jpg
- Licence: **CC BY 4.0**
- Cat. 1 and 4-like
- Shows paired panels for each ligand:
  - left, a "clean" view: no cartoon, only about 12 side chains plus the ligand;
  - right, the same site with the cartoon and a 2Fo–Fc mesh;
  - a small whole-protein thumbnail with a box marking the zoom;
  - a ChemDraw structure of each ligand.
- Style:
  - Renderer: PyMOL (stated).
  - Thick sticks (about 0.25 Å).
  - Ligand carbons orange. Residue carbons in the colour of their domain (olive green). The catalytic Asp518 and Asp616 picked out in magenta.
  - Waters as small red spheres. H-bonds as black dashes.
  - Labels in sans-serif, coloured like the residue.
  - White background. No outline, no ambient occlusion (AO), no depth of field (DOF).
  - Subsites −1 to +3 written on the sugar rings.
- Why it matters: the closest published analogue to our GH31 figure. Its best idea is the "thumbnail plus clean zoom plus 2D formula" layout.

### 02. GAA: isomaltose model, pocket surface and allosteric sites
- File: `02_activesite_surface_GAA_isomaltose_NatCommun2017_Fig4.jpg`
- Source: as 01, Fig. 4. Licence **CC BY 4.0**.
- Cat. 1 and 2
- Panel b is the useful one: an opaque surface coloured by domain (olive, magenta, slate, yellow), slightly transparent so the side chains show through. The disaccharide sits in sticks, with cyan/teal carbons, inside the pocket.
- Style:
  - PyMOL surface with `transparency` of about 0.2.
  - Ligand sticks thicker than the residue lines.
  - No labels inside the pocket.

### 03. GAA with α-Glc-cyclosulfamidates (Davies and Overkleeft labs)
- File: `03_activesite_GAA_cyclosulfamidate_JACS2022_Fig3.jpg`
- Source: Kok K., ... Artola M. *J. Am. Chem. Soc.* 144, 14819 (2022), doi:10.1021/jacs.2c05666. Fig. 3. https://pmc.ncbi.nlm.nih.gov/articles/PMC9389588/
- Image URL: https://pmc-oa-opendata.s3.amazonaws.com/PMC9389588.1/ja2c05666_0004.jpg
- Licence: **CC BY 4.0**
- Cat. 1
- Shows the "York look", the house style of Davies's group:
  - Pale salmon cartoon at about 50 % transparency as background only.
  - Residue carbons in salmon/pink sticks, slightly darker than the cartoon.
  - Ligand carbons in slate blue (left panel) or grey (right panel).
  - The density mesh is green for the ligand and blue for the residues.
  - A reference ligand from another PDB entry, drawn as thin grey lines.
  - Six to seven labelled residues. Catalytic roles in parentheses: "Asp518 (nuc.)", "Asp616 (a./b.)".
  - Black sans-serif labels. Thin grey dashed H-bonds. White background.
- Style:
  - CCP4mg/PyMOL look.
  - Thin sticks of about 0.15–0.2 Å.
  - Very low clutter. No outline.

### 04. CjAgd31B (GH31) with cyclophellitol cyclosulfate: unreacted versus covalent
- File: `04_activesite_CjAgd31B_ACSCentSci2017_Fig4.jpg`
- Source: Artola M. et al., *ACS Cent. Sci.* 3, 784 (2017), doi:10.1021/acscentsci.7b00214. Fig. 4. https://pmc.ncbi.nlm.nih.gov/articles/PMC5532717/
- Image URL: https://pmc-oa-opendata.s3.amazonaws.com/PMC5532717.1/oc-2017-002148_0004.jpg
- Licence: **ACS AuthorChoice** (free to read; non-commercial reuse with attribution)
- Cat. 1
- Shows the minimum active site: only the ligand plus the nucleophile and the acid/base, and nothing else. Pale mint cartoon, very transparent. Ligand and the two Asp in lilac-blue sticks. Each panel is captioned with the ring conformation (⁴C₁, ¹S₃).
- Style: three residues in all; the "less is more" extreme. Arrows lead from the unreacted to the reacted state.
- Why it matters: the model for a GH31 "catalytic pair" panel.

### 05. Taka-amylase A (GH13) with glycosylated cyclophellitols
- File: `05_activesite_TakaAmylase_JACS2021_Fig3.jpg`
- Source: Chen Y., Armstrong Z., ... Davies G. J., Overkleeft H. S. *J. Am. Chem. Soc.* 143, 2423 (2021), doi:10.1021/jacs.0c13059. Fig. 3. https://pmc.ncbi.nlm.nih.gov/articles/PMC7883350/
- Image URL: https://pmc-oa-opendata.s3.amazonaws.com/PMC7883350.1/ja0c13059_0003.jpg
- Licence: **CC BY-NC-ND 4.0** (style reference only)
- Cat. 1
- Shows a GH13 −1/+1 site in a 2 × 2 grid:
  - Trp83, His122, Asp206 (nucleophile), Glu230 (acid/base), Asp297, His296, Arg344, Asp340, Thr207, about 9 residues in all.
  - Green carbons for both protein and ligand.
  - Transparent green cartoon behind. Blue 2Fo–Fc mesh. Thick black dashes. Large black labels.
- Style: PyMOL look; thick sticks (about 0.25 Å); a single hue for everything.
- Why it matters: the direct GH13 analogue. The single hue is weak, because ligand and residues blur together; the better choice is to give the ligand its own colour.

### 06. Acarbose-preferred glucosidase Apg: surface-coloured cleft and subsite labels
- File: `06_cleft_activesite_Apg_acarbose_NatCommun2025_Fig3.jpg`
- Source: Huang et al., *Nat. Commun.* 16, 7839 (2025), doi:10.1038/s41467-025-62855-y. Fig. 3. https://pmc.ncbi.nlm.nih.gov/articles/PMC12373829/
- Image URL: https://pmc-oa-opendata.s3.amazonaws.com/PMC12373829.1/41467_2025_62855_Fig3_HTML.jpg
- Licence: **CC BY-NC-ND 4.0** (style reference only)
- Cat. 1 and 2
- Shows:
  - Panel A is a "clip diagram": the opaque surface is cut by the clipping plane, white and blue, and the acarbose lies in the lumen as cyan ball-and-stick with a grey mesh. Subsite numbers −3 … +1 sit beside each ring.
  - Panels B–D are the active site: pink residue carbons, catalytic residues in tan with orange labels (D448, E373, D336), and water spheres.
- Style:
  - Renderer: PyMOL (stated).
  - Ball-and-stick for the sugar, sticks for the residues.
  - Bold labels. Subsite numbers in black. Catalytic residues highlighted by label colour.
- Why it matters: "surface clipped through the cleft plus subsite numbers" is exactly what a starch-enzyme figure needs.

### 07. *R. bromii* amylosome GH13 enzymes: transparent surface over cartoon, ChimeraX
- File: `07_cleft_transparentsurface_amylosome_NatCommun2025_Fig3.jpg`
- Source: Wimmer et al., *Nat. Commun.* 16, 10763 (2025), doi:10.1038/s41467-025-65800-1. Fig. 3. https://pmc.ncbi.nlm.nih.gov/articles/PMC12663171/
- Image URL: https://pmc-oa-opendata.s3.amazonaws.com/PMC12663171.1/41467_2025_65800_Fig3_HTML.jpg
- Licence: **CC BY-NC-ND 4.0** (style reference only)
- Cat. 2 and 1
- Shows:
  - Panel d: five whole enzymes. Each is a pastel cartoon (teal for the GH13 AmyA domain, sand for the PulA domain, lilac for the CBMs) under a transparent surface of the same colour. Dashed boxes mark the zooms.
  - Panels e–g: the cleft with oligosaccharides docked or soaked in. The sugars are ball-and-stick: white carbons, red oxygens, small balls. Three to five labelled residues; the catalytic carbons are darker teal or brown. A dotted circle marks the −1 subsite.
- Style:
  - Renderer: ChimeraX (stated).
  - `lighting soft`-like ambient shading, gentle and with no hard shadows.
  - Surface transparency about 60–70 %. No silhouettes.
  - Thin labels in sans-serif. White background.
- Why it matters: the closest "modern ChimeraX" look for a GH13 cleft. **Recipe A imitates it.**

### 08. Human glycogen debranching enzyme (GH13 GT domain plus GH133): coloured surfaces, LigPlot+
- File: `08_cleft_surfaces_ligplot_hsGDE_NatCommun2025_Fig2.jpg`
- Source: Guan et al., *Nat. Commun.* 16, 5962 (2025), doi:10.1038/s41467-025-61077-6. Fig. 2. https://pmc.ncbi.nlm.nih.gov/articles/PMC12219210/
- Image URL: https://pmc-oa-opendata.s3.amazonaws.com/PMC12219210.1/41467_2025_61077_Fig2_HTML.jpg
- Licence: **CC BY-NC-ND 4.0** (style reference only)
- Cat. 2, 1 and **4**
- Shows:
  - e: an electrostatic surface (red to white to blue) with maltopentaose in green ball-and-stick down the cleft.
  - f: a hydrophobicity surface (cyan to white to gold, ChimeraX `mlp`) with the same sugar.
  - g: sticks over a cartoon, with Glc1–Glc5 labelled in blue.
  - h: a **LigPlot+** 2D diagram (green dashed H-bonds with distances; red "eyelashes" for hydrophobic contacts).
- Style:
  - Renderer: ChimeraX (stated), plus LigPlot+.
  - Colour-bar legends under the surfaces.
- Why it matters: shows how one sugar can be told in three surface colourings plus one 2D diagram. It is the template for category 4.

### 09. hsGDE catalytic residues against GH13 α-amylase and glucodextranase
- File: `09_catalytic_residues_hsGDE_NatCommun2025_Fig3.jpg`
- Source: as 08, Fig. 3. Licence **CC BY-NC-ND 4.0** (style reference only).
- Cat. 1 (superposition)
- Shows a whole-protein thumbnail with boxes, then six square insets. Each inset has a very transparent pastel cartoon, two or three catalytic acids in sticks, and two-colour labels (one per structure) for the superposed residues.
- Style: ChimeraX; thin sticks; strong cartoon transparency; square panels with thin black frames.
- Why it matters: the best published layout for comparing catalytic residues between GH13 and GH31 or another family.

### 10. Branching enzyme (GH13): rainbow surface with a long maltooligosaccharide over subsites −7…+5
- File: `10_cleft_subsites_branchingenzyme_NatCommun2024_Fig4.jpg`
- Source: Zong et al., *Nat. Commun.* 15, 8760 (2024), doi:10.1038/s41467-024-53018-6. Fig. 4. https://pmc.ncbi.nlm.nih.gov/articles/PMC11464733/
- Image URL: https://pmc-oa-opendata.s3.amazonaws.com/PMC11464733.1/41467_2024_53018_Fig4_HTML.jpg
- Licence: **CC BY-NC-ND 4.0** (style reference only)
- Cat. 2 and 1
- Shows:
  - a: an opaque surface in a rainbow from N to C terminus, with two superposed glucan chains in sticks (magenta, green) and subsite numbers −7 … +5.
  - b: the active site with H atoms shown, cyan residue carbons, green dashed H-bonds.
- Style: PyMOL (stated). The rainbow surface is loud. The subsite numbering along a long chain is the idea worth borrowing; the rainbow is not.

### 11. GH31 α-1,3-glucosidase LlGH31_u1 with nigero-oligosaccharides
- File: `11_activesite_GH31_nigerose_JBC2022_Fig6.jpg`
- Source: Ikegaya M. et al., *J. Biol. Chem.* 298, 101827 (2022), doi:10.1016/j.jbc.2022.101827. Fig. 6. https://pmc.ncbi.nlm.nih.gov/articles/PMC9061262/
- Image URL: https://pmc-oa-opendata.s3.amazonaws.com/PMC9061262.1/gr6.jpg
- Licence: **CC BY 4.0**
- Cat. 1 and 2
- Shows:
  - A: a row of ligands, each with its density mesh and subsite numbers.
  - B: a GH31 site with residue carbons coloured by the role of the residue (cyan for subsite −1, yellow for aromatic stacking, orange labels for the catalytic D394A and D341) and red distance labels.
  - C: a slate/white surface with Glc+1…+3.
- Style: Renderer UCSF Chimera stated (the image reads as PyMOL-like); thick sticks; bold black labels.

### 12. GH31 cleft surfaces: NtMGAM–acarbose against LlGH31_u1
- File: `12_cleft_surface_GH31_NtMGAM_JBC2022_Fig8.jpg`
- Source: as 11, Fig. 8. Licence **CC BY 4.0**.
- Cat. 2
- Shows a uniformly coloured opaque surface (yellow, or white with a slate loop) with the sugar in plain sticks. A dotted circle marks the pocket entrance that differs between the two enzymes.
- Style: one colour per surface; labels large and bold; the dotted circle used as an annotation.
- Why it matters: a GH31 MGAM/acarbose cleft, directly relevant to our GH31 example. Very readable, but the flat yellow looks dated.

### 13. GH1 β-glucosidase Td2F2: ChimeraX sticks and a silhouetted clipped surface
- File: `13_activesite_chimerax_GH1_JBC2026_Fig1.jpg`
- Source: Lin et al., *J. Biol. Chem.* 302, 111293 (2026), doi:10.1016/j.jbc.2026.111293. Fig. 1. https://pmc.ncbi.nlm.nih.gov/articles/PMC13015726/
- Image URL: https://pmc-oa-opendata.s3.amazonaws.com/PMC13015726.1/gr1.jpg
- Licence: **CC BY 4.0**
- Cat. 1 and 2
- Shows:
  - A, B, D: protein carbons tan (Chimera "tan"); disaccharide carbons lime green or cyan; subsites written *inside* the rings (−1′, +1′); magenta and blue dashed contacts with distances.
  - C: a white/grey surface clipped through the active site, with **black silhouette edges**, with residue sticks and the sugar showing through.
- Style: UCSF Chimera/ChimeraX (stated: Chimera, PyMOL, PLIP); thin sticks of about 0.15 Å; silhouettes on the surface panel.
- Why it matters: the clearest example of the silhouetted grey surface. **Recipe A, cleft variant.**

### 14. Porcine pancreatic α-amylase (GH13) with the peptide inhibitor helianthamide (Withers lab)
- File: `14_cleft_PPA_helianthamide_ACSCentSci2016_Fig2.jpg`
- Source: Tysoe C. et al., *ACS Cent. Sci.* 2, 154 (2016), doi:10.1021/acscentsci.5b00399. Fig. 2. https://pmc.ncbi.nlm.nih.gov/articles/PMC4819454/
- Image URL: https://pmc-oa-opendata.s3.amazonaws.com/PMC4819454.1/oc-2015-00399q_0004.jpg
- Licence: **ACS AuthorChoice** (style reference only)
- Cat. 2 and 1
- Shows:
  - a: a white opaque surface of the whole amylase with the inhibitor as CPK spheres (blue C, red O), the "where is it" view.
  - b, d: grey protein carbons against blue inhibitor carbons, with the background side chains faded as thin lines.
- Style: PyMOL (stated). A two-colour scheme (grey protein, one accent colour); background fading through depth cue and thinner lines.

### 15. Graphical abstract: GAA stabiliser
- File: `15_GA_GAA_cyclosulfamidate_JACS2022_TOC.jpg`
- Source: as 03, TOC graphic. Image URL https://pmc-oa-opendata.s3.amazonaws.com/PMC9389588.1/ja2c05666_0008.jpg. Licence **CC BY 4.0**.
- Cat. 3
- Shows:
  - left: a framed box with the 2D formula and a small dark-green cartoon of rhGAA;
  - right: a circular "magnifier" holding the active-site zoom from Fig. 3, linked by dotted lines to icons (tube, cells, zebrafish).
- Style: flat icons, a sans-serif type, one accent colour (magenta) for the disease name. Everything sits on white.
- Why it matters: shows how a crystal zoom becomes a TOC element through the magnifier circle.

### 16. Graphical abstract: the α-glucosidase conformational itinerary
- File: `16_GA_alphaglucosidase_itinerary_ACSCentSci2017_TOC.jpg`
- Source: as 04, TOC. Image URL https://pmc-oa-opendata.s3.amazonaws.com/PMC5532717.1/oc-2017-002148_0006.jpg. Licence **ACS AuthorChoice**.
- Cat. 3
- Shows a triptych in three rounded boxes: the crystal zoom (Michaelis complex) on the left, the 2D mechanism (⁴C₁ → ¹S₃) in the centre, the crystal zoom (covalent) on the right. The box frames are colour-coded (teal, peach, dark red) and joined by double arrows.
- Why it matters: the best template for a GA that puts structure and chemistry side by side.

### 17. Graphical abstract: GH31 activity-based probe workflow
- File: `17_GA_GH31_ABP_ACSCentSci2016_TOC.jpg`
- Source: Jiang J. et al., *ACS Cent. Sci.* 2, 351 (2016), doi:10.1021/acscentsci.6b00057. TOC. https://pmc.ncbi.nlm.nih.gov/articles/PMC4882745/. Image URL https://pmc-oa-opendata.s3.amazonaws.com/PMC4882745.1/oc-2016-000576_0007.jpg. Licence **ACS AuthorChoice**.
- Cat. 3
- Shows the enzyme as a stylised blue "Pac-Man" pictogram holding the 2D covalent adduct, then three outputs (gel, MS, crystal zoom).
- Style: a pictogram enzyme plus a real zoom in a small inset; light blue throughout.

### 18. Graphical abstract: whole α-amylase surface with a bound inhibitor
- File: `18_GA_PPA_surface_ACSCentSci2016_TOC.jpg`
- Source: as 14, TOC. Image URL https://pmc-oa-opendata.s3.amazonaws.com/PMC4819454.1/oc-2015-00399q_0001.jpg. Licence **ACS AuthorChoice**.
- Cat. 3
- Shows a white SES surface with the inhibitor in CPK spheres sitting in the cleft, next to the screening plot.
- Style: plain PyMOL surface with no AO. Clear, but flat; the ChimeraX `lighting soft` version (Recipe C) is the upgrade.

### 19. Graphical abstract: GH inhibitors and ABPs (review)
- File: `19_GA_GH_inhibitors_JACS2024_TOC.jpg`
- Source: Artola M., Aerts J. M. F. G., ..., Davies G. J., Overkleeft H. S. *J. Am. Chem. Soc.* 146, 24729 (2024), doi:10.1021/jacs.4c08840. TOC. https://pmc.ncbi.nlm.nih.gov/articles/PMC11403624/. Image URL https://pmc-oa-opendata.s3.amazonaws.com/PMC11403624.1/ja4c08840_0008.jpg. Licence **CC BY 4.0**.
- Cat. 3
- Shows a central mechanism box with four quadrants of BioRender-type icons; a rainbow-ish GH cartoon with CPK ligand in the top left.
- Style: teal headings in a large bold sans; light-blue rounded central panel.
- Note: the protein render here is dated (saturated purple and brown cartoon). Copy the layout, not the render.

### 20. Goodsell / Illustrate style: porcine pancreatic α-amylase (1PPI) with the acarbose-derived oligosaccharide
- File: `20_overview_Goodsell_amylase_1ppi_PDB101_MotM74.png` (converted from GIF)
- Source: D. S. Goodsell, RCSB PDB-101 Molecule of the Month #74 "Alpha-amylase" (Feb 2006), doi:10.2210/rcsb_pdb/mom_2006_2. https://pdb101.rcsb.org/motm/74
- Image URL: https://cdn.rcsb.org/pdb101/motm/74/74_1ppi.gif
- Licence: **CC BY 4.0**. Attribute "David S. Goodsell and the RCSB PDB" (PDB-101 image-download and how-to-cite pages).
- Cat. 2 and 3
- Shows every atom as a sphere: protein in light blue, each atom with its own pale highlight; the sugar in pale yellow lying in the cleft; thick black outlines round the molecule and thin outlines at depth breaks.
- Style:
  - Renderer: Illustrate (Goodsell's Fortran program).
  - Flat shading with no specular light, plus outlines.
  - Two colours only. White background.
- Why it matters: reads at stamp size, which makes it ideal for a graphical abstract or the "where is the cleft" thumbnail. **Recipe C imitates it.**

### 21. PDB-101 interactive-view screenshot (Jmol) of the same site: counter-example for the zoom
- File: `21_activesite_jmol_amylase_PDB101_MotM74.jpg`
- Source: PDB-101 MotM #74. Image URL https://cdn.rcsb.org/pdb101/motm/74/74-Alphaamylase-amylase.jpg. Licence **CC BY 4.0**.
- Cat. X
- Shows: black background, purple Cα trace, yellow ball-and-stick sugar, a huge green Ca²⁺ sphere, and the catalytic acids in grey/red sticks.
- Why it is a counter-example: the dark background, the backbone as a trace and the oversized ion are the "viewer screenshot" look the author dislikes.

### 22. Porcine MGAM N- and C-terminal domains (GH31) with acarbose-derived AC5: counter-example
- File: `22_counterexample_crowded_MGAM_JEIMC2026_Fig5.jpg`
- Source: Watanabe et al., *J. Enzyme Inhib. Med. Chem.* 41, 2612391 (2026), doi:10.1080/14756366.2025.2612391. Fig. 5. https://pmc.ncbi.nlm.nih.gov/articles/PMC12805851/. Image URL https://pmc-oa-opendata.s3.amazonaws.com/PMC12805851.1/IENZ_A_2612391_F0005_C.jpg. Licence **CC BY-NC 4.0**.
- Cat. X
- Shows: about 20 superposed residues in three colours, every one labelled twice (pig and human numbering), with no cartoon context.
- Why it is a counter-example: the same "all residues within 4 Å, all labelled" look as our current PyMOL renders. Labels overlap and nothing has priority.

### URL-only references (no downloadable licence stated, or software pages)
- **U1. ChimeraX gallery, "Calmodulin and target peptide":**
  - Image: https://www.rbvi.ucsf.edu/chimerax/gallery_samples/cam.png
  - Commands: https://www.rbvi.ucsf.edu/chimerax/gallery_samples/cam.cxc
  - What it does: `light soft` (with `light gentle` as the "less rustic" option), `graphics silhouettes true width 1`, a cartoon under a surface at 95 % transparency, pastel colours.
  - Licence not stated, so not downloaded.
- **U2. ChimeraX gallery, "Cyclodextrin pore" (a carbohydrate in a protein channel):**
  - Image: https://www.rbvi.ucsf.edu/chimerax/gallery_samples/xray.png
  - Commands: https://www.rbvi.ucsf.edu/chimerax/gallery_samples/xray.cxc
  - What it does: cyclodextrins as spheres (C steelblue, O firebrick); a white protein surface at 82 % transparency; `light gentle`; `cartoon style width 1.3`; `light depthcuestart .4 depthcueend .8`.
- **U3. RCSB Mol* 3D view of GAA–acarbose (PDB 5NN8):** https://www.rcsb.org/3d-view/5NN8. Mol* "illustrative" preset (outline plus flat shading plus AO), useful for quick previews.
- **U4. Illustrate (Goodsell):**
  - Web server: https://ccsb.scripps.edu/illustrate/
  - Source: https://github.com/ccsb-scripps/Illustrate (Fortran; it could be compiled here with gfortran).
- **U5. Molecular Nodes (Blender):** https://bradyajohnston.github.io/MolecularNodes/. Not installed (no Blender).
- **U6. CAZypedia GH13 and GH31 pages** have **no** structural figures (checked 2026-11-24). There is nothing to take from them.

---

## B. What the good examples have in common

1. **Few residues.** Three residues (example 04) up to about 10 (examples 01, 05). The nucleophile and the acid/base are always there, and always set apart by a different carbon colour or a different label colour. The Davies-group figures add the role in parentheses: "Asp518 (nuc.)", "Asp616 (a./b.)".
2. **The ligand has its own hue.** The sugar's carbons get a saturated accent (orange, gold, slate blue, cyan, lime). The protein's carbons get a muted hue (grey, tan, salmon, sage, or the domain colour). Weak figures use one hue for both (example 05).
3. **Protein context is faint.** The cartoon is drawn thin and pale, at 40–70 % transparency, or it is left out entirely (example 01, left). Background side chains are drawn as thin lines or faded by depth cue.
4. **Ligand sticks are thicker than residue sticks.** Roughly 0.20–0.25 Å against 0.15 Å. In ChimeraX papers the sugar is often ball-and-stick (example 07).
5. **Subsites are labelled −n … +n** beside or inside each glucose ring (examples 01, 06, 10, 11, 13). For a starch enzyme this is the single most informative annotation.
6. **H-bonds are thin dashes** in dark grey, with distances only when the text needs them.
7. **Labels** are sans-serif, one-letter code plus number (D357) or three-letter code (Asp357), placed off the atoms. Colour coding by role is used sparingly.
8. **Background is white**, with a light orthographic or near-orthographic view.
9. **Modern renders** (2024–2026, ChimeraX) add ambient occlusion (`lighting soft` or `gentle`) and, on surfaces, black silhouettes. Classic PyMOL renders (Nature Commun 2017, JACS 2021–22) use flat lighting with no outline. Both work because of points 1–4, not because of the lighting.
10. **Layout:** a small whole-protein thumbnail with a box, then the zoom (examples 01, 02, 07, 09). Graphical abstracts put the crystal zoom inside a circle or a rounded box next to a 2D formula (examples 15, 16).
11. **2D interaction diagrams** in these journals are LigPlot+ (example 08h) or PLIP / PoseView. Only example 08 uses one; the others draw ChemDraw formulas beside the 3D panel instead.

---

## C. Recipes reproducible with the local tools

All three recipes were test-run on 2026-11-24 on PDB 3WEO, sugar-beet GH31 α-glucosidase with
acarviosyl-maltooligosaccharide; the sugar is in chain D, residues AC1/GLC. The test files are
in the session scratchpad, outside the project. Commands are valid for ChimeraX 1.11.1 offscreen
and PyMOL 2.5.0. The camera view must still be chosen by hand for each enzyme (save it with
`view name v1` in ChimeraX, or `get_view` in PyMOL).

### Recipe A: "Nat Commun / ChimeraX clean" zoom and cleft (imitates 07, 09 and 13; colour logic of 01 and 03)

Active-site zoom (`.cxc`):
```
open XXXX.cif
delete solvent ; delete H ; delete :GOL,SO4,EDO,PEG,NAG,FUC     # strip additives
set bgColor white
graphics quality 5
lighting soft multiShadow 512          # ambient occlusion; 'lighting gentle' if soft looks smudgy
lighting depthCue false                # or keep depthCue to fade the background
graphics silhouettes true width 3 color black depthJump 0.01
hide atoms ; cartoon
cartoon style width 1.2 thickness 0.3
color protein #DCDCDC target c         # pale grey cartoon (or the domain colour at 30 % saturation)
transparency protein 50 target c       # optional: Davies-style faint cartoon
name lig  <sugar residues in subsites -2..+2>
name site <8-10 residues: nucleophile, acid/base, stacking Trp/Tyr/Phe, H-bonding Asp/Arg/His>
show lig atoms ; style lig stick ; size lig stickRadius 0.22     # or 'style lig ball' + 'size lig ballScale 0.22'
show site atoms ; style site stick ; size site stickRadius 0.16
color lig  & C #F4B942 target a        # gold sugar carbons
color site & C #6BAF92 target a        # sage residue carbons
color <nuc>,<acid/base> & C #C76BB5 target a   # catalytic pair in magenta
color (lig|site) & O #E03C31 target a
color (lig|site) & N #2E6DB4 target a
color site & S #E6C229 target a
hbonds lig restrict site reveal false color #404040 radius 0.05 dashes 8 name hb
label site residues text "{0.label_one_letter_code}{0.number}" height 0.8 color #222222 bgColor none
2dlabels text "−1" xpos .. ypos .. size 20      # subsite numbers (or add in Inkscape)
view lig ; clip near -9 far 9 position lig
save zoom.png width 2400 height 2000 supersample 3
```
Notes:
- `target a` is required on every atom colouring; without it the cartoon takes the element colours.
- Use `/D` (or the right chain) in the ligand spec: in mmCIF the oligosaccharide is a separate branched chain.
- Labels overlap at this size. The good figures place them by hand, so export without labels and add them in Inkscape at 7–8 pt Arial/Helvetica.

Cleft-overview variant (imitates 07 and 13):
```
surface /A
color /A #E8E8E8 target s              # light grey surface; or the domain colour, desaturated
transparency /A 65 target s            # see the cartoon and the sugar through it (07)
#  or opaque + clipped through the cleft (06A, 13C):
#  transparency /A 0 target s ; clip near -2 position lig
style lig ball ; size lig ballScale 0.25
color lig & C white target a           # 07 uses white-carbon ball-and-stick sugars
lighting soft ; graphics silhouettes true width 2
```

### Recipe B: "JACS / ACS Cent Sci / Nat Commun-2017 PyMOL" zoom (imitates 01, 03, 04 and 05)

Settings for a `.pml`; test-run with `/usr/bin/pymol -cq`, PyMOL 2.5:
```
bg_color white
set ray_opaque_background, 1
set orthoscopic, 1
set antialias, 2
set ray_trace_mode, 1          # 1 = colour + black outline; 0 = no outline (01/03 look); 3 = quantized "cartoonish"
set ray_trace_gain, 0.08       # outline thickness
set ray_trace_color, black
set ambient, 0.45
set direct, 0.55
set reflect, 0.25
set specular, 0.15
set shininess, 30
set ray_shadows, 0
set ambient_occlusion_mode, 1  # surface/stick AO in open-source PyMOL >= 1.8
set ambient_occlusion_scale, 15
set depth_cue, 0               # or 1 with fog_start 0.45 to fade the background
set valence, 0                 # no double-bond drawing
set cartoon_transparency, 0.45
set cartoon_tube_radius, 0.35
set cartoon_loop_radius, 0.25
set stick_radius, 0.16         # residues
set stick_radius, 0.22, lig    # sugar thicker
set dash_radius, 0.05
set dash_gap, 0.25
set dash_length, 0.12
set dash_color, grey40
set label_font_id, 5           # plain sans (7 = bold sans, too heavy)
set label_size, 16
set label_color, black
set float_labels, 1
set_color shell, [0.93,0.78,0.76]          # Davies salmon (03) - or [0.85,0.85,0.85] grey
set cartoon_color, shell                   # set AFTER atom colouring, else CA colour wins
set_color sugarC, [0.95,0.73,0.26]         # gold sugar
set_color siteC,  [0.42,0.69,0.57]         # sage residues
set_color catC,   [0.78,0.42,0.71]         # magenta catalytic pair (01)
color sugarC, lig and elem C
color siteC, site and elem C
color catC, <nuc+acid/base> and elem C
color red, (lig or site) and elem O
color slate, (lig or site) and elem N
distance hb, lig, site, mode=2
hide labels, hb
orient lig
zoom lig, 4
clip slab, 22
png zoom.png, width=2400, height=2000, dpi=300, ray=1
```
PyMOL 2.5 quirks found while testing:
- The usual `label ..., "%s%s" % (one_letter[resn], resi)` crashes PyMOL 2.5 under Python 3.12 (traceback inside the error handler).
- `cmd.label(..., space=)` does not exist in 2.5.
- What works is a short `python` block that loops over `cmd.get_model("site and name CA").atom` and calls `cmd.label(sel, repr("D357"))`. Alternatively use the flatpak PyMOL 3.1 (`flatpak run org.pymol.PyMOL -cq`).

For a density-mesh panel as in 01 and 03:
- `isomesh m, map, 1.0, lig, carve=1.6`
- `set mesh_width, 0.5`
- colour the mesh green for the ligand and slate for the residues.

### Recipe C: Goodsell / Illustrate-like overview for the cleft thumbnail and the graphical abstract (imitates 20; 18 upgraded)

ChimeraX:
```
open XXXX.cif
delete solvent ; delete H ; delete :GOL,SO4,NAG,FUC
set bgColor white
hide cartoons
show /A atoms ; style /A sphere                 # true Goodsell: atom spheres (or 'surface /A' for SES)
color /A #92B4E3 target a                      # Goodsell light blue; carbons slightly lighter:
color /A & C #A9C6EE target a
name lig <sugar> ; show lig atoms ; style lig sphere
color lig #F7C948 target a ; color lig & O #F09A3E target a   # warm yellow/orange sugar
lighting flat                                   # Illustrate look: flat + outlines
graphics silhouettes true width 3 color black depthJump 0.01
#  AO alternative (more 3D, as in the 2024-26 Nat Commun figures):
#  lighting soft ; lighting depthCue false
view ; turn ... (cleft facing the viewer)
save overview.png width 2000 height 2000 supersample 3
```
- Tested: flat plus silhouettes, and soft AO, both render offscreen. The soft-AO surface reads well even at thumbnail size.
- For the true Illustrate output, compile https://github.com/ccsb-scripps/Illustrate with gfortran and use its two-colour plus outline parameters. The ChimeraX version above is close enough for a GA.

### Which recipe for which figure
- **Active-site figure (GH13 and GH31 side by side):** Recipe A, zoom. Same camera logic and colours for both enzymes. Catalytic pair in magenta, labelled "D… (nuc.)" and "E…/D… (a./b.)". At most 8–10 residues. Subsite numbers −2…+2. If an interaction table is needed, add a LigPlot+ or PLIP 2D panel as in 08h. If the author prefers the classic Nat Commun 2017 / JACS look, use Recipe B with `ray_trace_mode 0`.
- **Cleft or whole-enzyme figure:** Recipe A, cleft variant: transparent or clipped grey surface with soft AO and silhouettes, and the oligosaccharide in ball-and-stick over the subsites. Add a Recipe C thumbnail with a box marking the zoom.
- **Graphical abstract:** Recipe C overview of each enzyme, plus a circular inset of the Recipe A zoom (layout of 15 and 16), plus a 2D formula of the α(1→4) substrate.
