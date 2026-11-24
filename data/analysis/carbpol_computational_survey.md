# Carbohydrate Polymers: how often does it publish purely computational papers? (survey 2026-11-24)

**Question.** Does *Carbohydrate Polymers* (Elsevier, ISSN 0144-8617, e-ISSN 1879-1344) take a purely computational paper based on a structural database? The case in hand is a census of PDB complexes of starch-degrading enzymes, covering subsite contacts and conservation, with published kinetics used as validation.

**What was read.** Classification uses **title and abstract only**. ScienceDirect full text returns HTTP 403 to automated fetch. None of the papers below are in Europe PMC full text, and OpenAlex gives no repository PDF for any of them except French et al. 2021, where the link is again ScienceDirect. No paper was read in full. The scope statement was read in full from the Elsevier shop page.

## 1. Data and method

| Source | Query | Records |
|---|---|---|
| Crossref `journals/0144-8617/works`, type journal-article, 2024-01-01 to 2026-11-24, no query | denominator (all item types) | 3966 |
| OpenAlex, `primary_location.source.issn:0144-8617`, from 2024-01-01 | articles / reviews / errata / paratext / retractions | 3006 / 445 / 139 / 69 / 39 |
| PubMed `"Carbohydr Polym"[ta] AND 2024/01/01:2026/11/24[dp]` | all records, all abstracts fetched | 3914 (3727 with abstract) |
| of which: research articles / reviews / errata, retractions, editorials | | **3242** / 489 / 183 |
| Crossref + OpenAlex keyword searches (28 terms: molecular dynamics, DFT, docking, in silico, bioinformatics, PDB, machine learning, glycoside hydrolase, amylase, CAZy, meta-analysis, ...) | candidate harvest | 1420 (only 270 had abstracts) |

Elsevier abstracts are mostly missing from Crossref and OpenAlex. PubMed carries them, so the full 2024–2026 output (3914 records) was downloaded from PubMed and screened. Keyword searching was not the screen.

1. **Computational screen.** The title or abstract mentions a computational method (MD, simulation, DFT/QM, docking, in silico, coarse-grained, force field, free energy, ML/neural network, database, PDB, meta-analysis, bioinformatics, phylogeny, Monte Carlo, ...). This gives **341 records**.
2. **Experimental screen.** Each of the 341 was scored for experimental cues: prepared, synthesized, extracted, measured, characterized, FTIR, XRD, NMR, SEM, rheology, in vitro/in vivo, assay, recombinant, ITC, films, hydrogels, and so on. All records with few cues (about 150) were read by hand from title and abstract. Every title carrying a computational keyword was also read.
3. **Classes.**
   - **A:** no synthesis, measurement or assay by the authors appears in the abstract.
   - **B:** computation plus the authors' own experiments, with computation a main part (a computational term in the title, three or more computational methods, or an explicit "experimental and simulation" design).
   - **C:** experimental work with one supporting modelling step, typically a docking or MD paragraph.

   The B/C split for the 2024–2026 window was made by this heuristic rule, not paper by paper. A was assigned by hand.
4. **Second search, 2010–2026.** PubMed searches for enzyme and structural-database papers (glycoside hydrolase, amylase, CAZy, CBM, PDB, Cambridge Structural Database, linkage torsions, energy maps, QM/MM, phylogeny).
5. **DOI check.** Every DOI listed below was resolved with `api.crossref.org/works/<doi>`. All 23 resolve to *Carbohydrate Polymers* with the volume and article number shown.

## 2. Counts, 2024-01-01 to 2026-11-24 (PubMed, abstracts)

| Class | n | Share of 3242 research articles |
|---|---|---|
| A: purely computational research article | **15** | **0.46 %** (≈ 1 in 215) |
| A: non-experimental, not original research (data-compiling systematic review; nomenclature proposal) | 2 | — |
| B: computation + own experiments, computation a main part | ≈ 73 | ≈ 2.3 % |
| C: experimental with minor modelling (docking/MD/DFT add-on) | ≈ 200 | ≈ 6.2 % |
| Reviews, retractions or errata that mention computation | 45 | — |
| False positives (e.g. "database matching" of NMR data, kinetic "model") | 9 | — |

The B/C figures are heuristic; read them as ±20 %.

By method, among B and C papers:

| Method | B | C |
|---|---|---|
| MD | 37 | 51 |
| DFT/QM | 9 | 40 |
| Docking | 11 | 39 |
| ML | 13 | 11 |

A papers by year:

| 2024 | 2025 | 2026 (to Nov) |
|---|---|---|
| 11 | 2 | 2 |

## 3. Class A papers, 2024–2026

All are research articles unless marked otherwise. "Validation" says what experimental data the abstract compares the model against.

| # | Authors | Year | Title | DOI | Method | Topic | Enzyme/protein–carbohydrate? | Validation |
|---|---|---|---|---|---|---|---|---|
| 1 | P. Peluso, R. Dallocchio, A. Dessì, A. Salgado, B. Chankvetadze, G.K.E. Scriba | 2024 | Molecular modeling study to unravel complexation of daclatasvir and its enantiomer by β-cyclodextrins. Computational analysis using quantum mechanics and molecular dynamics | 10.1016/j.carbpol.2024.122483 | QM + MD | Cyclodextrin–drug enantioseparation | No | CE enantioseparation results "as benchmark" (not stated as new work) |
| 2 | P. Kumar, R. Purohit | 2024 | Driving forces and large scale affinity calculations for piperine/γ-cyclodextrin complexes: mechanistic insights from umbrella sampling simulation and DFT calculations | 10.1016/j.carbpol.2024.122350 | 1 µs MD, umbrella sampling, DFT | Cyclodextrin inclusion | No | None stated |
| 3 | P. Kumar, V.K. Bhardwaj, R. Purohit | 2024 | Highly robust quantum mechanics and umbrella sampling studies on inclusion complexes of curcumin and β-cyclodextrin | 10.1016/j.carbpol.2023.121432 | MM/PBSA, umbrella sampling, ONIOM QM | Cyclodextrin inclusion | No | None stated |
| 4 | A. Anderson, R. García-Fandiño, Á. Piñeiro, M.S. O'Connor | 2024 | Unraveling the molecular dynamics of sugammadex-rocuronium complexation: a blueprint for cyclodextrin drug design | 10.1016/j.carbpol.2024.122018 | Docking, MD, metadynamics/PMF | Cyclodextrin drug | No | Literature binding constant |
| 5 | A. Yamaguchi, S. Arai, N. Arai | 2024 | Design strategy for blends of biodegradable polyester and thermoplastic starch based on a molecular dynamics study of the phase-separated interface | 10.1016/j.carbpol.2024.122005 | MD, tensile simulation | Starch/polyester blend interface | No | None stated |
| 6 | N. Barrios, J.G. Parra, R.A. Venditti, L. Pal | 2024 | Elucidation of temperature-induced water structuring on cellulose surfaces for environmental and energy sustainability | 10.1016/j.carbpol.2024.121799 | All-atom MD | Cellulose–water | No | None stated |
| 7 | A. Kolaříková, E. Kutálková, J. Hrnčiřík, M. Ingr | 2024 | Hyaluronan oligosaccharides form double-helical duplexes in water:1,4-dioxane mixed solvent | 10.1016/j.carbpol.2023.121632 | MD | Hyaluronan conformation | No | None stated |
| 8 | V.Á. Cabral, B. Govoni, H. Verli | 2024 | Unravelling carbohydrate binding module 21 (CBM21) dynamics of interaction with amylose | 10.1016/j.carbpol.2024.121792 | MD, metadynamics, MM and QM/MM, from a crystal structure | CBM21 (*Rhizopus oryzae*) binding amylose; questions the ¹S₃ Glc seen in the crystal | **Yes: the only protein–starch paper in class A** | Crystal structure |
| 9 | M.A. Bagherjeri, H. Monhemi, A.N.M.A. Haque, M. Naebe | 2024 | Molecular mechanism of cellulose dissolution in N-methyl morpholine-N-oxide: a molecular dynamics simulation study | 10.1016/j.carbpol.2023.121433 | µs MD | Cellulose dissolution | No | "Agreement with experimental observations" (literature) |
| 10 | W. Liu, A. Shomali, C. Zhang, B. Coasne, J. Carmeliet, D. Derome | 2024 | Nanostructure and interfacial mechanical properties of PEG/cellulose nanocomposites studied with molecular dynamics | 10.1016/j.carbpol.2024.122429 | MD, pull-out | Cellulose/PEG (wood conservation) | No | None stated |
| 11 | M. Másson | 2024 | The quantitative molecular weight-antimicrobial activity relationship for chitosan polymers, oligomers, and derivatives | 10.1016/j.carbpol.2024.122159 | **Meta-analysis of published data**: bilinear equation fitted to 29 MIC datasets from 7 studies (1984–2019) | Chitosan Mw versus antimicrobial activity | No | Published data only (single author) |
| 12 | Y. Hou, R.-Y. Dong | 2025 | Comparative analyses of thermal transport in cellulose at crystalline, paracrystalline and amorphous states | 10.1016/j.carbpol.2025.124063 | MD, phonon analysis | Cellulose thermal conductivity | No | None stated |
| 13 | A. Elahi, X. Yan, S. Chaudhuri | 2025 | Enhancing the OPLS-AA force field for cellulose Iβ: structural stability and surface functionalization capability with the CM5 charge model | 10.1016/j.carbpol.2025.123572 | Force-field development (OPLS-CM5), compared with CHARMM36 and GLYCAM06 | Cellulose Iβ | No | Literature unit cell (<1.5 % error) |
| 14 | S. Sang, L. Zhou, J. Ye, P. Chen, Z. Xiang | 2026 | Computational insights of noncovalent interactions in xylan hydrate crystal | 10.1016/j.carbpol.2025.124766 | Conformer search, DFT, MD, energy decomposition | Xylan hydrate crystal H positions | No | Published crystal structure |
| 15 | Z. Yu, X. Zheng, H. Yang, W. Tu, Y. Zhang | 2026 | Mechanical role of disordered cellulose in macrofibrils of softwood secondary cell walls | 10.1016/j.carbpol.2026.125484 | Reactive MD | Wood macrofibril mechanics | No | None stated |

**A, not original research:**

| # | Authors | Year | Title | DOI | Type and content |
|---|---|---|---|---|---|
| 16 | B.M. Guerreiro, J.C. Lima, J.C. Silva, F. Freitas | 2026 | Polysaccharides in cryopreservation: multidimensional systematic review of extremophilic traits and the role of selective pressure in structure-function relationships | 10.1016/j.carbpol.2025.124462 | Systematic review that compiles a literature database (XPOL-DB, 145 EPS) and analyses it. PubMed type: Systematic Review |
| 17 | A.M. Anderson, M.S. O'Connor, J. Pipkin, M. Malanga, T. Sohajda, T. Loftsson, L. Szente, R. García-Fandiño, Á. Piñeiro | 2025 | A comprehensive nomenclature system for cyclodextrins | 10.1016/j.carbpol.2025.123600 | Nomenclature proposal. PubMed type: Review |

**Borderline, classed B and not A:**

| Paper | DOI | Why not A |
|---|---|---|
| Richardson, Berti, Ravenscroft, Kuttel 2025, molecular modelling of GBS type II/III capsular polysaccharides | 10.1016/j.carbpol.2025.123964 | Modelling "corroborated by through-space NMR correlations"; the filter-retention data come from the vaccine process. The abstract does not say who acquired the NMR. |
| Borah et al. 2026, preferential solvation | 10.1016/j.carbpol.2026.125726 | "Simulations, theory, and experiments" |
| He et al. 2025, coarse-grained MD of HMSA droplets | 10.1016/j.carbpol.2025.124370 | Own QCM-D and interfacial-tension validation |
| Li et al. 2026, DPD of β-CD/Tween 60 | 10.1016/j.carbpol.2026.125630 | Own interfacial-tension measurements |

**What the A papers look like.**
- The 15 research papers are of three kinds:
  - 13 atomistic MD and/or QM studies on the physical chemistry of a named polysaccharide or cyclodextrin: cellulose (6), cyclodextrins (4), xylan, hyaluronan, starch blend.
  - 1 force-field paper.
  - 1 meta-analysis of published data (Másson).
- One of the 15 involves a protein–carbohydrate complex: CBM21 with amylose, MD and QM/MM from a crystal structure.
- None is a structural-database census, bioinformatics or conservation study, or a GH-family structure survey.
- Four of the 15 are cyclodextrin–guest complexation studies. Two of these (Kumar & Purohit) concern a single cyclodextrin type per complex series. The current scope bars "complexation of a drug with a single cyclodextrin", so that kind of paper is now riskier.
- The yearly count falls: 11 in 2024, 2 in 2025, 2 in 2026. This fits a tightened editorial stance, but small numbers do not prove it.

## 4. Enzyme–substrate structural, glycoside hydrolase and PDB-census papers in the journal, 2015–2026

| Paper | Year | DOI | Class | Relevance |
|---|---|---|---|---|
| French, Montgomery, Prevost, Edwards, Woods, "Comparison of cellooligosaccharide conformations in complexes with proteins with energy maps for cellobiose", vol. 264, 118004 | 2021 | 10.1016/j.carbpol.2021.118004 | **A, PDB census** | φ/ψ of cellooligosaccharide and lactose moieties in protein complexes from the PDB, compared with QM energy maps. **The only PDB-census paper found in the journal.** It is the closest precedent for both this project and the α(1→4) linkage work. |
| Cabral, Govoni, Verli, CBM21 with amylose (MD, QM/MM) | 2024 | 10.1016/j.carbpol.2024.121792 | A | Protein–starch simulation |
| Nam, French, Condon, Concha, "Segal crystallinity index revisited by the simulation of X-ray diffraction patterns", vol. 135, 1–9 | 2016 | 10.1016/j.carbpol.2015.08.035 | A (abstract not re-read here) | Simulation from crystal structures |
| Zinck, Perna, Hyldgaard, Hoff, Meyer, "Multi-analytical profiling of glycoside hydrolase 13 (GH13) α-amylases on potato starch", vol. 380, 125029 | 2026 | 10.1016/j.carbpol.2026.125029 | B/C | 20 GH13 α-amylases from several subfamilies, measured on starch. Subfamily predicted poorly; **active-site openness correlated with substrate preference**. The nearest topical precedent, though experimental. |
| de Araújo et al., GH9 BlCel9A (assays, X-ray, SAXS, MD) | 2024 | 10.1016/j.carbpol.2023.121739 | B | |
| Lindner, Bonin, Hellmann, Moerschbacher, "Three intertwining effects guide the mode of action of chitin deacetylase" (subsite capping, occupation, preference; includes an in silico analysis guide) | 2025 | 10.1016/j.carbpol.2024.122725 | Review | Subsite-level CE4 review |

Other enzyme papers found (2015–2026) are biochemical or crystallographic characterisations with own experiments. Examples, cited by PubMed ID only and not DOI-verified for this report:

| PMID | Paper |
|---|---|
| 38710568 | GH30 xylanase crystal structure |
| 34893255 | GH79 hyaluronate glycanohydrolase |
| 33910709 | GH48 Cel48B |
| 41475721 | GH28 polygalacturonase |
| 38616101 | GH174 fucanases |
| 37940306 | New GH family of endo-1,3-fucanase |
| 41320428 | Chitosanase residues |
| 39486976 | Xanthan-binding CBM |

No purely bioinformatic GH-family paper (sequence or structure conservation, CAZy survey) was found in the journal in 2015–2026.

## 5. Scope statement

Source: https://shop.elsevier.com/journals/carbohydrate-polymers/0144-8617, fetched 2026-11-24, HTTP 200. Quoted verbatim:

> "The role of the well-characterized carbohydrate polymer must be the major proportion of the work reported, not a peripheral topic. At least one named carbohydrate polymer must be cited and be the main focus of the paper and its title."

> "Topics of interest to the journal: structure-property relationships; analytical methods; chemical, enzymatic and physical modifications; biosynthesis; natural functions; interactions with other materials"

> "Topics not of interest to the journal: Bibliometric reviews; **Studies that involve only modelling without any comparison of model results with experimental data, either carried out by the authors or from the literature.** [...] studies of complexation of a drug with a single cyclodextrin [...] **production and isolation of enzymes which act on polysaccharides (studies on the mode of action of an enzyme on a polysaccharide are within the journal scope)** [...] **carbohydrate oligomers where the degree of polymerization is equal to or less than four** [...]"

The journal also states that "Nearly all scientific papers benefit from inclusion of a statement of hypothesis" and that the polysaccharide must be adequately characterised. The page names an open-access companion journal, *Carbohydrate Polymer Technologies and Applications*.

## 6. Conclusion

**Frequency (fact).** Purely computational research articles are rare in *Carbohydrate Polymers*: 15 of 3242 research articles in 2024 to November 2026, about 0.5 %, roughly 5 a year, with 2 a year in 2025 and 2026. Computation with the authors' own experiments is about 5 times as common (≈ 2 %). Experimental papers with a docking or MD add-on are about 13 times as common (≈ 6 %).

**Type (fact).** The accepted purely computational papers are:
- atomistic MD/QM studies of the physical chemistry of a named polysaccharide (cellulose, xylan, hyaluronan, starch blends) or of cyclodextrin complexes;
- one force-field paper;
- one meta-analysis of published data.

Only one involves a protein, CBM21 with amylose. The only PDB-census paper in the journal's 2015–2026 record is French et al. 2021 (cellooligosaccharide linkage conformations in protein complexes against QM maps). No bioinformatics or conservation census of a glycoside-hydrolase family was found.

**Reading against the scope (interpretation).** The written scope does not exclude a purely computational paper if the model is compared with experimental data, and literature data are explicitly allowed. A census validated against published kinetics therefore meets the stated rule. "Mode of action of an enzyme on a polysaccharide" is explicitly in scope. Three requirements bear directly on the project's paper:
1. Starch, amylose or amylopectin must be named in the title and be the main focus, not the enzymes as proteins.
2. The scope excludes studies of oligomers of DP ≤ 4. Most PDB ligands are short maltooligosaccharides, so the paper would need to frame itself around the polymer substrate and its subsite binding.
3. A hypothesis statement is expected.

The precedent is thin: one PDB-census paper in 12 years and no GH bioinformatics paper. A census of enzyme structures is therefore at the edge of what the journal usually publishes, not at its centre.
