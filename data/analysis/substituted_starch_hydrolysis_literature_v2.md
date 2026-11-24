# Enzymatic hydrolysis of substituted starches and oligosaccharides: which substituted positions each subsite tolerates

Literature collected 2026-11-24. Purpose: a reference set to check against the structural map of glucose
oxygens held by invariant residues at each subsite of GH13 (alpha-amylases), GH31 (alpha-glucosidases),
GH14 (beta-amylase) and GH15 (glucoamylase).

**Version 2 (2026-11-24, second search pass).** v1 is kept unchanged as `substituted_starch_hydrolysis_literature.md`.
What v2 changes:

- Now read in full (FT): van der Burgt et al. 1998 (§6.4, reducing-end analysis of alpha-amylase digests of methyl
  starch); Bock & Pedersen 1987 and 1988 (§8.5, deoxy and 6-substituted methyl beta-maltosides with glucoamylase);
  Huang et al. 2007, read as Chapter 4 of the Huang 2006 Wageningen thesis (§3.2, acetylated cowpea starch).
- New primary studies, FT: Okada, Ogawa & Usui 2000 (§8.10, 6-O-, 4-O- and 3-O-glycosylated maltooligosaccharides
  with porcine and human pancreatic alpha-amylase, cleavage frequencies at each bond); Omichi & Ikenaka 1988 (§8.11,
  6-deoxy-6-pyridylamino maltooligosaccharides); Nishio 2002 (§8.12, monodeoxy p-nitrophenyl alpha-glucosides with
  GH31 and GH13 alpha-glucosidases); Huijbrechts 2008 thesis (§2.11, allyloxy-hydroxypropyl starch).
- New primary studies, ABS: Mischnick 2001 (§6.5, methyl amylose, position of the methyl groups in the end
  residues of the digest products); Richardson, Cohen & Gorton 2001 (§2.8, hydroxypropyl); Wootton & Chaudhry 1981
  and Azemi & Wootton 1985 (§2.9–2.10, hydroxypropyl); Manelius 2002 and 2005 (§5.5–5.6, cationic); Schagerlöf 2005
  and 2006 (§9.2, cellulose-ether analogues).
- Still not read (closed access, no legal open copy found): Steeneken & Woortman 1994, Zhao 2012, Kavitha &
  BeMiller 1998, Richardson 2003 (full text), Manelius 2002/2005 (full text). See §13.1.
- Open-access PDFs downloaded in this pass are in `../Literatura_Referencia/` (outside the repository), §14.

## 0. How to read this report

**Subsite notation (Davies, Wilson & Henrissat 1997 convention).** The scissile alpha(1→4) glycosidic linkage lies
between subsite −1 (the glucosyl residue whose C1 is attacked, the glycone) and subsite +1 (the residue whose O4 leaves,
the aglycone). Subsites −2, −3 ... run toward the non-reducing end; +2, +3 ... toward the reducing end.

**Reading a product back onto subsites.** When an endo-acting enzyme cuts, the fragment on the non-reducing side
has its new reducing-end residue in subsite −1, and the fragment on the reducing side has its new non-reducing-end
residue in subsite +1. So, for a limit product of length n with the substituted residue at position k counted from
the reducing end (k = 1 is the reducing end):

- the cut that made the reducing end had the substituted residue at subsite −k (tolerated in a productive complex);
- the cut that made the non-reducing end had the substituted residue at subsite +(n − k + 1) (tolerated), provided
  that this end is not an original chain end;
- linkages *inside* a limit product were not cut under exhaustive conditions, which says the corresponding placements
  are not productive. For short products this is weaker evidence, because short oligosaccharides are poor substrates
  for many alpha-amylases even when unsubstituted.

Superscript notation of Hizukuri/Takeda: 6³-PG3 = maltotriose with a phosphate at C-6 of residue 3 counted from the
reducing end (= the non-reducing-end residue); 3²-PG4 = maltotetraose with a 3-phosphate on residue 2 from the
reducing end.

**Evidence labels.** FT = full text read (with page/figure); FT-thesis = read in the PhD thesis chapter that became the
paper; ABS = abstract only (from PubMed, Crossref or OpenAlex); NOT READ = identified, full text and abstract
unavailable here, cited only for completeness. All DOIs below were resolved through the Crossref API on 2026-11-24.

**Separation of fact and inference.** "Measured" = what the authors report. "Authors' inference" = their wording.
"Subsite reading (this report)" = my mapping of their products onto subsites with the rules above; it is not stated by
the authors unless said so.

---

## 1. Hydroxyethyl starch (HES)

### 1.1 Yoshida, Yamashita, Matsuo, Kishikawa 1973 — ABS
Yoshida M, Yamashita T, Matsuo J, Kishikawa T (1973) Enzymic degradation of hydroxyethyl starch. Part I. Influence of
the distribution of hydroxyethyl groups on the enzymic degradation of hydroxyethyl starch. *Starch/Stärke* 25(11):373–376.
doi:10.1002/star.19730251105

- Enzyme: pig pancreas alpha-amylase (GH13).
- Substrates: 2-O-hydroxyethyl starches (2-HES) and 6-O-hydroxyethyl starches (6-HES) at several MS/DS.
- Measured: rate of hydrolysis decreases with MS or DS for both.
- Key sentences (abstract, verbatim): "the more the unsubstituted anhydroglucose residues the specimen contained, the
  more rapidly the degradation took place." "2-O-hydroxyethylanhydroglucose residues played a more significant role on
  the resistance to α-amylase action than 6-O-hydroxyethyl-anhydroglucose residues. Thus, each sample of 2-HES has been
  more resistant to α-amylolytic attack than the corresponding specimen of 6-HES, regardless of the same content of
  unsubstituted anhydroglucose residues."
- Authors' inference: O2 substitution obstructs pancreatic alpha-amylase more than O6 substitution at equal content of
  unsubstituted units.

### 1.2 Yoshida & Kishikawa 1984 (Part II) — ABS
Yoshida M, Kishikawa T (1984) A study of hydroxyethyl starch. Part II. Degradation-sites of hydroxyethyl starch by pig
pancreas α-amylase. *Starch/Stärke* 36(5):167–169. doi:10.1002/star.19840360505

- Method: GLC of degradation sites.
- Key sentence (abstract, verbatim): "HES could be hydrolyzed by α-amylase only at the glucosidic bond of unsubstituted
  glucose residues, but not at the bond of 2-O-, 3-O-, and 6-O-hydroxyethylglucose residues."
- Subsite reading (this report): the abstract does not say whether "the bond of" a residue means its own C1 linkage
  (residue at −1) or both linkages flanking it (−1 and +1). The minimum certain reading is: a hydroxyethylated residue
  (any of O2, O3, O6) is not found at subsite −1 of a productive pancreatic alpha-amylase complex. The full text is
  needed to settle +1.

### 1.3 Yoshida, Minami, Kishikawa 1984 (Part III) — ABS
Yoshida M, Minami Y, Kishikawa T (1984) A study of hydroxyethyl starch. Part III. Comparison of metabolic fates between
2-O-hydroxyethyl starch and 6-O-hydroxyethyl starch in rabbits. *Starch/Stärke* 36(6):209–212.
doi:10.1002/star.19840360606

- Measured (in vivo, rabbit): 6-O-HES eliminated faster than 2-O-HES; 6-O-HES behaved like a 2-O-HES "which had a lower
  molar substitution of hydroxyethyl groups by 0.2 to 0.3".
- Authors' inference: "attributed to the substitution-site effect of hydroxyethyl groups on the enzymic degradation of
  hydroxyethyl starch by α-amylase."

### 1.4 Merkus, Mourits, De Galan, De Jong 1977 — ABS (context: where HES is substituted)
Merkus HG, Mourits JW, De Galan L, De Jong WA (1977) Substituent distribution in hydroxyethyl starch. *Starch/Stärke*
29(12):406–409. doi:10.1002/star.19770291203

- Measured: MS 0.03–1.9; "In the monosubstituted derivatives 70–80% of the hydroxyethyl groups is attached to the C2
  atom. The ratio of C3 to C6 monohydroxyethylated products changes from 2 to 0.2 when increasing the MS from 0.03 to
  1.9." At MS ≈ 2, ~50% of disubstituted products are 2-hydroxyethoxyethyl (chain extension on the new OH).

### 1.5 Treib et al. 1995 — ABS (clinical C2/C6 ratio)
Treib J, Haass A, Pindur G, Seyfert UT, Treib W, Grauer MT, Jung F, Wenzel E, Schimrigk K (1995) HES 200/0.5 is not
HES 200/0.5. Influence of the C2/C6 hydroxyethylation ratio of hydroxyethyl starch (HES) on hemorheology, coagulation
and elimination kinetics. *Thromb Haemost* 74(6):1452–1456. doi:10.1055/s-0038-1649964

- Two HES 200/0.5 with C2/C6 = 13.4 and 5.7, ten patients each, 10 days.
- Measured: after repeated dosing, intravascular mean MW 95 vs 84 kD; serum concentration rise 14.7 vs 8.6 mg/mL;
  PTT +30% vs +13% (p < 0.01).
- Authors' inference (verbatim): "The higher C2/C6 ratio of HES 200/0.5/13.4 slows down enzymatic degradation."

### 1.6 Treib, Baron, Grauer, Strauss 1999 — ABS (review)
Treib J, Baron JF, Grauer MT, Strauss RG (1999) An international view of hydroxyethyl starches. *Intensive Care Med*
25(3):258–268. doi:10.1007/s001340050833
- Verbatim: "The rate of metabolization depends mainly on the degree of hydroxyethyl substitution (ranging from 0.4 to
  0.7), and the C2/C6 ratio of hydroxyethylation. A high degree of substitution and a high C2/C6 ratio lead to a slow
  metabolization of HES".

### 1.7 Jungheinrich & Neff 2005 — ABS (review)
Jungheinrich C, Neff TA (2005) Pharmacokinetics of hydroxyethyl starch. *Clin Pharmacokinet* 44(7):681–699.
doi:10.2165/00003088-200544070-00002
- HES types with MS > 0.4 accumulate in plasma on repeated dosing; HES 130/0.4 does not; plasma clearance of HES 130/0.4
  "at least 20-fold higher than that for hetastarch" (HES 670/0.75).

**HES summary.** O2 hydroxyethylation obstructs pancreatic/serum alpha-amylase more than O6 at equal MS (in vitro,
Yoshida 1973; in vivo, Yoshida 1984 III; clinical, Treib 1995). Pancreatic alpha-amylase does not cleave at the
glucosidic bond of 2-, 3- or 6-O-hydroxyethyl glucose (Yoshida 1984 II). No subsite-resolved MS study of HES digests
was found with the tools available (see Gaps).

---

## 2. Hydroxypropyl starch (HP)

### 2.1 Leegwater & Luten 1971 — ABS
Leegwater DC, Luten JB (1971) A study on the in vitro digestibility of hydroxypropyl starches by pancreatin.
*Starch/Stärke* 23(12):430–432. doi:10.1002/star.19710231204
- Verbatim: "the digestibility decreased exponentially with increasing degree of substitution."

### 2.2 Leegwater 1972 — ABS
Leegwater DC (1972) A model for the structure of hydroxypropyl starches. *Starch/Stärke* 24(1):11–15.
doi:10.1002/star.19720240105
- Model: random distribution of HP groups; reducing power of pancreatin digests ∝ (1 − DS)^j.
- Fitted j = 5.3. Verbatim: "This result is interpreted as an indication that the pancreatic enzymes act only upon
  α(1→4) glucosidic chains containing five or more glucopyranose units."
- Worked meaning: for DS = 0.1, (1 − 0.1)^5.3 = 0.9^5.3 = 0.57, i.e. 57% of the reducing power of unmodified starch.
- Subsite reading (this report): a run of ~5 unsubstituted residues is a statistical estimate of the number of
  subsites that must be free of HP groups for porcine pancreatic alpha-amylase (whose binding region spans about five
  subsites, −3 to +2; see §8.2) — consistent with HP not being tolerated anywhere in the occupied binding region, but
  this is a model fit, not a positional measurement.

### 2.3 Conway & Hood 1976 — ABS
Conway RL, Hood LF (1976) Pancreatic alpha amylase hydrolysis products of modified and unmodified tapioca starches.
*Starch/Stärke* 28(10):341–343. doi:10.1002/star.19760281005
- Hog pancreatic alpha-amylase; hydroxypropyl distarch phosphate and hydroxypropyl starch gave "lesser amounts of the
  lower molecular weight oligosaccharides (G1—G4)" than unmodified starch.

### 2.4 Richardson, Nilsson, Bergquist, Gorton, Mischnick 2000 — ABS
Richardson S, Nilsson GS, Bergquist KE, Gorton L, Mischnick P (2000) Characterisation of the substituent distribution in
hydroxypropylated potato amylopectin starch. *Carbohydr Res* 328(3):365–373. doi:10.1016/s0008-6215(00)00103-8
- Substrates: HP amylopectin potato starch modified in granular slurry (HPPAPg) or in solution (HPPAPs).
- Enzymes: debranching enzymes; simultaneous alpha-amylase + amyloglucosidase (GH13 + GH15).
- Measured: no difference in O2/O3/O6 distribution between the two; debranching showed more substituents near branch
  points in HPPAPg; alpha-amylase + amyloglucosidase "liberated more unsubstituted glucose units" from HPPAPg.
- Authors' inference: more heterogeneous distribution in granular-modified starch. Method assumes that released free
  glucose is unsubstituted, i.e. that the enzymes do not release HP-glucose.

### 2.5 Kavitha & BeMiller 1998 — NOT READ
Kavitha R, BeMiller JN (1998) Characterization of hydroxypropylated potato starch. *Carbohydr Polym* 37(2):115–121.
doi:10.1016/s0144-8617(98)00056-3 (Crossref lists only the first author.)
- Cited by Bai (2013 thesis, p. 52) for "Inhibition on β-amylolysis has been reported for many chemically modified
  starches". Content not checked here.
- v2: still closed (OpenAlex: no open-access location). Second-hand content, from the Huang 2006 thesis (pp. 14–15,
  verbatim): Kavitha and BeMiller "investigated the isolated amylose fraction from hydroxypropylated potato starch with
  the combination of α-amylase and amyloglucosidase and found heterogeneous substitution within this component"; "the
  ß-limit value of isolated amylopectin fraction from hydroxypropylated potato starch was lower than the value of
  amylopectin from native starch and concluded that some substituents were present near the non-reducing ends of
  amylopectin branch chains." No subsite-level statement is reported second-hand.

### 2.6 Zhao, Schols, Chen, Jin, Buwalda, Gruppen 2012 — NOT READ
Zhao J, Schols HA, Chen Z, Jin Z, Buwalda P, Gruppen H (2012) Substituent distribution within cross-linked and
hydroxypropylated sweet potato starch and potato starch. *Food Chem* 133(4):1333–1340.
doi:10.1016/j.foodchem.2012.02.021
- Companion: Zhao J, Chen Z, Jin Z, de Waard P, Buwalda P, Gruppen H, Schols HA (2015) Level and position of
  substituents in cross-linked and hydroxypropylated sweet potato starches using NMR. *Carbohydr Polym* 131:424–431.
  doi:10.1016/j.carbpol.2015.06.005 — NOT READ. These are the Wageningen enzymatic-fingerprinting papers for HP; abstracts
  were not retrievable (Elsevier, not in PubMed).
- v2: still NOT READ. OpenAlex lists only the publisher and the research.wur.nl metadata record (no file); no open
  thesis containing these papers was found (the first author's PhD was not located in edepot.wur.nl).

### 2.7 Biliaderis 1982 — ABS (HP distarch phosphate, acetylated distarch phosphate)
Biliaderis CG (1982) Physical characteristics, enzymic digestibility and structure of chemically modified smooth pea and
waxy maize starches. *J Agric Food Chem* 30(5):925–930. doi:10.1021/jf00113a029
- Measured: debranching of AcDP-SP and HpDP-WM incomplete, "the debranched chains accounting for only 57-58% of the
  amylopectin present in the unmodified starches." beta-Amylolysis limits: WM 58, SP 66, HpDP-WM 40, AcDP-SP 34%.
- Authors' inference: acetylation "occurred exclusively in certain parts of the granule, whereas hydroxypropylation was
  more uniform".

### 2.8 Richardson, Cohen, Gorton 2001 — ABS (new in v2)
Richardson S, Cohen A, Gorton L (2001) High-performance anion-exchange chromatography–electrospray mass spectrometry for
investigation of the substituent distribution in hydroxypropylated potato amylopectin starch. *J Chromatogr A*
917(1–2):111–121. doi:10.1016/s0021-9673(01)00690-2 (abstract from PubMed 11403463)
- Method: enzymic hydrolysis of HP potato amylopectin, products identified on-line by HPAEC-PAD-ESI-MS, substituted
  and unsubstituted, without standards.
- Key sentence (abstract, verbatim): "the starch hydrolysing enzymes were hindered by the hydroxypropyl groups and
  preferentially cleaved glucosidic linkages between unsubstituted glucose units."
- Subsite reading (this report): the same rule as for the cationic ether (§5.1): an HP-glucose is disfavoured at −1
  and at +1. "Preferentially" means the abstract does not exclude a minor population of cuts next to an HP residue.

### 2.9 Wootton & Chaudhry 1981 (Part II) — ABS (new in v2)
Wootton M, Chaudhry MA (1981) In vitro digestion of hydroxypropyl derivatives of wheat starch. II. Effect of
substitution on the products of partial digestion by porcine pancreatic alpha-amylase. *Starch/Stärke* 33(5):168–170.
doi:10.1002/star.19810330507
- HP wheat starches, DS 0–0.17; partial digestion by PPA (GH13); digest separated into oligosaccharide and
  polysaccharide fractions.
- Measured (abstract, verbatim): "increasing DS caused a decrease in reducing value, an increase in average degree of
  polymerisation and a greater proportion of oligosaccharide in the hydrolysate. The oligosaccharide fraction contained
  most of the modifying groups in the total digest and had approximately double the hydroxypropyl content of the
  undigested starch."
- Reading (this report): PPA releases HP-bearing oligosaccharides, so HP residues are accepted at some subsites of
  productive PPA complexes (which ones is not resolved). Part I (doi:10.1002/star.19810330408, ABS): for gelatinised
  HP starch "the action pattern of the enzyme tended to favour non-repetitive attack, as degree of substitution
  increased."

### 2.10 Azemi & Wootton 1985 — ABS (new in v2)
Azemi BMNM, Wootton M (1985) Action pattern of porcine pancreatic alpha-amylase on hydroxypropyl derivatives of maize,
waxy maize and high amylose maize starches. *Starch/Stärke* 37(2):50–52. doi:10.1002/star.19850370204
- Blue value vs reducing value during PPA digestion. Verbatim: "Increasing substitution with hydroxypropyl groups led
  to increased random attack for both maize and high amylose maize starches. For waxy maize starches random attack
  increased up to a molar substitution of 0.07 and subsequently declined."
- Reading (this report): HP groups shorten the processive (multiple-attack) runs of PPA; this is a kinetic signature,
  not a positional one.

### 2.11 Huijbrechts 2008 (thesis, Chapters 2 and 5) — FT-thesis (new in v2; related ether: 1-allyloxy-2-hydroxypropyl, AHP)
Huijbrechts AML (2008) *Multifunctional starch derivatives: synthesis, characterization and properties*. PhD thesis,
Wageningen University. doi:10.18174/122072 (open access, edepot.wur.nl/122072)
- Substituent: 1-allyloxy-2-hydroxypropyl ether (from allyl glycidyl ether; a hydroxypropyl ether carrying an allyl
  ether on C3 of the propyl chain). Starches: amylose-extender maize (AEMS, 70% amylose) and waxy maize (WMS), both DS
  0.20 (thesis p. 42: "one AHP group per five glucose residues").
- Enzymes (thesis p. 31): pullulanase (B. licheniformis, Megazyme M2), alpha-amylase (B. subtilis, Fluka 10069),
  amyloglucosidase (Sigma A9268, stated "from Aspergillus oryzae"), used in sequence.
- Measured (pp. 40–42, MALDI-TOF-MS of the exhaustive digest): AHP-WMS gave resistant DP3 with one AHP (DS of the
  oligomer 0.33), DP4 and DP5 with one or two AHP, DP6 and DP7 with two or three AHP; AHP-AEMS gave DP3 to DP10,
  including DP3 with one and two AHP. Chapter 5 (p. 103): for AHP-WMS the resistant products ranged "from maltose (DP2)
  to DP9 with one to five allyl groups". Barley beta-amylase (GH14) released 89% of the maltose of native WMS from
  AHP-WMS (p. 102).
- Authors' inference (p. 42, verbatim): "the enzyme-resistant carbohydrate residues contain at least one AHP group and
  confirm that substitution of the polysaccharide chain prevents enzyme-induced degradation due to sterical hindrance."
- Subsite reading (this report): same logic as the acetylated DP3 of Chen 2004 (§3.1): a DP3 (and in Chapter 5 a DP2)
  carrying one AHP group survives glucoamylase, so at least one of the placements "AHP-residue at −1", "at +1" (and, for
  DP3, "at +2") is non-productive for GH15. Position of AHP inside the oligomers was not determined.

---

## 3. Acetylated starch

### 3.1 Chen, Schols, Voragen 2004 — FT-thesis
Chen Z, Schols HA, Voragen AGJ (2004) Differently sized granules from acetylated potato and sweet potato starches differ
in the acetyl substitution pattern of their amylose populations. *Carbohydr Polym* 56(2):219–226.
doi:10.1016/j.carbpol.2004.02.004
Read in: Chen Z (2003) *Physicochemical properties of sweet potato starches and their application in noodle products*.
PhD thesis, Wageningen University, Chapter 6, pp. 100–105. doi:10.18174/121460 (open access, edepot.wur.nl/121460).

- Substrates: amylose populations from acetylated potato and sweet potato starches; DS of amylose about 0.08 (potato),
  0.06 (sweet potato) — "on average, 1 acetyl group on every 13 and 17 glucose units".
- Enzymes: alpha-amylase (GH13), then amyloglucosidase (GH15). Products by HPSEC, HPAEC, MALDI-TOF-MS.
- Measured (p. 102): saponified (deacetylated) controls digested by alpha-amylase to DP ≤ 6, "whereas oligomers with a
  DP > 6 were found in the digests of all acetylated samples."
- Measured (p. 104–105): after amyloglucosidase, "The smallest enzyme resistant oligomer was trimer"; trimer about 28%,
  35% and 38% of the digests (potato, SuShu2, XuShu18). MALDI: "the enzyme resistant component of DP3 contains 1 acetyl
  group, DP4 contains 1 and 2 acetyl groups, and DP5, DP6 and DP7 contains 2 and 3 acetyl groups".
- Authors' inference (verbatim, p. 105): "the enzyme resistant residues contain at least 1 acetyl group which causes
  the hindrance for α-amylase and amyloglucosidase digestion." The position of the acetyl group inside the oligomers was
  not determined ("the precise location of acetyl with the oligomer will be investigated in future").
- Subsite reading (this report): a DP3 bearing one acetyl that survives exhaustive glucoamylase means glucoamylase
  could not remove its non-reducing-end glucose. With the acetyl on residue k (unknown), at least one of these
  placements is non-productive for GH15: acetylated residue at −1 (k = 3), at +1 (k = 2) or at +2 (k = 1). Which one
  is not resolved by these data.
- They cite Biliaderis (1982): acetyl substituents "exert a similar effect to hydroxyethyl groups which restrict
  α-amylase attack on adjacent α(1→4) glucosidic linkages" (thesis p. 102, paraphrasing Biliaderis).

### 3.2 Huang, Schols, Klaver, Jin, Voragen 2007 — FT-thesis (v2; was NOT READ in v1)
Huang J, Schols HA, Klaver R, Jin Z, Voragen AGJ (2007) Acetyl substitution patterns of amylose and amylopectin
populations in cowpea starch modified with acetic anhydride and vinyl acetate. *Carbohydr Polym* 67(4):542–550.
doi:10.1016/j.carbpol.2006.06.027
Read in: Huang J (2006) *Function-structure relationships of acetylated pea starches*. PhD thesis, Wageningen
University, Chapter 4, pp. 67–84 (same title as the paper). doi:10.18174/42177 (open access, edepot.wur.nl/42177).

- Enzymes (thesis p. 69): alpha-amylase B. subtilis (Fluka 10069, GH13); beta-amylase barley (Fluka 10100, GH14);
  pullulanase B. licheniformis (Megazyme M2, GH13); amyloglucosidase Sigma A9268 (stated "from Aspergillus oryzae",
  GH15). Products by HPSEC, HPAEC, MALDI-TOF-MS. Samples were compared with saponified (deacetylated) controls.
- Measured, Table 1 (thesis p. 72):

| sample | DS amylose | DS amylopectin | beta-limit amylose (%) | beta-limit amylopectin (%) |
|---|---|---|---|---|
| acetic anhydride (AA) | 0.099 | 0.029 | 40 | 55 |
| vinyl acetate (VA) | 0.092 | 0.039 | 37 | 47 |
| saponified | – | – | 86 | 63 |

- Measured, alpha-amylase alone (pp. 73–75, verbatim): "The largest unsubstituted oligomer was DP 8 for both AA AM and
  VA AM; the smallest substituted unit was DP4 with 1 acetyl group for AA AM, and DP3 with 1 acetyl group for VA AM.
  The highest substituted fractions were DP6 with 2 acetyl (DS 0.33), and DP4 with 2 acetyl (DS 0.5)". Amylopectin
  (p. 79): "the smallest substituted unit was DP6 with 1 acetyl group and DP5 with 1 acetyl group" (AA AP, VA AP).
- Measured, alpha-amylase + amyloglucosidase (p. 76, verbatim): "The fragments remaining were DP3 to DP9 with 1 to 4
  acetyl groups, and DP3 to DP8 with 1 to 3 acetyl groups for AA AM and VA AM, respectively. The highest substituted
  fragment was DP3 with 2 acetyl groups (DS 0.67)". Saponified amylose was converted to glucose only.
- Chapter 5 (acetylated yellow pea, p. 97): with alpha-amylase alone "the smallest substituted unit was DP4 with 1
  acetyl group for all four samples". General discussion (p. 116): "significant amounts of dimer were present in the
  digests of acetylated cowpea amylose and amylopectin after combined enzymatic hydrolysis" (substitution of the dimer
  not stated).
- Position of the acetyl group within the oligomers: not determined (p. 116: further study "to identify the glucosyl
  moiety within the oligomer substituted by the acetyl group ... and the precise location of the acetyl groups on the
  glucosyl unit").
- Subsite reading (this report): (i) a DP3 with one acetyl group released by the B. subtilis alpha-amylase alone (VA
  AM) was cut out of the chain by two productive cuts (unless one end is an original chain end); with the acetylated
  residue at position k from the reducing end (k = 1, 2, 3, unknown), the alpha-amylase accepted it at −k in one cut
  and at +(4 − k) in the other, so the acetylated residue was accepted somewhere within −3…+3 in both complexes; (ii)
  the AMG-resistant DP3 with one or two acetyl groups gives the same ambiguous GH15 constraint as Chen 2004 (§3.1).


---

## 4. Octenyl-succinic anhydride (OSA) starch

### 4.1 Bai, Kaufman, Wilson, Shi 2014 — ABS + FT-thesis
Bai Y, Kaufman RC, Wilson JD, Shi YC (2014) Position of modifying groups on starch chains of octenylsuccinic
anhydride-modified waxy maize starch. *Food Chem* 153:193–199. doi:10.1016/j.foodchem.2013.12.012
Read in: Bai Y (2013) *Preparation, structure and properties of octenylsuccinic anhydride modified starch*. PhD
thesis, Kansas State University, Chapter 4, pp. 51–55 (open access, krex.k-state.edu).

- Substrates: granular waxy maize OS starch, DS 0.018 (OS-S-L) and 0.092 (OS-S-H). OS is an ester mainly at O2 (and
  O3/O6) — position distribution not measured in this paper.
- Enzymes: amyloglucosidase (GH15), beta-amylase (GH14), isoamylase (GH13 debranching), successive beta-amylase.
- Measured (abstract + thesis Table 4.1 text):

| Quantity | native | DS 0.018 | DS 0.092 |
|---|---|---|---|
| amyloglucosidase hydrolysis (%) | ~100 (stated as complete for native) | 84.7 | 58.0 |
| beta-limit (%) | 55.9 | 52.8 | 34.4 |
| Mw of beta-limit dextrin | 1.08 × 10⁶ | 1.26 × 10⁶ | 7.04 × 10⁶ |
| branches resistant to isoamylase (%) | 0 | 13.3 | 26.7 |
| beta-amylolysis of debranched material (% maltose) | 100 | ~91 | ~70 |

- Key sentences (thesis p. 52, verbatim): "If the enzyme encounters an OS-substituted glucose units, action of
  amyloglucosidase would stop at the modifying group because the substitution groups interfere with the binding between
  enzyme and starch substrates". Abstract: "OS groups in the OS starch with a low DS were located on the repeat units
  near the branching points, whereas the OS substituents in the OS starch with a high DS occurred both near the
  branching points and the non-reducing ends."
- Authors' inference: OS blocks the exo-enzymes GH15 and GH14 and the debranching isoamylase; positions are inferred
  from where the enzymes stop. No subsite-level resolution.

### 4.2 Bai & Shi 2011 — FT-thesis (thesis Chapter 5)
Bai Y, Shi YC (2011) Structure and preparation of octenyl succinic esters of granular starch, microporous starch and
soluble maltodextrin. *Carbohydr Polym* 83(2):520–527. doi:10.1016/j.carbpol.2010.08.012
- OS maltodextrins: amyloglucosidase conversion 93.3% and 83.5% (granular route, low/high DS), 90.8% and 77.2%
  (solution route); beta-limits 39.3/37.5 and 36.8/34.9% (thesis pp. 73–75).
- Verbatim (thesis p. 74): "amylolysis clearly occurred at glucosidic bonds that were a few anhydroglucose units away
  from the substitution groups."

### 4.3 Whitney, Reuhs, Ovando Martinez, Simsek 2016 — ABS
Whitney K, Reuhs BL, Ovando Martinez M, Simsek S (2016) Analysis of octenylsuccinate rice and tapioca starches:
distribution of octenylsuccinic anhydride groups in starch granules. *Food Chem* 211:608–615.
doi:10.1016/j.foodchem.2016.05.096
- ¹H NMR: "OSA substitution was present only in amylose fractions of rice and tapioca starches" (3% OSA). Not an
  enzyme-subsite study; included because it conflicts with the amylopectin-branch-point location of Bai 2014 (different
  starches and DS).

---

## 5. Cationic starch

### 5.1 Richardson, Nilsson, Cohen, Momcilovic, Brinkmalm, Gorton 2003 — ABS
Richardson S, Nilsson G, Cohen A, Momcilovic D, Brinkmalm G, Gorton L (2003) Enzyme-aided investigation of the
substituent distribution in cationic potato amylopectin starch. *Anal Chem* 75(23):6499–6508. doi:10.1021/ac0300688
- Substrate: cationic (2-hydroxy-3-trimethylammonium propyl ether) potato amylopectin, modified in solution, slurry, dry.
- Enzymes: successive hydrolysis by "different enzymes" (alpha-amylase, amyloglucosidase, debranching; the abstract
  does not list them), products by ESI-MS and MALDI-MS.
- Key sentence (abstract, verbatim): "from the results obtained it was suggested that the enzymes cleave glucosidic
  linkages only between unsubstituted glucose units and, preferentially, linkages in sequences containing more than two
  adjacent unsubstituted units."
- Subsite reading (this report): "only between unsubstituted glucose units" = a cationic residue tolerated neither at
  −1 nor at +1; "preferentially ... more than two adjacent unsubstituted units" = −2 or +2 also mostly unsubstituted.

### 5.2 Wilke & Mischnick 1997 — ABS
Wilke O, Mischnick P (1997) Determination of the substitution pattern of cationic starch ethers. *Starch/Stärke*
49(11):453–458. doi:10.1002/star.19970491106
- O2/O3/O6 ratio by GLC-MS after permethylation and Hofmann elimination; near-random within the AGU.
- Verbatim: "Enzymatic degradation of the cationic starch ethers and measurement of the liberated glucose should give
  indication to the distribution along the polymer chain." (Proposal; no enzyme data in the abstract.)

### 5.3 Manelius, Nurmi, Bertoft 2000 (Cereal Chem) — ABS
Manelius R, Nurmi K, Bertoft E (2000) Enzymatic and acidic hydrolysis of cationized waxy maize starch granules.
*Cereal Chem* 77(3):345–353. doi:10.1094/cchem.2000.77.3.345
- Bacillus amyloliquefaciens alpha-amylase on granules of several DS. "The maximum rate of the enzymatic hydrolysis
  occurred in starches with intermediate DS. It appeared that the cationic substituents interfered with the binding to
  the active site of the enzyme at high levels of substitution." Branches resistant to isoamylase; beta-amylolysis
  resistant dextrins present.

### 5.4 Manelius, Buléon, Nurmi, Bertoft 2000 (Carbohydr Res) — ABS
Manelius R, Buléon A, Nurmi K, Bertoft E (2000) The substitution pattern in cationised and oxidised potato starch
granules. *Carbohydr Res* 329(3):621–633. doi:10.1016/s0008-6215(00)00206-8
- Verbatim: "The starches were partly resistant to the action of isoamylase and the successive beta-amylolysis,
  suggesting that substituents were found both close to the branches and near the nonreducing ends in the amylopectin
  component."
### 5.5 Manelius, Maaheimo, Nurmi, Bertoft 2002 — ABS (v2; was NOT READ in v1)
Manelius R, Maaheimo H, Nurmi K, Bertoft E (2002) Characterisation of fractions obtained by isoamylolysis and ion-exchange
chromatography of cationic waxy maize starch. *Starch/Stärke* 54(2):58–65.
doi:10.1002/1521-379X(200202)54:2<58::AID-STAR58>3.0.CO;2-1 (Crossref also holds a duplicate record,
doi:10.1002/1521-379x(200202)54:2<58::aid-star58>3.3.co;2-t). Abstract read on the VTT research portal
(cris.vtt.fi); full text closed.
- Isoamylase (GH13 debranching) hydrolysates of wet- and dry-cationised waxy maize starch, fractionated on
  CM-Sepharose; GPC and ¹H NMR.
- Measured (abstract, verbatim): bound dextrins "contained increased amounts, from 1.2 - 7.9, of cationic substituents
  per molecule and the degree of polymerisation increased with the density of substituents." "In the latter, the
  debranching was incomplete because of sterical hindrance by substituents at or close to the branch points. Most of
  the dextrins were partly hydrolysed by beta-amylase, but the more highly substituted fractions possessed also a
  population of [beta]-amylase resistant dextrins, suggesting substituents at the non-reducing ends."
- Reading (this report): for isoamylase a cationic group "at or close to" the branch residue blocks; for GH14 a
  cationic residue at the non-reducing end (subsite −2 or −1 of the first productive complex) blocks. Not resolved
  further.

### 5.6 Manelius, Nurmi, Bertoft 2005 — ABS (v2; was NOT READ in v1)
Manelius R, Nurmi K, Bertoft E (2005) Characterization of dextrins obtained by enzymatic treatment of cationic potato
starch. *Starch/Stärke* 57(7):291–300. doi:10.1002/star.200500364 (abstract from OpenAlex; full text closed)
- Enzymes: pullulanase, isoamylase, alpha-amylase; GPC; cation-exchange of the products.
- Measured (abstract, verbatim): "Isoamylase hydrolyzed the native potato starch sample readily, whereas hydrolysis with
  pullulanase was incomplete. Pullulanase hydrolyzed, however, cationic starch with higher DS (degree of substitution)
  more efficiently than isoamylase." "The high portion of uncharged dextrins after α-amylolysis suggested that the
  substitution pattern, on the molecular level, was non-random." "Possibly, the substituents influence substrate
  conformation and thereby alter the hydrolysis patterns."
- Reading (this report): the two GH13 debranching enzymes differ in tolerance of cationic groups near the branch point
  (pullulanase more tolerant than isoamylase at high DS). No subsite-level product structure.

---

## 6. Methylated starch

### 6.1 Steeneken & Woortman 1994 — NOT READ (cited second-hand)
Steeneken PAM, Woortman AJJ (1994) Substitution patterns in methylated starch as studied by enzymic degradation.
*Carbohydr Res* 258:207–221. doi:10.1016/0008-6215(94)84087-3
- No abstract exists in Crossref/Scilit. Second-hand (Bai 2013 thesis p. 74, citation unspecified in that sentence;
  also a US patent): "a minimum sequence length of two unsubstituted glucose residues is required for amylolysis to
  occur for methylated starch". This must be checked in the original before use.
- v2: still NOT READ. OpenAlex: closed; the only other location (TNO repository, resolver.tno.nl
  uuid:784a6c1d-56a3-4045-89c2-3a6e2a56440f) is a metadata record without a file. Second-hand summary in the Huang 2006
  thesis (p. 14, verbatim): Steeneken and Woortman "hydrolyzed methylated potato starch consecutively with α-amylase and
  glucoamylase. Larger amounts of glucose and highly substituted high molecular weight dextrins were obtained from the
  granular form than from the non-granular counterpart with similar DS. It is concluded that the granular form resulted
  in a heterogeneous distribution of the methyl groups within the starch polymers, whereas non-granular starch was
  homogeneously substituted." van der Burgt 1998 (§6.4) cites it for "blockwise substitution in granular methylated
  starches by GPC analysis after α-amylolysis followed by glucoamylolysis". Neither source repeats a subsite statement.

### 6.2 van der Burgt et al. 2000a — FT
van der Burgt YEM, Bergsma J, Bleeker IP, Mijland PJHC, van der Kerk-van Hoof A, Kamerling JP, Vliegenthart JFG (2000)
Distribution of methyl substituents in amylose and amylopectin from methylated potato starches. *Carbohydr Res*
325(3):183–191. doi:10.1016/s0008-6215(99)00329-8 (OA copy: dspace.library.uu.nl/handle/1874/5316)
- Granular potato starch methylated with dimethyl sulfate; glucoamylase (GH15) exhaustive digestion.
- Measured, Table 1 (p. 186), digestibility = glucose released by exhaustive glucoamylase / total glucose:

| sample | MS | gelatinised granule A (%) | amylopectin f1 (%) | amylose f2 (%) |
|---|---|---|---|---|
| P0 | 0 | 100 | 100 | 100 |
| P10 | 0.103 | 66 | 67 | 60 |
| P20 | 0.211 | 56 | 59 | 42 |
| P30 | 0.296 | 43 | 47 | 27 |

- Verbatim (p. 186): "The digestibility of the derivative decreases with increasing MS. Apparently, hydrolysis is
  hindered due to the presence of methyl substituents."
- Worked check: at MS 0.103, if glucoamylase stopped at the first methylated residue of every chain, and methyls were
  random, the released fraction would be far below 66%; the measured value therefore implies that a methyl group
  blocks only some placements (e.g. at −1 or +1, Lemieux et al. 1996, §8.4), not every residue behind it. (Inference,
  this report.)

### 6.3 van der Burgt et al. 2000b — FT
van der Burgt YEM, Bergsma J, Bleeker IP, Mijland PJHC, Kamerling JP, Vliegenthart JFG (2000) Substituent distribution
in highly branched dextrins from methylated starches. *Carbohydr Res* 327:423–429. doi:10.1016/s0008-6215(00)00062-8
(OA copy: dspace.library.uu.nl/handle/1874/5674)
- Enzymes: Bacillus subtilis alpha-amylase (GH13) exhaustive → alpha-limit dextrins (DP > 8, methanol-precipitated);
  then Bacillus cereus beta-amylase (GH14) → alpha,beta-limit dextrins.
- Measured (pp. 424–426): MS of alpha,beta-limit dextrins higher than MS of alpha-limit dextrins; "predominantly
  non-substituted maltose is released"; "in all a,b-limit dextrins the percentage of non-substituted glucose residues at
  the non-reducing end must be higher than 90%. Therefore, the average substitution level of terminal glucose residues
  is lower than 0.15 (MSterminal < 0.15)."
- Authors' inferences (verbatim): "the binding of the maltosyl groups that are subsequently cleaved will be sterically
  hindered by the presence of methyl substituents. Not only branching points stop the digestion, but also substituted
  glucose residues." (p. 424) "This implies that alpha- and beta-amylase both have a diminished ability for hydrolysis
  of the glycosidic linkages of saccharides in which the terminal non-reducing glucosyl group is methylated." (p. 426)
  "detailed information is needed on which partially methylated maltose fragments are tolerated by beta-amylase." (p. 428)
- Subsite reading (this report): beta-amylase releases maltose from subsites −2/−1; a methylated residue at −2 or −1
  blocks release. The finding that the non-reducing termini of the limit dextrins are > 90% unsubstituted means that
  the residue left at the new chain end (which sat in +1 at the last cut) is rarely methylated — i.e. cuts with a
  methylated residue at +1 are also disfavoured for GH14, and for the B. subtilis alpha-amylase that made the termini of
  the alpha-limit dextrins.

### 6.4 van der Burgt et al. 1998 — FT (v2; was NOT READ in v1)
van der Burgt YEM, Bergsma J, Bleeker IP, Mijland PJHC, van der Kerk-van Hoof A, Kamerling JP, Vliegenthart JFG (1998)
Distribution of methyl substituents over branched and linear regions in methylated starches. *Carbohydr Res*
312(4):201–208. doi:10.1016/s0008-6215(98)00238-9 (OA copy: dspace.library.uu.nl/handle/1874/5317; saved as
`Literatura_Referencia/vanderBurgt_1998_methylstarch_alphaamylase_limitdextrins.pdf`)

- **The problem.** Which glucose residues of a methylated starch can sit at subsite −1 when an alpha-amylase cuts it.
- **What the data contain.** Potato starch (P) and amylopectin potato starch (A) methylated in granular suspension with
  dimethyl sulfate, MS 0.04–0.29 (Table 1, p. 202). Monomer composition by methanolysis/GLC, Table 2 (p. 203), mol %:

| sample | MS (GLC) | Glc | 2-Me | 3-Me | 6-Me | 2,3-Me₂ | 2,6-Me₂ | 3,6-Me₂ |
|---|---|---|---|---|---|---|---|---|
| P10 | 0.103 | 90.1 | 6.3 | 1.9 | 1.3 | 0.2 | 0.2 | 0 |
| P20 | 0.211 | 80.7 | 12.2 | 3.3 | 2.0 | 0.9 | 0.7 | 0.2 |
| P30 | 0.296 | 73.7 | 16.4 | 4.2 | 2.4 | 1.7 | 1.3 | 0.3 |
| A20 | 0.191 | 82.5 | 11.0 | 3.1 | 1.8 | 0.8 | 0.6 | 0.2 |
| A30 | 0.293 | 74.0 | 16.2 | 4.2 | 2.3 | 1.8 | 1.2 | 0.3 |

  Verbatim (p. 203): "in all methylated starches the HO-2 group had by far the highest positional degree of
  substitution."
- **The method.** Gelatinised starch digested exhaustively with 0.1 mass% BAN 240L (Novo Nordisk), 21 h, 68.5 °C,
  pH 6.2 (p. 206); the paper calls it "α-amylase from Bacillus subtilis [2]" (BAN is Novo's trade name for a
  Bacillus amyloliquefaciens alpha-amylase — my note; GH13 either way). For P20, P30, A20, A30 the digest
  oligosaccharides were reduced with NaBD₄ (the reducing-end glucose becomes a glucitol), hydrolysed, and the
  monomers analysed by GLC-MS: methylated glucitols = methylated reducing ends; methylated glucononitriles = methylated
  internal or non-reducing-end residues (p. 204, p. 207).
- **What is read off.** Verbatim (p. 204): "It was found that methyl substituents were present mainly in the
  glucononitriles; only a trace of 6-O-methylglucitol was found. It is therefore concluded that α-amylolysis hardly
  occurs at methylated glucose sites; only 6-O-methylation is tolerated." Fig. 1: the average DP of the digests rises
  with MS; "methylation of hydroxyl groups decreases the ability for binding of α-amylase at subsites, thus limiting the
  number of (1→4) linkages that can be cleaved."
- **The checks (worked, this report).** If subsite −1 were indifferent to methylation, the methylated fraction of the
  new reducing ends would equal the methylated fraction of all residues. For P30 that is 100 − 73.7 = 26.3% of
  residues, of which 2-O-methylated (2-Me + 2,3-Me₂ + 2,6-Me₂) = 16.4 + 1.7 + 1.3 = 19.4% and 3-O-methylated
  (3-Me + 2,3-Me₂ + 3,6-Me₂) = 4.2 + 1.7 + 0.3 = 6.2%. Measured: no 2-O- or 3-O-methylglucitol, a trace of
  6-O-methylglucitol (6-Me = 2.4% of residues). The reading does not depend on the enzyme being fully exhaustive.
- **What it does not include.** Only the reducing ends (−1 side) were analysed; the non-reducing ends (+1 side) of the
  products were not. "Trace" is not quantified. The digest includes branched dextrins whose reducing ends were also
  made by the alpha-amylase.
- **Subsite reading (this report).** For this GH13 liquefying alpha-amylase: 2-O-methyl and 3-O-methyl glucose are not
  found at subsite −1 of productive complexes; 6-O-methyl glucose is found there rarely (the only methylated residue
  seen at −1).
- Other results: MS of branched regions 60–140% higher than linear regions (P series) and 140–260% higher (A series)
  (Table 3, p. 205); FABMS of trimers/tetramers showed heterogeneous substitution, more so in linear regions (Table 4).

### 6.5 Mischnick 2001 — ABS (new in v2)
Mischnick P (2001) Specificity of microbial α-amylase and amyloglucosidase hydrolysis of methyl amyloses.
*Starch/Stärke* 53(3–4):110–120. doi:10.1002/1521-379x(200104)53:3/4<110::aid-star110>3.0.co;2-y (abstract from
OpenAlex; full text closed)
- Enzymes: alpha-amylase from B. licheniformis (GH13) and amyloglucosidase from A. niger (GH15), exhaustive digestion
  of methyl amyloses made under homogeneous and heterogeneous conditions. Products analysed by FAB- and MALDI-MS after
  reduction and O-methylation/O-ethylation, "differentiating between non-reducing terminal glucosyl residues, 1→4-linked
  inner glucosyl units and reducing glucose end groups".
- Key sentences (abstract, verbatim): "the reducing glucose end group was usually unsubstituted, while the non-reducing
  glucosyl residues could be 2-, 6- or even 2,6-di-O-substituted. In contrast, all 3-O-methyl groups were located in
  the inner 1,4-linked glucosyl units, indicating inhibition of both enzymes."
- Second-hand from the full text (Huang 2006 thesis, p. 116, verbatim): Mischnick "found that the position of
  O-methylation was an important structural feature for the action of the α-amylase and amyloglucosidase, showing
  inhibition in the following order: 3-O-methylation > 2-O-methylation > 6-O-methylation."
- Subsite reading (this report): (i) the reducing end of a product was at −1 when it was made: methylated residues are
  "usually" not at −1 of the B. licheniformis alpha-amylase (agrees with §6.4 for the other Bacillus enzyme). (ii) The
  non-reducing terminal residue of a product was at +1 of the last productive complex that made that end (alpha-amylase
  cut, or glucoamylase removal of the residue before it); since AMG then stopped, it is also a residue AMG could not
  accept at −1. So 2-O-, 6-O- and 2,6-di-O-methyl glucose are accepted at +1 (of GH13 and/or GH15 — the abstract does
  not separate the two enzymes) and are not accepted at −1 of GH15 (agrees with Lemieux 1996, §8.4). (iii) 3-O-methyl
  glucose was never terminal or reducing: not accepted at −1 or at +1 of either enzyme (agrees with the GH15 requirement
  for O3 at +1, Bock & Pedersen 1987, §8.5).

---

## 7. Phosphorylated starch (natural C-6 and C-3 phosphate, potato) — the best subsite-resolved data set

### 7.1 Takeda 1987 (award review of the Kagoshima work) — FT (Japanese text, English abstract, Fig. 1)
Takeda Y (1987) Actions of various amylases in the vicinity of phosphate groups of potato starch and location of the
phosphate groups in the starch. *J Jpn Soc Starch Sci* 34(3):225–233. doi:10.5458/jag1972.34.225 (OA, J-STAGE)

- Substrate: potato starch (563 ppm bound P; C-3 : C-6 phosphate ≈ 1 : 3, p. 226) for endo-enzymes; phosphodextrin
  (debranched potato amylopectin, DP ≈ 40, ~1 phosphate per molecule, non-reducing ends unphosphorylated) for exo-enzymes.
- Fig. 1 (p. 227), limit products, non-reducing end on the left (Ø = reducing glucose):

| enzyme (family) | product with 3-phosphate | product with 6-phosphate |
|---|---|---|
| sweet-potato beta-amylase (GH14) | (stops ≥ 1 glucosyl before the 3-P residue; text) | chains ending P6-G-G… or G-P6-G… (odd/even) |
| Aspergillus niger glucoamylase (GH15) | G-(P3)G-G-G… (one Glc left on the non-reducing side) | (P6)G-G-G… (6-P residue left terminal) |
| Bacillus subtilis alpha-amylase (GH13) | 3²-PG4 = G-G-(P3)G-Ø | 6²-PG3 = G-(P6)G-Ø |
| porcine pancreatic alpha-amylase (GH13) | 3²-PG4 = G-G-(P3)G-Ø | 6³-PG3 = (P6)G-G-Ø |
| Aspergillus oryzae alpha-amylase (GH13) | 3³-PG4 = G-(P3)G-G-Ø | 6³-PG3 = (P6)G-G-Ø |

  (The English abstract lists A. oryzae as giving "3²-PG4 and 6³-PG3"; Fig. 1 and the 1986 abstract, §7.4, both show
  the 3-phosphate on the second residue from the non-reducing end, i.e. 3³-PG4. Fig. 1 is followed here.)
- Key sentences (English abstract, verbatim): "Sweet potato β-amylase and Aspergillus niger glucoamylase cleaved the
  (1→4)-α-D-linkage of the 6-phosphorylated residue at the non-reducing side but not that of the 3-phosphorylated one,
  leaving at least one glucosyl residue attached to the 3-phosphorylated one." "The phosphate group at C-6 was as
  obstructive as the branch linkage of amylopectin, but less obstructive than the phosphate group at C-3. Their
  inhibitory behavior seemed to reflect a bulk rather than an electronegative effect."
- Text p. 227 (Japanese, my translation): for 6-P chains beta-amylase stops at bond Ā when the number of glucose
  residues on the non-reducing side is even and at bond B̄ when odd; glucoamylase degrades 6-P chains up to bond Ā
  (the bond on the non-reducing side of the phosphorylated residue) and stops at bond B̄ for 3-P chains, releasing one
  glucose fewer; "the 3-phosphate is a larger obstacle than the 6-phosphate by one glucose residue".

**Subsite reading (this report), from Fig. 1 and the rules of §0:**

| enzyme | 6-P tolerated at | 6-P not productive at | 3-P tolerated at | 3-P not productive at |
|---|---|---|---|---|
| A. niger glucoamylase (GH15) | +1 | −1 | — (+2 and beyond, implied) | +1 (and −1) |
| sweet-potato beta-amylase (GH14) | +1 | −1, −2 | — | +1 (and −1, −2) |
| B. subtilis alpha-amylase (GH13) | −2, +2 | −1, +1 | −2, +3 | −1, +1, +2 |
| porcine pancreatic alpha-amylase (GH13) | −3, +1 | −1, −2 | −2, +3 | −1, +1, +2 |
| A. oryzae alpha-amylase (GH13) | −3, +1 | −1, −2 | −3, +2 | −2, −1, +1 |

  Caveats: (i) "not productive" comes from linkages that remain uncut in short limit products; (ii) a product end could
  be an original chain end, which removes the "+" entries for that product; (iii) +n entries assume the enzyme has at
  least n aglycone subsites. The striking points are consistent with the authors' text: no enzyme puts a phosphorylated
  residue at −1; the 3-phosphate is excluded from more subsites than the 6-phosphate; and the three GH13 enzymes differ
  in where they accept the 6-P (PPA and TAKA at +1 and −3; B. subtilis at ±2) and the 3-P (TAKA at +2/−3; the others at
  −2/+3).

### 7.2 Takeda, Hizukuri, Ozono, Suetake 1983 — ABS
Takeda Y, Hizukuri S, Ozono Y, Suetake M (1983) Actions of porcine pancreatic and Bacillus subtilis α-amylases and
Aspergillus niger glucoamylase on phosphorylated (1→4)-α-D-glucan. *Biochim Biophys Acta* 749(3):302–311.
doi:10.1016/0167-4838(83)90240-6
- Primary paper for the PPA and B. subtilis rows above. Verbatim: "the two alpha-amylases cleave different sites in the
  vicinity of the phosphate group at C-6 of the glucosyl residue but the same site in the vicinity of that at C-3."
  Glucoamylase cleaves on the non-reducing side of the 6-phosphoryl residue "but not that of the 3-phosphoryl residue,
  leaving a glucosyl residue attached to the 3-phosphoryl residue. These results indicate that the phosphate group at
  C-3 is more obstructive than that at C-6 for the actions of these three amylases."

### 7.3 Abe, Takeda, Hizukuri 1982 — ABS
Abe J, Takeda Y, Hizukuri S (1982) Action of glucoamylase from Aspergillus niger on phosphorylated substrate. *Biochim
Biophys Acta* 703(1):26–33. doi:10.1016/0167-4838(82)90006-1
- Glucoamylase purified free of alpha-amylase and phosphatase. Verbatim: "The glucoamylase action was prevented by the
  esterified phosphate groups of the substrate. Thus, the extensive action of the glucoamylase on potato starch exposed
  the 6-posphorylglucosyl residue of the starch at the non-reducing terminal and large molecular weight limit dextrins
  remained. The concomitant action of the phosphatase was necessary for the complete degradation of the starch."
- Subsite reading: 6-phosphate at −1 not tolerated by GH15 (A. niger).

### 7.4 Takeda & Hizukuri 1986 — ABS
Takeda Y, Hizukuri S (1986) Actions of Aspergillus oryzae alpha-amylase, potato phosphorylase, and rabbit muscle
phosphorylase a and b on phosphorylated (1→4)-α-D-glucan. *Carbohydr Res* 153(2):295–307.
doi:10.1016/s0008-6215(00)90271-4
- A. oryzae alpha-amylase gave 6³-phosphorylmaltotriose and a maltotetraose with the 3-phosphate on the second residue
  from the non-reducing end. Verbatim: "the enzyme hydrolyses the same linkages in the vicinity of the 6-phosphorylated
  residue as porcine-pancreatic alpha-amylase, but hydrolyses different linkages in the vicinity of the
  3-phosphorylated residue when compared with B. subtilis and pancreatic alpha-amylases." Phosphorylases (GT35) "were
  unable to by-pass the phosphorylated D-glucosyl residue".

### 7.5 Takeda & Hizukuri 1981; 1982 — NOT READ (content via §7.1)
Takeda Y, Hizukuri S (1981) Re-examination of the action of sweet-potato beta-amylase on phosphorylated (1→4)-α-D-glucan.
*Carbohydr Res* 89(1):174–178. doi:10.1016/s0008-6215(00)85242-8
Takeda Y, Hizukuri S (1982) Location of phosphate groups in potato amylopectin. *Carbohydr Res* 102(1):321–327.
doi:10.1016/s0008-6215(00)88077-5

### 7.6 Kamasaka et al. 1997 — ABS + p. 238 read (scanned)
Kamasaka H, To-o K, Kusaka K, Kuriki T, Kometani T, Hayashi H, Okada S (1997) The structures of phosphoryl
oligosaccharides prepared from potato starch. *Biosci Biotechnol Biochem* 61(2):238–244. doi:10.1271/bbb.61.238
- Enzymes: bacterial liquefying alpha-amylase (BLA, GH13), glucoamylase (GA, GH15), pullulanase; then bacterial
  saccharifying alpha-amylase (BSA) after GA for analysis.
- Measured (abstract, verbatim): PO-1 "was made up of 3-phosphoryl oligosaccharides (3³-phosphoryl maltotetraose and
  3⁴-phosphoryl maltopentaose) and 6-phosphoryl oligosaccharides (6³-phosphoryl maltotriose, 6²-phosphoryl maltotriose,
  6³-phosphoryl maltotetraose, and 6⁴-phosphoryl maltopentaose)." "A small amount of 6²-phosphoryl maltose was also
  detected after treating with GA alone, indicating that 6²-phosphoryl maltotriose existed in the PO-1 fraction."
- Subsite reading (this report): every 3-P product carries exactly one glucose on the non-reducing side of the 3-P
  residue (glucoamylase limit, as in §7.1), and two or three residues on the reducing side — i.e. the liquefying
  alpha-amylase cut with the 3-P residue at −3 or −4, never at −2 in this fraction. GA cleaving G-(P6)G-G to (P6)G-G
  confirms 6-P tolerated at +1 by GH15.

### 7.7 Hansen et al. 2009 — ABS (conformation of the linkage next to a phosphate)
Hansen PI, Spraul M, Dvortsak P, Larsen FH, Blennow A, Motawia MS, Engelsen SB (2009) Starch phosphorylation — maltosidic
restrains upon 3′- and 6′-phosphorylation investigated by chemical synthesis, molecular dynamics and NMR spectroscopy.
*Biopolymers* 91(3):179–193. doi:10.1002/bip.21111
- Verbatim: "When maltose is phosphorylated in the 3′-position, MD and NMR show that the glucosidic space is seriously
  restricted to one narrow potential energy well which is strongly offset from the global potential energy well of
  maltose and almost 50 degrees from the Phi angle of the alpha-maltose crystal structure."
- Relevance: one physical reason why the 3-phosphate obstructs more subsites than the 6-phosphate (§7.1) besides direct
  contact with the enzyme: it changes the φ/ψ of the adjacent alpha(1→4) linkage.

### 7.8 Wikman et al. 2011 — ABS (location context)
Wikman J, Larsen FH, Motawia MS, Blennow A, Bertoft E (2011) Phosphate esters in amylopectin clusters of potato tuber
starch. *Int J Biol Macromol* 48(4):639–649. doi:10.1016/j.ijbiomac.2011.02.005

---

## 8. Specifically modified maltooligosaccharides and glycosides used to map subsites

### 8.1 Omichi, Hase, Ikenaka 1992 — ABS
Omichi K, Hase S, Ikenaka T (1992) Studies on the active site of human α-amylases: examination of the third subsite S3′
of the aglycone-binding site by control of substrate binding mode. *J Biochem* 111(1):4–7.
doi:10.1093/oxfordjournals.jbchem.a123716
- Substrates: IG-G-G-G-AG-G, IG-G-G-G-AG-M, IG-G-G-G-AG-φ (IG = 6-deoxy-6-iodo-Glc at the non-reducing end; AG =
  6-amino-6-deoxy-Glc; M methyl; φ phenyl).
- Measured: human salivary and pancreatic alpha-amylases (GH13) gave IG-G-G + G-AG-G (and analogues). Km for the G and M
  substrates equal, twice that of the φ substrate.
- Subsite reading (from the products, this report): IG at −3, cut between −1 and +1, AG at +2. So a 6-iodo-6-deoxy
  residue is accepted at −3 and a 6-amino-6-deoxy residue at +2. Authors: S3′ "has no affinity for the glucose residue
  and is not a subsite but a hydrophobic environment."

### 8.2 Nagamine, Omichi, Ikenaka 1988 — ABS
Nagamine Y, Omichi K, Ikenaka T (1988) A comparison of the modes of actions of human salivary and pancreatic
α-amylases on modified maltooligosaccharides. *J Biochem* 104(4):667–670. doi:10.1093/oxfordjournals.jbchem.a122529
- Phenyl alpha-maltopentaoside/-tetraoside with iodo, amino or carboxyl at the first or penultimate residue from the
  non-reducing end. Verbatim: results "suggested the presence of five subsites (S3, S2, S1, S1′, and S2′) and a
  hydrophobic amino acid residue at subsite S3 in the active site of HPA".

### 8.3 Omichi, Shiosaki, Matsubara, Ikenaka 1989 — ABS
Omichi K, Shiosaki K, Matsubara K, Ikenaka T (1989) Actions of three human α-amylases expressed in yeast on modified
substrates and the amino acid residues causing their different actions. *J Biochem* 106(4):646–650.
doi:10.1093/oxfordjournals.jbchem.a122910
- Substrates G5φ with the non-reducing-end C6 converted to CH2NH2 (AG5φ), COOH (CG5φ) or CH2I (IG5φ). AG5φ and CG5φ gave
  more Gφ and less G2φ; IG5φ more G2φ. Verbatim: the residues responsible "are located in subsite S3, and interact with
  the CH2OH residue of the non-reducing-end glucose residue of G5 phi." A C6 modification at −3 shifts the binding
  register rather than blocking hydrolysis.
- Also: Omichi K, Hase S, Ikenaka T (1990) Modified malto-oligosaccharides as inhibitors of human alpha-amylases.
  *Carbohydr Res* 208:312–316. doi:10.1016/0008-6215(90)80118-m — no abstract in PubMed; NOT READ.

### 8.4 Lemieux, Spohr, Bach, Cameron, Palcic, Frandsen, Stoffer, Svensson 1996 — ABS
Lemieux RU, Spohr U, Bach M, Cameron DR, Palcic MM, Frandsen TP, Stoffer BB, Svensson B (1996) Chemical mapping of the
active site of the glucoamylase of Aspergillus niger. *Can J Chem* 74(3):319–335. doi:10.1139/v96-036
- Substrates: methyl alpha-isomaltoside with each of the seven OH groups either deoxygenated or O-methylated. Primed (′)
  = non-reducing glucosyl (subsite −1); unprimed = reducing unit (subsite +1). Enzyme GH15.
- Verbatim: "both deoxygenation and O-methylation of OH-4 (reducing unit), OH-4′, or OH-6′ strongly hindered
  hydrolysis, whereas the kinetics were virtually unaffected by either the substitutions at OH-2 or structural changes
  at C-1. The substitutions at OH-3 caused increases of 2.1 and 1.9 kcal/mol in the ΔΔG≠. In contrast, whereas
  deoxygenation of either OH-2′ or OH-3′ caused much smaller (0.96 and 0.52 kcal/mol) increases in ΔΔG≠, the
  mono-O-methylations resulted in severe steric hindrance to the formation of the activated complex."
- Subsite reading: at −1, O-methyl at O2, O3, O4 or O6 not tolerated; at +1, O-methyl at O2 tolerated, at O3 costs
  1.9 kcal/mol, at O4 strongly hindering. Caveat: these +1 numbers are for the alpha(1→6)-linked isomaltoside, whose
  reducing unit sits in +1 with a different orientation from the alpha(1→4) maltoside; the authors argue the docked
  model represents both, but transfer to maltose-type substrates is an assumption.
- Relation used for ΔΔG≠: ΔΔG≠ = −RT ln[(kcat/KM)analogue / (kcat/KM)parent]. Worked: 1.9 kcal/mol at 298 K
  (RT = 0.593 kcal/mol) means (kcat/KM) reduced by exp(1.9/0.593) = exp(3.20) ≈ 25-fold.

### 8.5 Bock & Pedersen 1987 (Part I) — FT (v2; was ABS in v1)
Bock K, Pedersen H (1987) The substrate specificity of the enzyme amyloglucosidase (AMG). Part I. Deoxy derivatives.
*Acta Chem Scand B* 41:617–628. doi:10.3891/acta.chem.scand.41b-0617
Source: the journal's free back-issue archive (actachemscand.org/pdf/acta_vol_41b_p0617-0628.pdf); the site did not
answer from here on 2026-11-24, so the PDF was taken from the Internet Archive capture of that same publisher URL
(2023-11-19). Saved as `Literatura_Referencia/Bock_1987_AMG_deoxy_maltosides.pdf` (scanned; read from page images).

- Notation: methyl β-maltoside; primed = non-reducing glucosyl unit (subsite −1); unprimed = reducing unit carrying the
  β-OMe (subsite +1). Enzyme: AMG (A. niger glucoamylase, GH15), gift from Novo; 0.1 M acetate pH 4.5, 27 °C.
- Table 3 (p. 622), qualitative (identity of each compound from Table 1 ¹H NMR, my reading of the deoxy positions):

| compound | modification | subsite | substrate? |
|---|---|---|---|
| 1a | 6-deoxy | +1 | yes |
| 2a | 6′-deoxy | −1 | no |
| 3a | 4-O-α-Glc-1,5-anhydro-D-glucitol (no anomeric OMe, "1-deoxy") | +1 | yes |
| 5a | 3-deoxy | +1 | no |
| 9a | 2′-deoxy | −1 | yes |
| 11a | 3′-deoxy | −1 | yes |
| 14a | 4′-deoxy | −1 | no |
| 15a | 6,6′-dideoxy | −1, +1 | "very poor substrate on prolonged exposure" |
| 16a | 2-deoxy analogue of 3a (1,5-anhydro-2-deoxy-glucitol unit) | +1 | yes |

- Table 4 (p. 622), ¹H-NMR progress-curve kinetics, 50 µg enzyme per experiment:

| compound | S₀ (mM) | K_M (mM) | V_M (µmol min⁻¹) | k₂/K_M (M⁻¹ s⁻¹) |
|---|---|---|---|---|
| methyl β-maltoside | 5.98 | 0.33 | 0.13 | 9.7 × 10³ |
| methyl β-maltoside | 31.2 | 0.52 | 0.30 | 1.4 × 10⁴ |
| 1a (6-deoxy, +1) | 7.23 | 0.36 | 0.14 | 9.0 × 10³ |
| 3a (1,5-anhydroglucitol at +1) | 17.8 | 5.2 | 0.89 | 4.2 × 10³ |
| 9a (2′-deoxy, −1) | 7.27 | 3.6 | 0.39 | 2.5 × 10³ |
| 16a (2-deoxy anhydroglucitol at +1) | 22.3 | 9.1 | 0.47 | 1.2 × 10³ |

- Key sentence (abstract, verbatim): "three hydroxy groups (3, 4′ and 6′) are essential for the compounds to act as
  substrate for the enzyme." Text (p. 622): removal of the other hydroxy groups "generally gives substrates which are
  hydrolysed faster than methyl maltoside, due to weaker association to the enzyme, but which have transition states of
  comparable energy."
- Worked example: 2′-deoxy vs parent (S₀ ≈ 6–7 mM): k₂/K_M ratio = 2.5 × 10³ / 9.7 × 10³ = 0.26, i.e.
  ΔΔG = −RT ln 0.26 = 0.593 × 1.35 = 0.80 kcal/mol at 298 K (the experiment was at 300 K; the difference is < 1%).
- Subsite reading: −1 requires O4′ and O6′; O2′ and O3′ dispensable. +1 requires O3; O6 and O2 (in the anhydroglucitol
  series) and the anomeric substituent dispensable (with a 16-fold K_M penalty when the anomeric OMe is absent:
  5.2 / 0.33).

### 8.5b Bock & Pedersen 1988 (Part II, 6-substituted maltosides) — FT (v2; was NOT READ in v1)
Bock K, Pedersen H (1988) The substrate specificity of the enzyme amyloglucosidase (AMG). Part II. 6-Substituted maltose
derivatives. *Acta Chem Scand B* 42:75–85. doi:10.3891/acta.chem.scand.42b-0075
Source: same publisher archive, Internet Archive capture (2023-12-02) of actachemscand.org/pdf/acta_vol_42b_p0075-0085.pdf.
Saved as `Literatura_Referencia/Bock_1988_AMG_6-substituted_maltosides.pdf` (scanned; read from page images).

- **The problem.** Whether the 6-position of the residue at subsite +1, and of the residue at −1, can carry something
  other than OH in a glucoamylase substrate.
- **What the data contain.** Methyl β-maltosides modified at C6 of the reducing unit (+1) or at C6′ (−1), structures in
  the schemes on pp. 76–77; qualitative substrate test (Table 1, p. 78); initial rates and half-lives by NMR
  (Table 2, p. 78; 27 °C, pH 4.3, 5.5 × 10⁻¹⁰ M AMG).
- Table 1 (p. 78):

| compound | modification | subsite | substrate? |
|---|---|---|---|
| 3b | 1,6-anhydro (reducing unit) | +1 | no |
| 5b | 6-chloro-6-deoxy | +1 | yes |
| 6b | 6-bromo-6-deoxy | +1 | yes |
| 7b | 6-deoxy-6-iodo | +1 | yes |
| 8b | 6-azido-6-deoxy | +1 | yes |
| 9 | 6-amino-6-deoxy | +1 | no |
| 10 | 6-acetamido-6-deoxy | +1 | yes |
| 11b | 6-deoxy-xylo-hex-5-enopyranoside (5,6-ene) | +1 | yes |
| 12b | 6-deoxy-6-fluoro | +1 | yes |
| 13b | uronic acid (6-COOH) | +1 | no |
| 14b | uronic acid methyl ester (6-COOMe) | +1 | yes |
| 15b | 6′-bromo-6′-deoxy | −1 | no |
| 16b | 3′,6′-anhydro | −1 | no |
| 17b | 6′-deoxy-6′-fluoro | −1 | no |

- Table 2 (p. 78), initial rate v₀ and half-life (pairs below the first block were measured in the same tube, i.e. in
  competition):

| substrate | conc. (mM) | v₀ (mM min⁻¹) | t½ (min) |
|---|---|---|---|
| methyl 6-chloro-6-deoxy-β-maltoside | 7.2 | 0.19 | 254 |
| methyl 6-deoxy-β-maltoside | 8.3 | 4.07 | 12 |
| methyl β-maltoside | 7.3 | 1.85 | 27 |
| methyl 3′-deoxy-β-maltoside | 8.3 | 4.16 | 13 |
| 6-chloro (+ maltoside 3.8 mM) | 3.6 | 0.20 (maltoside 0.32) | 215 (127) |
| 3′-deoxy (+ maltoside 3.8 mM) | 4.2 | 2.79 (maltoside 1.11) | 19 (42) |
| 6-fluoro (+ maltoside 3.8 mM) | 4.0 | 0.45 (maltoside 0.47) | 106 (102) |
| 3-deoxy (+ maltoside 3.8 mM) | 4.2 | 0.00 (maltoside 4.90) | – (9) |

- **What is read off (authors, verbatim, p. 78).** "the enzyme can accept a wide range of modifications in the
  6-position provided that the substituents are not charged groups. If this is the case, as in the 6-amino and
  6-carboxy derivatives, then the enzyme cannot hydrolyse the compounds, indicating that the active site accommodates
  the substrates with the 6-position in a hydrophobic environment." The 6-chloro and 6-fluoro compounds "are bound very
  tightly to the active site of the enzyme since they show substantial inhibition of the reaction of the normal
  substrate. This is not the case for the 6′-bromo or 6′-fluoro compounds, which are neither substrates nor inhibitors".
- **The checks (worked, this report).** Competition: with maltoside alone at 3.8 mM, v₀ = 1.11–4.90 mM min⁻¹ in the
  other pairs; with 3.6 mM 6-chloro present the maltoside v₀ drops to 0.32 and its t½ rises to 127 min, which is the
  inhibition the authors report. The 6-fluoro compound and the maltoside are cleaved at the same rate (0.45 vs 0.47).
- **What it does not include.** One enzyme (A. niger AMG); disaccharide substrates only (no +2 subsite occupied);
  qualitative test for most compounds; the 3-deoxy row shows the 3-deoxy compound is not cleaved, not that it binds.
- **Subsite reading.** +1 of GH15: C6 substituents Cl, Br, I, N₃, F, NHAc, COOMe and the 5,6-ene are accepted; a
  charged C6 (CH₂NH₃⁺, COO⁻) and the 1,6-anhydro ring are not. −1: every C6′ change tested (F, Br, H, 3′,6′-anhydro)
  abolishes hydrolysis. Compare §7.1: in the polymer, A. niger glucoamylase does accept a 6-phosphate (an anionic
  ester on O6) at +1; the 6-carboxylate is rejected in the disaccharide. The two observations differ in the
  substituent (O6 retained and extended vs CH₂OH replaced by COO⁻) and in the substrate length (polymer with +2, +3
  occupied vs a methyl maltoside); which factor explains the difference is not tested in either paper.

### 8.6 Sierks, Bock, Refn, Svensson 1992 — ABS
Sierks MR, Bock K, Refn S, Svensson B (1992) Active site similarities of glucose dehydrogenase, glucose oxidase, and
glucoamylase probed by deoxygenated substrates. *Biochemistry* 31(37):8972–8977. doi:10.1021/bi00152a038
- Verbatim: activity "found to depend heavily on the 4- and 6-OH groups (4′- and 6′-OH in maltose), where changes in
  binding energies from 10 to 18 kJ/mol suggested strong hydrogen bonds"; "The 3-OH (3′-OH in maltose) was involved in
  weaker interactions, while the 2-OH (2′-OH in maltose) had a very small if any role in transition-state binding."

### 8.7 Palcic, Skrydstrup, Bock, Le, Lemieux 1993 — ABS
Palcic MM, Skrydstrup T, Bock K, Le N, Lemieux RU (1993) Substrate recognition by amyloglucosidase: evaluation of
conformationally biased isomaltosides. *Carbohydr Res* 250(1):87–92. doi:10.1016/0008-6215(93)84157-2
- Methyl beta-maltoside hydrolysed 30–50 times faster than methyl alpha-isomaltoside; OH-6′, OH-4′, OH-4 in key polar
  contacts. 6R-C-methyl isomaltoside Km 0.9 mM vs 24.5 mM parent vs 90 mM 6S; kcat 0.85 ± 0.20 s⁻¹ for all.

### 8.8 Ogawa et al. 2004 (GH31) — FT
Ogawa M, Nishio T, Hakamata W, Matsuishi Y, Hoshino S, Kondo A, Kitagawa M, Kawachi R, Oku T (2004) Substrate hydroxyl
groups are involved in the ionization of catalytic carboxyl groups of Aspergillus niger α-glucosidase. *J Appl
Glycosci* 51(1):9–14. doi:10.5458/jag.51.9 (OA)
- Enzyme: A. niger alpha-glucosidase, stated as GH31 (p. 10).
- Verbatim (p. 9–10): chosen "because of its strong hydrolytic activity against 2- and 3-deoxy glucosides modified at
  the glycon portion, although it possesses no activity toward 4- and 6-deoxy glucosides" (their ref. 8).
- Subsite reading: at −1 of this GH31 enzyme, O4 and O6 are required, O2 and O3 are not required for turnover (O2/O3
  affect the pKa of the catalytic carboxylates, pH optimum 4.4 → 5.6 on deoxygenation). No data here on bulky
  substituents at O2/O3.

### 8.9 Kandra and co-workers — ABS
Kandra L, Gyémánt G, Farkas E, Lipták A (1997) Action pattern of porcine pancreatic alpha-amylase on three different
series of β-maltooligosaccharide glycosides. *Carbohydr Res* 298(3):237–242. doi:10.1016/s0008-6215(96)00310-2
- 4,6-O-benzylidene at the non-reducing end: "a favourable interaction between the Bnl group and the subsites (-3) and
  (-5) but an unfavourable one with subsite (-4)". Also: "the glucopyranose residue could be replaced by the aglycon
  group" (2-chloro-4-nitrophenyl/4-nitrophenyl accepted at aglycone subsites).
Kandra L, Gyémánt G (2000) Examination of the active sites of human salivary α-amylase (HSA). *Carbohydr Res*
329(3):579–585. doi:10.1016/s0008-6215(00)00221-4
- At least six subsites (−4 … +2); benzylidene "favourable" at −3, "unfavourable" at −4; "subsite (-4) must be more
  hydrophylic than hydrophobic."
Kandra L, Gyémánt G, Remenyik J, Hovánszki G, Lipták A (2002) Action pattern and subsite mapping of Bacillus
licheniformis α-amylase (BLA) with modified maltooligosaccharide substrates. *FEBS Lett* 518(1-3):79–82.
doi:10.1016/s0014-5793(02)02649-2
- "The binding region of BLA is composed of five glycone, three aglycone-binding sites and a 'barrier' subsite."

### 8.10 Okada, Ogawa, Usui 2000 — FT (new in v2)
Okada M, Ogawa K, Usui T (2000) Substrate specificity of human and bovine pancreatic α-amylases using end-blocked
maltooligosaccharides. *J Appl Glycosci* 47(1):35–44. doi:10.5458/jag.47.35 (OA, J-STAGE; saved as
`Literatura_Referencia/Okada_2000_pancreatic_amylase_end-blocked_maltooligosaccharides.pdf`). The enzymes actually used
are porcine pancreatic (PPA) and human pancreatic (HPA) alpha-amylase (GH13), as the figures and tables state.

- **The problem.** Which subsites of PPA and HPA accept a glucose residue that carries a sugar on O6, on O3 or on O4.
- **What the data contain.** p-Nitrophenyl α-maltooligosides GnP (n = 3–5 Glc, P = p-nitrophenyl at the reducing end)
  whose non-reducing-end Glc carries β-D-Gal on O4 (Gal⁴GnP), β-D-Gal on O6 (Gal⁶GnP) or β-D-GlcNAc on O3 (NGnP).
  Cleavage frequency of each bond in the initial stage (≤ 15% hydrolysis), HPLC of the p-nitrophenyl products (Figs. 1
  and 2, p. 41); k₀ and K_m (Tables 3 and 4, pp. 41–42).
- **Subsite convention.** The authors use five subsites 1–5 with the catalytic site between 3 and 4 (Fig. 3, p. 42); in
  the Davies notation 1, 2, 3 = −3, −2, −1 and 4, 5 = +1, +2. In Figs. 1–2 every product is drawn against the
  catalytic line, so the subsite of the substituted residue is read directly.
- **Cleavage frequencies (Figs. 1, 2), with the subsite of the substituted Glc (this report's translation):**

| substrate | PPA: substituted Glc at (frequency) | HPA: substituted Glc at (frequency) |
|---|---|---|
| Gal⁶G5P | −3 (0.61), −4 (0.39) | −4 (0.83), −3 (0.17) |
| Gal⁶G4P | −3 (0.96), −4 (0.03), −2 (0.01) | −3 (0.91), −4 (0.06), −2 (0.03) |
| Gal⁶G3P | not a substrate | not a substrate |
| NG5P (3-O-GlcNAc) | −3 (0.96), −4 (0.04) | −3 (0.96), −4 (0.04) |
| NG4P | −3 (1.0) | −3 (1.0) |
| NG3P | hydrolysed "very slowly", GP : phenol = 0.40 : 0.60 (text p. 40) | same tendency |
| Gal⁴G4P (Gal as a fifth residue) | Gal at −3 (0.97), at −4 (0.03) | Gal at −3 (0.88), −4 (0.12) |
| Gal⁴G3P | Gal at −3 (1.0) | Gal at −3 (1.0) |

  NG3P: GP release means the 3-O-GlcNAc-Glc sat at −2; phenol release means it sat at −3 with p-nitrophenyl at +1.
  For Gal⁴GnP the Gal itself occupies a subsite (it is the terminal residue), so its subsite is given.
- **Kinetics (Tables 3, 4):**

| substrate | PPA k₀/K_m (s⁻¹ mM⁻¹) | HPA k₀/K_m (s⁻¹ mM⁻¹) |
|---|---|---|
| G6P / G5P / G4P / G3P | 6500 / 4400 / 1300 / 12 | 7600 / 3500 / 1400 / 11 |
| Gal⁶G5P | 1300 (K_m 0.27 mM, k₀ 352 s⁻¹) | 530 (K_m 0.76, k₀ 405) |
| Gal⁶G4P | 54 (K_m 1.60, k₀ 86) | 7 (K_m 1.80, k₀ 13) |
| NG5P / NG4P / NG3P | 3700 / 870 / 2 | 3900 / 730 / 2 |
| Gal⁴G5P / Gal⁴G4P / Gal⁴G3P | 5900 / 3500 / 820 | 6900 / 3300 / 480 |

  Worked: Gal⁶G4P vs G4P, PPA 54 / 1300 = 0.042, HPA 7 / 1400 = 0.005 (the authors: "1.2 and 0.2% of those of G5P",
  the comparison they make with the parent of equal total length); ΔΔG = −RT ln(0.042) = 0.593 × 3.17 = 1.9 kcal/mol
  (PPA) and 0.593 × 5.30 = 3.1 kcal/mol (HPA) at 298 K, taking G4P as reference.
- **Authors' inference (verbatim, p. 42).** "subsite 2 was, in contrast to subsite 1, strict for binding the glycosyl
  group at the nonreducing end since it did not interact with sugar residues other than glucose." "the 6-O-substituted
  Gal group of Gal⁶G5P acts as a barrier for binding with enzyme." "The introduction of the 3-O-substituted GlcNAc group
  at the terminal position of NG5P does not appreciably influence the enzyme" action.
- **What it does not include.** Only the non-reducing-end residue was substituted, so +1/+2 tolerance of a substituted
  Glc is not tested; the substituent is a whole hexose, bulkier than methyl, hydroxypropyl or phosphate.
- **Subsite reading (this report).** Mammalian pancreatic GH13: a 6-O-glycosylated Glc is accepted at −3 and −4,
  almost never at −2 (≤ 0.03), never at −1; a 3-O-glycosylated Glc is accepted at −3 and −4, at −2 only in a very slow
  reaction (k₀/K_m = 2); a galactose (an epimer, as a terminal residue) is accepted at −3 but not at −2. For HPA, Gal⁶
  pushes binding outward (−4 preferred in Gal⁶G5P), while PPA prefers −3.

### 8.11 Omichi & Ikenaka 1988 — FT (new in v2; Japanese text, English abstract)
Omichi K, Ikenaka T (1988) Differential assays of human pancreatic and salivary α-amylases with fluorogenic
maltooligosaccharides as substrates. *J Jpn Soc Starch Sci* 35(2):165–171. doi:10.5458/jag1972.35.165 (OA, J-STAGE;
saved as `Literatura_Referencia/Omichi_1988_amylase_6-pyridylamino_maltooligosaccharides.pdf`)
- Substrate FG6R = FG-G-G-G-G-GR, where FG = 6-deoxy-6-[(2-pyridyl)amino]-D-Glc at the non-reducing end and GR =
  glucitol at the reducing end; prepared from 6-aldehydo amylose by reductive amination, glucoamylase trimming,
  partial B. subtilis liquefying alpha-amylase digestion and a second glucoamylase trimming (p. 166–167).
- Measured (p. 167, Fig. 3, HPLC of early digests): human pancreatic (HPA) and salivary (HSA) alpha-amylase (GH13) give
  FG-G-G (FG3) and FG-G-G-G (FG4); the bond diagram reads HPA 57% / 43% and HSA 83% / 17% (read left to right over
  the bonds giving FG3 and FG4 — my reading of the diagram; the text says only that the FG3/FG4 ratio differs between the
  two enzymes and uses it to assay them).
- Measured (p. 166, my translation of the Japanese): when the modified glucose residue reaches the non-reducing end,
  glucoamylase does not act further; FG6R and FG5P "are not acted on by exo-type glucoamylase or α-glucosidase because
  the non-reducing-end glucose residue is modified with the 2-pyridylamino group". (Source of the α-glucosidase not
  restated in the pages read.)
- Subsite reading (this report): a 6-deoxy-6-pyridylamino Glc is accepted at −3 (FG3) and −4 (FG4) by HPA and HSA;
  it is not accepted at −1 by GH15 glucoamylase (and the α-glucosidase used).

### 8.12 Nishio 2002 (award review) — FT (new in v2; Japanese text, English abstract and tables)
Nishio T (2002) Effect of sugar hydroxyl groups on activities and specificities of α-glycosidase and lipase. *J Appl
Glycosci* 49(1):45–55. doi:10.5458/jag.49.45 (OA, J-STAGE; saved as
`Literatura_Referencia/Nishio_2002_alpha-glucosidase_deoxy_substrates.pdf`)
- Substrates: p-nitrophenyl α-D-glucopyranoside (1) and its 2-, 3-, 4-, 6-monodeoxy derivatives (2–5); all ⁴C₁ by
  NMR/MM2 (p. 46). Relative rate = activity on the deoxy substrate / activity on 1 × 100 (Table 1, p. 47):

| α-glucosidase source | 2-deoxy | 3-deoxy | 4-deoxy | 6-deoxy |
|---|---|---|---|---|
| Saccharomyces cerevisiae | – | – | – | – |
| Bacillus stearothermophilus | – | – | – | – |
| honeybee I | – | – | – | – |
| honeybee III | – | – | – | – |
| rice | 175 | – | – | – |
| sugar beet | 244 | – | – | – |
| flint corn | 231 | 3.7 | – | – |
| Aspergillus niger | 259 | 11.9 | – | – |

  (– = not hydrolysed.) The text (p. 47, Japanese) relates the two behaviours to the two α-glucosidase families
  ("family I and II" by sequence); the enzymes in the first four rows are GH13 and those of rice, sugar beet, flint corn
  and A. niger are GH31 (standard correspondence of family I/II to GH13/GH31 — my note).
- Subsite reading: −1 of the four GH31 enzymes: O4 and O6 required; O2 not required (2-deoxy faster than the parent,
  175–259%); O3 required by rice and sugar beet enzymes, partly dispensable for flint corn (3.7%) and A. niger (11.9%).
  −1 of the four GH13 α-glucosidases: all four OH groups required. Only deoxygenation was tested, no bulky
  substituents, and only the glycone (−1).

---

## 9. Carboxymethyl starch

### 9.1 Steingaß, Adden, Mischnick 2024 — ABS (review)
Steingaß F, Adden A, Mischnick P (2024) Analysis of the substituent distribution in carboxymethyl-1,4-glucans on
different structural levels — an ongoing challenge. *Polysaccharides* 5(3):332–357. doi:10.3390/polysaccharides5030022
- Verbatim: "Due to the ionic character and acid functionality, all the attempts of the direct depolymerization of CMG
  and further sample preparation suffered from bias, side reactions, and multiple ion formation in MS. Finally, the
  transformation of CMGs to the corresponding hydroxyethylglucans (HEGs) by the reduction of the esterified carboxy
  groups with LiAlH4 opened the window for quantitative oligomer MS analysis."
- No primary paper on enzymatic (amylase) degradation of carboxymethyl starch with positional read-out was found (see
  Gaps). The full text of this review (MDPI) could not be fetched (HTTP 403).
- v2: MDPI still returned HTTP 403 to this host; still ABS.

### 9.2 Cellulose-ether analogues with endoglucanases — ABS (new in v2; context only, not starch enzymes)
Schagerlöf H, Richardson S, Momcilovic D, Brinkmalm G, Wittgren B, Tjerneld F (2005, online; issue 2006) Characterization of chemical
substitution of hydroxypropyl cellulose using enzymatic degradation. *Biomacromolecules* 7(1):80–85.
doi:10.1021/bm050430n
Schagerlöf H, Johansson M, Richardson S, Brinkmalm G, Wittgren B, Tjerneld F (2006) Substituent distribution and
clouding behavior of hydroxypropyl methyl cellulose analyzed using enzymatic degradation. *Biomacromolecules*
7(12):3474–3481. doi:10.1021/bm0604799
- Verbatim (both abstracts): "The different endoglucanases showed varying degradation capability ... depending on
  structure of the active site." MALDI-TOF-MS of short products; no subsite-level statement in the abstracts.
- Relevance: same enzymatic-fingerprinting logic (Lund group, as Richardson 2003) on β(1→4) glucans; included only to
  mark that the cellulose literature also stops at "unsubstituted regions are degraded".

---

## 10. Oxidised starch (supplementary)
Zhu Q, Bertoft E (1997) Enzymic analysis of the structure of oxidized potato starches. *Int J Biol Macromol*
21(1-2):131–135. doi:10.1016/s0141-8130(97)00052-4 — ABS. Hypochlorite-oxidised starch: lower beta-limit; "After
debranching, approximately 30% of the chains in the HO-sample and approximately 20% in the PO-sample remained partly
resistant to successive beta-amylolysis."

---

## 11. Summary table

Status column: FT = full text read; FT-thesis = read as the thesis chapter; ABS = abstract only; NOT READ = neither.
Rows marked **new** were added in v2; rows marked **upd** changed status or content since v1. All DOIs resolved on the
Crossref API on 2026-11-24.

| substituent | position | enzyme (family) | tolerated at | not tolerated / obstructive at | evidence type | status | reference |
|---|---|---|---|---|---|---|---|
| phosphate | O6 | A. niger glucoamylase (GH15) | +1 | −1 | limit products | ABS; FT (review Fig. 1) | Abe 1982 doi:10.1016/0167-4838(82)90006-1; Takeda 1983 doi:10.1016/0167-4838(83)90240-6; Takeda 1987 doi:10.5458/jag1972.34.225 |
| phosphate | O3 | A. niger glucoamylase (GH15) | (+2) | +1, −1 | limit products | ABS; FT (review) | Takeda 1983; Takeda 1987; Kamasaka 1997 doi:10.1271/bbb.61.238 |
| phosphate | O6 | sweet-potato beta-amylase (GH14) | +1 | −1, −2 | limit products | FT (review Fig. 1) | Takeda 1987; Takeda & Hizukuri 1981 doi:10.1016/s0008-6215(00)85242-8 |
| phosphate | O3 | sweet-potato beta-amylase (GH14) | — | +1 (and −1, −2) | limit products | FT (review) | Takeda 1987 |
| phosphate | O6 | B. subtilis alpha-amylase (GH13) | −2, +2 | −1, +1 | product structure | ABS; FT (review) | Takeda 1983 |
| phosphate | O6 | porcine pancreatic and A. oryzae alpha-amylase (GH13) | −3, +1 | −1, −2 | product structure | ABS; FT (review) | Takeda 1983; Takeda & Hizukuri 1986 doi:10.1016/s0008-6215(00)90271-4 |
| phosphate | O3 | B. subtilis and porcine pancreatic alpha-amylase (GH13) | −2, +3 | −1, +1, +2 | product structure | ABS; FT (review) | Takeda 1983 |
| phosphate | O3 | A. oryzae alpha-amylase (GH13) | −3, +2 | −2, −1, +1 | product structure | ABS; FT (review) | Takeda & Hizukuri 1986; Takeda 1987 |
| phosphate | O3 | bacterial liquefying alpha-amylase (GH13) | −3, −4 | (−2 not seen) | product structure | ABS | Kamasaka 1997 |
| **new** methyl | O2, O3 (O6 trace) | Bacillus liquefying alpha-amylase, BAN 240L (GH13) | −1: 6-O-Me (trace only) | −1: 2-O-Me, 3-O-Me (none found at reducing ends; 19.4% and 6.2% of all residues at MS 0.30) | reducing-end GLC-MS of exhaustive digests | FT | van der Burgt 1998 doi:10.1016/s0008-6215(98)00238-9 |
| **new** methyl | O2, O3, O6 | B. licheniformis alpha-amylase (GH13) + A. niger AMG (GH15) | +1: 2-O-, 6-O-, 2,6-di-O-Me | −1 ("usually unsubstituted"); 3-O-Me at neither −1 nor +1; inhibition 3-O > 2-O > 6-O | MS of end residues of digest products | ABS (ordering second-hand from full text via Huang 2006 thesis) | Mischnick 2001 doi:10.1002/1521-379x(200104)53:3/4<110::aid-star110>3.0.co;2-y |
| methyl | mixed | glucoamylase (GH15) | partial digestion: 66/56/43% at MS 0.10/0.21/0.30 | — | exhaustive digestion | FT | van der Burgt 2000a doi:10.1016/s0008-6215(99)00329-8 |
| methyl | mixed | B. subtilis alpha-amylase, B. cereus beta-amylase | — | GH14 −2/−1 (maltosyl unit) and +1; GH13 +1 (terminal residues of limit dextrins > 90% unsubstituted) | limit dextrins | FT | van der Burgt 2000b doi:10.1016/s0008-6215(00)00062-8 |
| methyl (O-Me) | each OH of methyl isomaltoside | A. niger glucoamylase (GH15) | +1: O2 | −1: O2′, O3′, O4′, O6′; +1: O4 strongly, O3 (1.9 kcal/mol) | kinetics of synthetic analogues | ABS | Lemieux 1996 doi:10.1139/v96-036 |
| **upd** deoxy | each OH of methyl β-maltoside | A. niger glucoamylase (GH15) | −1: 2′-deoxy (K_M 3.6 vs 0.33 mM), 3′-deoxy; +1: 6-deoxy (k₂/K_M 9.0 × 10³ vs 9.7 × 10³), 1,5-anhydroglucitol (K_M 5.2 mM), its 2-deoxy | −1: 4′-deoxy, 6′-deoxy; +1: 3-deoxy (v₀ = 0) | qualitative test + NMR kinetics | FT | Bock & Pedersen 1987 doi:10.3891/acta.chem.scand.41b-0617; Sierks 1992 doi:10.1021/bi00152a038 (ABS) |
| **new** C6 substituents (Cl, Br, I, N₃, F, NHAc, COOMe, 5,6-ene) | C6 of +1; C6′ of −1 | A. niger glucoamylase (GH15) | +1: all neutral C6 groups listed (6-F as fast as parent; 6-Cl ~10-fold slower, tight binding) | +1: CH₂NH₃⁺, COO⁻, 1,6-anhydro; −1: 6′-F, 6′-Br, 6′-H, 3′,6′-anhydro | qualitative test + NMR rates | FT | Bock & Pedersen 1988 doi:10.3891/acta.chem.scand.42b-0075 |
| deoxy | O2, O3, O4, O6 of glycone | A. niger alpha-glucosidase (GH31) | −1: 2-deoxy, 3-deoxy hydrolysed | −1: 4-deoxy, 6-deoxy not hydrolysed | kinetics | FT | Ogawa 2004 doi:10.5458/jag.51.9 |
| **new** deoxy | O2, O3, O4, O6 of glycone | GH31 α-glucosidases (rice, sugar beet, flint corn, A. niger) | −1: 2-deoxy (175–259% of parent); 3-deoxy only flint corn 3.7%, A. niger 11.9% | −1: 4-deoxy, 6-deoxy (all four); 3-deoxy (rice, sugar beet) | relative rates, pNP glucosides | FT | Nishio 2002 doi:10.5458/jag.49.45 |
| **new** deoxy | O2, O3, O4, O6 of glycone | GH13 α-glucosidases (yeast, B. stearothermophilus, honeybee I, III) | — | −1: every monodeoxy | relative rates | FT | Nishio 2002 |
| **new** β-Gal on O6 (6-O-glycosyl) | O6 of non-reducing-end Glc | PPA (GH13) | −3 (0.61 in Gal⁶G5P; 0.96 in Gal⁶G4P), −4 (0.39; 0.03) | −2 (0.01), −1 (never); k₀/K_m 54 vs 1300 for G4P | cleavage frequencies + kinetics | FT | Okada 2000 doi:10.5458/jag.47.35 |
| **new** β-Gal on O6 | O6 of non-reducing-end Glc | HPA (GH13) | −4 (0.83 in Gal⁶G5P; 0.06), −3 (0.17; 0.91) | −2 (0.03), −1 (never); k₀/K_m 7 vs 1400 for G4P | cleavage frequencies + kinetics | FT | Okada 2000 |
| **new** β-GlcNAc on O3 | O3 of non-reducing-end Glc | PPA, HPA (GH13) | −3 (0.96–1.0), −4 (0.04) | −2 (only in NG3P, k₀/K_m = 2 vs 12 for G3P) | cleavage frequencies + kinetics | FT | Okada 2000 |
| **new** 6-deoxy-6-(2-pyridylamino) | C6 of non-reducing-end Glc | human pancreatic and salivary alpha-amylase (GH13) | −3, −4 (FG3:FG4 57:43 HPA, 83:17 HSA, diagram reading) | GH15 glucoamylase and α-glucosidase: −1 | product HPLC | FT | Omichi & Ikenaka 1988 doi:10.5458/jag1972.35.165 |
| 6-iodo-6-deoxy, 6-amino-6-deoxy | O6 replaced | human salivary/pancreatic alpha-amylase (GH13) | 6-I at −3; 6-NH₂ at +2 | C6 changes at −3 shift binding register | product structure | ABS | Omichi 1992 doi:10.1093/oxfordjournals.jbchem.a123716; Omichi 1989 doi:10.1093/oxfordjournals.jbchem.a122910 |
| 4,6-O-benzylidene | O4 + O6 of non-reducing end | PPA, HSA (GH13) | −3 (favourable), −5 (PPA) | −4 (unfavourable) | product patterns | ABS | Kandra 1997 doi:10.1016/s0008-6215(96)00310-2; Kandra 2000 doi:10.1016/s0008-6215(00)00221-4 |
| hydroxyethyl | O2, O3, O6 | pig pancreatic alpha-amylase (GH13) | — | −1 (at least); O2 > O6 in obstruction | degradation sites by GLC; rates | ABS | Yoshida 1984 II doi:10.1002/star.19840360505; Yoshida 1973 doi:10.1002/star.19730251105 |
| hydroxyethyl | C2/C6 ratio | serum alpha-amylase (in vivo) | — | higher C2/C6 → slower degradation | clinical PK | ABS | Treib 1995 doi:10.1055/s-0038-1649964 |
| hydroxypropyl | mixed | pancreatin (GH13 + others) | — | ~5 consecutive unsubstituted units needed | kinetic statistics | ABS | Leegwater 1972 doi:10.1002/star.19720240105 |
| hydroxypropyl | mixed | alpha-amylase + amyloglucosidase | free glucose released only from unsubstituted runs | — | enzymatic + GLC-MS | ABS | Richardson 2000 doi:10.1016/s0008-6215(00)00103-8 |
| **new** hydroxypropyl | mixed | starch-hydrolysing enzymes (alpha-amylase, AMG) | — | −1/+1 ("preferentially cleaved glucosidic linkages between unsubstituted glucose units") | HPAEC-ESI-MS of products | ABS | Richardson 2001 doi:10.1016/s0021-9673(01)00690-2 |
| **new** hydroxypropyl | mixed | PPA (GH13), partial digests | HP residues inside released oligosaccharides (oligosaccharide fraction ≈ 2 × HP content of starch) | — | fractionation + HP assay | ABS | Wootton & Chaudhry 1981 doi:10.1002/star.19810330507 |
| **new** allyloxy-hydroxypropyl ether | not resolved | pullulanase + B. subtilis alpha-amylase + AMG | — | DP3 (and DP2) with one AHP resist GH15 (−1, +1 or +2) | MALDI-TOF-MS | FT-thesis | Huijbrechts 2008 thesis doi:10.18174/122072 |
| cationic ether | mixed (O2, O3, O6 near random) | alpha-amylase, amyloglucosidase, debranching | −2/+2 occasionally | −1, +1 (cleavage only between unsubstituted units; prefers > 2 adjacent unsubstituted) | ESI/MALDI-MS of digests | ABS | Richardson 2003 doi:10.1021/ac0300688 |
| **new** cationic ether | mixed | isoamylase (GH13), beta-amylase (GH14) | — | isoamylase: substituent at or near the branch residue; GH14: substituent at non-reducing end | ion exchange + GPC + NMR | ABS | Manelius 2002 doi:10.1002/1521-379X(200202)54:2<58::AID-STAR58>3.0.CO;2-1 |
| **new** cationic ether | mixed | pullulanase vs isoamylase (GH13) | pullulanase more tolerant at high DS | isoamylase less tolerant | GPC, ion exchange | ABS | Manelius 2005 doi:10.1002/star.200500364 |
| acetyl | O2/O3/O6 (not resolved) | alpha-amylase + amyloglucosidase | — | DP3 with one acetyl resists GH15 (−1, +1 or +2) | HPAEC + MALDI | FT-thesis | Chen 2004 doi:10.1016/j.carbpol.2004.02.004 |
| **upd** acetyl | not resolved | B. subtilis alpha-amylase (GH13); AMG (GH15); barley beta-amylase (GH14) | GH13: DP3 + 1 Ac released → acetylated residue accepted within −3…+3 | GH15: DP3 + 1–2 Ac resistant; GH14: beta-limit 40/37% vs 86% (amylose) | MALDI-TOF-MS, beta-limits | FT-thesis | Huang 2007 doi:10.1016/j.carbpol.2006.06.027 |
| octenyl succinate | mainly O2 (not measured) | amyloglucosidase, beta-amylase, isoamylase | cleavage "a few anhydroglucose units away" | exo-enzymes stop at OS residue | limits, GPC, HPAEC | ABS + FT-thesis | Bai 2014 doi:10.1016/j.foodchem.2013.12.012; Bai & Shi 2011 doi:10.1016/j.carbpol.2010.08.012 |

---

## 12. Overall picture (revised in v2)

Facts first (with source), then my reading, in separate sentences.

1. **Subsite −1.** Facts: no phosphate at −1 for GH13, GH14, GH15 (Takeda 1983, 1987; Abe 1982); no hydroxyethyl at
   −1 of PPA (Yoshida 1984 II); no 2-O- or 3-O-methyl glucose at the reducing ends of a Bacillus alpha-amylase digest,
   6-O-methyl only as a trace (van der Burgt 1998); reducing ends of methyl amylose digests "usually unsubstituted"
   (Mischnick 2001); GH15 needs O4′ and O6′, and every C6′ change tested abolishes hydrolysis (Bock & Pedersen 1987,
   1988); GH31 α-glucosidases need O4 and O6 and four GH13 α-glucosidases need all four OH (Nishio 2002); no 6-O- or
   3-O-glycosylated Glc at −1 of PPA or HPA (Okada 2000). Reading: −1 accepts deoxygenation at O2 (and, weakly, O3) but
   no added group anywhere, in all four families; the one exception seen is a trace of 6-O-methyl glucose at −1 of a
   GH13 liquefying amylase.

2. **Subsite +1.** Facts: GH15 accepts 6-phosphate (Takeda), every neutral C6 replacement tested (Bock & Pedersen 1988),
   2-O-, 6-O- and 2,6-di-O-methyl (Mischnick 2001, GH13 and/or GH15) and O2 methylation (Lemieux 1996); GH15 rejects
   3-deoxy (Bock & Pedersen 1987), 3-phosphate (Takeda), 3-O-methyl (Mischnick 2001), and a charged C6 in the
   disaccharide (Bock & Pedersen 1988). Reading: +1 of GH15 is open at O2 and O6 (neutral or, in the polymer, the
   6-phosphate ester) and closed at O3. For GH13, a 6-phosphate is accepted at +1 by PPA and TAKA-amylase and not by B.
   subtilis amylase (Takeda), and cationic/HP residues are disfavoured at +1 (Richardson 2001, 2003).

3. **Subsite −2 of GH13.** Facts: PPA and HPA accept a 6-O-glycosylated Glc at −2 in ≤ 3% of cuts and a
   3-O-glycosylated Glc only in a reaction with k₀/K_m = 2 (Okada 2000); the authors call subsite 2 (= −2) "strict";
   6-phosphate not at −2 for PPA, but 3-phosphate is accepted at −2 by PPA and B. subtilis amylase (Takeda 1983).
   Reading: −2 of mammalian GH13 tolerates a small O3 substituent (phosphate) but not a hexose on O3 or O6, nor a
   phosphate on O6.

4. **Subsites −3, −4 (and −5) of GH13.** Facts: accepted there: 6-O-Gal, 3-O-GlcNAc (Okada 2000), 6-deoxy-6-pyridylamino
   (Omichi 1988), 6-iodo (Omichi 1992), 6-phosphate (Takeda: PPA, TAKA at −3), 4,6-O-benzylidene (Kandra: favourable at
   −3, unfavourable at −4). Reading: from −3 outward the GH13 binding cleft accommodates even a hexose substituent; the
   enzymes differ in which of −3/−4 they prefer (HPA puts Gal⁶ at −4 in Gal⁶G5P, PPA at −3).

5. **Ranking by hydroxyl.** Facts: C-3 phosphate excluded from more subsites than C-6 phosphate (Takeda 1983); 2-O-
   hydroxyethyl obstructs more than 6-O (Yoshida 1973); methyl: 3-O > 2-O > 6-O (Mischnick 2001, second-hand ranking;
   abstract: 3-O-Me never in end residues); 6-O-Me is the only methyl seen at −1 (van der Burgt 1998). Reading: across
   substituent types the order of obstruction is O3 > O2 > O6.

6. **Random polymer derivatives.** Facts: products of HP, cationic, acetyl, AHP and OSA starches are explained by "cuts
   between unsubstituted units" (Richardson 2001, 2003) and enzyme-resistant DP2–DP3 oligomers carrying one substituent
   (Chen 2004, Huang 2007, Huijbrechts 2008). Reading: these data agree with points 1–2 but cannot place a substituent
   at a single subsite, because the position of the substituent inside the oligomer was never determined.

7. **Points of tension to keep in view (facts).** (i) Methyl at +1: van der Burgt 2000b found > 90% unsubstituted
   non-reducing ends on α,β-limit dextrins (GH14 last), while Mischnick 2001 found 2-O-, 6-O-methylated non-reducing
   termini after GH13 + GH15. Different enzymes, and different questions (fraction vs presence); not a contradiction
   until the full text of Mischnick 2001 gives numbers. (ii) Charged C6 at +1 of GH15: 6-phosphate ester accepted in
   the polymer (Takeda), 6-carboxylate rejected in a disaccharide (Bock & Pedersen 1988) — see §8.5b.

## 13. Gaps (revised in v2)

- **Subsite-resolved product structure for HES, HP, acetyl, OSA, carboxymethyl, AHP and cationic starch** still not
  found: these studies report extent, DP/DS of resistant oligomers or "unsubstituted runs", never the residue that carries
  the substituent. The positional data sets are: phosphate (Takeda/Hizukuri, polymer), methyl reducing ends (van der
  Burgt 1998, FT) and methyl end residues (Mischnick 2001, ABS only), and synthetic oligosaccharides (Okada 2000, Omichi,
  Kandra; Bock & Pedersen, Lemieux for GH15).
- **+1 and +2 of GH13 with a defined substituent**: only 6-phosphate (Takeda), 6-amino at +2 (Omichi 1992, ABS) and
  methyl at +1 (Mischnick 2001, ABS, enzymes not separated). Okada 2000 and Omichi 1988 substituted only the
  non-reducing end.
- **GH31 beyond −1**: only deoxy data at −1 (Ogawa 2004; Nishio 2002); nothing on bulky groups at −1 or on any group at
  +1/+2. Omichi 1988 states only that the α-glucosidase used does not act on a 6-pyridylamino non-reducing end.
- **GH14 with defined O-alkyl positions**: only phosphate (O3, O6), statistical methyl data, beta-limits of OSA,
  acetyl, AHP, cationic starch.
- **Hydroxyethyl at +1**: Yoshida 1984 II full text still unread.
- **Mischnick 2001 full text** is now the most valuable single unread paper (numbers behind "usually unsubstituted",
  2-/6-/2,6-methylated termini, and the 3-O > 2-O > 6-O ranking).

### 13.1 Remaining unavailable (closed access; no legal open copy found on 2026-11-24)

| paper | DOI | what was checked |
|---|---|---|
| Steeneken & Woortman 1994, Carbohydr Res 258:207 | 10.1016/0008-6215(94)84087-3 | OpenAlex closed; TNO record without file; Utrecht/WUR repositories; no thesis |
| Mischnick 2001, Starch/Stärke 53:110 | 10.1002/1521-379x(200104)53:3/4<110::aid-star110>3.0.co;2-y | OpenAlex closed (abstract read) |
| Zhao et al. 2012, Food Chem 133:1333 (and Zhao 2015, Carbohydr Polym 131:424) | 10.1016/j.foodchem.2012.02.021; 10.1016/j.carbpol.2015.06.005 | OpenAlex closed; research.wur.nl metadata only; no thesis in edepot |
| Kavitha & BeMiller 1998, Carbohydr Polym 37:115 | 10.1016/s0144-8617(98)00056-3 | OpenAlex closed; content known second-hand only |
| Richardson et al. 2003, Anal Chem 75:6499 | 10.1021/ac0300688 | OpenAlex closed; Lund LUP record 128984 and Richardson's thesis record (LUP 41995) have no file |
| Richardson, Cohen, Gorton 2001, J Chromatogr A 917:111 | 10.1016/s0021-9673(01)00690-2 | closed (abstract read, PubMed) |
| Manelius et al. 2002, Starch/Stärke 54:58; Manelius et al. 2005, Starch/Stärke 57:291 | see §5.5–5.6 | Wiley closed (abstracts read); Manelius 2000 Åbo Akademi thesis not online |
| Yoshida & Kishikawa 1984 (HES Part II) | 10.1002/star.19840360505 | Wiley closed |
| Wootton & Chaudhry 1981 I/II; Azemi & Wootton 1985 | 10.1002/star.19810330408; 10.1002/star.19810330507; 10.1002/star.19850370204 | Wiley closed (abstracts read) |
| Richardson & Gorton 2003 review, Anal Chim Acta 497:27; Mischnick & Momcilovic 2010 review, Adv Carbohydr Chem Biochem 64:117 | 10.1016/j.aca.2003.08.005; 10.1016/s0065-2318(10)64004-8 | closed, no abstract |
| Steingaß, Adden, Mischnick 2024 review (OA, gold) | 10.3390/polysaccharides5030022 | open access, but MDPI returns HTTP 403 to this host; the author can open it in a browser |

Tool limits in this pass: Scholar Gateway (Wiley) and Consensus free quotas exhausted; Firecrawl out of credits;
actachemscand.org unreachable (the Acta Chem Scand PDFs were taken from the Internet Archive capture of the publisher's
own free archive URL). No Sci-Hub or other unlicensed source was used.

## 14. Sources read in full or in thesis form

v1:
- Takeda 1987, J-STAGE PDF (Fig. 1, p. 227; English abstract p. 225).
- van der Burgt 2000a (Utrecht repository PDF), Table 1 p. 186.
- van der Burgt 2000b (Utrecht repository PDF), pp. 423–428.
- Ogawa et al. 2004 (J-STAGE PDF), pp. 9–10.
- Chen Z 2003 PhD thesis, Wageningen, Chapter 6 (= Chen, Schols, Voragen 2004).
- Bai Y 2013 PhD thesis, Kansas State, Chapters 4–5 (= Bai et al. 2014; Bai & Shi 2011).
- Kamasaka et al. 1997, p. 238 (scanned page; abstract and introduction).

v2 (PDFs saved in `<LOCAL_DISK>/Documents/Investigacion/Literatura_Referencia/`):
- van der Burgt et al. 1998, all 8 pages — `vanderBurgt_1998_methylstarch_alphaamylase_limitdextrins.pdf` (Utrecht
  repository, hdl 1874/5317).
- Bock & Pedersen 1987, Tables 1, 3, 4 and text pp. 617–623 — `Bock_1987_AMG_deoxy_maltosides.pdf`.
- Bock & Pedersen 1988, pp. 75–78, 81–82 — `Bock_1988_AMG_6-substituted_maltosides.pdf`.
- Okada, Ogawa, Usui 2000, pp. 35–44 (Figs. 1–3, Tables 3–4) — `Okada_2000_pancreatic_amylase_end-blocked_maltooligosaccharides.pdf`.
- Omichi & Ikenaka 1988, pp. 165–167 — `Omichi_1988_amylase_6-pyridylamino_maltooligosaccharides.pdf`.
- Nishio 2002, abstract and Table 1 (p. 47) — `Nishio_2002_alpha-glucosidase_deoxy_substrates.pdf`.
- Huang J 2006 PhD thesis, Wageningen, Chapters 1 (pp. 14–15), 4 (= Huang et al. 2007), 5 (p. 97), 6 (p. 116) —
  `Huang_2006_thesis_acetylated_pea_cowpea_starch_WUR.pdf`.
- Huijbrechts AML 2008 PhD thesis, Wageningen, Chapter 2 (pp. 31, 40–43) and Chapter 5 (pp. 102–104) —
  `Huijbrechts_2008_thesis_allyloxyhydroxypropyl_starch_WUR.pdf`.
All other entries: abstract only, or not read, as marked.
