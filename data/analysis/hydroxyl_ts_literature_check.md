# Literature check: sugar-hydroxyl contributions to transition-state vs ground-state binding at subsite -1

Date: 2026-11-23. Read-only literature verification (no project file changed except this one).

Sources and access. Metadata and abstracts came from the Europe PMC REST API (`/search?query=EXT_ID:<pmid>&resultType=core`).
Full text was read for one paper only: Frandsen et al. 2000, Plant Physiol (PMC59001, Table II). The ACS (Biochemistry,
Acc. Chem. Res.) and Wiley (Eur. J. Biochem.) pages returned HTTP 403, and the Europe PMC PDF of McCarter et al. 1992 could
not be downloaded. Wherever this file says "abstract", the numbers come from the abstract only and the full tables were
**not** checked.

Definitions used below. In this kind of study the parent substrate and its deoxy (OH → H) or deoxyfluoro (OH → F) analogue
are compared by kinetics:
- ΔΔG‡ = −RT ln[(kcat/Km)analogue / (kcat/Km)parent]. This is the loss of **transition-state** binding (relative to free
  E + S) when that hydroxyl is removed.
- The ground-state (Michaelis complex) contribution comes from Km or Ki: ΔΔG_GS = −RT ln[Ks,parent / Ks,analogue].
- RT at 25 °C = 2.479 kJ/mol. A factor of 10 in kcat/Km is therefore 5.7 kJ/mol.

---

## Q1. Which papers measured per-hydroxyl TS versus ground-state contributions with deoxy/deoxyfluoro substrates?

**Answer.** Your recollection is correct. Namchuk & Withers (1995) is the paper that separates ground-state from
transition-state contributions for each glycone hydroxyl of a retaining glycosidase (Agrobacterium β-glucosidase, GH1).
The same Withers group did the same for E. coli β-galactosidase (GH2, McCarter et al. 1992), and earlier for two
non-glycosidase glucosyl-transfer enzymes (glycogen phosphorylase and phosphoglucomutase). Zechel & Withers (2000,
Acc. Chem. Res.) is a real review of glycosidase mechanisms that covers non-covalent interactions. Only its abstract could
be read, so whether it restates the per-hydroxyl numbers was **not verified**.

### 1a. Namchuk & Withers 1995: Agrobacterium sp. β-glucosidase (GH1, retaining)
Namchuk, M. N.; Withers, S. G. Mechanism of *Agrobacterium* β-glucosidase: kinetic analysis of the role of noncovalent
enzyme/substrate interactions. *Biochemistry* **1995**, 34, 16194–16202. DOI 10.1021/bi00049a035. PMID 8519777.

Verbatim (abstract):
> "Contributions of each hydroxyl group to binding of these substrates at the ground state are relatively weak
> (interaction energies of 3.3 kJ/mol or smaller) but are much greater at the two transition states (glycosylation and
> deglycosylation). The strongest transition state interactions were at the 2 position (at least 18 and 22 kJ/mol for
> glycosylation and deglycosylation, respectively) with the interactions at the 3 and 6 positions contributing at least
> another 9 kJ/mol of binding energy at both transition states. The interaction at the 4 position is less crucial to
> transition state binding but important for stabilization of the glycosyl-enzyme intermediate."

| Hydroxyl | Ground state (kJ/mol) | TS, glycosylation (kJ/mol) | TS, deglycosylation (kJ/mol) |
|---|---|---|---|
| 2-OH | ≤ 3.3 | ≥ 18 | ≥ 22 |
| 3-OH | ≤ 3.3 | ≥ 9 | ≥ 9 |
| 4-OH | ≤ 3.3 | "less crucial" (no number in abstract) | important for the glycosyl-enzyme intermediate |
| 6-OH | ≤ 3.3 | ≥ 9 | ≥ 9 |

The values are lower bounds ("at least") because each OH → H or OH → F replacement also changes the intrinsic reactivity of
the sugar. The authors separate this electronic effect by comparing with the spontaneous hydrolysis rates. The individual
per-hydroxyl table values were **not read** (ACS paywall).

### 1b. McCarter, Adam & Withers 1992: E. coli (lacZ) β-galactosidase (GH2, retaining)
McCarter, J. D.; Adam, M. J.; Withers, S. G. Binding energy and catalysis. Fluorinated and deoxygenated glycosides as
mechanistic probes of *Escherichia coli* (lacZ) β-galactosidase. *Biochem. J.* **1992**, 286, 721–727.
DOI 10.1042/bj2860721. PMID 1417731. PMC1132963.

Verbatim (abstract):
> "rates found to be two to nine orders of magnitude lower than that for the parent compound. These large rate
> reductions result primarily from the loss of transition-state binding interactions due to the replacement of sugar
> hydroxy groups, and such interactions are estimated to contribute at least 16.7 kJ (4 kcal).mol-1 to binding at the
> 3, 4 and 6 positions and more than 33.5 kJ (8 kcal).mol-1 at the 2 position."

TS: ≥ 16.7 kJ/mol each at 3-, 4- and 6-OH, and > 33.5 kJ/mol at 2-OH. The abstract gives **no** ground-state numbers; the
full paper was not readable here.

### 1c. Related Withers studies on non-glycosidase glucosyl-transfer enzymes (useful as a contrast)
- Street, I. P.; Rupitz, K.; Withers, S. G. Fluorinated and deoxygenated substrates as probes of transition-state structure
  in glycogen phosphorylase. *Biochemistry* **1989**, 28, 1581–1587. DOI 10.1021/bi00430a024. PMID 2719921. Abstract:
  Vm/Km lower by "10(2) to 10(5) times", which gives "estimates of the strengths of hydrogen bonds to individual sugar
  hydroxyls at the transition state". The abstract gives no per-hydroxyl numbers.
- Percival, M. D.; Withers, S. G. Binding energy and catalysis: deoxyfluoro sugars as probes of hydrogen bonding in
  phosphoglucomutase. *Biochemistry* **1992**, 31, 498–505. DOI 10.1021/bi00117a028. PMID 1531025. Verbatim:
  > "Interactions at each hydroxyl are found to strengthen only slightly upon progressing from the ground state to the
  > transition state in contrast to that seen with glycogen phosphorylase [...] where transition-state interactions became
  > much stronger. This is in accord with the mechanisms for these two enzymes where no distortion of the sugar ring occurs
  > for phosphoglucomutase, whereas considerable distortion is expected for glycogen phosphorylase."

  This is the control case. The rise in hydroxyl interactions from ground state to TS is linked to the sugar ring being
  distorted and carrying oxocarbenium character at the TS. It is not a general property of every hydrogen bond to a sugar.

### 1d. Reviews
- Zechel, D. L.; Withers, S. G. Glycosidase mechanisms: anatomy of a finely tuned catalyst. *Acc. Chem. Res.* **2000**,
  33, 11–18. DOI 10.1021/ar970172+. PMID 10639071. The abstract mentions "the roles of noncovalent interactions" in TS
  stabilisation. Whether the body quotes the Namchuk numbers was **not verified**.
- Not read, but listed as a reference of Frandsen 2002: McCarter, J. D.; Withers, S. G. Mechanisms of enzymatic glycoside
  hydrolysis. *Curr. Opin. Struct. Biol.* 1994 (PMID 7712292).

---

## Q2. Has the same measurement been done on α-glycosidases (GH13, GH31, glucoamylase, CGTase)?

**Answer. Yes, for the transition state, and mainly by the Svensson / Sierks / Frandsen / Palcic / Bock group.** It has been
done for glucoamylase (GH15, inverting), for three GH13 α-glucosidases, and for a GH31 α-glucosidase (barley). There is one
important difference from Namchuk & Withers: these studies report ΔΔG‡ from kcat/Km (or Vmax/Km) only. **They do not
separate the ground-state part.** The one partial exception is the early glucoamylase work by Bock & Pedersen 1987, which
measured Km and Vm (numbers not read). I found **no** deoxy-substrate study of α-amylase proper or of CGTase that reports
energies. The only related study is qualitative (Evers et al. 1997, below).

### 2a. GH13: three yeast/bacterial α-glucosidases and oligo-1,6-glucosidase. Closest to the manuscript
Frandsen, T. P.; Palcic, M. M.; Svensson, B. Substrate recognition by three family 13 yeast α-glucosidases.
*Eur. J. Biochem.* **2002**, 269, 728–734. DOI 10.1046/j.0014-2956.2001.02714.x. PMID 11856334.

Verbatim (abstract):
> "The transition-state stabilization energy, DeltaDeltaG, contributed by the individual OH-groups was calculated from
> the activities for the parent and the deoxy analogs [...] This analysis of the energetics gave DeltaDeltaG values for
> all three enzymes ranging from 16.1 to 24.0 kJ.mol-1 for OH-2', -3', -4', and -6', i.e. the OH-groups of the
> nonreducing sugar ring. These OH-groups interact with enzyme via charged hydrogen bonds. In contrast, OH-2 and -3 of the
> reducing sugar contribute to transition-state stabilization, by 5.8 and 4.1 kJ.mol-1, respectively, suggesting that
> these groups participate in neutral hydrogen bonds. The OH-4 group is found to be unimportant in this respect [...]"

The substrate is methyl α-isomaltoside. The non-reducing ring (primed) sits at subsite −1 and the reducing ring at +1.

| Ring | Hydroxyls | ΔΔG‡ (kJ/mol) |
|---|---|---|
| −1 (non-reducing) | 2′, 3′, 4′, 6′ | 16.1–24.0 (range over the three enzymes) |
| +1 (reducing) | OH-2, OH-3 | 5.8, 4.1 (the abstract does not say which enzyme) |
| +1 | OH-4 | unimportant |

The per-enzyme table was **not read** (Wiley 403). The abstract says the results are compared with earlier data for
glucoamylase and for "a retaining alpha-glucosidase of glycoside hydrolase family 31". From the reference list, that GH31
enzyme is very likely the barley enzyme in 2b (inference, not read in the 2002 text).

### 2b. GH31: barley malt high-pI α-glucosidase. Full table read
Frandsen, T. P.; Lok, F.; Mirgorodskaya, E.; Roepstorff, P.; Svensson, B. Purification, enzymatic characterization, and
nucleotide sequence of a high-isoelectric-point α-glucosidase from barley malt. *Plant Physiol.* **2000**, 123, 275–286.
DOI 10.1104/pp.123.1.275. PMID 10806244. PMC59001.

Table II (verbatim values; Vmax/Km and ΔΔG‡ in kJ/mol). Vmax and Km of the analogues were not measured separately; the
methods say "Because these were sparse, we did not determine Vmax and Km":

| Substrate | Vmax/Km | ΔΔG‡ (kJ/mol) | Subsite of the removed OH |
|---|---|---|---|
| Methyl β-maltoside (parent) | 9.51 × 10⁻³ | — | — |
| Methyl 1-deoxy-β-maltoside | 1.05 × 10⁻² | −0.3 | +1 |
| Methyl 1,2-dideoxy-β-maltoside | 1.17 × 10⁻² | −0.3 | +1 |
| Methyl 3-deoxy-β-maltoside | 6.73 × 10⁻⁴ | 7.0 | +1 |
| 2′-Deoxy maltose | 2.83 × 10⁻⁴ | 9.3 | −1 |
| 3′-Deoxy maltose | 1.66 × 10⁻⁴ | 10.7 | −1 |
| 4′-Deoxy maltose | 7.53 × 10⁻⁶ | 18.9 | −1 |
| 6′-Deoxy maltose | 8.09 × 10⁻⁶ | 18.7 | −1 |

Verbatim: "ΔΔG‡ values of approximately 19 kJ mol−1 were correlated with the elimination of OH-4′ and -6′"; "The OH-2′ and
-3′ at the non-reducing, and the OH-3 at the reducing end ring each contributed 7 to 11 kJ mol−1". The abstract adds:
"all OH-groups in the non-reducing ring and OH-3 in the reducing ring of maltose formed important hydrogen bonds to the
enzyme in the transition state complex."

Worked check, 4′-OH: ratio = 7.53 × 10⁻⁶ / 9.51 × 10⁻³ = 7.92 × 10⁻⁴; ln = −7.14; ΔΔG‡ = 2.479 × 7.14 = 17.7 kJ/mol
at 25 °C. The paper gives 18.9. That value is reproduced with RT at about 37 °C (2.58 × 7.14 = 18.4) or with a slightly
different parent value. The assay temperature was not checked; the difference is under 1.5 kJ/mol.

### 2c. Glucoamylase (GH15, inverting; exo-acting on starch)
- Bock, K.; Pedersen, H. The substrate specificity of the enzyme amyloglucosidase (AMG). Part I. Deoxy derivatives.
  *Acta Chem. Scand. B* **1987**, 41, 617–628. DOI 10.3891/acta.chem.scand.41b-0617. PMID 3122479. All eight
  monodeoxy methyl β-maltosides were made. Verbatim: "three hydroxy groups (3, 4' and 6') are essential for the compounds to
  act as substrate for the enzyme." KM and VM were measured by ¹H NMR, so this paper could in principle separate ground
  state from TS. The numbers were **not read**.
- Sierks, M. R.; Bock, K.; Refn, S.; Svensson, B. Active site similarities of glucose dehydrogenase, glucose oxidase, and
  glucoamylase probed by deoxygenated substrates. *Biochemistry* **1992**, 31, 8972–8977. DOI 10.1021/bi00152a038.
  PMID 1390684. Verbatim: "The activity of all three enzymes was found to depend heavily on the 4- and 6-OH groups (4'- and
  6'-OH in maltose), where changes in binding energies from 10 to 18 kJ/mol suggested strong hydrogen bonds [...] The 3-OH
  (3'-OH in maltose) was involved in weaker interactions, while the 2-OH (2'-OH in maltose) had a very small if any role
  in transition-state binding."
- Sierks, M. R.; Svensson, B. Kinetic identification of a hydrogen bonding pair in the glucoamylase-maltose transition state
  complex. *Protein Eng.* **1992**, 5, 185–188. DOI 10.1093/protein/5.2.185. PMID 1350675. This combines deoxy substrates
  with mutants (a double-mutant-cycle type of analysis) and pairs Glu180 with maltose 2-OH.
- Frandsen, T. P.; Stoffer, B. B.; Palcic, M. M.; Hof, S.; Svensson, B. Structure and energetics of the
  glucoamylase-isomaltose transition-state complex [...]. *J. Mol. Biol.* **1996**, 263, 79–89. DOI 10.1006/jmbi.1996.0557.
  PMID 8890914. Abstract: OH-4′, OH-6′ and OH-4 are critical; kcat reduced "10(3) to 10(5)-fold" for deoxy and O-methyl
  analogues; with a conformationally biased isomaltoside the 4′- and 6′-deoxy analogues show "full delta delta G++ 18 to 20
  kJ mol-1".
- Sierks, M. R.; Svensson, B. Energetic and mechanistic studies of glucoamylase using molecular recognition of maltose OH
  groups coupled with site-directed mutagenesis. *Biochemistry* **2000**, 39, 8585–8592. DOI 10.1021/bi9925219.
  PMID 10913265. Ten active-site mutants were each profiled with the deoxy-maltose series (see Q3). Verbatim conclusion:
  "Binding interactions between the enzyme and the glucosyl group in subsite -1, particularly with the 4'- and 6'-OH
  groups, play an important role in substrate binding, while subsite 1 interactions may play a more important role in
  product release."
- Palcic, M. M.; Skrydstrup, T.; Bock, K.; Le, N.; Lemieux, R. U. Substrate recognition by amyloglucosidase: evaluation of
  conformationally biased isomaltosides. *Carbohydr. Res.* **1993**, 250, 87–92. DOI 10.1016/0008-6215(93)84157-2.
  PMID 8143293. "OH-6', OH-4', and OH-4 [...] are involved in key polar interactions with the enzyme".

### 2d. α-Amylase / CGTase: qualitative only
Evers, B.; Petríček, M.; Thiem, J. Specificity of amylases and cyclodextrin-glucanotransferase in reactions with
2-deoxy-maltooligosaccharides. *Carbohydr. Res.* **1997**, 300, 153–159. DOI 10.1016/s0008-6215(97)00040-2. PMID 9203340.
Verbatim: "Porcine pancreatic alpha-amylase did not tolerate the deoxygenated substrate, nor were cyclization experiments
with cyclodextrin-glucanotransferase (CGT) successful." Thermomonospora curvata α-amylase, β-amylase and amyloglucosidase
did cleave the 2-deoxy oligomers. No energies are reported. No per-hydroxyl energy study of α-amylase or CGTase was found
in Europe PMC or web searches. That is a negative search result, not proof that none exists.

### Pattern across α-glucosidases (my reading, not a quoted claim)
In GH13 (Frandsen 2002), GH31 (Frandsen 2000) and GH15 (Sierks 1992, 2000; Frandsen 1996), the hydroxyls of the subsite −1
glucose (2′, 3′, 4′, 6′) carry roughly 10–24 kJ/mol each of TS stabilisation. The +1 ring carries 0–7 kJ/mol. 4′-OH and
6′-OH are the largest in every case (about 18–20 kJ/mol). This matches the Withers β-glycosidase picture, apart from which
hydroxyl dominates: 2-OH in the β-enzymes, 4′/6′-OH in the α-glucosidases. Only the β-glycosidase studies (Namchuk 1995 and
probably McCarter 1992) give a measured ground-state counterpart. **For α-glycosidases, the statement "much more than to
ground-state binding" is not supported by a published measurement that I could find.** Only the TS side is measured.

---

## Q3. Polar (hydrogen-bonding) versus apolar/stacking contacts at −1, and tolerance of apolar substitution

**Answer.** There is partial precedent, mainly from glucoamylase. I found **no paper stating that apolar subsite residues
are interchangeable while polar ones are invariant** across GH13/GH31. One large PDB survey goes against part of the
manuscript's framing (Hudson et al. 2015, below).

1. **Glucoamylase, hydroxyl profiles of mutants (closest precedent).** Sierks & Svensson 2000 (above). Verbatim:
   > "Replacing the active site acidic groups, Asp55, Glu180, or Asp309, with the corresponding amides, and the neutral
   > Trp178 with the basic Arg, all substantially reduced the binding energy contribution of the 4'- and 6'-OH groups of
   > maltose at subsite -1 [...] the substitutions Tyr116 --> Ala, Ser119 --> Tyr, or Trp120 --> Phe also did not
   > substantially alter the energy profiles even though Trp120 has a critical role in directing conformational changes
   > necessary for activity. Since the mutations at Trp120 and Asp176 reduced k(cat) values by 50- and 12-fold,
   > respectively, a large effect on k(cat) is not necessarily accompanied by changes in hydroxyl group binding energy
   > contributions."

   Removing charged hydrogen-bond partners wipes out the subsite −1 hydroxyl energies. Conservative aromatic changes (Trp→Phe,
   Tyr→Ala) leave the hydroxyl profile intact, although kcat falls. This is GH15, not GH13/GH31, and the substitutions are
   aromatic, not Ile/Val/Leu/Met.
2. **Glucoamylase Arg54 (hydrogen bond to 4′-OH).** Frandsen, T. P.; Christensen, T.; Stoffer, B.; Lehmbeck, J.; Dupont,
   C.; Honzatko, R. B.; Svensson, B. *Biochemistry* **1995**, 34, 10162–10169. DOI 10.1021/bi00032a009. PMID 7640270.
   Verbatim: "The two mutants of Arg54 display losses in transition-state stabilization of 16-21 kJ mol-1 in the
   hydrolysis of different maltooligodextrins, which originate from a [(1.2-1.8) x 10(3)]-fold reduction in kcat and
   changes in Km ranging from 25% to 300% of the wild-type values." Here the hydrogen-bond loss goes almost entirely into
   kcat (the TS), with Km nearly unchanged. This is direct mutant-side support for "the H-bond acts at the TS, not in the
   ground state".
3. **GH13 invariant histidines.** Søgaard, M.; Kadziola, A.; Haser, R.; Svensson, B. *J. Biol. Chem.* **1993**, 268,
   22480–22484. DOI 10.1016/s0021-9258(18)41554-2. PMID 7901200. Barley α-amylase H93N and H290N: kcat/Km 0.3% and 1.2%
   of wild type. Verbatim: "these histidines participate in transition state stabilization but not directly in catalysis.
   This finding agrees with H bonding to OH groups of the valienamine ring of acarbose". MacGregor, E. A.; Janeček, Š.;
   Svensson, B. *Biochim. Biophys. Acta* **2001**, 1546, 1–20. DOI 10.1016/s0167-4838(00)00302-2. Verbatim: "most enzymes
   of the family also contain two histidine residues critical for transition state stabilisation."
4. **Charged versus neutral hydrogen bonds (general principle these papers use).** Fersht, A. R. et al. Hydrogen bonding and
   biological specificity analysed by protein engineering. *Nature* **1985**, 314, 235–238. DOI 10.1038/314235a0.
   Verbatim: deleting a partner to leave "an unpaired, uncharged hydrogen-bond donor or acceptor weakens binding energy by
   only 0.5-1.5 kcal mol-1. But the presence of an unpaired and charged donor or acceptor weakens binding by a further
   approximately 3 kcal mol-1." Frandsen 2002 uses this to call the 16–24 kJ/mol −1 contacts "charged hydrogen bonds".
5. **Stacking / CH–π contacts are weaker and orientationally permissive.**
   - Asensio, J. L.; Ardá, A.; Cañada, F. J.; Jiménez-Barbero, J. Carbohydrate–aromatic interactions. *Acc. Chem. Res.*
     **2013**, 46, 946–954. DOI 10.1021/ar300024d. Verbatim: "Experimental values measured in water are somewhat smaller,
     approximately 1.5 kcal/mol for each interaction between a monosaccharide and an aromatic ring." and "Despite their
     relatively modest influence on the stability of carbohydrate/protein complexes, the aromatic platforms play a major
     role in determining the specificity".
   - Keys, A. M.; Kastner, D. W.; Kiessling, L. L.; Kulik, H. J. The energetic landscape of CH–π interactions in
     protein–carbohydrate binding. *Chem. Sci.* **2025**, 16, 1746–1761. DOI 10.1039/d4sc06246a. Verbatim: "numerous
     orientations allow for highly favorable interaction strengths." Also Keys et al., *J. Biol. Chem.* **2025**, 301,
     110379. DOI 10.1016/j.jbc.2025.110379. Verbatim: "the favored orientation of a CH-π stacking interaction is
     controlled by the location of hydrogen bonds in the protein binding site." Both studies are on galactose-binding
     proteins, not GH13/GH31. They support "hydrogen bonds fix the geometry; stacking is permissive".
6. **Caution: aliphatic residues in sugar pockets.** Hudson, K. L.; Bartlett, G. J.; Diehl, R. C.; Agirre, J.; Gallagher,
   T.; Kiessling, L. L.; Woolfson, D. N. Carbohydrate–aromatic interactions in proteins. *J. Am. Chem. Soc.* **2015**, 137,
   15152–15160. DOI 10.1021/jacs.5b08424. Verbatim: "Within the carbohydrate-binding pockets, aliphatic hydrophobic residues
   are disfavored, whereas aromatic side chains are enriched. The greatest preference is for tryptophan with an increased
   prevalence of 9-fold." Across all carbohydrate-binding proteins, apolar packing against sugars is mainly aromatic. A
   finding that Ile/Val/Leu/Met are "always present" at GH13/GH31 −1 sites therefore needs to be set against this survey,
   and the manuscript's own counts should say whether Trp/Tyr/Phe were classed as polar or as stacking.
7. **Side-chain restriction as extra TS stabilisation (relevant to 6-OH).** Quirke, J. C. K.; Crich, D. Glycoside
   hydrolases restrict the side chain conformation of their substrates to gain additional transition state stabilization.
   *J. Am. Chem. Soc.* **2020**, 142, 16965–16973. DOI 10.1021/jacs.0c05592. This analyses 513 GH–pyranoside structures.
   Verbatim: "most glucosidases and β-mannosidases preferentially bind their substrates in the most reactive
   gauche,gauche (gg) conformation, thereby maximizing stabilization of the corresponding oxocarbenium ion-like transition
   state".

---

## What could not be verified
- The per-hydroxyl tables of Namchuk & Withers 1995, McCarter et al. 1992, Sierks et al. 1992, Sierks & Svensson 2000,
  Frandsen et al. 1996 and Frandsen et al. 2002 (publisher 403 / PDF not retrievable). Only abstract numbers are given
  above for these.
- Whether Zechel & Withers 2000 restates the Namchuk numbers.
- The ground-state numbers in Bock & Pedersen 1987 (Km, Vm measured).
- The identity of the GH31 enzyme compared in Frandsen 2002 (inferred to be barley, from the reference list).
- Any deoxy-substrate energy study of α-amylase proper or CGTase (none found).
