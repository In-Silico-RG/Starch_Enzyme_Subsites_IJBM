# GH13 vs GH31: clan membership and common ancestry — literature check

Date: 2026-11-23 (23:31 -05). Purpose: verify the uncited draft sentence that GH13 and GH31 "do not share a
common ancestor for this activity — they are unrelated retaining folds that converged on the same chemistry".
No manuscript file was edited.

How the sources were reached: CAZy pages (cazy.org/Glycoside-Hydrolases.html, cazy.org/help, cazy.org/GH13.html)
and CAZypedia wiki source (`index.php?title=...&action=raw`) were fetched with curl on 2026-11-23; abstracts from
PubMed E-utilities; full text of Naumoff 2005 from PMC (open access). Quotes below are copied from those
downloads.

---

## Q1. Clan of GH13 and of GH31

**Answer (from sources).** GH13 is in clan GH-H with GH70 and GH77. GH31 is in clan GH-D with GH27 and GH36.
The recollection in the brief is correct.

**CAZy clan table** (cazy.org/Glycoside-Hydrolases.html, page dated "Last update: 2026-11-06" per WebFetch):

| Clan | Fold | Families (as listed) |
|---|---|---|
| GH-D | (β/α)8 | 27, 31, 36 |
| GH-H | (β/α)8 | 13, 70, 77 |
| GH-L | (α/α)6 | 15, 65, 125, 178 |
| GH-T | (β/α)7 | 57, 119 |

The other (β/α)8 clans in the same table: GH-A, GH-K (18, 20, 85), GH-R (29, 107).

**CAZy GH13 family page** (cazy.org/GH13.html): "Mechanism Retaining a Clan GH-H 3D Structure Status ( β / α ) 8
barrel Catalytic Nucleophile/Base Asp (experimental) Catalytic Proton Donor Glu (experimental)". (The CAZy GH31,
GH14, GH15, GH57 family pages timed out and were not read.)

**CAZypedia GH13** (Svensson B, Janecek S; page "Glycoside Hydrolase Family 13", infobox "Clan: GH-H"):
> "The α-amylase family represents a clan GH-H of three glycoside hydrolase families GH13, GH70 and GH77
> [MacGregor2001], and should be distinguished from the second smaller α-amylase family GH57 [Janecek2005]. A
> remote homology to the family GH31 has also been discussed [Janecek2007]."

**CAZypedia GH31** (Zhang R, Williams S; page "Glycoside Hydrolase Family 31", infobox "Clan: GH-D"):
> "Among these, the crystallographic study of the Sulfolobus solfataricus α-glucosidase (MalA) is notable as the
> first experimental report to highlight the 3-D structural relationship of GH31 enzymes to those of GH27 and
> GH36 [Ernst2006]; these three families now compose clan GH-D."
> "Enzymes of family GH31 are retaining α-glycosidases." "These structures reveal a common (β/α)8 barrel
> catalytic domain."

**Ernst et al. 2006** abstract: "Structural comparisons with other GH families suggest that the GH31 enzymes
belong to clan GH-D."

## Q2. What "different clan" means formally

**CAZy Help** (cazy.org/help):
> "Clans are groups of families sharing a fold and catalytic machinery. They represent a higher level in the
> hierarchy of CAZy classification and are only defined for Glycoside Hydrolases at present. [...] Attribution of
> GH families to clans usually follows structural determination."

**CAZy GH classification page** (cazy.org/Glycoside-Hydrolases.html):
> "Because the fold of proteins is better conserved than their sequences, some of the families can be grouped in
> 'clans' : (i) when new sequences are found to be related to more than one family, (ii) when the sensitivity of
> sequence comparison methods is increased or (iii) when structural determinations demonstrate the resemblance
> between members of different families [3]"

**CAZypedia "Glycoside hydrolases"** page (citing Henrissat & Bairoch 1996):
> "A 'clan' is a group of families that possess significant similarity in their tertiary structure, catalytic
> residues and mechanism. Families within clans are thought to have a common evolutionary ancestry."

**Janecek, Svensson & MacGregor 2007** abstract, on why GH31 is not in GH-H:
> "Although both the alpha-amylase super-family, i.e. the glycoside hydrolase (GH) clan GH-H (the GH families 13,
> 70 and 77), and family GH31 share some characteristics, their different catalytic machinery prevents
> classification of GH31 in clan GH-H."

**What the sources say.** Clan membership requires a shared fold *and* shared catalytic machinery. Being in the
same clan implies (is "thought to" reflect) common ancestry. None of the definitions states the converse: being
in different clans is not defined as absence of common ancestry. For GH13/GH31 specifically, the stated reason
for separate clans is different catalytic machinery, not absence of detectable relationship.

## Q3. Published arguments for distant common ancestry of GH13 (GH-H) and GH31 (GH-D)

**For common ancestry (four sources found; none found arguing the opposite):**

1. **Janecek S, Svensson B, MacGregor EA (2007)**, abstract:
   > "A significant but remote evolutionary relatedness is, however, proposed for clan GH-H with GH31. A sequence
   > alignment, based on the idea that residues equivalent in the primordial catalytic GH-H/GH31
   > (beta/alpha)(8)-barrel may not be found in the present-day GH-H and GH31 structures at strictly equivalent
   > positions, shows remote sequence homologies covering beta3, beta4, beta7 and beta8 of the GH-H and GH31
   > (beta/alpha)(8)-barrels. Structure comparison of GH13 alpha-amylase and GH31 alpha-xylosidase guided alignment
   > of GH-H and GH31 members for construction of evolutionary trees. The closest sequence relationship displayed
   > by GH31 is to GH77 of clan GH-H."
   This paper is cited for exactly this point on the CAZypedia GH13 page (quoted in Q1).

2. **Rigden DJ (2002)**, title and abstract:
   Title: "Iterative database searches demonstrate that glycoside hydrolase families 27, 31, 36 and 66 share a
   common evolutionary origin with family 13."
   > "Here we identify triose phosphate isomerase (TIM) barrels in GH families 27, 31, 36 and 66. We further show
   > that iterated sequence database searches provide evidence for their sharing a common evolutionary origin with
   > GH family 13. The catalytic, nucleophilic residue common to all these families is thereby determined and
   > candidate catalytic proton donors identified within each family."

3. **Naumoff DG (2005)**, full text (Background):
   > "The majority of known glycosidases with the α-galactosidase activity [EC3.2.1.22] belong to families GH27
   > and GH36, that form clan GH-D [2, 20]. This clan and family GH31 compose the α-galactosidase superfamily
   > [21-24]. This superfamily has a distant relationship with clan GH-H [25, 26], which we have proposed to name
   > the α-glucosidase superfamily [24]. Both superfamilies contain proteins sharing the same enzymatic mechanism
   > (retention), a similar (β/α)8-barrel fold of the catalytic domain [2], and use substrates only with the axial
   > orientation of the glycosidic bond [4]."
   Abstract conclusion: "These data suggest a common evolutionary origin of glycosidases representing different
   families and clans."
   Ref. 25 is Henrissat B (1998) Biochem Soc Trans 26:153-156 (a CAZy founder); ref. 26 is Rigden 2002.

4. **Henrissat B (1998)** "Glycosidase families", Biochem Soc Trans 26:153-156 — cited by Naumoff 2005 for the
   GH-D/GH-H distant relationship. **Not read**: no abstract in PubMed, no open full text. Its content is known
   here only through Naumoff's citation.

**Against (sources stating GH13 and GH31 are unrelated / independently evolved):** none found. Searches: PubMed
(Janecek S[au] AND GH31; "common evolutionary origin glycoside hydrolase family 13 family 31"), Europe PMC
("GH-H" AND "GH31" AND evolution/homology/ancestor, 14 hits screened by title). The only "unrelated" statements
found concern other pairs:
- GH13 vs GH57 (both retaining α-amylases, different folds and machinery). Janecek & Gabrisko 2016 abstract: "at
  least two different amino acid sequences utilizing two different catalytic machineries have evolved to execute
  the same α-amylolytic specificity. The two have been classified [...] in the glycoside hydrolase (GH) families
  GH13 and GH57." CAZypedia GH57 page: "the independence of the family GH57 with regard to GH-H clan is at present
  based not only on differences in the catalytic domain, but more importantly, on the differences in the catalytic
  machineries and conserved sequence regions".
- GH13 vs GH14 β-amylase. Janecek 1994 title: "...versus the barrel of beta-amylase: evolutionary distance is a
  reflection of unrelated sequences."

**Broader context, not specific to GH13/GH31.** Nagano, Orengo & Thornton 2002 (TIM-barrel review): the barrels'
"sequences are so diverse that even the most powerful approaches find few relationships, yet whose active sites
all cluster at one end of the barrel." (Full text not read; GH13/GH31 placement in their 21 superfamilies not
checked.)

## Q4. Clans of GH14, GH15, GH57

From the CAZy clan table (cazy.org/Glycoside-Hydrolases.html, 2026-11-23):
- **GH14** (β-amylase): not listed in any clan (no clan). CAZypedia GH14 page is unassigned (empty).
- **GH15** (glucoamylase): clan **GH-L**, (α/α)6, with GH65, GH125, GH178. CAZypedia GH15 infobox "Clan: GH-L".
- **GH57**: clan **GH-T**, (β/α)7, with GH119. CAZypedia GH57 infobox "Clan: GH-T". (Janecek & Gabrisko 2016
  abstract still described GH57 as "only been predicted to maybe define a future clan with the family GH119";
  the current CAZy table lists GH-T.)

---

## Interpretation (separate from the sources)

- GH13 and GH31 are in different CAZy clans because their catalytic machinery differs, per Janecek et al. 2007.
  CAZy's definition makes a shared clan sufficient evidence of common ancestry; it does not make a different
  clan evidence of separate origin.
- Three published analyses (Rigden 2002, Naumoff 2005, Janecek et al. 2007; plus Henrissat 1998 as cited) argue
  for a distant common ancestor of GH-H and GH31/GH-D. The 2007 paper is by the GH13 CAZypedia curators and is
  cited on that page. No source was found arguing that GH13 and GH31 are unrelated.
- Therefore the draft sentence "do not share a common ancestor ... unrelated retaining folds that converged" is
  contradicted by the literature as written. A defensible version: GH13 and GH31 sit in different clans (GH-H
  and GH-D) with different catalytic machinery; at most a remote common ancestry of their (β/α)8 barrels has been
  proposed (Janecek et al. 2007; Rigden 2002). A shared contact pattern at subsite −1 could then be either
  inheritance from that remote ancestor or convergence; the clan difference alone does not decide between them.
- The clean "unrelated, convergent" example in amylolytic enzymes is GH13 vs GH57 (different fold, (β/α)8 vs
  (β/α)7, and different machinery, same retaining α-amylase specificity), and GH13 vs GH14/GH15 (different
  mechanism and/or fold).

## Sources not reached

- Henrissat & Bairoch 1996 (Biochem J 316:695-696; PMC1217404): no machine-readable text returned (Europe PMC
  full-text endpoint empty). Its clan definition is known here only through the CAZypedia page that cites it.
- Henrissat 1998 (Biochem Soc Trans): no abstract, no open text.
- Lombard et al. 2014 and Drula et al. 2022: abstracts read; neither abstract defines clans or discusses
  GH13/GH31. Full texts not read.
- CAZypedia clan pages (Clan_GH-H, Clan_GH-D, Clans_of_Glycoside_Hydrolases): HTTP 404 at those URLs.
- CAZy family pages for GH31, GH14, GH15, GH57: timed out.

## References

- CAZy. Glycoside Hydrolase family classification. http://www.cazy.org/Glycoside-Hydrolases.html (accessed 2026-11-23).
- CAZy. Classifications in CAZy: short description. http://www.cazy.org/help (accessed 2026-11-23).
- CAZy. Glycoside Hydrolase Family 13. http://www.cazy.org/GH13.html (accessed 2026-11-23).
- CAZypedia Consortium. Glycoside Hydrolase Family 13 (authors Svensson B, Janecek S). https://www.cazypedia.org/index.php/Glycoside_Hydrolase_Family_13 (accessed 2026-11-23).
- CAZypedia Consortium. Glycoside Hydrolase Family 31 (authors Zhang R, Williams S). https://www.cazypedia.org/index.php/Glycoside_Hydrolase_Family_31 (accessed 2026-11-23).
- CAZypedia Consortium. Glycoside hydrolases. https://www.cazypedia.org/index.php/Glycoside_hydrolases (accessed 2026-11-23).
- CAZypedia Consortium. Glycoside Hydrolase Family 15; Family 57; Family 14 (accessed 2026-11-23).
- Henrissat B, Bairoch A (1996). Updating the sequence-based classification of glycosyl hydrolases. Biochem J 316:695-696. doi:10.1042/bj3160695. (Not read.)
- Henrissat B (1998). Glycosidase families. Biochem Soc Trans 26:153-156. doi:10.1042/bst0260153. (Not read.)
- Rigden DJ (2002). Iterative database searches demonstrate that glycoside hydrolase families 27, 31, 36 and 66 share a common evolutionary origin with family 13. FEBS Lett 523:17-22. doi:10.1016/S0014-5793(02)02879-X.
- Nagano N, Orengo CA, Thornton JM (2002). One fold with many functions: the evolutionary relationships between TIM barrel families based on their sequences, structures and functions. J Mol Biol 321:741-765. doi:10.1016/S0022-2836(02)00649-6.
- Naumoff DG (2005). GH97 is a new family of glycoside hydrolases, which is related to the alpha-galactosidase superfamily. BMC Genomics 6:112. doi:10.1186/1471-2164-6-112.
- Ernst HA, Lo Leggio L, Willemoës M, Leonard G, Blum P, Larsen S (2006). Structure of the Sulfolobus solfataricus alpha-glucosidase: implications for domain conservation and substrate recognition in GH31. J Mol Biol 358:1106-1124. doi:10.1016/j.jmb.2006.02.056.
- Janecek S, Svensson B, MacGregor EA (2007). A remote but significant sequence homology between glycoside hydrolase clan GH-H and family GH31. FEBS Lett 581:1261-1268. doi:10.1016/j.febslet.2007.02.036.
- MacGregor EA, Janecek S, Svensson B (2001). Relationship of sequence and structure to specificity in the alpha-amylase family of enzymes. Biochim Biophys Acta 1546:1-20. doi:10.1016/S0167-4838(00)00302-2.
- Janecek S (1994). Parallel beta/alpha-barrels of alpha-amylase, cyclodextrin glycosyltransferase and oligo-1,6-glucosidase versus the barrel of beta-amylase: evolutionary distance is a reflection of unrelated sequences. FEBS Lett 353:119-123. doi:10.1016/0014-5793(94)01019-6.
- Janecek S, Svensson B, MacGregor EA (2014). α-Amylase: an enzyme specificity found in various families of glycoside hydrolases. Cell Mol Life Sci 71:1149-1170. doi:10.1007/s00018-013-1388-z.
- Janecek S, Gabrisko M (2016). Remarkable evolutionary relatedness among the enzymes and proteins from the α-amylase family. Cell Mol Life Sci 73:2707-2725. doi:10.1007/s00018-016-2246-6.
- Lombard V, Golaconda Ramulu H, Drula E, Coutinho PM, Henrissat B (2014). The carbohydrate-active enzymes database (CAZy) in 2013. Nucleic Acids Res 42:D490-D495. doi:10.1093/nar/gkt1178.
- Drula E, Garron ML, Dogan S, Lombard V, Henrissat B, Terrapon N (2022). The carbohydrate-active enzyme database: functions and literature. Nucleic Acids Res 50:D571-D577. doi:10.1093/nar/gkab1045.
