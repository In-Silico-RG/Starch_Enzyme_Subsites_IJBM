# Style guide for the IJBM draft, measured from published papers (2026-11-24)

## Sources read

IJBM research articles, open access, full text (versions of record unless marked otherwise):

- **[C23]** Christensen SJ, Madsen MS, Zinck SS, Hedberg C, Sørensen OB, Svensson B, Meyer AS. Enzymatic potato
  starch modification and structure-function analysis of six diverse GH77 4-α-glucanotransferases. *Int J Biol
  Macromol* 224 (2023) 105–114. doi:10.1016/j.ijbiomac.2022.10.107 (DTU Orbit, publisher PDF). Starch modification
  plus enzyme structure.
- **[C24]** Caseiro C, McGregor NGS, Alves VD, Carvalho AL, Romão MJ, Davies GJ, Fontes CMGA, Bule P. Family GH157
  enzyme exhibits broad linkage tolerance and a dual endo/exo-β-glucanase activity on β-glucans. *Int J Biol
  Macromol* 282 (2024) 137402. doi:10.1016/j.ijbiomac.2024.137402 (White Rose, publisher PDF). Enzyme–carbohydrate
  structure (crystal complexes, subsites).
- **[D24]** Ding J, Li X, Jin Z, Abou Hachem M, Bai Y. Efficient glycosylation of polyphenols via dynamic
  complexation of cyclodextrin and synchronous coupling reaction with cyclodextrin glycosyltransferase in water.
  *Int J Biol Macromol* 280 (2024) 136065. doi:10.1016/j.ijbiomac.2024.136065 (DTU Orbit, journal pre-proof).
  Starch-active GH13 enzyme. This copy is the only one of the three that shows the highlights.

Papers by A. F. Combariza. Only one full text could be read: **[CSC09]** Combariza, Sastre, Corma, *Diffusion
Fundamentals* 11 (2009) 26, 2 pp. Preprints.org blocked automated download of the 2019 preprint
(10.20944/preprints201902.0124), and the publisher copies of the other papers are closed. For **[AFC-abs]** I read
the public abstracts (Europe PMC) of: J Phys Chem B 2005 (microwave-driven zeolites), PCCP 2009 (H2 in MOF-5),
PCCP 2012 (MOF confinement), Chem Soc Rev 2013 (silica zeolite force fields) and Arch Insect Biochem Physiol 2023
(dioscorin–trypsin docking/MD). The author-related observations below therefore cover abstracts only, plus one
2-page paper.

Our draft: `manuscript/ijbm_main.tex` as it stood on 2026-11-24 (after the revert in 4fb64fc7). The prose
counts leave out figures, tables, citations and math.

---

## 1. Section lengths (words in running text)

| Section | [C23] | [D24] | [C24] | **Ours** |
|---|---|---|---|---|
| Abstract | ~165, 1 para | ~245, 1 para | ~270, 1 para | **273, 1 para** |
| Introduction | ~465, 2 paras | ~700, 6 paras | ~830, 5 paras | **775, 5 paras** |
| Methods | ~1000, 10 subsections | ~1150, 11 subsections | ~3000, 15 subsections (at end; crystallography) | **1780, 5 subsections** |
| Results & discussion | ~2300, 5 subsections | ~2200, 6 subsections | ~4200, 7 subsections | **2870, 5 subsections** |
| Conclusion(s) | ~330, 2 paras | ~150, 1 para | ~190, 1 para | **243, 1 para** (+150 Limitations) |
| Figures / tables | 7 / 3 | 7 / 4 | 7 / 3 | **4 / 3** |

Totals are similar. The difference is in how long the units are. Methods paragraphs in the references run about
100–200 words, one per subsection ([C23] 2.4 "Combined temperature-pH optima", ~120 words). Ours run 230–430
words (§2.2, 353; §2.4, 431). Results paragraphs in the references run about 100–350 words: in [C23] most are
120–200, and the longest in [D24] (§3.1) is ~350. Ours run up to 460 words (the literature-comparison paragraph,
lines 435–458).

## 2. How the Introduction is built

- **[C23]**, 2 paragraphs. P1 (~330 words) moves from use (modified starch, clean label) to the enzyme class, its
  CAZy families, what is known, and the gap in one sentence: "There are currently no reports available providing a
  systematic comparison…". P2 (~135 words) gives the aim as a hypothesis ("The present study was undertaken to test
  the hypothesis that…") and one sentence on what was done.
- **[D24]**, 6 paragraphs: importance → strategy (glycosylation) → enzyme classes in 1 paragraph → current fix and
  its problem → proposed system, ending on the gap ("However, the factors limiting … remain unclear") → aim ("The aim
  was to…").
- **[C24]**, 4 paragraphs plus a one-sentence close: biological and industrial relevance → enzyme activities, CAZy
  families and folds in one paragraph, ending on the gap ("very limited or no structural and mechanistic
  information") → what is known of the family, then "Aiming to fill this gap…" with a 6-sentence preview of the
  findings → "Overall, this study provides…".
- Pattern: context, then the enzyme classification in **one paragraph at most**, then the gap in one sentence, then
  the aim (in two of three papers, a preview of the results). No dataset sizes and no caveats appear in the
  Introduction.

## 3. How Results/Discussion paragraphs are built

- **Topic sentence first.** It says what was tested or what was found, with a pointer to the figure or table:
  "Substrate specificity of LaGH157 and BcGH157 was evaluated with natural polysaccharides (Table 1)." [C24 §2.1];
  "Tt4αGT demonstrated superior potato starch modifying efficacy." [C23 §3.4]; "The diversity of glycosylated
  products hampered the quantification of yield…" [D24 §3.2].
- **Order inside a paragraph:** result with numbers → comparison with the literature ("in good agreement with
  previous reports [31,34,36]" [C23]) → one sentence of interpretation, marked as such ("This could indicate…", "We
  speculate…" [C23]).
- **Number density.** Numerals per 100 words of Results prose, leaving out reference, figure and family labels:
  [D24] 5.7, [C24] 12 (the §2.1–2.2 sample includes a list of specific activities), **ours 15.4**. Our two densest
  paragraphs carry 64 and 68 numerals (lines 403–420 and 435–458).
- **Statistics.** Data are given as mean ± SD with n stated once ("triplicate… mean ± standard deviation… one-way
  ANOVA… P < 0.05", [D24 §2.11]). Table entries carry Tukey superscript letters ([C23] Tables 1–3, footnote "at p <
  0.05 (one way ANOVA… Tukey-Kramer)"). Fit quality is a single number ("R² ≥ 0.90", [C23 §3.2]). Each comparison
  gets one test and one threshold. Sensitivity checks, if any, are left out of the prose. **Ours:** every contrast
  carries three p values (free, within family, within location), and the paragraph at lines 312–330 lists eight
  settings where the result holds and five where it weakens.

## 4. Bullet lists

None of the three papers has a bullet or numbered list in the body (0 bullet characters in the extracted text).
[D24] lists its two candidate mechanisms inline as "(1) k1 > k2 … ; (2) k1 < k2 …". Our (i)–(iii) hypotheses in
§3.4 follow this inline style. AFC also rejected bullet lists on 2026-11-24 (LOG 15:57).

## 5. Captions

- **Figures.** Each caption opens with a short title phrase, then describes the panels: "Fig. 1. Mechanism, domain
  architecture and phylogeny of GH77 4αGTs. (A)…" [C23]; "Fig. 4. The 3D structure of LaGH157. a, c. Ribbon…" [C24].
  Lengths: [C23] about 30–115 words (median ~45); [D24] 15–60; [C24] 60–185 (median ~100; the panels are structural).
  The captions hold no methods and no results. **Ours:** 37–113 words, the right range. However, Fig. 1 and Fig. 2
  open directly with "(a)" and have no title phrase.
- **Tables.** The title is one line of 4–20 words: "Rate of product formation." [C23 Table 1]; "Kinetic parameters
  of LaGH157 and BcGH157 on pNP-Lam2 and MLGs." [C24 Table 2]. Definitions, the test used and ND go in footnotes
  under the table. **Ours:** table captions of 84, 112 and 84 words, with the method (permutation design, scoring
  rule, baselines) inside the caption.

## 6. Methods: level of detail and what goes to the SI

- Each subsection names the operation and gives the parameters (concentrations, cut-offs, software with version)
  in ~100–200 words. Equations are numbered and every symbol is defined underneath ([C23] eq. 1 substrate inhibition,
  eq. 2 amylose fraction).
- The SI holds the routine evidence: SDS-PAGE, sequence identity matrix and progress curves [C23 Figs S1–S3]; the
  HPLC programme [D24 Table S1]; primers and constructs [C24 Suppl. Table 3].
- Figures are not described in Methods. Ours describes Fig. 2 at lines 234–236.
- [C24] puts Methods after the Conclusion (Results first). IJBM accepts both orders.

## 7. Abstract

All three are single unstructured paragraphs of 165–270 words. They open with 1–2 sentences of context, then "Here
we report" / "In this study", then the main results with few numbers ([C23] 1 number; [C24] 4; [D24] 7, all
yields or rate constants), and close on one implication sentence. None reports a p value or a subset count. The
author's own abstracts follow the same pattern ([AFC-abs]: 157–245 words). They open on the problem ("Despite the
sustained use of forcefield methodologies… few reviews… are available" [Chem Soc Rev 2013]), state the finding
plainly ("It was found that the confinement enhances…" [PCCP 2012]), and close with the consequence. **Ours:** 273
words, ~30 numerals, two p values and subset tallies ("GH14 and GH15 4 of 5 each", "18 of the 21 scored cells…").

## 8. Highlights

[D24] has 4 highlights. Each is a single claim of about 50–85 characters (the Elsevier limit is 85 characters,
spaces included) and contains no numbers or codes: "Dynamic complexation of CD promotes efficient polyphenol
glycosylation by CGTase." / "The rate-limiting step (k1 > k2) of the system is elucidated for the first time."
**Ours** (`submission/ijbm_highlights.txt`): 5 highlights, two with tallies ("11 of 14…", "12 of 14 'not tolerated'
cells agree; a subsite-only rule gets 10 of 12"), which a reader cannot interpret without the paper.

## 9. How limitations are stated

None of the three papers has a Limitations section. Limits appear once, inline, where they apply, in one sentence:
"It is outside the scope of the present study to map preferred reaction paths and subsite affinities." [C23 §3.3];
"The role of the DUF in Bt4αGT and Ll4αGT was not identified." [C23 §3.5]; "…which remains to be explored." [C24];
"…will be further explored in subsequent studies" [D24 Conclusion]. Unfinished work is never mentioned. **Ours:** a
Limitations subsection (150 words, which mentions unfinished MD, line 608), plus caveats repeated in the Abstract,
Introduction (lines 153–155), Results (lines 338–340, 440, 448) and Conclusions.

## 10. Tone

- Passive voice dominates. "We" is used to report a decision or an interpretation ("We speculate that…", "We found
  no significant correlation…" [C23]; "Here we report…" [C24]). The author uses "we" in the same way ("We have
  performed equilibrium and nonequilibrium molecular dynamics…", J Phys Chem B 2005; "We suggest a two-step
  mechanism…").
- Hedges are short and attached to an interpretation, not to a result: "likely", "may indicate", "perhaps", "We
  speculate" [C23]; "It is speculated that…" [D24].
- The references use evaluative adverbs ("Remarkably", "Interestingly", "Surprisingly" [C23]; "for the first time",
  "highest… reported so far" [D24]). This is recorded as an observation and **not** recommended for our draft (AFC
  rule: no bragging).
- Terms. The references use field-standard abbreviations, collected in one footnote ([C23] "Abbreviations: …").
  Ours coins about 12 codes (H/Ar/Al, m, m_map, INV/VAR/MC/EXP/BUR, F/s/N, O6*) plus working words ("cell", "held",
  "fixed").

---

## Do / don't for editing the draft

AFC reverted a broad readability pass on 2026-11-24 (commit 4fb64fc7; LOG 16:00). Make edits one place at a time
and agree each with AFC first.

**Do**
1. Open every Results paragraph with what was tested or found and the figure/table it rests on; give numbers after
   that.
2. Keep paragraphs to ≤ 250 words in Results and ≤ 200 in Methods. Split at the change of question, not
   mid-argument.
3. Report one test per comparison in the text. Move the alternative permutation schemes and sensitivity settings to
   Table 1 / SI and point to them in one sentence.
4. Start each figure caption with a title phrase. Keep table titles to one line and put definitions and the test in
   table footnotes.
5. Write each limitation once, where it applies, in one sentence.
6. Keep the Abstract to context → what was done → 3–5 key numbers → one implication.
7. Write highlights as claims in words (≤ 85 characters), without tallies.
8. Collect the coined codes in one place (a short abbreviations footnote or a definitions table) and use field
   terms where they exist.

**Don't**
1. Don't use bullet lists in the body, or captions that are paragraphs (AFC 2026-11-24).
2. Don't put dataset sizes, caveats or a figure walk-through in the aim paragraph of the Introduction.
3. Don't describe figures in Methods.
4. Don't mention unfinished work (the MD in §S3) in the paper.
5. Don't repeat the same score ("12 of 14 … 10 of 12") in Abstract, Results, Limitations, Conclusions and
   Highlights.
6. Don't add evaluative adverbs or claims of priority.

## The 10 places in our draft that deviate most

| # | Where (ijbm_main.tex) | Deviation | Norm (source) |
|---|---|---|---|
| 1 | Abstract, l. 64–79 | 273 words, ~30 numerals, 2 p values, subset tallies, a caveat ("18 of the 21 scored cells rest partly on phosphorylated glucan") | 165–270 words, 1–7 numbers, no p, closes on an implication [C23][C24][D24]; [AFC-abs] |
| 2 | Introduction ¶2, l. 114–129 | 239 words of CAZy taxonomy (clans, folds, both mechanisms, subfamily counts 35/49) | classification in 1 paragraph at most, ending on the gap [C24 ¶2][D24 ¶3] |
| 3 | Introduction ¶5, l. 148–158 | Aim mixed with dataset counts (637, 79, 18, 114), a caveat ("is a hypothesis…") and a Fig. 1 walk-through | aim/hypothesis plus what was done, ~135 words [C23 ¶2]; aim plus findings preview [C24] |
| 4 | Methods §2.2–2.3, l. 192–237 | Paragraphs of 353 and 305 words; five statistical designs in one paragraph; Fig. 2 described in Methods (l. 234–236) | one operation per ~100–200-word subsection; a statistics subsection of 3 lines [D24 §2.11] |
| 5 | Methods §2.4, l. 238–278 | 431 + 277 words; five state codes and three aggregation levels defined in running text | equations/definitions set out with symbols defined [C23 eq. 1–2]; procedural detail to SI [D24 Table S1] |
| 6 | Results §3.1 ¶2, l. 312–330 | ~340 words, ~37 numerals, three p values per contrast, lists of settings where it holds/fails | one test, mean ± SD, P < 0.05 [D24]; Tukey letters in tables [C23] |
| 7 | Results §3.1 ¶3, l. 332–345 | Conclusion of the subsection followed by self-limiting sentences and subfamily make-up of GH31 groups | interpretation in one sentence, limits once, inline [C23 §3.3, §3.5] |
| 8 | Results §3.2, l. 435–458 | Longest paragraph (460 words, 64 numerals): agreement, baselines, permutation, subsets, source details; defensive aside ("which their scoring rule makes easy") | ≤ ~350 words per paragraph [D24 §3.1]; split by question |
| 9 | Captions: Fig. 1 (l. 163), Fig. 2 (l. 350), Tables 1–3 (l. 361, 480, 545) | Figure captions without a title phrase; table captions of 84–112 words that contain method and scoring rules | title phrase then panels [C23][C24]; table title 4–20 words plus footnotes [C23 Table 1][C24 Table 2] |
| 10 | §3.5 Limitations (l. 605–614) + Conclusions (l. 616–630) + Highlights | Separate Limitations section mentioning unfinished MD; Conclusions and Highlights repeat the 12/14 vs 10/12 tally | limits inline, no unfinished work [C23][C24][D24]; highlights as claims without tallies [D24] |

Leftover from the 2026-11-24 rename: "trivial map" still appears in Table 2 (l. 497, "Trivial map (−1/+1
blocked)"), in Limitations (l. 612, "a trivial subsite-only map") and in Conclusions (l. 624, "a trivial map"). The
rest of the text now says "subsite-only rule".
