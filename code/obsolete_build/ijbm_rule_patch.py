"""Author decision (via coordinator, 2026-11-24): replace the old rule material of Results with
resultados/analysis/44_rule_test/{rule_section,tab_rule,tab_tests}.tex and make the kept paragraphs
consistent with the new statement. Used by build_ijbm.py.
# WARNING (2026-11-24): obsolete, kept for the record. Running it overwrites the manuscript. Do not run."""
import re

A44 = "<PROJECT_ROOT>/resultados/analysis/44_rule_test/"


def sub1(t, old, new):
    assert t.count(old) == 1, ("not unique/not found", old[:70], t.count(old))
    return t.replace(old, new)


def tables_single_spaced(t):
    return t.replace("\\centering\\small", "\\centering\\renewcommand{\\baselinestretch}{1}\\small")


def build_results(L, conv):
    # old lines: 120-169 (GH13 cut, one ring away, GH31 incl. Fig 2 sentence), 171-182 Fig 2,
    # 184-223 old Table 1, 225-231 "The rule"; kept: 233-329
    assert L(120, 120).startswith("\\paragraph{At the cut in GH13.}")
    assert L(225, 225).startswith("\\paragraph{The rule.}")
    assert L(233, 233).startswith("\\paragraph{The rule is not carried")
    rule = open(A44 + "rule_section.tex").read()
    rule = "\n".join(l for l in rule.split("\n") if not l.startswith("% ----"))
    rule = re.sub(r"\\paragraph\{([^}]*?)\.\}", r"\\paragraph{\1}", rule)  # elsarticle adds the run-in period
    tab_rule = open(A44 + "tab_rule.tex").read()
    tab_rule = sub1(tab_rule, "\\begin{table}[h]\n\\centering\\small", "\\begin{table}[p]\n\\centering\\renewcommand{\\baselinestretch}{1}\\scriptsize\\setlength{\\tabcolsep}{3pt}")
    tab_rule = sub1(tab_rule, "p{4.6cm}", "p{4.6cm}")
    tab_tests = open(A44 + "tab_tests.tex").read()
    tab_tests = sub1(tab_tests, "\\centering\\small", "\\centering\\renewcommand{\\baselinestretch}{1}\\scriptsize\\setlength{\\tabcolsep}{3pt}")
    fig2 = L(171, 182)  # figure environment kept unchanged (to be replaced later)
    fig2_sentence = L(168, 169)  # "Figure~\ref{fig:classes} shows the three kinds ..." (was the end of the GH31 paragraph)
    assert fig2_sentence.startswith("Figure~\\ref{fig:classes}")
    kept = conv(L(233, 329))

    # --- consistency edits in kept paragraphs ---
    kept = sub1(kept,
        "In GH13, the 24 entries whose labelled units are glucose only contact the six\n$-1$ positions",
        "In GH13, the 24 entries whose labelled units are glucose only contact the six\npositions contacted at $-1$ in at least 80\\% of the entries (Asp179, Asp289, His288, Arg177,\nHis92, Tyr51)")
    kept = sub1(kept,
        "glucose-only entries (5, from two proteins) contact all nine in 5 of 5.",
        "glucose-only entries (5, from two proteins) contact all nine such positions (Asp327, His600,\nAsp443, Asp542, Trp441, Trp539, Phe575, Ile328, Ile364) in 5 of 5.")
    kept = sub1(kept,
        "the third group\nof the rule cannot be tested there.",
        "neither the variability\nof their aliphatic residues nor that of their outer-subsite ($-2$) hydrogen-bonders can be\ntested there.")
    kept = sub1(kept,
        "Read together with Table~\\ref{tab:rule}, the\nresidues fixed across proteins are those whose contacts carry the large, measured\ntransition-state energy; aromatic contacts, smaller and permissive in orientation, keep\ntheir class but not their identity; aliphatic contacts, with no specific role, change.",
        "Read together with Table~\\ref{tab:rule}, the\nresidues fixed in every lineage, the side-chain hydrogen-bonders of the rings at subsites\n$-1$ and $+1$, are those whose contacts carry the large, measured transition-state energy;\naromatic contacts, smaller and permissive in orientation, keep an aromatic ring more often\nthan their identity; aliphatic contacts, with no specific role, change, as do contacts of\nevery kind at the outer subsites.")
    kept = sub1(kept,
        "the barley $\\alpha$-amylase\nhistidines above lose 99\\%",
        "the barley $\\alpha$-amylase~1\nhistidines His93 and His290 (His92 and His288 of the GH13 reference) lose 99\\%")
    kept = sub1(kept,
        "One of the seven non-catalytic GH31 positions of\nTable~\\ref{tab:rule} has a measured substitution:",
        "One of the non-catalytic GH31 positions of\nTable~\\ref{tab:rule} has a measured substitution:")
    kept = sub1(kept,
        "For the other six (their equivalents across GH31 proteins fixed by structural\nalignment~\\citep{Tan2010})",
        "For His600, Trp441, Trp539, Phe575, Ile328 and Ile364 (their equivalents across GH31 proteins\nfixed by structural alignment~\\citep{Tan2010})")
    kept = sub1(kept,
        "The full\nrule rests on 19 positions in GH13 and GH31;",
        "The test\nrests on 36 positions in GH13 and GH31 (Table~\\ref{tab:rule});")

    return (rule.strip() + "\n\n" + tab_rule.strip() + "\n\n" + tab_tests.strip() + "\n\n"
            + fig2_sentence + "\n\n" + fig2 + "\n\n" + kept)


def build_conclusion(concl):
    concl = sub1(concl,
        "every protein is decided by the kind of contact it makes. Side chains that hydrogen-bond\nthe sugar are fixed, including the catalytic machinery and residues one ring away from the\nbond; aromatic rings stacked against it stay aromatic; aliphatic side chains packed around\nit change from protein to protein.",
        "every lineage is decided jointly by where it touches the sugar and by the kind of contact it\nmakes. Side chains that hydrogen-bond the rings at subsites $-1$ and $+1$ are fixed,\nincluding the catalytic machinery and residues one ring away from the bond; aromatic rings\nstacked against them are kept as aromatic rings; aliphatic side chains packed around them,\nand residues of every kind at the outer subsites, change from lineage to lineage.")
    return concl
