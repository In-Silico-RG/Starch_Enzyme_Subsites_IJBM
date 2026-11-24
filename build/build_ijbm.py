"""Build manuscript_ijbm/ijbm_main.tex and carbpol_SI.tex from manuscript_short sources.
Body text is copied line-for-line; only markup changes (\\cite->\\citep, \\paragraph->\\subsection,
achemso front matter -> elsarticle front matter)."""
import re, sys
root = "/media/aldo/Aldo/Documents/Investigacion/Amilosa/Proyecto_F_Enzimas_Amiloliticas/"
src = open(root + "manuscript_short/grip_finding.tex").read().split("\n")
L = lambda a, b: "\n".join(src[a - 1:b])  # inclusive, 1-indexed


def conv(t):
    t = t.replace("\\cite{", "\\citep{")
    # \paragraph{Heading.} -> \subsection{Heading}
    t = re.sub(r"\\paragraph\{([^}]*?)\.?\}\s*", lambda m: "\\subsection{" + m.group(1) + "}\n", t)
    return t


def para_body(a, b):
    t = L(a, b)
    return conv(re.sub(r"^\\paragraph\{[^}]*\}\s*", "", t))

assert src[88].startswith("\\paragraph{Database and superposition.}")
assert src[102].startswith("\\paragraph{Position along the chain.}")
assert src[111].startswith("\\paragraph{Kind of contact.}")
assert src[119].startswith("\\paragraph{At the cut in GH13.}")
assert src[330].startswith("\\section{Conclusion}")
assert src[22].startswith("DRAFT ABSTRACT")

abstract = L(23, 40)
intro = conv(L(45, 85))
m1 = para_body(89, 101)
m2 = para_body(103, 110)
m3 = para_body(112, 118)
import sys; sys.path.insert(0, "/tmp/claude-1000/-media-aldo-Aldo-Documents-Investigacion-Amilosa/73ce7e7c-dcde-4c2d-b271-84f9a35a4992/scratchpad")
import ijbm_rule_patch as RP
results = RP.build_results(L, conv)
concl = RP.build_conclusion(conv(L(333, 341)))
ack = L(344, 346)
title = "What Fixes a Residue at the Catalytic Subsite of Starch-Degrading Glycoside\nHydrolases Is the Kind of Contact It Makes with the Sugar"

main = r"""%% ijbm_main.tex -- International Journal of Biological Macromolecules (Elsevier) submission format of
%% manuscript_short/grip_finding.tex. Text from the source except the rule section (44_rule_test, author decision 2026-11-24)
%% and the consistency edits listed in ../manuscript_carbpol/JOURNAL_RULES.md; journal rules there.
%% Build: pdflatex, bibtex, pdflatex x2.
\documentclass[preprint,12pt]{elsarticle}

\usepackage{amsmath,amssymb}
\usepackage{booktabs}
\usepackage{array}
\usepackage{siunitx}
\usepackage{graphicx}
\usepackage{lineno}
\usepackage{setspace}% double spacing: IJBM length guidance is in double-spaced pages (unofficial source)
\usepackage[hidelinks]{hyperref}
\graphicspath{{figures/}}
\sisetup{group-separator={\,},group-minimum-digits=4}
\newcommand{\pend}[1]{\textbf{[PENDIENTE: #1]}}

\journal{International Journal of Biological Macromolecules}

\begin{document}

\begin{frontmatter}

\title{""" + title + r"""}

\author{Aldo F. Combariza\corref{cor1}}
\ead{aldo.combariza@unisucre.edu.co}
\cortext[cor1]{Corresponding author.}
\affiliation{organization={Grupo de Investigaci\'{o}n {\itshape in silico},
  Departamento de Biolog\'{\i}a y Qu\'{\i}mica, Universidad de Sucre},
  addressline={Carrera 28 No.\ 5-267},
  city={Sincelejo},
  state={Sucre},
  country={Colombia}}

\begin{abstract}
""" + abstract + r"""
\end{abstract}

%% Graphical abstract: Elsevier specification 1328 x 531 px minimum (5:2, w x h), 300 dpi,
%% TIFF/EPS/PDF/MS Office; uploaded as a separate file in Editorial Manager.
\begin{graphicalabstract}
\fbox{\parbox[c][5.2cm][c]{13cm}{\centering \pend{graphical abstract, to be prepared by the
author; 1328\,$\times$\,531 px minimum, 300 dpi}}}
\end{graphicalabstract}

%% Highlights are in ijbm_highlights.tex (uploaded as a separate Word file).

\begin{keyword}
\pend{keyword 1} \sep \pend{keyword 2} \sep \pend{keyword 3} \sep \pend{keyword 4}
\sep \pend{keyword 5} \sep \pend{keyword 6}
\end{keyword}

\end{frontmatter}

\linenumbers
\doublespacing

\section{Introduction}

""" + intro + r"""

\section{Materials and methods}

\subsection{Database and superposition}
""" + m1 + r"""

\subsection{Position along the chain}
""" + m2 + r"""

\subsection{Kind of contact}
""" + m3 + r"""

\section{Results and discussion}

""" + results + r"""

\section{Conclusions}

""" + concl + r"""

\section*{CRediT authorship contribution statement}
\textbf{Aldo F. Combariza} (sole author): \pend{CRediT roles to be selected by the author from
the 14 CRediT roles}.

\section*{Declaration of competing interest}
\pend{statement of competing interests, to be made by the author}

\section*{Declaration of generative AI and AI-assisted technologies in the manuscript preparation process}
\pend{declaration to be written by the author (Elsevier template: tool used, purpose, extent of
oversight)}

\section*{Funding}
\pend{funding sources, or a statement that the work received no specific grant}

\section*{Data availability}
\pend{repository and DOI not yet decided}

\section*{Acknowledgements}
""" + ack + r"""

\section*{Appendix A. Supplementary material}
Supplementary material (Sections S1--S4, Tables S1--S3) accompanies this article
(ijbm\_SI.pdf).

\bibliographystyle{elsarticle-num}
\bibliography{ijbm}

\end{document}
"""
open(root + "manuscript_ijbm/ijbm_main.tex", "w").write(main)

# ---------------- SI ----------------
si = open(root + "manuscript_short/grip_finding_SI.tex").read()
i0 = si.index("\\begin{document}") + len("\\begin{document}")
i1 = si.index("\\bibliography{grip_finding}")
body = si[i0:i1].replace("\\cite{", "\\citep{")
SI = r"""%% ijbm_SI.tex -- Supplementary material (Appendix A) for ijbm_main.tex,
%% International Journal of Biological Macromolecules format. Text copied unchanged from manuscript_short/grip_finding_SI.tex.
\documentclass[preprint,12pt]{elsarticle}
\usepackage{amsmath,amssymb,booktabs,graphicx}
\usepackage[hidelinks]{hyperref}
\graphicspath{{figures/}}
\renewcommand{\thesection}{S\arabic{section}}
\renewcommand{\thetable}{S\arabic{table}}
\renewcommand{\theequation}{S\arabic{equation}}
\renewcommand{\thefigure}{S\arabic{figure}}
\journal{International Journal of Biological Macromolecules}
\begin{document}
\begin{frontmatter}
\title{Supplementary material for: """ + title + r"""}
\author{Aldo F. Combariza}
\ead{aldo.combariza@unisucre.edu.co}
\affiliation{organization={Grupo de Investigaci\'{o}n {\itshape in silico},
  Departamento de Biolog\'{\i}a y Qu\'{\i}mica, Universidad de Sucre},
  addressline={Carrera 28 No.\ 5-267},
  city={Sincelejo},
  state={Sucre},
  country={Colombia}}
\end{frontmatter}
""" + body + r"""
\bibliographystyle{elsarticle-num}
\bibliography{ijbm}
\end{document}
"""
open(root + "manuscript_ijbm/ijbm_SI.tex", "w").write(SI)
print("ok")
