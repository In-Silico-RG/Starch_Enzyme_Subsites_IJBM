"""Part 1: which CAZy activities break amylose or amylopectin.

Input
  data/cazy_20261111_activities.tsv  (from cazy_family_pages.py: one row per family x activity,
                                            with the bond cut written as residue -1, reacting bond, residue +1 ...)
  data/cazy_20261111_families.tsv    (mechanism, clan, fold, catalytic residues per family)
  cazy_data.zip                            (sequence counts per family)
Selection, applied to the bond actually cut (residue -1, reacting bond, residue +1):
  A  hydrolysis of an alpha-D-Glcp-(1->4)-D-Glcp bond inside an alpha(1->4) glucan chain
     (residue +1 is aDGlcp followed by a further 1,4 bond), or of an alpha(1->6) branch bond
     carried by an alpha(1->4) chain (1,6 cut with a 1,4 bond on residue +1 or on residue -1/-2)
  B  same bonds broken by another reaction: phosphorolysis (Pi), elimination (lyase), oxidation (AA13, "1o,4")
  C  oligo-1,6-glucosidase EC 3.2.1.10 (1,6 cut next to a 1,X bond, written generic);
     alpha-glucosidase EC 3.2.1.20, which CAZy writes with a generic aglycone (1,X alcohol):
     starch action depends on the enzyme, to be checked
  excluded (with reason): transferases (T), isomerases (M), cuts in disaccharides only, trehalose,
     all-alpha(1->6) dextran chains, mixed 1,3 glucans
GT35: its activity is not on an activity page; taken from its family page ("Activities in Family":
  alpha-1,4-glucan phosphorylase, EC 2.4.1.1, retaining), added by hand below.
"""
import csv
import zipfile
from collections import Counter

R = "../data/"
acts = list(csv.DictReader(open(R + "cazy_20261111_activities.tsv"), delimiter="\t"))
fams = {r["family"]: r for r in csv.DictReader(open(R + "cazy_20261111_families.tsv"), delimiter="\t")}

acts.append({"family": "GT35", "activity_id": "page", "EC": "2.4.1.1",
             "name": "alpha-1,4-glucan phosphorylase", "res_m2": "", "bond_m1": "", "res_m1": "[aDGlcp",
             "reacting_bond": "1,4", "res_p1": "aDGlcp", "bond_p1": "1,4", "res_p2": "aDGlcp", "reactant": "Pi"})


def classify(a):
    m1, b, p1, bp1, p2 = a["res_m1"], a["reacting_bond"], a["res_p1"], a["bond_p1"], a["res_p2"]
    kind = a["activity_id"][:1]
    glc_m1 = "aDGlcp" in m1
    if a["EC"] == "3.2.1.20":
        return "C", "alpha-glucosidase, generic aglycone in CAZy"
    if (a["reacting_bond"] == "1,6" and a["bond_p1"] == "1,X" and a["activity_id"].startswith("H")
            and "aDGlcp" in a["res_m1"] and "aDGlcp" in a["res_p1"]):
        return "C", "alpha(1->6) cut next to an unspecified (1,X) bond: limit dextrins / isomaltose"
    if not glc_m1:
        return None, "not an alpha-Glc bond"
    if kind in ("T",):
        return "X", "transferase (rearranges, does not break the chain)"
    if kind == "M":
        return "X", "isomerase"
    if p1.endswith("]") or not p1:
        return "X", "cut in a disaccharide/small substrate only"
    if b in ("1,1",) or "1,1" in bp1:
        return "X", "trehalose-type substrate"
    chain14 = bp1 == "1,4" or a["bond_m1"] == "1,4" or "1,4" in p1
    if b == "1,4" or b == "1o,4":
        if "Glcp" not in p1 or bp1 == "1,3":
            return "X", "not an alpha(1->4) glucan"
        if b == "1o,4":
            return "B", "oxidative cleavage"
        if kind == "L":
            return "B", "elimination (lyase)"
        if a["reactant"] in ("aPi", "Pi", "ePi"):
            return "B", "phosphorolysis"
        if kind == "H":
            return "A", "hydrolysis of alpha(1->4)"
    if b == "1,6":
        if not chain14:
            return "X", "all-alpha(1->6) chain (dextran)"
        if kind == "H":
            return "A", "hydrolysis of alpha(1->6) branch"
    return "X", "other"


def endo_exo(a):
    return "exo (non-reducing end)" if a["res_m1"].startswith("[") or a["res_m2"].startswith("[") else "endo"


# sequence counts per family (whole family, all subfamilies)
counts = Counter()
with zipfile.ZipFile("../cazy_data.zip") as z, z.open(z.namelist()[0]) as f:
    for line in f:
        counts[line.split(b"\t", 1)[0].decode().split("_")[0]] += 1

rows, excluded = [], []
for a in acts:
    g, why = classify(a)
    if g is None:
        continue
    rec = {"group": g, "reason": why, "family": a["family"], "EC": a["EC"], "name": a["name"],
           "bond": f'{a["res_m1"]} -({a["reacting_bond"]})- {a["res_p1"]}', "attack": endo_exo(a),
           "reactant": a["reactant"]}
    (excluded if g == "X" else rows).append(rec)

cols = ["group", "family", "EC", "name", "attack", "bond", "reactant", "reason",
        "mechanism", "clan", "fold", "nucleophile_base", "proton_donor", "n_sequences"]
with open(R + "part1_starch_breaking_activities.tsv", "w") as f:
    f.write("\t".join(cols) + "\n")
    for r in sorted(rows, key=lambda r: (r["group"], r["family"], r["EC"])):
        fi = fams.get(r["family"], {})
        f.write("\t".join([r["group"], r["family"], r["EC"], r["name"], r["attack"], r["bond"], r["reactant"],
                           r["reason"], fi.get("Mechanism", ""), fi.get("Clan", ""),
                           fi.get("3D Structure Status", ""), fi.get("Catalytic Nucleophile/Base", ""),
                           fi.get("Catalytic Proton Donor", ""), str(counts[r["family"]])]) + "\n")
with open(R + "part1_excluded_activities.tsv", "w") as f:
    f.write("family\tEC\tname\tbond\treason\n")
    for r in sorted(excluded, key=lambda r: (r["reason"], r["family"])):
        f.write("\t".join([r["family"], r["EC"], r["name"], r["bond"], r["reason"]]) + "\n")
print(len(rows), "kept;", len(excluded), "excluded")
