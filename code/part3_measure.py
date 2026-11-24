"""Part 3, steps 2-3: geometry of the catalytic site in one reference structure per family (AFC 2026-11-11).

Data: PDB coordinate files (RCSB, legacy .pdb format), downloaded once to ../data/part3/pdb/.
Catalytic residues: M-CSA (data/part3/mcsa_starch_entries.json), which gives them in the author
numbering of a reference PDB entry of the same protein; for families without an M-CSA entry, CAZypedia
residues are used for the structure CAZypedia names. The residue name found at each number is checked
(a mismatch, e.g. ASN where ASP is expected, marks an inactive mutant and is reported, not hidden).

Measured, per structure:
  d_CC   distance between the two catalytic carboxyl carbons (Asp CG, Glu CD), in Å
  d_OO   shortest distance between an oxygen of one carboxylate and an oxygen of the other, in Å
  For the bound ligand: the sugar at subsite -1 is the carbohydrate residue whose C1 is closest to the
  midpoint of the two catalytic carboxyl carbons; the glycosidic atom X is the O or N of another residue
  within 1.75 Å of that C1; the +1 sugar is the residue that owns X.
  d(nuc/base carboxyl C -> C1), d(acid carboxyl O -> X)
  ring puckering of the -1 sugar (Cremer-Pople Q, theta, phi; atoms O5 C1 C2 C3 C4 C5)
  phi = O5-C1-X-C4' and psi = C1-X-C4'-C5' of the scissile linkage
Output: ../data/part3/part3_site_geometry.tsv
"""
import math
import urllib.request
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent.parent / "data/part3"
PDBDIR = OUT / "pdb"

# family, label, PDB, chain, (nucleophile-or-base resname, number), (acid resname, number), source, mechanism
SITES = [
    ("GH13", "alpha-amylase, barley AMY2 (M-CSA ref)", "1AMY", "A", ("ASP", 179), ("GLU", 204), "M-CSA 397 (auth numbering 1amy)", "retaining"),
    ("GH13", "maltogenic alpha-amylase Novamyl", "1QHO", "A", ("ASP", 228), ("GLU", 256), "M-CSA 905 (auth numbering 1qho)", "retaining"),
    ("GH13", "neopullulanase T. vulgaris (M-CSA ref)", "1BVZ", "A", ("ASP", 325), ("GLU", 354), "M-CSA 433 (auth numbering 1bvz)", "retaining"),
    ("GH13", "neopullulanase T. vulgaris, glucan bound", "2D2O", "A", ("ASP", 325), ("GLU", 354), "M-CSA 433 numbering, same protein", "retaining"),
    ("GH13", "isoamylase P. amyloderamosa", "1BF2", "A", ("ASP", 375), ("GLU", 435), "M-CSA 421 (auth numbering 1bf2)", "retaining"),
    ("GH13", "oligo-1,6-glucosidase B. cereus", "1UOK", "A", ("ASP", 199), ("GLU", 255), "M-CSA 285 (auth numbering 1uok)", "retaining"),
    ("GH14", "beta-amylase soybean (M-CSA ref)", "1BYA", "A", ("GLU", 380), ("GLU", 186), "M-CSA 436 (auth numbering 1bya)", "inverting"),
    ("GH14", "beta-amylase soybean, maltotetraose bound", "1WDR", "A", ("GLU", 380), ("GLU", 186), "M-CSA 436 numbering, same protein", "inverting"),
    ("GH15", "glucoamylase A. awamori + acarbose", "1AGM", "A", ("GLU", 400), ("GLU", 179), "M-CSA 393 (auth numbering 1agm)", "inverting"),
    ("GH31", "maltase-glucoamylase human + acarbose", "2QMJ", "A", ("ASP", 443), ("ASP", 542), "M-CSA 976 (auth numbering 2qmj)", "retaining"),
    ("GH57", "4-alpha-glucanotransferase T. litoralis, no carbohydrate ligand (CA, TRS only; family reference)", "1K1X", "A", ("GLU", 123), ("ASP", 214), "CAZypedia GH57", "retaining"),
]
CARB = {"GLC", "BGC", "MAL", "AGL", "AC1", "ACI", "GLD", "G6D", "BMA", "MAN", "DAF", "ACR", "QPS", "AAO",
        "Z9N", "ADH", "HSD", "NDG", "GAL", "XYP", "BGC", "4GL", "TRE"}


def fetch(pdb):
    PDBDIR.mkdir(parents=True, exist_ok=True)
    p = PDBDIR / f"{pdb}.pdb"
    if not p.exists():
        with urllib.request.urlopen(f"https://files.rcsb.org/download/{pdb}.pdb", timeout=120) as r:
            p.write_bytes(r.read())
    return p


def parse(path):
    atoms = []
    for line in open(path):
        if line.startswith(("ATOM", "HETATM")):
            alt = line[16]
            if alt not in " A":
                continue
            atoms.append(dict(rec=line[:6].strip(), name=line[12:16].strip(), resn=line[17:20].strip(),
                              chain=line[21], resi=int(line[22:26]), x=np.array([float(line[30:38]),
                              float(line[38:46]), float(line[46:54])]), el=line[76:78].strip() or line[12:14].strip()[0]))
        if line.startswith("ENDMDL"):
            break
    return atoms


def carboxyl(atoms, chain, resn, resi):
    res = [a for a in atoms if a["chain"] == chain and a["resi"] == resi and a["rec"] == "ATOM"]
    if not res:
        return None, None, None, "missing"
    found = res[0]["resn"]
    table = {"ASP": ("CG", ("OD1", "OD2")), "GLU": ("CD", ("OE1", "OE2")),
             "ASN": ("CG", ("OD1", "ND2")), "GLN": ("CD", ("OE1", "NE2"))}  # amide = inactive-mutant proxy
    cname, onames = table.get(found, (None, None))
    if cname is None:
        return None, None, found, f"{found} found, {resn} expected (mutant?)"
    c = next(a["x"] for a in res if a["name"] == cname)
    os_ = [a["x"] for a in res if a["name"] in onames]
    return c, os_, found, "ok" if found == resn else f"MUTANT {found} found, {resn} expected (amide used as proxy)"


def dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    return math.degrees(math.atan2(np.dot(np.cross(b1n, v), w), np.dot(v, w)))


def cremer_pople(ring):
    """ring: 6 coordinates in the order O5 C1 C2 C3 C4 C5. Returns Q (Å), theta, phi (deg)."""
    r = np.array(ring) - np.mean(ring, axis=0)
    j = np.arange(6)
    r1 = (r * np.sin(2 * np.pi * j / 6)[:, None]).sum(0)
    r2 = (r * np.cos(2 * np.pi * j / 6)[:, None]).sum(0)
    n = np.cross(r1, r2); n /= np.linalg.norm(n)
    z = r @ n
    q2c = math.sqrt(2 / 6) * (z * np.cos(4 * np.pi * j / 6)).sum()
    q2s = -math.sqrt(2 / 6) * (z * np.sin(4 * np.pi * j / 6)).sum()
    q3 = (1 / math.sqrt(6)) * (z * (-1) ** j).sum()
    q2 = math.hypot(q2c, q2s)
    Q = math.sqrt(q2 ** 2 + q3 ** 2)
    theta = math.degrees(math.atan2(q2, q3))
    phi = math.degrees(math.atan2(q2s, q2c)) % 360
    return Q, theta, phi


def pucker_name(theta):
    if theta < 30:
        return "4C1 chair"
    if theta > 150:
        return "1C4 chair"
    if 60 <= theta <= 120:
        return "boat / skew-boat"
    return "half-chair / envelope"


def check_against_mcsa():
    """Every residue number taken from M-CSA must equal an auth_resid M-CSA gives for that reference PDB."""
    import json
    E = json.load(open(OUT / "mcsa_starch_entries.json"))
    known = {(c["pdb_id"].upper(), c["code"].upper(), c["auth_resid"]) for e in E for r in e["residues"]
             for c in r["residue_chains"]}
    for fam, label, pdb, ch, nb, ac, src, mech in SITES:
        if src.startswith("M-CSA") and "auth numbering" in src:
            for resn, num in (nb, ac):
                assert (pdb, resn, num) in known, f"{pdb} {resn}{num} not in M-CSA"
    print("residue numbers checked against M-CSA auth_resid")


def main():
    check_against_mcsa()
    rows = []
    for fam, label, pdb, ch, nb, ac, src, mech in SITES:
        atoms = parse(fetch(pdb))
        c1, o1, f1, s1 = carboxyl(atoms, ch, *nb)
        c2, o2, f2, s2 = carboxyl(atoms, ch, *ac)
        rec = dict(family=fam, structure=label, pdb=pdb, mechanism=mech, residue_source=src,
                   nuc_or_base=f"{nb[0]}{nb[1]}:{s1}", acid=f"{ac[0]}{ac[1]}:{s2}")
        if c1 is None or c2 is None:
            rows.append(rec); print(pdb, "residues not found", s1, s2); continue
        rec["d_CC_A"] = round(float(np.linalg.norm(c1 - c2)), 2)
        rec["d_OO_min_A"] = round(min(float(np.linalg.norm(a - b)) for a in o1 for b in o2), 2)
        mid = (c1 + c2) / 2
        het = [a for a in atoms if a["rec"] == "HETATM" and a["resn"] != "HOH"]
        groups = {}
        for a in het:
            groups.setdefault((a["chain"], a["resn"], a["resi"]), []).append(a)
        rings = []  # (key, role->coord, kind)
        for key, g in groups.items():
            A = {a["name"]: a["x"] for a in g}
            if key[1] == "AC1" and all(k in A for k in ("C1B", "C2B", "C3B", "C4A", "C5B", "C7B")):
                # acarbose deposited as one residue: valienamine ring C7B-C1B-C2B-C3B-C4A-C5B, linked by N4A to C4
                rings.append((key, {"O5": A["C7B"], "C1": A["C1B"], "C2": A["C2B"], "C3": A["C3B"], "C4": A["C4A"],
                                    "C5": A["C5B"], "X": A["N4A"], "C4p": A["C4"], "C5p": A["C5"]},
                              "valienamine (acarbose), ring O5 position is carbon C7"))
            elif all(k in A for k in ("O5", "C1", "C2", "C3", "C4", "C5")):
                rings.append((key, {k: A[k] for k in ("O5", "C1", "C2", "C3", "C4", "C5")}, "pyranose"))
        best = None
        for key, R, kind in rings:
            d = float(np.linalg.norm(R["C1"] - mid))
            if d < 9 and (best is None or d < best[0]):
                best = (d, key, R, kind)
        if best is None:
            rec["ligand_at_-1"] = "none within 9 Å of the site"
        else:
            _, key, R, kind = best
            # partner X bonded to C1 from another residue, unless the acarbose residue carries it internally
            if "X" not in R:
                for k2, g2 in groups.items():
                    if k2 == key:
                        continue
                    for a in g2:
                        if a["el"] in ("O", "N") and float(np.linalg.norm(a["x"] - R["C1"])) < 1.75:
                            B = {b["name"]: b["x"] for b in g2}
                            R = dict(R, X=a["x"], Xname=f"{k2[1]} {k2[0]}{k2[2]} {a['name']}")
                            if "C4" in B and float(np.linalg.norm(B["C4"] - a["x"])) < 1.75:
                                R["C4p"], R["C5p"] = B["C4"], B.get("C5")
                            if a["el"] == "N":
                                kind = "valienamine-type ring deposited as " + key[1] + " (C1 bonded to N4 of the next unit)"
            else:
                R["Xname"] = f"{key[1]} {key[0]}{key[2]} N4A"
            rec["ligand_at_-1"] = f"{key[1]} {key[0]}{key[2]}: {kind}"
            rec["d_nuc/base_C_to_C1_A"] = round(float(np.linalg.norm(c1 - R["C1"])), 2)
            Q, th, ph = cremer_pople([R[k] for k in ("O5", "C1", "C2", "C3", "C4", "C5")])
            rec["Q_A"], rec["theta_deg"], rec["phiCP_deg"], rec["pucker_at_-1"] = round(Q, 3), round(th, 1), round(ph, 1), pucker_name(th)
            if "X" in R:
                rec["plus1"] = R["Xname"]
                rec["d_acid_O_to_X_A"] = round(min(float(np.linalg.norm(o - R["X"])) for o in o2), 2)
                if "C4p" in R:
                    rec["phi_deg"] = round(dihedral(R["O5"], R["C1"], R["X"], R["C4p"]), 1)
                    if R.get("C5p") is not None:
                        rec["psi_deg"] = round(dihedral(R["C1"], R["X"], R["C4p"], R["C5p"]), 1)
                else:
                    rec["linkage_note"] = "X not bonded to C4 of the next unit (not a (1-4) linkage)"
            else:
                rec["plus1"] = "no residue bonded to C1 (product or single sugar at -1)"
        rows.append(rec)
        print(pdb, {k: v for k, v in rec.items() if k not in ("structure", "residue_source")})
    cols = ["family", "structure", "pdb", "mechanism", "residue_source", "nuc_or_base", "acid", "d_CC_A", "d_OO_min_A",
            "ligand_at_-1", "d_nuc/base_C_to_C1_A", "plus1", "d_acid_O_to_X_A", "Q_A", "theta_deg", "phiCP_deg",
            "pucker_at_-1", "phi_deg", "psi_deg", "linkage_note"]
    with open(OUT / "part3_site_geometry.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


if __name__ == "__main__":
    main()
