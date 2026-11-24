"""Part 2: the 3D structures of the starch-breaking enzymes (AFC 2026-11-11).

Problem: for every family/activity of the Part 1 table, how many structures exist, how good they are,
which carry a carbohydrate ligand, and from which organisms.
Data: CAZy structure pages https://www.cazy.org/<FAM>_structure.html — one row per characterized
protein with a structure: protein name, EC number(s), organism, GenBank, UniProt, and for each PDB entry
(with chains) the carbohydrate ligands bound and the resolution in Å (blank for NMR).
Method: fetch each page once (cached), parse the protein rows and their nested PDB rows; keep the rows
whose EC is one of the Part 1 activities of that family (in_scope), the others are kept and flagged.
Outputs:
  ../data/part2_structures_all.tsv        one row per (protein, PDB entry)
  ../data/part2_summary_family_ec.tsv     per family x EC: proteins, PDB entries, best resolution,
                                                entries with a carbohydrate ligand, organisms
Usage: python3 part2_structures.py CACHE_DIR
"""
import csv
import html
import re
import sys
import time
import urllib.request
from collections import defaultdict
from pathlib import Path

R = Path(__file__).resolve().parent.parent / "resultados"


def fetch(name, cache):
    p = cache / f"{name}.html"
    if not p.exists():
        for a in range(5):
            try:
                req = urllib.request.Request(f"https://www.cazy.org/{name}.html", headers={"User-Agent": "Mozilla/5.0"})
                with urllib.request.urlopen(req, timeout=90) as r:
                    p.write_bytes(r.read())
                break
            except Exception:
                time.sleep(5 * (a + 1))
        else:
            raise RuntimeError(f"could not fetch {name}")
        time.sleep(1.0)
    return p.read_text(errors="replace")


def clean(x):
    x = re.sub(r"<br\s*/?>", "; ", x, flags=re.I)
    return re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", " ", x))).strip(" ;")


def parse(page):
    """Yield (protein, ECs, organism, uniprot, pdb, chains, ligands, resolution)."""
    rows = page.split('<tr valign="top" onmouseover=')
    for chunk in rows[1:]:
        # protein-level cells come before the nested table
        head, _, nested = chunk.partition("<table width='100%' border='0'>")
        cells = re.findall(r"<td[^>]*>(.*?)</td>", head, flags=re.S)
        if len(cells) < 5:
            continue
        protein = clean(cells[0])
        ecs = list(dict.fromkeys(re.findall(r"\d+\.\d+\.\d+\.(?:\d+|-)", clean(cells[1]))))
        organism = clean(cells[2])
        uniprot = clean(cells[4])
        for tr in re.findall(r'<tr valign="top">(.*?)</tr>', nested, flags=re.S):
            tds = re.findall(r"<td[^>]*>(.*?)</td>", tr, flags=re.S)
            if len(tds) < 3:
                continue
            m = re.search(r"structureId=([0-9A-Za-z]{4})", tds[0])
            if not m:
                continue
            ch = re.search(r"\[([^\]]*)\]", clean(tds[0]))  # link text may be wrapped in <b>
            yield protein, ecs, organism, uniprot, m.group(1).upper(), ch.group(1) if ch else "", clean(tds[1]), clean(tds[2])


def main(cache_dir):
    cache = Path(cache_dir)
    p1 = list(csv.DictReader(open(R / "part1_starch_breaking_activities.tsv"), delimiter="\t"))
    scope = defaultdict(set)
    for r in p1:
        scope[r["family"]].add(r["EC"])
    out = []
    for fam in sorted(scope, key=lambda f: (f[:2], int(re.sub(r"\D", "", f)))):
        page = fetch(f"{fam}_structure", cache)
        n = 0
        for protein, ecs, org, uni, pdb, chains, lig, res in parse(page):
            ins = any(e in scope[fam] for e in ecs)
            out.append(dict(family=fam, in_scope="yes" if ins else "no", EC=";".join(ecs), protein=protein,
                            organism=org, uniprot=uni, pdb=pdb, chains=chains, carbohydrate_ligands=lig,
                            resolution=res))
            n += 1
        print(fam, n, "PDB rows", flush=True)
    cols = ["family", "in_scope", "EC", "protein", "organism", "uniprot", "pdb", "chains",
            "carbohydrate_ligands", "resolution"]
    with open(R / "part2_structures_all.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in out:
            f.write("\t".join(r[c].replace("\t", " ") for c in cols) + "\n")
    # summary per family x EC (in scope only); a PDB entry counted once per family x EC
    agg = defaultdict(lambda: {"proteins": set(), "pdb": set(), "lig": set(), "res": [], "org": set()})
    for r in out:
        if r["in_scope"] != "yes":
            continue
        for e in r["EC"].split(";"):
            if e not in scope[r["family"]]:
                continue
            a = agg[(r["family"], e)]
            a["proteins"].add((r["protein"], r["organism"]))
            a["pdb"].add(r["pdb"])
            if r["carbohydrate_ligands"]:
                a["lig"].add(r["pdb"])
            try:
                a["res"].append(float(r["resolution"]))
            except ValueError:
                pass
            a["org"].add(r["organism"])
    names = {(r["family"], r["EC"]): r["name"] for r in p1}
    groups = {(r["family"], r["EC"]): r["group"] for r in p1}
    with open(R / "part2_summary_family_ec.tsv", "w") as f:
        f.write("group\tfamily\tEC\tactivity\tproteins\tpdb_entries\tpdb_with_carbohydrate\tbest_resolution_A\t"
                "median_resolution_A\torganisms\n")
        for (fam, e), a in sorted(agg.items(), key=lambda kv: (groups.get(kv[0], ""), kv[0])):
            res = sorted(a["res"])
            med = res[len(res) // 2] if res else ""
            f.write(f"{groups.get((fam, e), '')}\t{fam}\t{e}\t{names.get((fam, e), '')}\t{len(a['proteins'])}\t"
                    f"{len(a['pdb'])}\t{len(a['lig'])}\t{res[0] if res else ''}\t{med}\t"
                    f"{'; '.join(sorted(a['org']))}\n")
    print(len(out), "rows written")


if __name__ == "__main__":
    main(sys.argv[1])
