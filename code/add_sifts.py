"""Add SIFTS mappings to the experimental database (AFC 2026-11-23).

Problem: residue numbers differ between depositions of the same protein and between proteins; to compare sites
we need, for every deposited protein residue, which UniProt protein and which UniProt residue it is, and the
domains (Pfam, InterPro, CATH) and EC numbers of every chain.
Data:
  - SIFTS residue-level files, https://ftp.ebi.ac.uk/pub/databases/msd/sifts/split_xml/<mid>/<pdb>.xml.gz: for each
    residue of each chain, the author residue number and insertion code (PDBresnum), the UniProt accession and
    residue number, and whether the residue is observed;
  - PDBe API https://www.ebi.ac.uk/pdbe/api/mappings/<pdb>: Pfam, InterPro, CATH, EC per chain, with ranges.
  Raw files cached in ../data/db/sifts_raw/ (not in git).
Method:
  1. table sifts_residue: one row per observed residue mapped to UniProt (pdb, chain, author number, insertion
     code, residue name, UniProt accession, UniProt number, UniProt one-letter residue);
  2. table sifts_domain: one row per (pdb, chain, database, accession, range) for Pfam, InterPro, CATH, EC;
  3. table contact_uniprot: every protein residue of table contact with its UniProt accession and number (NULL
     when SIFTS maps none).
Output: tables sifts_residue, sifts_domain, contact_uniprot in ../data/db/enzimas_experimental.sqlite.
"""
import gzip
import json
import re
import sqlite3
import time
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

R = Path(__file__).resolve().parent.parent / "resultados"
DB = R / "db" / "enzimas_experimental.sqlite"
RAW = R / "db" / "sifts_raw"
NS = "{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}"


def get(url, f):
    if f.exists() and f.stat().st_size > 0:
        return True
    for a in range(5):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                f.write_bytes(r.read())
            return True
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return False
            time.sleep(5 * (a + 1))
        except Exception:  # noqa: BLE001
            time.sleep(5 * (a + 1))
    return False


def residues(p):
    f = RAW / f"{p}.xml.gz"
    if not get(f"https://ftp.ebi.ac.uk/pub/databases/msd/sifts/split_xml/{p[1:3].lower()}/{p.lower()}.xml.gz", f):
        return None
    out = []
    root = ET.fromstring(gzip.open(f).read())
    for res in root.iter(NS + "residue"):
        pdb_x = unp = None; observed = True
        for x in res:
            if x.tag == NS + "crossRefDb":
                if x.get("dbSource") == "PDB":
                    pdb_x = x
                elif x.get("dbSource") == "UniProt" and unp is None:
                    unp = x
            elif x.tag == NS + "residueDetail" and x.text == "Not_Observed":
                observed = False
        if not observed or pdb_x is None or unp is None or pdb_x.get("dbResNum") in (None, "null"):
            continue
        m = re.fullmatch(r"(-?\d+)([A-Za-z]*)", pdb_x.get("dbResNum"))
        if not m:
            continue
        n, ins = int(m.group(1)), m.group(2)
        out.append((p, pdb_x.get("dbChainId"), n, ins, pdb_x.get("dbResName"), unp.get("dbAccessionId"),
                    int(unp.get("dbResNum")), unp.get("dbResName")))
    return out


def domains(p):
    f = RAW / f"{p}_mappings.json"
    if not get(f"https://www.ebi.ac.uk/pdbe/api/mappings/{p.lower()}", f):
        return []
    d = json.loads(f.read_text()).get(p.lower()) or {}
    out = []
    for dbname in ("Pfam", "InterPro", "CATH", "EC"):
        for acc, v in (d.get(dbname) or {}).items():
            for m in v.get("mappings", []):
                st, en = m.get("start") or {}, m.get("end") or {}
                out.append((p, m.get("chain_id"), dbname, acc, v.get("identifier") or v.get("name"),
                            st.get("author_residue_number"), en.get("author_residue_number"), m.get("domain")))
    return out


def main():
    RAW.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(DB)
    pdbs = [r[0] for r in con.execute("SELECT pdb FROM entry ORDER BY pdb")]
    res_rows, dom_rows, failed = [], [], []
    for n, p in enumerate(pdbs, 1):
        r = residues(p)
        if r is None:
            failed.append(p)
        else:
            res_rows += r
        dom_rows += domains(p)
        if n % 100 == 0:
            print(n, "/", len(pdbs), len(res_rows), "residues", flush=True)
    for t in ("sifts_residue", "sifts_domain", "contact_uniprot"):
        con.execute(f"DROP TABLE IF EXISTS {t}")
    con.execute("CREATE TABLE sifts_residue (pdb, chain, auth_seq, ins_code, resname, uniprot, uniprot_resnum, uniprot_res)")
    con.executemany("INSERT INTO sifts_residue VALUES (?,?,?,?,?,?,?,?)", res_rows)
    con.execute("CREATE INDEX i_sres ON sifts_residue(pdb, chain, auth_seq, ins_code)")
    con.execute("CREATE TABLE sifts_domain (pdb, chain, db, accession, identifier, auth_start, auth_end, domain)")
    con.executemany("INSERT INTO sifts_domain VALUES (?,?,?,?,?,?,?,?)", dom_rows)
    con.execute("""CREATE TABLE contact_uniprot AS
        SELECT c.rowid AS contact_rowid, c.instance_id, i.pdb, c.chain, c.seq, c.ins_code, c.comp,
               s.uniprot, s.uniprot_resnum
        FROM contact c JOIN instance i USING(instance_id)
        LEFT JOIN sifts_residue s ON s.pdb = i.pdb AND s.chain = c.chain AND s.auth_seq = c.seq AND s.ins_code = c.ins_code
        WHERE c.partner_kind = 'protein'""")
    con.commit()
    q = lambda s: con.execute(s).fetchall()
    summ = dict(entries_with_residue_mapping=q("SELECT COUNT(DISTINCT pdb) FROM sifts_residue")[0][0],
                residues_mapped=len(res_rows), uniprot_accessions=q("SELECT COUNT(DISTINCT uniprot) FROM sifts_residue")[0][0],
                domains_by_db=q("SELECT db, COUNT(*), COUNT(DISTINCT accession) FROM sifts_domain GROUP BY 1"),
                protein_contacts=q("SELECT COUNT(*), SUM(uniprot IS NOT NULL) FROM contact_uniprot")[0],
                resname_mismatch=q("""SELECT COUNT(*) FROM contact_uniprot cu JOIN sifts_residue s
                    ON s.pdb=cu.pdb AND s.chain=cu.chain AND s.auth_seq=cu.seq AND s.ins_code=cu.ins_code
                    WHERE s.resname != cu.comp""")[0][0], failed=failed)
    print(json.dumps(summ, indent=1))
    (R / "db" / "sifts_summary.json").write_text(json.dumps(summ, indent=1))
    con.close()


if __name__ == "__main__":
    main()
