"""Shared helpers for supporting calculations 1-3 (laptop)."""
import os, string, sqlite3, gemmi

P = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
CIF = os.path.join(P, 'data', 'census', 'cif')
DB = os.path.join(P, 'data', 'db', 'enzimas_experimental.sqlite')
CALC = os.path.join(P, 'data', 'calc')


def db():
    return sqlite3.connect('file:%s?mode=ro' % DB, uri=True)


def load(pdb):
    """Model 1, no hydrogens, first alternate conformation kept."""
    st = gemmi.read_structure(os.path.join(CIF, pdb + '.cif'))
    st.setup_entities()
    st.remove_hydrogens()
    st.remove_alternative_conformations()
    for i in range(len(st) - 1, 0, -1):
        del st[i]
    return st


def single_letter_chains(st):
    """Rename chains with >1-character names so PDB format holds them.
    Returns {new_name: original_name}."""
    used = {ch.name for ch in st[0] if len(ch.name) == 1}
    pool = [c for c in string.ascii_uppercase + string.ascii_lowercase + string.digits if c not in used]
    mapping = {}
    for ch in st[0]:
        if len(ch.name) == 1:
            mapping.setdefault(ch.name, ch.name)
    for ch in st[0]:
        if len(ch.name) > 1:
            orig = ch.name
            new = next((k for k, v in mapping.items() if v == orig), None)
            if new is None:
                new = pool.pop(0)
                mapping[new] = orig
            ch.name = new
    return mapping


def is_protein_residue(res):
    info = gemmi.find_tabulated_residue(res.name)
    return res.het_flag == 'A' and info is not None and info.is_amino_acid()
