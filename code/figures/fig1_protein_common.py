"""Shared definitions for the Figure 1 candidates (protein overview, IJBM manuscript).

Every atom drawn comes from a deposited mmCIF (data/census/cif, model 1, first alternate conformation,
no hydrogens). Nothing is drawn by hand: the cartoon, surface and sticks are PyMOL representations of those
coordinates; the scatter panels plot ring centroids computed from the same files and moved into the family
frame with the stored superpositions (data/calc/03_frames/superposition.tsv).

Sugar units, their rings and their subsites: data/analysis/33_subsites/units.tsv (one row per ring;
reference atom C1/C1<letter>; subsite label), rings found with ring_units() of code/analysis/a33_subsites.py
(the same function that produced the subsite table, so a unit here is exactly a unit there).

Run with ~/miniconda3/bin/python (gemmi 0.7.5, numpy, networkx, matplotlib); renderer: PyMOL (pymol -cq).
"""
import sys
import csv
from pathlib import Path

import numpy as np
import gemmi

P = Path(__file__).resolve().parents[2]
CIF = P / 'data/census/cif'
OUT = P / 'figures/candidates/fig1_protein'
WORK = OUT / 'work'
WORK.mkdir(parents=True, exist_ok=True)
UNITS33 = P / 'data/analysis/33_subsites/units.tsv'
UNITS44 = P / 'data/analysis/44_rule_test/units_used.tsv'
SUP = P / 'data/calc/03_frames/superposition.tsv'

sys.path.insert(0, str(P / 'code/analysis'))
from a33_subsites import ring_units  # noqa: E402

# Style of the approved Figure 2 (code/figures/fig2_v3_common.py / fig2_v3_struct.py): sugar gold with the
# subsite numbers written on the rings, catalytic residues magenta, faint grey cartoon, grey surface,
# Liberation Sans, ChimeraX soft lighting + black silhouettes on white.
SUBSITES = ['-3', '-2', '-1', '+1', '+2', '+3']
SUGAR_C = '#F4B942'            # fig2_v3 SUGAR_C
SUGAR_EDGE = '#B87800'         # fig2_v3 darker orange text colour, used for point edges / ring outlines
CTX_C = '#C8C8C8'              # sugar units outside -3..+3 (fig2_v3 CTX_C)
OUTSIDE_C = CTX_C
O_COL = '#E03C31'              # fig2_v3 O_COL
N_COL = '#2E6DB4'              # fig2_v3 N_COL
CAT_C = '#B0308F'              # fig2_v3 CAT_USED (catalytic carbons)
CARTOON = '#DCDCDC'            # fig2_v3 CARTOON
SURF_C = '#D5D9DE'             # fig2_v3 thumbnail surface
# domains: the grey of the cartoon / surface with a faint tint, so the domains can be told apart
DOM_COL = {'A': '#D5D9DE', 'B': '#E6DCC8', 'C': '#D3E0CF', 'N': '#E3D3DE', 'C2': '#E4E4E4'}
FONT = 'Liberation Sans'
FONT_FILE = '/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf'
FS = 7.0          # base (fig2_v3)
FS_RES = 7.5      # residue labels (fig2_v3 FS_RES)
FS_SUB = 8.0      # subsite numbers, bold (fig2_v3 FS_SUB)
FS_LET = 9.0      # panel letters, bold
FS_TITLE = 7.5    # panel titles


def sub_label(s):
    """'-1' -> '−1' (true minus sign) for printing."""
    return s.replace('-', '−')


def hexrgb(h):
    h = h.lstrip('#')
    return [round(int(h[i:i + 2], 16) / 255, 4) for i in (0, 2, 4)]


def load(pdb):
    st = gemmi.read_structure(str(CIF / f'{pdb}.cif'))
    st.setup_entities()
    st.remove_hydrogens()
    st.remove_alternative_conformations()
    while len(st) > 1:
        del st[1]
    return st


def units_of(pdb):
    """Rows of 33_subsites/units.tsv for one entry."""
    return [r for r in csv.DictReader(open(UNITS33), delimiter='\t') if r['pdb'] == pdb]


def find_res(model, chain, seq, name=None):
    for ch in model:
        if ch.name == chain:
            for r in ch:
                if r.seqid.num == int(seq) and (name is None or r.name == name):
                    return r
    raise KeyError((chain, seq, name))


def unit_atoms(model, row):
    """Atoms (name -> xyz) of one sugar unit (one ring) and its ring atoms, from a units.tsv row.
    Atoms shared between rings by graph distance (a33 a2u) go to every ring at that distance."""
    res = find_res(model, row['chain'], row['seq'], row['comp'])
    units, a2u, _, _ = ring_units(res)
    k = next(i for i, u in enumerate(units) if u['ref'] == row['ref_atom'])
    pos = {a.name: np.array(a.pos.tolist()) for a in res}
    atoms = {n: pos[n] for n, ks in a2u.items() if k in ks}
    ring = {n: pos[n] for n in units[k]['atoms']}
    return atoms, ring


def superpositions():
    """(family, pdb, chain) -> (R, t) with x_frame = R x + t."""
    out = {}
    for r in csv.DictReader(open(SUP), delimiter='\t'):
        if not r['status'].startswith('accepted'):
            continue
        v = [float(x) for x in r['R11..R33_t1..t3'].split()]
        out[(r['family'], r['pdb'], r['chain'])] = (np.array(v[:9]).reshape(3, 3), np.array(v[9:]))
    return out
