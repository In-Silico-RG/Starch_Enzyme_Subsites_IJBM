"""Chemistry panel of the Figure 1 candidates: (i) an alpha-1,4-linked maltohexaose segment spanning subsites
-3..+3 (non-reducing end at left, reducing end at right); (ii) one glucose unit of the chain with R at O2, O3 and
O6 and the four common starch substituents.

Every 2D structure is an RDKit depiction (rdDepictor, CoordGen) of a SMILES built here; no coordinates are
placed by hand. The chain is assembled with Chem.molzip from one alpha-D-glucopyranose template
(OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O; CIP C1 S, C2 R, C3 S, C4 S, C5 R checked at run time), so every
stereocentre is the template's. Unit k's C1 is joined to O4 of unit k+1: the glycosidic oxygen between -1 and +1
is O4 of the +1 unit (the scissile bond is C1(-1)-O4(+1)).
Substituents (attachment point *): hydroxypropyl *CC(O)C; acetyl *C(C)=O; octenyl succinyl (half-ester of
2-octenylsuccinic anhydride, ester at the ring-distal carbonyl as written by the author)
*C(=O)CC(C/C=C/CCCCC)C(=O)O; carboxymethyl *CC(=O)O.

Images are drawn by MolDraw2DCairo with Liberation Sans at the pixel size that matches the target point
size of the figure; atom drawing coordinates are returned so matplotlib can place the subsite marks.
"""
import math
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDepictor, AllChem
from rdkit.Chem.Draw import rdMolDraw2D

import fig1_protein_common as C

FONT = C.FONT_FILE      # Liberation Sans, as the rest of the figure (fig2_v3)
GLC = 'OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O'
# template with C1 -> [*:a] (O1 removed) and O4-H -> O4-[*:b]; atom order O6 C6 C5 O5 C1 * C2 O2 C3 O3 C4 O4 *
NAMES = ['O6', 'C6', 'C5', 'O5', 'C1', None, 'C2', 'O2', 'C3', 'O3', 'C4', 'O4', None]
SUBST = [('hydroxypropyl', '*CC(O)C'), ('acetyl', '*C(C)=O'),
         ('octenyl succinyl', '*C(=O)CC(C/C=C/CCCCC)C(=O)O'), ('carboxymethyl', '*CC(=O)O')]


def check_template():
    m = Chem.MolFromSmiles(GLC)
    from rdkit.Chem import rdCIPLabeler
    rdCIPLabeler.AssignCIPLabels(m)
    order = ['O6', 'C6', 'C5', 'O5', 'C1', 'O1', 'C2', 'O2', 'C3', 'O3', 'C4', 'O4']
    lab = {order[i]: a.GetProp('_CIPCode') for i, a in enumerate(m.GetAtoms()) if a.HasProp('_CIPCode')}
    assert lab == {'C5': 'R', 'C1': 'S', 'C2': 'R', 'C3': 'S', 'C4': 'S'}, lab
    return lab


def unit_smiles(k, n):
    """alpha-D-Glcp unit k of n (k = 0 non-reducing end). C1 bonded to [*:k+1] (next unit's O4) unless reducing
    end (free alpha O1); O4 bonded to [*:k] (previous unit's C1) unless non-reducing end."""
    c1 = f'[*:{k + 1}]' if k < n - 1 else 'O'
    o4 = f'O[*:{k}]' if k > 0 else 'O'
    return f'OC[C@H]1O[C@H]({c1})[C@H](O)[C@@H](O)[C@@H]1{o4}'


def build_chain(n=6):
    frags = []
    for k in range(n):
        m = Chem.MolFromSmiles(unit_smiles(k, n))
        names = ['O6', 'C6', 'C5', 'O5', 'C1', 'X1', 'C2', 'O2', 'C3', 'O3', 'C4', 'O4', 'X4']
        idx = 0
        for a in m.GetAtoms():
            nm = names[idx]
            if k == n - 1 and nm == 'X1':
                nm = 'O1'
            if k == 0 and nm == 'O4' and idx == 11:
                pass
            a.SetProp('nm', nm)
            a.SetIntProp('unit', k)
            idx += 1
            if k == 0 and idx == 12:      # non-reducing end has no X4
                break
        frags.append(m)
    combo = frags[0]
    for f in frags[1:]:
        combo = Chem.CombineMols(combo, f)
    mol = Chem.molzip(combo)
    Chem.SanitizeMol(mol)
    return mol


def depict(mol, straighten_units=None):
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(mol)
    if straighten_units:
        conf = mol.GetConformer()
        X = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y] for i in range(mol.GetNumAtoms())])
        cen = np.array([X[[a.GetIdx() for a in mol.GetAtoms() if a.HasProp('unit') and a.GetIntProp('unit') == k
                            and a.GetProp('nm') in ('C1', 'C2', 'C3', 'C4', 'C5', 'O5')]].mean(0)
                        for k in range(straighten_units)])
        v = np.polyfit(cen[:, 0], cen[:, 1], 1)
        d = cen[-1] - cen[0]
        ang = math.atan2(d[1], d[0])
        R = np.array([[math.cos(-ang), -math.sin(-ang)], [math.sin(-ang), math.cos(-ang)]])
        X = (X - X.mean(0)) @ R.T
        for i in range(mol.GetNumAtoms()):
            conf.SetAtomPosition(i, (float(X[i, 0]), float(X[i, 1]), 0.0))
    return mol


def draw(mol, w, h, font_px, png, bond_px=2.0, padding=0.04, atom_labels=None, highlight=None):
    d = rdMolDraw2D.MolDraw2DCairo(w, h)
    o = d.drawOptions()
    o.fontFile = FONT
    o.fixedFontSize = int(font_px)
    o.minFontSize = int(font_px)
    o.maxFontSize = int(font_px)
    o.bondLineWidth = bond_px
    o.scaleBondWidth = False
    o.padding = padding
    o.clearBackground = True       # white; a transparent Cairo PNG gives hollow strokes after resampling
    o.annotationFontScale = 1.0
    o.additionalAtomLabelPadding = 0.08
    o.useBWAtomPalette()
    o.setAtomPalette({-1: (0, 0, 0), 0: (0, 0, 0), 1: (0, 0, 0), 6: (0, 0, 0), 7: tuple(C.hexrgb(C.N_COL)),
                      8: tuple(C.hexrgb(C.O_COL))})
    if atom_labels:
        for i, lab in atom_labels.items():
            o.atomLabels[i] = lab
    rdMolDraw2D.PrepareMolForDrawing  # noqa (kept explicit: no wedge recomputation beyond RDKit defaults)
    if highlight:
        d.DrawMolecule(mol, highlightAtoms=list(highlight), highlightAtomColors=highlight,
                       highlightBonds=[])
    else:
        d.DrawMolecule(mol)
    d.FinishDrawing()
    Path(png).write_bytes(d.GetDrawingText())
    xy = {}
    for a in mol.GetAtoms():
        p = d.GetDrawCoords(a.GetIdx())
        xy[a.GetIdx()] = (p.x, p.y)
    return xy


def chain_panel(w, h, font_px, png, bond_px):
    check_template()
    mol = depict(build_chain(6), straighten_units=6)
    for a in mol.GetAtoms():               # O2, O3, O6 named on unit -2 (k = 1) as RDKit atom notes
        if a.GetIntProp('unit') == 1 and a.GetProp('nm') in ('O2', 'O3', 'O6'):
            a.SetProp('atomNote', a.GetProp('nm'))
    xy = draw(mol, w, h, font_px, png, bond_px)
    info = dict(units={}, atoms={})
    for a in mol.GetAtoms():
        k = a.GetIntProp('unit'); nm = a.GetProp('nm')
        info['atoms'][(k, nm)] = xy[a.GetIdx()]
    for k in range(6):
        ring = [info['atoms'][(k, n)] for n in ('C1', 'C2', 'C3', 'C4', 'C5', 'O5')]
        info['units'][C.SUBSITES[k]] = dict(ring=ring, cen=tuple(np.mean(ring, 0)))
    # scissile bond: C1 of -1 (unit 2) to O4 of +1 (unit 3)
    info['cut'] = (info['atoms'][(2, 'C1')], info['atoms'][(3, 'O4')])
    return info


def glucose_R_panel(w, h, font_px, png, bond_px):
    """one glucose unit of the chain, chain continuing at C1 and O4, R at O2, O3, O6"""
    smi = '[*:6]OC[C@H]1O[C@H]([*:1])[C@H](O[*:2])[C@@H](O[*:3])[C@@H]1O[*:4]'
    m = Chem.MolFromSmiles(smi)
    labels = {}
    names = {}
    for a in m.GetAtoms():
        mp = a.GetAtomMapNum()
        if mp:
            a.SetAtomMapNum(0)
            if mp in (2, 3, 6):
                labels[a.GetIdx()] = 'R'
                names[a.GetIdx()] = f'R{mp}'
            else:
                labels[a.GetIdx()] = ' '   # chain continues (bond stub)
                names[a.GetIdx()] = 'chain'
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(m)
    xy = draw(m, w, h, font_px, png, bond_px, atom_labels=labels)
    # positions of the ring O2/O3/O6 (neighbours of the R atoms)
    opos = {}
    for i, nm in names.items():
        if nm.startswith('R'):
            o = m.GetAtomWithIdx(i).GetNeighbors()[0].GetIdx()
            opos['O' + nm[1:]] = xy[o]
            opos[nm] = xy[i]
    return opos


def substituent(smi, w, h, font_px, png, bond_px):
    m = Chem.MolFromSmiles(smi)
    labels = {a.GetIdx(): 'R' if False else 'O' for a in m.GetAtoms() if a.GetAtomicNum() == 0}
    # the attachment atom is the ring oxygen O2/O3/O6: draw it as "O" bonded to the substituent, with a stub
    em = Chem.RWMol(m)
    for i in list(labels):
        em.GetAtomWithIdx(i).SetAtomicNum(8)
        stub = em.AddAtom(Chem.Atom(0))
        em.AddBond(i, stub, Chem.BondType.SINGLE)
    m2 = em.GetMol()
    Chem.SanitizeMol(m2)
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(m2)
    lab = {a.GetIdx(): ' ' for a in m2.GetAtoms() if a.GetAtomicNum() == 0}
    xy = draw(m2, w, h, font_px, png, bond_px, atom_labels=lab)
    return xy


if __name__ == '__main__':
    print(check_template())
    mol = build_chain(6)
    print(Chem.MolToSmiles(mol))
    print(Chem.rdMolDescriptors.CalcMolFormula(mol))
