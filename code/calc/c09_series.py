#!/usr/bin/env python3
"""Calculation 9, congeneric series: MM/GBSA of the iminosugar inhibitors of human lysosomal alpha-glucosidase
(GAA, UniProt P10253, ChEMBL target CHEMBL2608) whose complexes are deposited and whose potency against GAA is in
ChEMBL (tables chembl_pdb_ligand, chembl_activity of the database, read-only).

Usage (GUANE, Amber env with RDKit/gemmi/propka): c09_series.py PDB COMP     e.g. c09_series.py 5NN5 NOJ

The problem: one number per inhibitor that ranks binding, from the deposited pose, to set against pChEMBL.
Method (one structure per complex, no MD: the 07 docking poses do not exist and no node is free for MD of four
~900-residue complexes):
 1. chain A of the entry, standard residues; Ca/Na/Mg/Zn within 3.5 A of it; the inhibitor (COMP); waters, other
    ligands and N-glycans removed. Protonation as in c08_build.py (PROPKA pH 7, metal-bound carboxylates charged,
    His tautomer from reduce, disulfides SG-SG < 2.5 A).
 2. Inhibitor: deposited heavy atoms, bond orders from the component SMILES, ring nitrogen protonated (+1, the
    iminium paired with the catalytic carboxylates), hydrogens by RDKit; GAFF2 types, AM1-BCC charges (no geometry
    optimisation).
 3. ff14SB + GAFF2, mbondi2 radii; minimisation in GB (igb 5, 0.15 M salt, cut-off 16 A, Born radii from atoms
    within 15 A), 300 steepest-descent + 700 conjugate-gradient steps (the first trial with no cut-off and 2000 steps
    ran at ~9 s per step on one legacy core), protein heavy atoms restrained to the
    deposited positions with 5 kcal/mol/A^2 except residues within 6 A of the inhibitor (free), inhibitor free.
 4. MMPBSA.py on the minimised complex: GB igb 5, saltcon 0.150; per-residue decomposition idecomp 1 for residues
    within 6 A of the inhibitor.
Output: data/calc/09_mmgbsa/series/<PDB>_<COMP>/ (FINAL_RESULTS_MMPBSA.dat, FINAL_DECOMP_MMPBSA.dat, logs).
"""
import os, sys, re, subprocess
import numpy as np
import gemmi
from rdkit import Chem
from rdkit.Chem import AllChem
sys.path.insert(0, os.path.dirname(__file__))
from common import db, CIF, CALC

PDB, COMP = sys.argv[1], sys.argv[2]
ENV = os.environ.get('AMBER_ENV', os.path.expanduser('~/miniconda3/envs/sc_modeling'))
E = dict(os.environ, AMBERHOME=ENV, PATH=f'{ENV}/bin:' + os.environ['PATH'])
OUT = os.path.join(CALC, '09_mmgbsa', 'series', f'{PDB}_{COMP}')
os.makedirs(OUT, exist_ok=True)
os.chdir(OUT)
log = open('build.log', 'w', buffering=1)


def say(*a):
    s = ' '.join(map(str, a)); print(s); log.write(s + '\n')


def run(cmd):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True, env=E)
    log.write(f'$ {cmd}\n{r.stdout[-2000:]}{r.stderr[-2000:]}\n')
    if r.returncode:
        raise SystemExit(f'failed: {cmd}\n{r.stdout[-2000:]}{r.stderr[-2000:]}')
    return r.stdout


c = db()
smi = c.execute('select smiles from component where comp_id=?', (COMP,)).fetchone()[0]
st = gemmi.read_structure(os.path.join(CIF, PDB + '.cif'))
st.setup_entities(); st.remove_hydrogens(); st.remove_alternative_conformations()
model = st[0]
AA = set(gemmi.expand_one_letter_sequence('ACDEFGHIKLMNPQRSTVWY', gemmi.ResidueKind.AA))
lig = [r for ch in model for r in ch if r.name == COMP and ch.name == 'A']
assert len(lig) >= 1, 'ligand not found in chain A'
lig = lig[0]
prot = [r for r in model['A'] if r.name in AA and r.entity_type == gemmi.EntityType.Polymer]
LX = np.array([a.pos.tolist() for a in lig])
PX = np.array([a.pos.tolist() for r in prot for a in r])
ions = [r for ch in model for r in ch if r.name in ('CA', 'NA', 'MG', 'ZN') and len(r) == 1
        and np.min(np.linalg.norm(PX - np.array(r[0].pos.tolist()), axis=1)) < 3.5]
say(PDB, COMP, 'protein residues', len(prot), '| ions', [(i.name, i.seqid.num) for i in ions], '| ligand atoms', len(lig))

# protonation (same rules as c08_build.py)
metal = [i[0].pos for i in ions]
ss = []
cys = [(k, r) for k, r in enumerate(prot) if r.name == 'CYS' and r.find_atom('SG', '*')]
for i in range(len(cys)):
    for j in range(i + 1, len(cys)):
        if cys[i][1]['SG'][0].pos.dist(cys[j][1]['SG'][0].pos) < 2.5:
            ss.append((cys[i][0], cys[j][0]))
mcarb = {k for k, r in enumerate(prot) if r.name in ('ASP', 'GLU') and any(a.name in ('OD1', 'OD2', 'OE1', 'OE2') and any(a.pos.dist(m) < 3.0 for m in metal) for a in r)}


def pdbline(rec, n, name, resn, ch, resi, xyz, el):
    nm = name if len(name) == 4 else ' ' + name.ljust(3)
    return '%-6s%5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n' % (rec, n % 100000, nm, resn[:3], ch, resi % 10000, *xyz, el)


with open('propka_in.pdb', 'w') as f:
    n = 0
    for k, r in enumerate(prot):
        for a in r:
            n += 1; f.write(pdbline('ATOM', n, a.name, r.name, 'A', k + 1, a.pos.tolist(), a.element.name))
    for a in lig:
        n += 1; f.write(pdbline('HETATM', n, a.name, 'LIG', 'L', 1, a.pos.tolist(), a.element.name))
    f.write('END\n')
subprocess.run([sys.executable, '-m', 'propka', '--quiet', 'propka_in.pdb'], check=True, capture_output=True)
pka = {}
for line in open('propka_in.pka'):
    m = re.match(r'\s+(ASP|GLU|HIS|LYS)\s+(\d+)\s+A\s+(-?\d+\.\d+)\s+', line)
    if m:
        pka.setdefault(int(m.group(2)), float(m.group(3)))
with open('reduce_in.pdb', 'w') as f:
    for k, r in enumerate(prot):
        for a in r:
            f.write(pdbline('ATOM', 1, a.name, r.name, 'A', k + 1, a.pos.tolist(), a.element.name))
# reduce: from the env if present; AmberTools 26 on GUANE has no reduce, so reduce_out.pdb is made on the laptop
# (AmberTools 24 reduce, same input written by this script) and placed in the output folder beforehand.
if os.path.exists(f'{ENV}/bin/reduce'):
    red = subprocess.run([f'{ENV}/bin/reduce', '-build', '-quiet', 'reduce_in.pdb'], capture_output=True, text=True, env=E).stdout
else:
    red = open('reduce_out.pdb').read()
assert 'HIS' in red or not any(r.name == 'HIS' for r in prot), 'no reduce output'
hisH = {}
for line in red.splitlines():
    if line.startswith('ATOM') and line[17:20] == 'HIS':
        hisH.setdefault(int(line[22:26]), set()).add(line[12:16].strip())
names, changes = [], []
for k, r in enumerate(prot):
    nm, i = r.name, k + 1
    if nm == 'HIS':
        h = hisH.get(i, set()); nm = 'HID' if ('HD1' in h and 'HE2' not in h) else 'HIE'
        if pka.get(i, 0) > 7.0:
            nm = 'HIP'
    elif nm == 'ASP' and pka.get(i, 0) > 7.0 and k not in mcarb:
        nm = 'ASH'
    elif nm == 'GLU' and pka.get(i, 0) > 7.0 and k not in mcarb:
        nm = 'GLH'
    elif nm == 'LYS' and pka.get(i, 99) < 7.0:
        nm = 'LYN'
    elif nm == 'CYS' and any(k in p for p in ss):
        nm = 'CYX'
    names.append(nm)
    if nm not in (r.name, 'HID', 'HIE', 'CYX'):
        changes.append(f'{r.name}{r.seqid.num}->{nm} ({pka.get(i)})')
say('protonation changes (pH 7):', changes, '| disulfides', len(ss))

# ligand
blk = ''.join(pdbline('HETATM', k + 1, a.name, 'LIG', 'L', 1, a.pos.tolist(), a.element.name) for k, a in enumerate(lig)) + 'END\n'
m = Chem.MolFromPDBBlock(blk, removeHs=True, proximityBonding=True)
tmpl = Chem.MolFromSmiles(smi)
m = AllChem.AssignBondOrdersFromTemplate(tmpl, m)
Chem.AssignStereochemistryFrom3D(m)
ringN = [a for a in m.GetAtoms() if a.GetSymbol() == 'N' and a.IsInRing()]
assert len(ringN) == 1, 'expected one ring N'
ringN[0].SetFormalCharge(1); ringN[0].SetNoImplicit(False); Chem.SanitizeMol(m)
mh = Chem.AddHs(m, addCoords=True)
pos = mh.GetConformer().GetPositions()
hn = 0
with open('lig_in.mol2', 'w') as f:
    bonds = list(mh.GetBonds())
    f.write(f'@<TRIPOS>MOLECULE\nLIG\n{mh.GetNumAtoms()} {len(bonds)} 1 0 0\nSMALL\nNO_CHARGES\n\n@<TRIPOS>ATOM\n')
    for a in mh.GetAtoms():
        if a.GetSymbol() == 'H':
            hn += 1; nm = f'H{hn}'
        else:
            nm = a.GetPDBResidueInfo().GetName().strip()
        p = pos[a.GetIdx()]
        f.write(f'{a.GetIdx() + 1:6d} {nm:<5s} {p[0]:10.4f}{p[1]:10.4f}{p[2]:10.4f} {a.GetSymbol():<5s} 1 LIG  0.0000\n')
    f.write('@<TRIPOS>BOND\n')
    for k, b in enumerate(bonds, 1):
        bt = {Chem.BondType.SINGLE: '1', Chem.BondType.DOUBLE: '2', Chem.BondType.AROMATIC: 'ar'}[b.GetBondType()]
        f.write(f'{k:6d} {b.GetBeginAtomIdx() + 1:5d} {b.GetEndAtomIdx() + 1:5d} {bt}\n')
say('ligand', Chem.MolToSmiles(Chem.RemoveHs(mh)), 'charge', Chem.GetFormalCharge(mh))
run('antechamber -i lig_in.mol2 -fi mol2 -o lig.mol2 -fo mol2 -c bcc -nc 1 -at gaff2 -rn LIG '
    '-ek "qm_theory=\'AM1\', maxcyc=0, scfconv=1.d-10" -pf y -dr no')
run('parmchk2 -i lig.mol2 -f mol2 -o lig.frcmod -s 2')

# complex PDB
lines, n, resi = [], 0, 0
for k, r in enumerate(prot):
    resi += 1
    for a in r:
        if a.name == 'OXT':
            continue
        n += 1; lines.append(pdbline('ATOM', n, a.name, names[k], 'A', resi, a.pos.tolist(), a.element.name))
    if k + 1 < len(prot) and prot[k + 1].seqid.num != r.seqid.num + 1:
        lines.append('TER\n')                       # chain break: separate termini
lines.append('TER\n')
ion_res = []
for r in ions:
    resi += 1; n += 1; ion_res.append(resi)
    lines.append(pdbline('HETATM', n, r.name, r.name, 'A', resi, r[0].pos.tolist(), r[0].element.name)); lines.append('TER\n')
resi += 1; lig_resi = resi
open('complex_in.pdb', 'w').writelines(lines + ['END\n'])
bcmd = '\n'.join(f'bond c.{a + 1}.SG c.{b + 1}.SG' for a, b in ss)
leap = f"""source leaprc.protein.ff14SB
source leaprc.gaff2
source leaprc.water.tip3p
set default PBradii mbondi2
loadamberparams lig.frcmod
LIG = loadmol2 lig.mol2
r = loadpdb complex_in.pdb
{bcmd.replace('c.', 'r.')}
c = combine {{ r LIG }}
saveamberparm r receptor.prmtop receptor.inpcrd
saveamberparm LIG ligand.prmtop ligand.inpcrd
saveamberparm c complex.prmtop complex.inpcrd
savepdb c complex_start.pdb
quit
"""
open('tleap.in', 'w').write(leap)
o = run('tleap -f tleap.in')
assert 'Exiting LEaP: Errors = 0' in o, o[-2000:]
gaps = o.count('TER')
say('tleap ok')

# minimisation in GB; restraints on protein heavy atoms farther than 6 A from the inhibitor
cpx = gemmi.read_structure('complex_start.pdb')
LXc = np.array([a.pos.tolist() for ch in cpx[0] for r in ch if r.name == 'LIG' for a in r])
free = sorted({r.seqid.num for ch in cpx[0] for r in ch if r.name != 'LIG'
               and min(np.linalg.norm(LXc - np.array(a.pos.tolist()), axis=1).min() for a in r) < 6.0})
nres = sum(1 for ch in cpx[0] for r in ch)
restr = [x for x in range(1, nres + 1) if x not in free and x not in (nres,)]
def ranges(v):
    out, s = [], v[0]
    for a, b in zip(v, v[1:] + [None]):
        if b != a + 1:
            out.append(f'{s}-{a}' if s != a else f'{s}'); s = b
    return ','.join(out)
mask = f':{ranges(restr)}&!@H='
open('min.in', 'w').write(f"""GB minimisation, protein heavy atoms > 6 A from the inhibitor restrained
 &cntrl
  imin=1, maxcyc=1000, ncyc=300, ntb=0, igb=5, saltcon=0.150, cut=16.0, rgbmax=15.0, ntpr=100,
  ntr=1, restraint_wt=5.0, restraintmask='{mask}',
 /
""")
run('sander -O -i min.in -o min.out -p complex.prmtop -c complex.inpcrd -r min.rst7 -ref complex.inpcrd')
run('cpptraj -p complex.prmtop -y min.rst7 -x min.nc')
pres = ranges(sorted(set(free + [nres])))
open('mmgbsa.in', 'w').write(f"""MM/GBSA, single minimised structure
&general
  keep_files=0, verbose=1,
/
&gb
  igb=5, saltcon=0.150,
/
&decomp
  idecomp=1, dec_verbose=0, print_res="{pres}",
/
""")
run('MMPBSA.py -O -i mmgbsa.in -o FINAL_RESULTS_MMPBSA.dat -do FINAL_DECOMP_MMPBSA.dat -cp complex.prmtop -rp receptor.prmtop -lp ligand.prmtop -y min.nc')
txt = open('FINAL_RESULTS_MMPBSA.dat').read()
dg = re.search(r'DELTA TOTAL\s+(-?\d+\.\d+)', txt)
say('free residues (within 6 A):', len(free), '| DELTA TOTAL (kcal/mol):', dg.group(1) if dg else None)
