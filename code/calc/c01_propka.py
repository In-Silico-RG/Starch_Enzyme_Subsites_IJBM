"""Calculation 1: PROPKA 3 on one entry. Usage: c01_propka.py PDB [PDB ...]
Input: model 1 of the cached mmCIF, hydrogens and waters removed, first altloc kept,
protein and hetero groups (ligands) as deposited. Output per entry in
data/calc/01_propka/runs/<PDB>/ : <PDB>.pdb (input), <PDB>.pka (PROPKA), chainmap.tsv."""
import os, sys, subprocess
sys.path.insert(0, os.path.dirname(__file__))
from common import load, single_letter_chains, CALC

OUT = os.path.join(CALC, '01_propka', 'runs')


def run(pdb):
    d = os.path.join(OUT, pdb)
    os.makedirs(d, exist_ok=True)
    if os.path.exists(os.path.join(d, pdb + '.pka')):
        return pdb, 'done'
    st = load(pdb)
    st.remove_waters()
    st.remove_empty_chains()
    cmap = single_letter_chains(st)
    with open(os.path.join(d, 'chainmap.tsv'), 'w') as f:
        f.write('propka_chain\tauth_chain\n')
        for k, v in cmap.items():
            f.write('%s\t%s\n' % (k, v))
    st.write_pdb(os.path.join(d, pdb + '.pdb'))
    r = subprocess.run([sys.executable, '-m', 'propka', '--quiet', pdb + '.pdb'], cwd=d,
                       capture_output=True, text=True)
    with open(os.path.join(d, 'propka.log'), 'w') as f:
        f.write(r.stdout + r.stderr)
    ok = os.path.exists(os.path.join(d, pdb + '.pka'))
    return pdb, 'ok' if ok else 'FAILED rc=%d' % r.returncode


import gemmi
if __name__ == '__main__':
    for p in sys.argv[1:]:
        print(*run(p), flush=True)
