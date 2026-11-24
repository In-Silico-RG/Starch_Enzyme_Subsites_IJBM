"""Calculation 1, collection: PROPKA summaries of all entries into
  data/calc/01_propka/propka_all.tsv        every titratable group
  data/calc/01_propka/propka_near_ligand.tsv groups of residues within 4.0 A of a sugar/'other'
                                                   instance (database contact, partner_kind='protein')
UniProt numbers from sifts_residue; catalytic flags from 03_frames/catalytic_residues_mapped.tsv
(reference catalytic residues of the family carried to each chain by the sequence alignment of
calculation 3; flagged when the mapped residue has the reference type and its functional atoms
lie within 3.0 A of the reference ones after superposition)."""
import os, sys, csv, collections
sys.path.insert(0, os.path.dirname(__file__))
from common import CALC, db

D = os.path.join(CALC, '01_propka')
con = db()
fam = dict(con.execute("select pdb, family from entry"))
sifts = {}
for p, ch, s, ins, un, unum in con.execute("select pdb, chain, auth_seq, ins_code, uniprot, uniprot_resnum from sifts_residue"):
    sifts[(p, ch, int(s), ins or '')] = (un, unum)
cat = {}
cf = os.path.join(CALC, '03_frames', 'catalytic_residues_mapped.tsv')
for r in csv.DictReader(open(cf), delimiter='\t'):
    if r['role'] == 'site' or not r['auth_seq']:
        continue
    ok = r['type_check'] == 'same type' and float(r['functional_atom_distance_to_reference_A']) <= 3.0
    cat[(r['pdb'], r['chain'], int(r['auth_seq']), r['ins_code'])] = (
        r['role'] + ('' if ok else '?'), r['reference_residue'], r['functional_atom_distance_to_reference_A'])
near = collections.defaultdict(list)
for iid, ch, s, ins, comp, dmin, p, lcomp, lch, lseq, cls in con.execute(
        """select c.instance_id, c.chain, c.seq, c.ins_code, c.comp, c.min_distance_A, i.pdb, i.comp_id,
           i.auth_asym, i.auth_seq, i.ligand_class from contact c join instance i on c.instance_id=i.instance_id
           where i.ligand_class in ('sugar','other') and c.partner_kind='protein'"""):
    near[(p, ch, int(s), ins or '')].append((dmin, '%s:%s%d(%s,%.2f)' % (lch, lcomp, lseq, cls, dmin)))
H = ['pdb', 'family', 'chain', 'resname', 'auth_seq', 'uniprot', 'uniprot_resnum', 'group_type', 'pKa',
     'model_pKa', 'shift', 'catalytic_role', 'catalytic_reference_residue', 'catalytic_superposition_distance_A',
     'min_distance_to_sugar_other_A', 'sugar_other_ligands_within_4A']
rows, missing = [], []
hetidx = {}
for p in sorted(fam):
    f = os.path.join(D, 'runs', p, p + '.pka')
    if not os.path.exists(f):
        missing.append(p); continue
    cmap = dict(l.rstrip('\n').split('\t') for l in list(open(os.path.join(D, 'runs', p, 'chainmap.tsv')))[1:])
    on = False
    for line in open(f):
        if line.startswith('SUMMARY OF THIS PREDICTION'):
            on = True; continue
        if on and line.startswith('-----'):
            break
        if not on:
            continue
        t = line.split()
        if len(t) < 5 or t[0] == 'Group':
            continue
        if t[1].lstrip('-').isdigit():
            resname, num, ch, pka, mpka = t[0], int(t[1]), t[2], float(t[3]), float(t[4])
            gtype = resname
        else:
            # ligand group: PROPKA names it by component, atom and chain only; the residue number is
            # read from the input PDB (several numbers joined by ';' when the atom name is not unique)
            resname, atom, ch, pka, mpka = t[0], t[1], t[2], float(t[3]), float(t[4])
            if hetidx.get(p) is None:
                hetidx[p] = collections.defaultdict(list)
                for l in open(os.path.join(D, 'runs', p, p + '.pdb')):
                    if l.startswith('HETATM'):
                        hetidx[p][(l[17:20].strip(), l[21], l[12:16].strip())].append(l[22:27].strip())
            nums = hetidx[p].get((resname, ch, atom), [])
            num = ';'.join(nums)
            gtype = 'ligand atom %s (%s)' % (atom, t[5] if len(t) > 5 else '')
        ch = cmap.get(ch, ch)
        key = (p, ch, int(num) if isinstance(num, int) else num, '')
        un = sifts.get(key, ('', ''))
        c = cat.get(key, ('', '', ''))
        nl = sorted(near.get(key, []))
        rows.append([p, fam[p], ch, resname, num, un[0], un[1], gtype, '%.2f' % pka, '%.2f' % mpka, '%.2f' % (pka - mpka),
                     c[0], c[1], c[2], '%.2f' % nl[0][0] if nl else '', ';'.join(x[1] for x in nl)])
with open(os.path.join(D, 'propka_all.tsv'), 'w') as f:
    f.write('\t'.join(H) + '\n')
    for r in rows:
        f.write('\t'.join(map(str, r)) + '\n')
with open(os.path.join(D, 'propka_near_ligand.tsv'), 'w') as f:
    f.write('\t'.join(H) + '\n')
    for r in rows:
        if r[14] != '' and r[3] not in ('N+', 'C-') and not r[7].startswith('ligand'):
            f.write('\t'.join(map(str, r)) + '\n')
n_near = sum(1 for r in rows if r[14] != '' and r[3] not in ('N+', 'C-') and not r[7].startswith('ligand'))
n_cat = sum(1 for r in rows if r[11])
print('entries with pka', len(fam) - len(missing), 'missing', missing)
print('groups', len(rows), 'near ligand', n_near, 'catalytic-flagged', n_cat,
      'with uniprot', sum(1 for r in rows if r[6] != ''))
