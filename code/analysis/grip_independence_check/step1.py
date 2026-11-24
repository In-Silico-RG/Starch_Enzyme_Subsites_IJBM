"""Select -1 entries exactly as a33_subsites.py, assign one protein per entry, dump chain sequences."""
import sqlite3, json, sys
import pandas as pd, gemmi
P=__import__('os').path.abspath(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '..', '..', '..')) + '/'
S=sys.argv[1]
con=sqlite3.connect(f'file:{P}data/db/enzimas_experimental.sqlite?mode=ro',uri=True)
entry=pd.read_sql('select pdb,family,organism,protein,uniprot,resolution_A from entry',con).set_index('pdb')
u=pd.read_csv(P+'data/analysis/33_subsites/units.tsv',sep='\t',dtype={'subsite':str,'frame_chain':str})
u['subsite']=u.subsite.fillna('')
rep=u[(u.subsite=='-1')&(u.rscc>=0.80)&u.in_tables&u.family.isin(['GH13','GH31'])]
rows=[]
sifts=pd.read_sql('select distinct pdb,chain,uniprot from sifts_residue',con)
for (f,p),g in rep.groupby(['family','pdb']):
    fc=g.frame_chain.iloc[0]
    ss=sorted(set(sifts[(sifts.pdb==p)&(sifts.chain==fc)].uniprot.dropna()))
    eu=entry.loc[p,'uniprot'] or ''
    if len(ss)==1: acc,src=ss[0],'SIFTS frame chain'
    elif eu.strip() and ';' not in eu: acc,src=eu.strip(),'entry table'
    else: acc,src='','none'
    st=gemmi.read_structure(P+f'data/census/cif/{p}.cif'); st.setup_entities()
    seq=''
    for ent in st.entities:
        if fc in ent.subchains or any(sc.startswith(fc) for sc in ent.subchains):
            pass
    # chain sequence: from entity whose subchain belongs to frame chain
    ch=st[0][fc]
    sc=ch.get_polymer().subchain_id()
    ent=st.get_entity_of(ch.get_polymer())
    seq=gemmi.one_letter_code(ent.full_sequence) if ent and ent.full_sequence else gemmi.one_letter_code([r.name for r in ch.get_polymer()])
    rows.append(dict(family=f,pdb=p,frame_chain=fc,acc=acc,acc_source=src,entry_uniprot=eu,sifts_frame_chain=';'.join(ss),
                     organism=entry.loc[p,'organism'],protein=entry.loc[p,'protein'],resolution_A=entry.loc[p,'resolution_A'],
                     chain_seq=seq.upper().replace('-','')))
d=pd.DataFrame(rows); d.to_csv(S+'/entries.tsv',sep='\t',index=False)
print(d.groupby('family').size()); print(d[d.acc_source!='SIFTS frame chain'][['pdb','acc','acc_source','entry_uniprot','organism']])
print(d.groupby('family').acc.nunique())
