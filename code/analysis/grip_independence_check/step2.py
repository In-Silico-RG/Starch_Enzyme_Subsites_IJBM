import json, sys, itertools
import numpy as np, pandas as pd
from Bio import Align
from Bio.Align import substitution_matrices
P=__import__('os').path.abspath(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '..', '..', '..')) + '/'
S=sys.argv[1]
d=pd.read_csv(S+'/entries.tsv',sep='\t',keep_default_na=False)
TAX={'Alicyclobacillus sp. 18711':'Bacteria; Bacillota','Bacteroides ovatus ATCC 8483':'Bacteria; Bacteroidota',
     'Sutcliffiella halmapala':'Bacteria; Bacillota','Alkalimonas sp. NCh-2':'Bacteria; Pseudomonadota',
     'Blastobotrys adeninivorans CBS 8244':'Eukaryota; Fungi; Ascomycota'}
def uni(acc):
    j=json.load(open(f'{S}/dl/{acc}.json'))
    if j.get('entryType','').startswith('Inactive'): return None
    lin=[x['scientificName'] for x in j['lineage']] if 'lineage' in j else []
    rank={x['scientificName']:x.get('rank') for x in j.get('lineage',[])}
    sk=lin[0] if lin else ''
    king=[x for x in lin if rank.get(x)=='kingdom']; phy=[x for x in lin if rank.get(x)=='phylum']
    name=j['proteinDescription'].get('recommendedName',j['proteinDescription'].get('submissionNames',[{}])[0]).get('fullName',{}).get('value','') if 'proteinDescription' in j else ''
    return dict(seq=j['sequence']['value'],org=j['organism']['scientificName'],superkingdom=sk,
                kingdom=(king[0] if king else ''),phylum=(phy[0] if phy else ''),uname=name)
# protein key + sequence per entry
keys,seqs,meta={},{},{}
for r in d.itertuples():
    info=uni(r.acc) if r.acc else None
    if info:
        k=r.acc; seqs[k]=info['seq']; meta[k]=dict(source='UniProt '+r.acc,**{x:info[x] for x in ['org','superkingdom','kingdom','phylum','uname']})
    else:
        k='PDB:'+r.pdb; seqs[k]=r.chain_seq
        t=TAX[r.organism].split('; ')
        meta[k]=dict(source=f'chain {r.frame_chain} of {r.pdb} mmCIF'+(f' (UniProt {r.acc} inactive)' if r.acc else ' (no UniProt mapping)'),
                     org=r.organism,superkingdom=t[0],kingdom=(t[1] if len(t)>2 else ''),phylum=t[-1],uname=r.protein)
    keys[r.pdb]=k
d['protein_key']=d.pdb.map(keys)
al=Align.PairwiseAligner(); al.mode='global'; al.substitution_matrix=substitution_matrices.load('BLOSUM62')
al.open_gap_score=-10; al.extend_gap_score=-0.5; al.end_gap_score=0
def ident(a,b):
    aln=al.align(a,b)[0]
    s1,s2=aln[0],aln[1]
    cols=[(x,y) for x,y in zip(s1,s2) if x!='-' and y!='-']
    nid=sum(x==y for x,y in cols)
    return nid/len(cols), nid/min(len(a),len(b)), len(cols)
def single_link(names,M,thr):
    par={n:n for n in names}
    def f(x):
        while par[x]!=x: par[x]=par[par[x]]; x=par[x]
        return x
    for i,j in itertools.combinations(range(len(names)),2):
        if M[i,j]>=thr: par[f(names[i])]=f(names[j])
    roots={}; return {n:roots.setdefault(f(n),len(roots)+1) for n in names}
out={}
for fam in ['GH13','GH31']:
    dd=d[d.family==fam]
    names=sorted(dd.protein_key.unique())
    # merge unaccessioned entries with 100% identical chain sequences
    n=len(names); M=np.eye(n); Ms=np.eye(n); L=np.zeros((n,n),int)
    for i,j in itertools.combinations(range(n),2):
        a,b,c=ident(seqs[names[i]],seqs[names[j]]); M[i,j]=M[j,i]=a; Ms[i,j]=Ms[j,i]=b; L[i,j]=L[j,i]=c
    # chain-sequence check: entry-level
    ents=list(dd.pdb); cs=dict(zip(dd.pdb,dd.chain_seq)); m=len(ents); E=np.eye(m)
    for i,j in itertools.combinations(range(m),2):
        if keys[ents[i]]==keys[ents[j]] and cs[ents[i]]==cs[ents[j]]: E[i,j]=E[j,i]=1; continue
        E[i,j]=E[j,i]=ident(cs[ents[i]],cs[ents[j]])[0]
    out[fam]=dict(names=names,M=M.tolist(),Ms=Ms.tolist(),L=L.tolist(),ents=ents,E=E.tolist())
json.dump(dict(out=out,meta=meta,keys=keys,seqlen={k:len(v) for k,v in seqs.items()}),open(S+'/ident.json','w'))
d.to_csv(S+'/entries2.tsv',sep='\t',index=False)
print('done')
