import json, urllib.request
exec(open('rcsb_q.py').read().split("subs=")[0])
comps={}
for g in ['glucopyranose','glucopyranoside','glucopyranosyl','glucitol','maltose','maltosyl']:
    for ph in ['2-O','3-O','6-O','4-O','O-methyl','O-acetyl','2-hydroxypropyl','carboxymethyl']:
        body={"query":{"type":"group","logical_operator":"and","nodes":[T("chem_comp.name","contains_words",g),T("chem_comp.name","contains_phrase",ph)]},
              "return_type":"mol_definition","request_options":{"return_all_hits":True}}
        try: j=q(body)
        except Exception as e: print('ERR',g,ph,e); continue
        for x in j.get('result_set',[]): comps.setdefault(x['identifier'],set()).add(ph)
old=json.load(open('rcsb_hits.json'))['comps']
allc=sorted(set(comps)|set(old))
ids=[c for c in allc if not c.startswith('PRD')]
print('comps',len(ids))
pf=["PF00128","PF01055","PF01373","PF00723","PF02806","PF16863","PF21365","PF10566"]
fam_node={"type":"group","logical_operator":"or","nodes":[T("rcsb_polymer_entity_annotation.annotation_id","exact_match",p,'text') for p in pf]}
res={}
for c in ids:
    chem_node={"type":"group","logical_operator":"or","nodes":[T("rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id","exact_match",c,'text'),
                T("rcsb_branched_entity_container_identifiers.chem_comp_monomers","exact_match",c,'text')]}
    body={"query":{"type":"group","logical_operator":"and","nodes":[fam_node,chem_node]},"return_type":"entry","request_options":{"return_all_hits":True}}
    try: j=q(body)
    except Exception as e: print('ERR',c,e); continue
    e=[x['identifier'] for x in j.get('result_set',[])]
    if e: res[c]=e
for c,e in sorted(res.items()): print(c, e)
json.dump(res,open('rcsb_hits2.json','w'),indent=1)
