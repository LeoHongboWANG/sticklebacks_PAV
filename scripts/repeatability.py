#!/usr/bin/env python3
"""Reconstructed ecotype-label permutation; historical random seed is unknown."""
import argparse
import csv
import gzip
import json
import math
from pathlib import Path
import numpy as np


def opened(path):
    return gzip.open(path,'rt',newline='') if str(path).endswith('.gz') else open(path,newline='')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--counts',required=True,help='Long CSV with GeneID, Population, EcotypeSimple, PopulationSize, AbsentCount')
    p.add_argument('--outdir',required=True)
    p.add_argument('--exclude',action='append',default=None,help='Default SWE-LUN; may be repeated')
    p.add_argument('--keep-genes',help='Optional final matrix GeneID list, one ID per line; NOT supplied in this package')
    p.add_argument('--permutations',type=int,default=2000)
    p.add_argument('--seed',type=int,default=20261026)
    p.add_argument('--min-fw',type=int,help='Default ceil(number of freshwater populations / 2)')
    a=p.parse_args()
    if a.permutations<0: p.error('permutations must be nonnegative')
    exclude=set(a.exclude or ['SWE-LUN']); keep=set(Path(a.keep_genes).read_text().split()) if a.keep_genes else None
    if keep is not None and not keep: raise ValueError('Empty final GeneID list')
    rows={}; pops={}; sizes={}; raw_genes=set()
    with opened(a.counts) as f:
        for row in csv.DictReader(f):
            gene=row['GeneID'];raw_genes.add(gene)
            if row['Population'] in exclude or (keep is not None and gene not in keep):continue
            pop=row['Population'];eco=row['EcotypeSimple'];n=int(row['PopulationSize']);absent=int(row['AbsentCount'])
            if eco not in {'Marine','Freshwater'} or n<1 or not 0<=absent<=n:raise ValueError('Invalid population row')
            if pop in pops and (pops[pop]!=eco or sizes[pop]!=n):raise ValueError('Inconsistent population metadata')
            pops[pop]=eco;sizes[pop]=n
            if pop in rows.setdefault(gene,{}):raise ValueError('Duplicate gene/population')
            rows[gene][pop]=int(2*absent>=n)
    if keep is not None and not keep <= raw_genes:raise ValueError('Requested final genes missing from counts')
    if not rows:raise ValueError('No genes retained')
    genes=sorted(rows); populations=sorted(pops)
    if any(set(rows[g])!=set(populations) for g in genes):raise ValueError('Incomplete gene-by-population matrix')
    x=np.array([[rows[g][s] for s in populations] for g in genes],dtype=np.uint8)
    fw=np.array([pops[s]=='Freshwater' for s in populations]);nfw=int(fw.sum());nmar=len(fw)-nfw
    if not nfw or not nmar:raise ValueError('Both ecotypes are required')
    k=a.min_fw if a.min_fw is not None else math.ceil(nfw/2)
    if not 1<=k<=nfw:raise ValueError('Invalid freshwater threshold')
    fw_count=x[:,fw].sum(1);mar_count=x[:,~fw].sum(1);selected=(fw_count>=k)&(mar_count<=1)
    observed=int(selected.sum());out=Path(a.outdir);out.mkdir(parents=True,exist_ok=True)
    with (out/'gene_screen.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['GeneID','FWPopGE50','MarinePopGE50','Selected'])
        w.writerows((g,int(f),int(m),int(s)) for g,f,m,s in zip(genes,fw_count,mar_count,selected))
    (out/'selected_gene_ids.txt').write_text(''.join(g+'\n' for g,s in zip(genes,selected) if s))
    # Only genes with k..nfw+1 total qualifying populations can meet the screen under any labeling.
    totals=x.sum(1);potential=x[(totals>=k)&(totals<=nfw+1)]
    if len(potential):patterns,weights=np.unique(potential,axis=0,return_counts=True)
    else:patterns=np.empty((0,len(populations)),dtype=np.uint8);weights=np.empty(0,dtype=int)
    ptot=patterns.sum(1);rng=np.random.default_rng(a.seed);null=[]
    for start in range(0,a.permutations,100):
        masks=np.zeros((len(populations),min(100,a.permutations-start)),dtype=np.uint8)
        for j in range(masks.shape[1]):masks[rng.choice(len(populations),nfw,replace=False),j]=1
        f=patterns@masks;m=ptot[:,None]-f
        null.extend((((f>=k)&(m<=1))*weights[:,None]).sum(0).astype(int).tolist())
    with (out/'new_permutation_distribution.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['Permutation','Count']);w.writerows(enumerate(null,1))
    report={'raw_genes':len(raw_genes),'screened_genes':len(genes),'populations':len(populations),'freshwater_populations':nfw,'marine_populations':nmar,'excluded':sorted(exclude),'freshwater_minimum':k,'marine_maximum':1,'observed':observed,'permutations':len(null),'seed':a.seed,'historical_seed_recovered':False,'final_gene_filter_applied':keep is not None,'matrix_version_status':'final GeneID list supplied' if keep is not None else 'local unfiltered gene universe; reconcile with final manuscript matrix'}
    if null:report.update(null_mean=float(np.mean(null)),p_corrected=(1+sum(v>=observed for v in null))/(len(null)+1))
    (out/'summary.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))


if __name__=='__main__':main()
