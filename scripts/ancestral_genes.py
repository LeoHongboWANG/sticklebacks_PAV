#!/usr/bin/env python3
"""Reconstructed ascertainment from binary PAV data; not the historical script."""
import argparse
import csv
import gzip
import json
from pathlib import Path
import numpy as np


def opened(path):
    return gzip.open(path, 'rt', newline='') if str(path).endswith('.gz') else open(path, newline='')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--outgroup-pav', required=True, help='TSV, ID then binary sample columns; .gz accepted')
    p.add_argument('--outdir', required=True)
    p.add_argument('--pav', help='Optional full binary PAV TSV')
    p.add_argument('--samples', help='Optional TSV with SampleID, Population (e.g. lost.txt)')
    a = p.parse_args()
    if bool(a.pav) != bool(a.samples):
        p.error('--pav and --samples must be supplied together')
    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    ancestor = set(); seen = set()
    with opened(a.outgroup_pav) as f, (out/'ancestral_presence.csv').open('w', newline='') as g:
        r = csv.reader(f, delimiter='\t'); h = next(r); n = len(h)-1
        if n < 1 or len(set(h[1:])) != n: raise ValueError('Invalid or duplicate sample columns')
        w = csv.writer(g); w.writerow(['GeneID','OutgroupPresent','OutgroupN','PresenceFraction','AncestralGT90'])
        for row in r:
            if len(row) != n+1 or not set(row[1:]) <= {'0','1'}: raise ValueError('Invalid binary PAV row')
            gene = row[0]
            if gene in seen: raise ValueError('Duplicate GeneID: '+gene)
            seen.add(gene); present = row[1:].count('1'); keep = 10*present > 9*n
            if keep: ancestor.add(gene)
            w.writerow([gene,present,n,present/n,int(keep)])
    (out/'ancestral_gene_ids.txt').write_text(''.join(x+'\n' for x in sorted(ancestor)))
    report = {'outgroup_samples':n,'input_genes':len(seen),'ancestral_genes':len(ancestor),'criterion':'pooled presence > 0.90; strict inequality'}
    if a.pav:
        with opened(a.samples) as f: samples = list(csv.DictReader(f, delimiter='\t'))
        ids = [x['SampleID'] for x in samples]
        if len(set(ids)) != len(ids): raise ValueError('Duplicate SampleID')
        totals = np.zeros(len(ids), dtype=np.int64); found = set()
        with opened(a.pav) as f:
            r = csv.reader(f, delimiter='\t'); h = next(r)
            if len(set(h)) != len(h): raise ValueError('Duplicate PAV columns')
            idx = [h.index(s) for s in ids]
            for row in r:
                if row[0] not in ancestor: continue
                if row[0] in found: raise ValueError('Duplicate ancestral GeneID')
                found.add(row[0]); vals = [row[i] for i in idx]
                if not set(vals) <= {'0','1'}: raise ValueError('Nonbinary PAV')
                totals += np.array([v == '0' for v in vals])
        if found != ancestor: raise ValueError('Ancestral genes missing from full matrix')
        with (out/'sample_ancestral_loss.csv').open('w', newline='') as f:
            w=csv.writer(f); w.writerow(['SampleID','Population','RecomputedLostGenes','SavedLostGenes','Difference'])
            differences=[]
            for x,t in zip(samples,totals):
                saved=x.get('Lost_gene_count',''); diff=int(t)-int(saved) if saved else ''
                w.writerow([x['SampleID'],x['Population'],int(t),saved,diff]); differences.append(diff)
        report.update(samples=len(ids), sample_loss_disagreements=sum(d != 0 for d in differences if d != ''))
    (out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report))


if __name__ == '__main__': main()
