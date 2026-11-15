"""Rebuild the released Ensembl-113 exon-union table from its public GTF."""
import argparse
from collections import defaultdict
import gzip
from pathlib import Path
import re
import urllib.request

import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
URL='https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.chr.gtf.gz'


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--gtf',type=Path);args=ap.parse_args()
    path=args.gtf or ROOT/'Data/Homo_sapiens.GRCh38.113.chr.gtf.gz'
    if not path.exists():
        path.parent.mkdir(parents=True,exist_ok=True)
        urllib.request.urlretrieve(URL,path)
    intervals=defaultdict(list);names={}
    gidpat=re.compile(r'gene_id "([^"]+)"');namepat=re.compile(r'gene_name "([^"]+)"')
    with gzip.open(path,'rt') as handle:
        for line in handle:
            if line.startswith('#'):continue
            fields=line.rstrip('\n').split('\t')
            if len(fields)!=9 or fields[2]!='exon':continue
            gid,name=gidpat.search(fields[8]),namepat.search(fields[8])
            if not gid or not name:continue
            key=(gid.group(1),fields[0]);names[key]=name.group(1)
            intervals[key].append((int(fields[3]),int(fields[4])))
    rows=[]
    for key,values in intervals.items():
        length=0;lo=hi=-1;n_exons=0
        for start,end in sorted(values):
            if start>hi+1:
                if hi>=0:length+=hi-lo+1
                lo,hi=start,end;n_exons+=1
            else:hi=max(hi,end)
        length+=hi-lo+1
        rows.append(dict(gene_id=key[0],chromosome=key[1],gene=names[key],
                         gene_length=length,n_exon_union_intervals=n_exons))
    table=pd.DataFrame(rows).sort_values(['gene','gene_length','gene_id'],ascending=[True,False,True]).drop_duplicates('gene')
    destination=ROOT/'Provenance/ensembl113_exon_union_lengths.tsv'
    destination.parent.mkdir(parents=True,exist_ok=True)
    table.to_csv(destination,sep='\t',index=False)
    print(destination,flush=True)


if __name__=='__main__':main()
