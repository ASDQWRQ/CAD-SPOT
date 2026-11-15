"""Download only public analysis inputs, with explicit repository identifiers."""
import argparse
import json
from pathlib import Path
import tarfile
import urllib.request

ROOT=Path(__file__).resolve().parents[1]
EXTERNAL={
 'GSM6177599':'NYU_BRCA0_Vis','GSM6177601':'NYU_BRCA1_Vis','GSM6177603':'NYU_BRCA2_Vis',
 'GSM6177607':'NYU_GIST1_Vis','GSM6177609':'NYU_GIST2_Vis','GSM6177612':'NYU_LIHC1_Vis',
 'GSM6177614':'NYU_OVCA1_Vis','GSM6177617':'NYU_OVCA3_Vis','GSM6177618':'NYU_PDAC1_Vis',
 'GSM6177623':'NYU_UCEC3_Vis'}


def download(url,path):
    if path.exists() and path.stat().st_size>0:
        print('Present',path,flush=True)
        return
    path.parent.mkdir(parents=True,exist_ok=True)
    partial=path.with_name(path.name+'.partial')
    print('Download',url,flush=True)
    with urllib.request.urlopen(url,timeout=180) as source,partial.open('wb') as out:
        while True:
            chunk=source.read(1024*1024)
            if not chunk:break
            out.write(chunk)
    partial.replace(path)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--datasets',nargs='+',choices=['pdac','external','annotation'],default=['pdac','annotation']);args=ap.parse_args()
    if 'pdac' in args.datasets:
        path=ROOT/'Data/GSE274103_RAW.tar'
        download('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/GSE274103/suppl/GSE274103_RAW.tar',path)
        target=ROOT/'Data/GSE274103/raw';target.mkdir(parents=True,exist_ok=True)
        with tarfile.open(path) as f:f.extractall(target,filter='data')
        download('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/GSE274103/soft/GSE274103_family.soft.gz',ROOT/'Data/GSE274103_family.soft.gz')
    if 'annotation' in args.datasets:
        download('https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.chr.gtf.gz',ROOT/'Data/Homo_sapiens.GRCh38.113.chr.gtf.gz')
    if 'external' in args.datasets:
        for gsm,sample in EXTERNAL.items():
            filename=f'{gsm}_{sample}_processed_filtered_feature_bc_matrix.h5'
            download(f'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6177nnn/{gsm}/suppl/{filename}',ROOT/'Data/GSE203612/raw'/filename)
