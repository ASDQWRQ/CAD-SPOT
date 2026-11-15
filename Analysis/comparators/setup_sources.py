"""Restore the recorded official source commits and minimal compatibility diffs."""
import argparse
import json
from pathlib import Path
import subprocess

ROOT=Path(__file__).resolve().parent


def run(git,args,folder=None,patch=None,check=True):
    command=[git]
    if folder:command+=['-C',str(folder)]
    result=subprocess.run(command+args,input=patch.encode('utf-8') if patch else None,capture_output=True)
    if check and result.returncode:
        raise RuntimeError(result.stderr.decode('utf-8',errors='replace'))
    return result


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--git',default='git');ap.add_argument('--vendor-dir',type=Path,default=ROOT/'vendor');args=ap.parse_args()
    sources=json.loads((ROOT/'source_versions.json').read_text(encoding='utf-8'))
    args.vendor_dir.mkdir(parents=True,exist_ok=True)
    for name,record in sources.items():
        folder=args.vendor_dir/name
        if not folder.exists():
            folder.mkdir()
            run(args.git,['init'],folder)
            run(args.git,['config','core.autocrlf','false'],folder)
            run(args.git,['remote','add','origin',record['url']],folder)
            run(args.git,['fetch','--depth','1','origin',record['commit']],folder)
            run(args.git,['checkout','--detach','FETCH_HEAD'],folder)
        current=run(args.git,['rev-parse','HEAD'],folder).stdout.decode().strip()
        if current!=record['commit']:
            raise ValueError(f'{folder}: existing checkout has another commit; use a new vendor directory')
        patch=record.get('local_patch','')
        if patch:
            already=run(args.git,['apply','--reverse','--check','-'],folder,patch,check=False)
            if already.returncode:
                run(args.git,['apply','--check','-'],folder,patch)
                run(args.git,['apply','-'],folder,patch)
        asset=record.get('runtime_asset')
        if asset and not (folder/asset).exists():raise FileNotFoundError(folder/asset)
        print(name,current,'source and recorded patch verified',flush=True)


if __name__=='__main__':main()
