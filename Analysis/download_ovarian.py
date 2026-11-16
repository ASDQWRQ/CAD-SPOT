"""Retrieve GEO files with verified byte ranges when long streams stall."""
import argparse
import http.client
import json
from pathlib import Path
import time
import urllib.request
import tarfile

ROOT = Path(__file__).resolve().parents[1]
URL = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE211nnn/GSE211956/suppl/GSE211956_RAW.tar'


def download_ranges(url, destination, chunk=1024*1024):
    partial = destination.with_suffix(destination.suffix + '.ranges.partial')
    start = partial.stat().st_size if partial.exists() else 0
    total = None
    with partial.open('ab') as output:
        while total is None or start < total:
            end = start + chunk - 1 if total is None else min(start + chunk - 1, total - 1)
            for attempt in range(4):
                try:
                    req = urllib.request.Request(url + f'?cadspot_range={start}_{end}_{attempt}',
                        headers={'Range': f'bytes={start}-{end}', 'User-Agent': 'CAD-SPOT-data-retrieval/2.2'})
                    with urllib.request.urlopen(req, timeout=40) as response:
                        crange = response.headers.get('Content-Range', '')
                        prefix = f'bytes {start}-'
                        if response.status != 206 or not crange.startswith(prefix):
                            raise ValueError(f'Expected exact range response, got {response.status}/{crange}')
                        bounds, denominator = crange.split('/')
                        stop = int(bounds.split('-')[-1])
                        total = int(denominator)
                        try:
                            payload = response.read()
                        except http.client.IncompleteRead as truncated:
                            payload = truncated.partial
                        if not payload or len(payload) > stop-start+1 or stop != min(end, total-1):
                            raise ValueError('Empty or incorrectly ranged GEO response')
                        # A closed stream may still contain a valid range prefix.
                        # Append only received bytes and request its missing tail next.
                        if len(payload) < stop-start+1:
                            print(f'Resuming a short range after {len(payload)} verified-prefix bytes', flush=True)
                    break
                except Exception:
                    if attempt == 3:
                        raise
                    time.sleep(2*(attempt+1))
            output.write(payload)
            output.flush()
            start += len(payload)
            print(f'GEO download {start}/{total} bytes ({100*start/total:.1f}%)', flush=True)
    if partial.stat().st_size != total:
        raise ValueError('Incomplete archive')
    with tarfile.open(partial) as archive:
        members = archive.getmembers()
        if not members:
            raise ValueError('Empty GEO archive')
    partial.replace(destination)
    return total


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--folder', type=Path, default=ROOT/'Data/GSE211956')
    args = p.parse_args()
    args.folder.mkdir(parents=True, exist_ok=True)
    path = args.folder/'GSE211956_RAW.tar'
    if path.exists():
        try:
            with tarfile.open(path) as archive:
                archive.getmembers()
        except (tarfile.TarError, EOFError):
            path.rename(path.with_suffix('.tar.failed-stream.partial'))
    if not path.exists():
        download_ranges(URL, path)
    raw = args.folder/'raw'
    raw.mkdir(exist_ok=True)
    with tarfile.open(path) as archive:
        archive.extractall(raw, filter='data')
        names = [m.name for m in archive.getmembers()]
    (args.folder/'archive_members.json').write_text(json.dumps(names, indent=2)+'\n')
    print('GEO archive extracted:', len(names), 'members', flush=True)


if __name__ == '__main__':
    main()
