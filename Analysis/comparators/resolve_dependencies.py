"""Complete staged dependencies without upgrading or copying shared numeric/CUDA packages."""
from __future__ import annotations
import argparse
import importlib.metadata as metadata
from pathlib import Path
import subprocess
import sys
import common
from packaging.requirements import Requirement
from packaging.utils import canonicalize_name


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--apply', action='store_true')
    args = parser.parse_args()
    protected = {'numpy', 'pandas', 'scipy', 'scikit-learn', 'torch', 'sympy', 'mpmath'}
    missing, conflicts = set(), []
    for dist in metadata.distributions(path=[str(common.DEPS)]):
        for requirement in dist.requires or []:
            r = Requirement(requirement)
            if r.marker and not r.marker.evaluate({'extra': ''}):
                continue
            name = canonicalize_name(r.name)
            try:
                version = metadata.version(r.name)
            except metadata.PackageNotFoundError:
                if name not in protected:
                    missing.add(str(r).split(';')[0].strip())
            else:
                if r.specifier and version not in r.specifier:
                    conflicts.append(dict(parent=dist.metadata['Name'], requirement=str(r), installed=version, protected=name in protected))
    common.write_json(common.ROOT / 'dependency_resolution.json', dict(missing=sorted(missing), conflicts=conflicts))
    print('Missing:', sorted(missing), flush=True)
    print('Conflicts:', conflicts, flush=True)
    if args.apply and missing:
        subprocess.run([sys.executable, '-m', 'pip', 'install', '--disable-pip-version-check', '--no-deps',
                        '--target', str(common.DEPS), *sorted(missing)], check=True)


if __name__ == '__main__':
    main()
