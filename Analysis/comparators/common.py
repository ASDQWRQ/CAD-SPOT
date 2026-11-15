"""Shared, candidate-blind preparation and auditable method exports."""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import random
import subprocess
import sys
import time
import traceback

ROOT = Path(__file__).resolve().parent
DEPS = Path(os.environ.get('COMPARATOR_DEPS', str(ROOT.parents[1] / 'work/comparator_deps')))
if DEPS.exists():
    sys.path.insert(0, str(DEPS))

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree

SAMPLES = [f'GSM{8443449+i}_PDAC-p{i+1}' for i in range(5)]


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, default=str, allow_nan=False), encoding='utf-8')


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def bool_labels(series: pd.Series) -> np.ndarray:
    normalized = series.astype(str).str.lower()
    if not normalized.isin(['true', 'false', '1', '0']).all():
        raise ValueError('Candidate labels must be explicit True/False or 1/0')
    return normalized.isin(['true', '1']).to_numpy()


def seed_all(seed: int) -> None:
    import torch
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.benchmark = False


def graph(coords: np.ndarray) -> sparse.csr_matrix:
    """Radius-limited first-ring graph, not unrestricted kNN across tissue gaps."""
    tree = cKDTree(coords)
    pitch = float(np.median(tree.query(coords, k=2)[0][:, 1]))
    pairs = tree.query_pairs(1.25 * pitch, output_type='ndarray')
    if not len(pairs):
        raise ValueError('No spatial edges; inspect coordinate units')
    row = np.concatenate([pairs[:, 0], pairs[:, 1]])
    col = np.concatenate([pairs[:, 1], pairs[:, 0]])
    return sparse.csr_matrix((np.ones(len(row), dtype=np.float32), (row, col)), shape=(len(coords), len(coords)))


def log_normalize(counts: sparse.spmatrix, library_size: np.ndarray) -> sparse.csr_matrix:
    x = sparse.diags(1e4 / library_size) @ counts.astype(np.float32)
    x = x.tocsr()
    x.data = np.log1p(x.data)
    return x


class RunLog:
    def __init__(self, path: Path, args: object, repository: str):
        self.path = path
        self.start = time.monotonic()
        import torch
        packages = {}
        for name in ['numpy', 'scipy', 'pandas', 'scikit-learn', 'torch', 'anndata', 'scanpy', 'squidpy', 'pyro-ppl']:
            try:
                packages[name] = importlib.metadata.version(name)
            except importlib.metadata.PackageNotFoundError:
                packages[name] = 'not-installed'
        source = ROOT / 'vendor' / repository
        try:
            commit = subprocess.run(['git', '-C', str(source), 'rev-parse', 'HEAD'], capture_output=True, text=True)
            diff = subprocess.run(['git', '-C', str(source), 'diff'], capture_output=True, text=True)
            commit_text = commit.stdout.strip() if commit.returncode == 0 else ''
            diff_text = diff.stdout if diff.returncode == 0 else ''
        except OSError:
            commit_text, diff_text = '', ''
        source_provenance = 'live_git'
        if not commit_text:
            frozen = json.loads((ROOT / 'source_versions.json').read_text(encoding='utf-8'))[repository]
            commit_text, diff_text = frozen['commit'], frozen['local_patch']
            source_provenance = 'vendored_source_versions_manifest'
        self.data = dict(status='running', arguments=vars(args), python=sys.version, platform=platform.platform(),
                         packages=packages, repository=repository, commit=commit_text, source_diff=diff_text, source_provenance=source_provenance,
                         cuda_available=torch.cuda.is_available(), utc_started=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()))
        if torch.cuda.is_available():
            self.data['gpu'] = torch.cuda.get_device_name(0)
            self.data['gpu_total_mib'] = torch.cuda.get_device_properties(0).total_memory / 2**20
        write_json(path, self.data)

    def finish(self, status: str, error: str | None = None) -> None:
        self.data.update(status=status, wall_seconds=time.monotonic()-self.start)
        if error:
            self.data['error'] = error
        write_json(self.path, self.data)


def export_scores(out: Path, adata: object, scores: np.ndarray, method: str, seed: int,
                  probability: bool = False, smoke: bool = False) -> np.ndarray:
    out.mkdir(parents=True, exist_ok=True)
    z = np.asarray(scores, dtype=float)
    if z.shape[0] != adata.n_obs or z.ndim != 2 or not np.isfinite(z).all():
        raise ValueError('Invalid score dimensions or nonfinite values')
    if np.min(z) < -1e-6:
        raise ValueError('Negative final embeddings cannot be called nonnegative compositions')
    native = z.copy()
    z = np.maximum(z, 0)
    sums = z.sum(axis=1)
    valid = sums > 0
    weights = z[valid] / sums[valid, None]
    names = [f'factor_{i+1:02d}' for i in range(z.shape[1])]
    base = adata.obs[['sample_id', 'barcode']].reset_index(drop=True).copy()
    base['method'] = method
    base['seed'] = seed
    base['smoke_test'] = smoke
    base['score_kind'] = 'posterior_topic_proportion' if probability else 'row_normalized_nonnegative_embedding'
    base['valid_factor_vector'] = valid
    base['assignment_status'] = np.where(valid, 'assigned', 'unassigned_zero_vector')
    base['dominant_domain'] = None
    base.loc[valid, 'dominant_domain'] = np.asarray(names)[weights.argmax(axis=1)]
    base.to_csv(out / 'all_spot_assignments.tsv.gz', sep='\t', index=False)
    result = pd.concat([base.loc[valid].reset_index(drop=True), pd.DataFrame(weights, columns=names)], axis=1)
    result.to_csv(out / 'spot_factors.tsv.gz', sep='\t', index=False)
    pd.concat([base[['sample_id', 'barcode']], pd.DataFrame(native, columns=names)], axis=1).to_csv(
        out / 'raw_factor_scores.tsv.gz', sep='\t', index=False)
    write_json(out / 'export_metadata.json', dict(method=method, seed=seed, smoke_test=smoke,
               n_spots=int(valid.sum()), n_input_spots=adata.n_obs, n_factors=z.shape[1], factor_labels_aligned_across_methods=False,
               posterior_proportions=probability, zero_sum_rows=int((~valid).sum()), pseudocount_added=False))
    if not valid.any():
        raise ValueError('No valid final factor vectors; raw vectors retained but no comparator evidence')
    return valid


def execute_logged(log: RunLog, function: object) -> None:
    try:
        function()
    except Exception:
        log.finish('failed', traceback.format_exc())
        raise
    else:
        log.finish('completed')
