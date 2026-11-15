param(
    [string]$Python = 'python',
    [string]$Target = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot '../../work/comparator_deps'))
)
$ErrorActionPreference = 'Stop'
New-Item -ItemType Directory -Force -Path $Target | Out-Null
# Core numeric packages and torch deliberately stay in the shared runtime.
$packages = @(
    'anndata', 'h5py', 'matplotlib', 'scanpy', 'squidpy', 'seaborn',
    'pyro-ppl', 'pyro-api', 'torchinfo', 'gseapy', 'discotoolkit',
    'awkward', 'awkward-cpp', 'ortools', 'kneed',
    'umap-learn', 'pynndescent', 'numba', 'llvmlite', 'natsort',
    'legacy-api-wrap', 'session-info2', 'zarr', 'numcodecs', 'donfig',
    'google-crc32c', 'array-api-compat', 'pillow', 'contourpy', 'cycler',
    'fonttools', 'kiwisolver', 'pyparsing', 'packaging', 'tqdm',
    'scikit-image', 'opencv-python', 'imagecodecs', 'scikit-network',
    'igraph', 'leidenalg', 'adjustText', 'setuptools<81',
    'requests', 'urllib3', 'charset-normalizer', 'certifi', 'idna',
    'networkx', 'joblib', 'threadpoolctl', 'texttable', 'patsy', 'statsmodels',
    'absl-py', 'protobuf', 'immutabledict', 'typing-extensions', 'fsspec',
    'scikit-misc', 'pydantic', 'pydantic-core', 'pydantic-settings',
    'python-dotenv', 'annotated-types', 'typing-inspection', 'PyYAML', 'click'
)
& $Python -m pip install --disable-pip-version-check --target $Target --no-deps @packages
if ($LASTEXITCODE -ne 0) { throw 'Dependency install failed; shared environment unchanged.' }
Write-Output "Dependencies staged in $Target. Use PYTHONPATH; import checks remain required."
