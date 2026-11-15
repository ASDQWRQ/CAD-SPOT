"""Compatibility entry point; the current implementation is in Software/."""
from pathlib import Path
import runpy
import sys

if __name__ == '__main__':
    folder=Path(__file__).resolve().parent/'Software'
    sys.path.insert(0,str(folder))
    runpy.run_path(str(folder/'cadspot.py'),run_name='__main__')
