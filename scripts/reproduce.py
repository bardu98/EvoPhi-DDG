"""Recompute the paper artifacts from frozen predictions and pathway data."""
from pathlib import Path
import argparse
import contextlib
import importlib
import os
import sys
import time

ROOT = Path(__file__).resolve().parent.parent
os.environ.setdefault('MPLCONFIGDIR', str(ROOT / '.cache' / 'matplotlib'))
sys.path.insert(0, str(ROOT / 'analysis'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--figures-only', action='store_true', help='Redraw using existing statistics')
    args = parser.parse_args()
    stages = ['refresh_figures'] if args.figures_only else [
        'evophi_physics', 'analysis', 'classification', 'ensemble_ablation',
        's669_multisite_paper', 'benchmark_summary', 'refresh_figures']
    logs = ROOT / '.cache'
    logs.mkdir(exist_ok=True)
    for name in stages:
        start = time.monotonic()
        print(f'{name}: running', flush=True)
        with (logs / f'{name}.log').open('w') as log, contextlib.redirect_stdout(log):
            importlib.import_module(name).main()
        print(f'{name}: complete ({time.monotonic()-start:.1f} s)', flush=True)
    from export_results import main as export
    export()


if __name__ == '__main__':
    main()
