"""Check input hashes, benchmark labels and regenerated numerical results."""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'analysis'))
from label_policy import apply_labels
from analysis import score


def compare(actual, expected, path):
    if isinstance(expected, dict):
        if actual.keys() != expected.keys():
            raise AssertionError(f'{path}: keys differ')
        for key in expected:
            compare(actual[key], expected[key], f'{path}.{key}')
    elif isinstance(expected, list):
        if len(actual) != len(expected):
            raise AssertionError(f'{path}: lengths differ')
        for i, (a, e) in enumerate(zip(actual, expected)):
            compare(a, e, f'{path}[{i}]')
    elif isinstance(expected, (int, float)) and not isinstance(expected, bool):
        np.testing.assert_allclose(actual, expected, rtol=1e-9, atol=1e-10, err_msg=path)
    elif actual != expected:
        raise AssertionError(f'{path}: {actual!r} != {expected!r}')


def main():
    manifest = json.loads((ROOT / 'data/checksums.json').read_text())
    for name, digest in manifest.items():
        actual = hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
        if actual != digest:
            raise AssertionError(f'{name}: SHA256 mismatch')
    stats = json.loads((ROOT / 'reference/statistics.json').read_text())
    for dataset, filename, n in [('S669','s669',669),('S461','s461',461),('S96','s96',96),('PTmul-NR','ptmul_nr',82)]:
        frame = pd.read_csv(ROOT / f'data/{filename}.csv')
        assert len(frame) == n
        labels = apply_labels(frame, dataset)
        np.testing.assert_array_equal(frame.y, labels.y)
        for key, values in [('J',frame.J),('M',frame.M),('E',(frame.J+frame.M)/2)]:
            compare(score(frame.y,values), stats[dataset][key], f'{dataset}.{key}')
    for path in sorted((ROOT / 'reference').glob('*.json')):
        compare(json.loads((ROOT / 'generated' / path.name).read_text()),
                json.loads(path.read_text()), path.name)
    print('Verified input hashes, all benchmark labels, primary metrics and reference statistics.')


if __name__ == '__main__':
    main()
