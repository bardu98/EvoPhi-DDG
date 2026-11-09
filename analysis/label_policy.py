"""Apply the user-selected JanusDDG experimental labels to all manuscript scores.

The frozen table preserves all benchmark records, including repeated S669
measurements. Keys are internal variant codes, not ambiguous protein IDs.
"""
from pathlib import Path

import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent.parent
POLICY='Experimental labels distributed with JanusDDG, applied identically to every predictor.'


def apply_labels(frame,dataset):
    reference=pd.read_csv(HERE/'data/experimental_labels.csv')
    reference=reference.loc[reference.dataset.eq(dataset)]
    if reference.empty or not reference.code.is_unique or not frame.code.is_unique:
        raise ValueError(f'{dataset}: missing or nonunique variant codes in label mapping')
    if set(frame.code)!=set(reference.code):
        raise ValueError(f'{dataset}: benchmark membership differs from the frozen label reference')
    result=frame.copy()
    result['y']=reference.set_index('code').loc[result.code,'y'].to_numpy()
    if not np.isfinite(result.y).all():
        raise ValueError(f'{dataset}: nonfinite experimental labels')
    result.attrs['label_policy']=POLICY
    return result
