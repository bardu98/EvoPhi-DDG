"""Validate frozen EvoPhi-DDG scores against the independent head-ablation export."""
from pathlib import Path
import numpy as np
import pandas as pd
HERE=Path(__file__).resolve().parent.parent
ROOT=HERE

def apply_model(frame,dataset):
    result=frame.copy()
    reference=pd.read_csv(HERE/'data/model_controls'/f'{dataset}_predictions.csv',index_col=0)
    if not result.code.is_unique or set(result.code)!=set(reference.index):
        raise ValueError(f'{dataset}: model-control membership mismatch')
    reference=reference.loc[result.code]
    np.testing.assert_allclose(result.M_anti,reference.original_AS,rtol=0,atol=5e-5)
    np.testing.assert_allclose(result.M,reference.zero_hidden_bias_AS,rtol=0,atol=1e-12)
    result['M']=reference.zero_hidden_bias_AS.to_numpy()
    result.attrs['primary_msa_model']='EvoPhi-DDG'
    return result

def physical_edges():
    return np.load(HERE/'data/physics/evophi_edges.npy',allow_pickle=False)
