"""Refresh manuscript figures from cached data, preserving all numerical results."""
import json
from pathlib import Path

import pandas as pd

import analysis
import classification
import ensemble_ablation
import s669_multisite_paper
import bias_comparison

HERE=Path(__file__).resolve().parent.parent


def read_json(filename):
    return json.loads((HERE/'generated'/filename).read_text())


def main():
    analysis.figures(pd.read_csv(HERE/'data/s669.csv'),pd.read_csv(HERE/'data/ptmul_nr.csv'),
        pd.read_csv(HERE/'data/transitivity.csv'),read_json('statistics.json'))
    stats=read_json('classification_statistics.json')['datasets'];bundles={}
    columns={'JanusDDG':'janus_score','MSAesm original':'msa_original_score',
        'MSAesm anti':'msa_anti_score','EvoPhi-DDG':'evophi_score','Ensemble':'ensemble_score'}
    for name,filename in [('S669','s669'),('S461','s461'),('S96','s96'),('PTmul-NR','ptmul_nr')]:
        data=pd.read_csv(HERE/f'data/classification_{filename}.csv')
        models={model:data[column].to_numpy() for model,column in columns.items() if model in stats[name]['metrics']}
        bundles[name]=(stats[name],data,data.class_stabilizing.to_numpy(dtype=bool),models)
    classification.make_figure(bundles,stats)
    ensemble_ablation.figure(read_json('ensemble_ablation_statistics.json')['datasets'],
        pd.read_csv(HERE/'data/ensemble_ablation_path_summary.csv'))
    s669_multisite_paper.figures()
    bias_comparison.main()
    classification.make_figure({'S96':bundles['S96']},stats,filename='classification_s96')
    print('Refreshed all figures (PDF and PNG) from cached EvoPhi-DDG data.')


if __name__=='__main__':main()
