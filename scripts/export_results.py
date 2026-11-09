"""Export unrounded, tidy numerical results and readable benchmark tables."""
from pathlib import Path
import json
import shutil
import sys
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/'analysis'))
from analysis import score


def main():
    output = ROOT/'results'
    output.mkdir(exist_ok=True)
    records=[]
    for dataset, filename in [('S669','s669'),('S461','s461'),('S96','s96'),('PTmul-NR','ptmul_nr')]:
        f=pd.read_csv(ROOT/f'data/{filename}.csv')
        models=[('JanusDDG',f.J),('EvoPhi-DDG',f.M),('Ensemble',(f.J+f.M)/2),
                ('MSAesm_ddG antisymmetrized',f.M_anti)]
        if 'M_original' in f:
            models.append(('MSAesm_ddG original',f.M_original))
        if dataset=='PTmul-NR':
            models.extend([('EvoPhi-DDG sequential',f.M_sequence),('EvoPhi-DDG all orders',f.M_orders)])
        for name,pred in models:
            records.append(dict(dataset=dataset,model=name,**score(f.y,pred)))
    regression=pd.DataFrame(records)
    regression.to_csv(output/'regression.csv',index=False)
    for source,target in [('data/transitivity.csv','same_site_transitivity.csv'),
                          ('data/s669_multisite_summary.csv','multisite_transitivity.csv'),
                          ('generated/classification_metrics.csv','classification.csv')]:
        shutil.copy2(ROOT/source,output/target)
    stats=json.loads((ROOT/'generated/statistics.json').read_text())
    bootstrap=[]
    for dataset,values in stats.items():
        for comparator,key in [('JanusDDG','J'),('EvoPhi-DDG','M')]:
            row=values['bootstrap']['ensemble_minus_'+key]
            bootstrap.append(dict(dataset=dataset,comparator=comparator,delta_pearson=row['delta'],
                                  ci95_low=row['ci95'][0],ci95_high=row['ci95'][1]))
    pd.DataFrame(bootstrap).to_csv(output/'bootstrap.csv',index=False)
    ablation=json.loads((ROOT/'generated/ensemble_ablation_statistics.json').read_text())
    rows=[]
    for dataset,values in ablation['datasets'].items():
        for name in ['E_original','E_whole_only','E_anti','E_phi']:
            if name in values:
                v=values[name]
                rows.append(dict(dataset=dataset,construction=name,
                    **{k:v[k] for k in ['n','pearson','spearman','mae','rmse','reversal_mae','reversal_rmse','reversal_max']}))
    pd.DataFrame(rows).to_csv(output/'ensemble_controls.csv',index=False)
    literature=pd.read_csv(ROOT/'data/literature_comparisons.csv')
    shutil.copy2(ROOT/'data/literature_comparisons.csv',output/'literature_comparisons.csv')
    latex=[r'\begin{tabular}{llrrr}\toprule',r'Method & Source & Pearson & MAE & RMSE\\\midrule']
    for row in literature[literature.dataset.eq('S669')].itertuples():
        method=row.method.replace('_',r'\_')
        source=r'\citep{'+row.bibtex_keys.replace(';',',')+'}'
        values=['---' if pd.isna(value) else f'{value:.2f}' for value in [row.pearson,row.mae,row.rmse]]
        latex.append(' & '.join([method,source]+values)+r' \\')
    latex.append(r'\midrule')
    f=pd.read_csv(ROOT/'data/s669.csv')
    models=[('JanusDDG',f.J),('MSAesm_ddG, current inputs',f.M_original),
            ('MSAesm_ddG, antisymmetrized',f.M_anti),('EvoPhi-DDG',f.M),
            ('Mean with original MSAesm_ddG',(f.J+f.M_original)/2),
            ('Fixed equal-weight ensemble',(f.J+f.M)/2)]
    from table_formatting import performance_rows
    values=[[score(f.y,p)[k] for k in ['pearson','mae','rmse']] for _,p in models]
    formatted=performance_rows([name.replace('_',r'\_') for name,_ in models],values,['.2f']*3,['max','min','min'])
    for row in formatted:
        label,cells=row.split(' & ',1)
        latex.append(label+' & This study & '+cells)
    latex.append(r'\bottomrule\end{tabular}')
    (ROOT/'generated/literature_rows.tex').write_text('\n'.join(latex)+'\n')
    lines=['# Benchmark results','','All methods use the frozen experimental labels distributed with JanusDDG. '
           'The ensemble is the equal-weight mean of JanusDDG and EvoPhi-DDG. '
           'Errors are in kcal/mol. Values below use two decimals; CSV and JSON files retain full precision.','',
           '| Benchmark | Variants | Model | Pearson | Spearman | MAE | RMSE |',
           '| --- | ---: | --- | ---: | ---: | ---: | ---: |']
    for r in regression[regression.model.isin(['JanusDDG','EvoPhi-DDG','Ensemble'])].itertuples():
        lines.append(f'| {r.dataset} | {r.n} | {r.model} | {r.pearson:.2f} | {r.spearman:.2f} | {r.mae:.2f} | {r.rmse:.2f} |')
    lines += ['','S461 is a curation sensitivity analysis of S669, rather than an independent test. '
              'S96 also shares biological variants with S669; see the overlap audit in `generated/`. '
              'The ensemble and EvoPhi-DDG both round to Pearson 0.66 on PTmul-NR; the ensemble does not exceed '
              'EvoPhi-DDG in unrounded Pearson on the full set.','',
              '## Numerical source files','',
              '- `results/regression.csv`: individual models, ensemble and multiple-substitution sensitivity analyses.',
              '- `results/classification.csv`: individual classifiers and ensemble after the experimental exclusion band.',
              '- `results/bootstrap.csv`: paired differences in Pearson and protein-bootstrap intervals.',
              '- `results/same_site_transitivity.csv`: same-site closure summaries.',
              '- `results/multisite_transitivity.csv`: S669 multisite closure summaries.',
              '- `results/ensemble_controls.csv`: original-head, projected-head and zero-bias controls.',
              '- `generated/`: LaTeX tables and complete statistical outputs, including additional intervals.','']
    (ROOT/'docs/results.md').write_text('\n'.join(lines))
    print('Exported numerical results to results/ and docs/results.md.')


if __name__=='__main__':
    main()
