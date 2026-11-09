"""Compare fixed Janus ensembles before/after MSA antisymmetrization.

Reuses cached original checkpoint predictions and identical S669 pathways.
No inference, alignment search, fitting, calibration or threshold selection.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr,spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from classification import score as classification_score
from table_formatting import best_cells,bold,performance_rows
from figure_formatting import label as display_label, format_axes
from label_policy import apply_labels, POLICY
from evophi_policy import apply_model,physical_edges

HERE=Path(__file__).resolve().parent.parent
DATASETS={'S669':'s669.csv','S461':'s461.csv','S96':'s96.csv','PTmul-NR':'ptmul_nr.csv'}
COLORS={'E_original':'#8c668e','E_anti':'#b86829','E_phi':'#168274','JanusDDG':'#315da8'}


def score(y,p,reverse):
    y=np.asarray(y);p=np.asarray(p);defect=p+np.asarray(reverse);mask=abs(y)>.5
    c=classification_score(y[mask]>.5,p[mask])
    return dict(n=len(y),pearson=float(pearsonr(y,p).statistic),spearman=float(spearmanr(y,p).statistic),
        mae=float(abs(p-y).mean()),rmse=float(np.sqrt(np.mean((p-y)**2))),
        reversal_mae=float(abs(defect).mean()),reversal_rmse=float(np.sqrt(np.mean(defect**2))),
        reversal_max=float(abs(defect).max()),classification=c)


def bootstrap(frame):
    groups=list(frame.groupby('group',sort=False).indices.values());rng=np.random.default_rng(42)
    y=frame.y.to_numpy();o=frame.E_original.to_numpy();a=frame.E_anti.to_numpy();values=[]
    def metrics(idx):
        def corr(p):
            yy=y[idx]-y[idx].mean();pp=p[idx]-p[idx].mean();return np.dot(yy,pp)/np.sqrt(np.dot(yy,yy)*np.dot(pp,pp))
        return [corr(a)-corr(o),np.mean(abs(o[idx]-y[idx]))-np.mean(abs(a[idx]-y[idx])),
            np.sqrt(np.mean((o[idx]-y[idx])**2))-np.sqrt(np.mean((a[idx]-y[idx])**2))]
    for _ in range(5000):
        idx=np.concatenate([groups[k] for k in rng.integers(len(groups),size=len(groups))])
        values.append(metrics(idx))
    point=metrics(np.arange(len(frame)));bounds=np.quantile(values,[.025,.975],axis=0)
    return dict(n_groups=len(groups),n_resamples=5000,**{name:dict(value=float(point[i]),ci95=bounds[:,i].tolist())
        for i,name in enumerate(['pearson_anti_minus_original','mae_gain_anti','rmse_gain_anti'])})


def paths():
    cache=HERE/'data/physics'
    f=pd.read_csv(cache/'variants.csv');site=f.site_index.to_numpy();aa='ACDEFGHIKLMNPQRSTVWY'
    wt=np.array([aa.index(r.wt_seq[int(r.pos)]) for r in f.itertuples()]);mut=np.array([aa.index(r.mut_seq[int(r.pos)]) for r in f.itertuples()])
    j=np.load(cache/'janus_edges.npy').astype(float)
    original=np.load(cache/'original_msa_edges.npy').astype(float)
    anti=np.load(cache/'msa_edges.npy').astype(float)
    assert np.max(abs((original-original.transpose(0,2,1))/2-anti))<2e-5
    draws=np.load(cache/'paths.npz');records=[]
    models={'JanusDDG':j,'E_original':(j+original)/2,'E_anti':(j+anti)/2,'E_phi':(j+physical_edges())/2}
    for name,edges in models.items():
        direct=edges[site,wt,mut]
        for n in [1,3,5,7,9]:
            for repeat in range(10):
                chain=np.column_stack([wt,draws[str(n)][repeat],mut])
                total=edges[site[:,None],chain[:,:-1],chain[:,1:]].sum(axis=1);error=total-direct
                records.append(dict(model=name,n_intermediates=n,repeat=repeat,
                    pearson=float(pearsonr(direct,total).statistic),mae=float(abs(error).mean()),
                    rmse=float(np.sqrt(np.mean(error**2))),bias=float(error.mean())))
    repeats=pd.DataFrame(records);repeats.to_csv(HERE/'data/ensemble_ablation_path_repeats.csv',index=False)
    summary=repeats.groupby(['model','n_intermediates'],sort=False).agg(pearson_mean=('pearson','mean'),
        pearson_std=('pearson','std'),mae_mean=('mae','mean'),rmse_mean=('rmse','mean'),bias_mean=('bias','mean')).reset_index()
    summary.to_csv(HERE/'data/ensemble_ablation_path_summary.csv',index=False)
    return summary


def table(stats):
    rows=[r'\begin{tabular}{lrrrrr}\toprule Fixed predictor & Pearson & MAE & RMSE & Reversal MAE & ROC-AUC\\\midrule']
    for name,result in stats.items():
        rows.append(r'\multicolumn{6}{l}{\textit{'+name+r'}}\\')
        labels=[r'JanusDDG + original MSAesm\_ddG',r'JanusDDG + antisymmetric MSAesm\_ddG','JanusDDG + EvoPhi-DDG']
        keys=['E_original','E_anti','E_phi']
        if name=='PTmul-NR':
            labels=['Janus + directional WT star','Janus + whole-projected raw star','Janus + component/whole projection','JanusDDG + EvoPhi-DDG']
            keys=['E_original','E_whole_only','E_anti','E_phi']
        values=[[result[k][m] if m!='roc_auc' else result[k]['classification']['roc_auc']
            for m in ['pearson','mae','rmse','reversal_mae','roc_auc']] for k in keys]
        rows.extend(performance_rows(labels,values,['.4f','.3f','.3f','.4f','.4f'],['max','min','min','min','max']))
        rows.append(r'\addlinespace')
    rows.append(r'\bottomrule\end{tabular}')
    (HERE/'generated/ensemble_ablation_rows.tex').write_text('\n'.join(rows)+'\n')


def figure(stats,physical):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':9,
        'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'savefig.bbox':'tight'})
    fig,axes=plt.subplots(1,3,figsize=(7.2,2.9),layout='constrained')
    x=np.arange(4)
    for offset,key,label in [(-.24,'E_original',display_label('Mean original MSA')),(0.,'E_anti',display_label('Mean anti')),(.24,'E_phi','Mean (EvoPhi-DDG)')]:
        axes[0].bar(x+offset,[stats[d][key]['pearson'] for d in DATASETS],width=.22,color=COLORS[key],label=label)
    axes[0].set(xticks=x,xticklabels=['S669','S461','S96','PTmul'],ylim=(.5,.80),ylabel='Pearson vs experiment',title='a  Prediction accuracy')
    axes[0].tick_params(axis='x',labelsize=7)
    for key,label in [('JanusDDG','JanusDDG'),('E_original',display_label('Mean original MSA')),('E_anti',display_label('Mean anti')),('E_phi','Mean (EvoPhi-DDG)')]:
        f=physical[physical.model.eq(key)]
        axes[1].errorbar(f.n_intermediates,f.pearson_mean,yerr=f.pearson_std,color=COLORS[key],label=label,marker='o',ms=3,capsize=2,lw=1)
        axes[2].plot(f.n_intermediates,f.rmse_mean,color=COLORS[key],marker='o',ms=3,lw=1)
    axes[1].set(xlabel='Intermediate residues',ylabel='Pearson: direct vs path',title='b  Same-site transitivity')
    axes[2].set(xlabel='Intermediate residues',ylabel='Closure RMSE (kcal/mol)',title='c  Absolute pathway error')
    for ax in axes[1:]:ax.set_xticks([1,3,5,7,9]);ax.grid(alpha=.15)
    handles,labels=axes[1].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,1.15),ncol=2,frameon=False,fontsize=7)
    format_axes(fig,skip=[(ax,'x') for ax in axes])
    fig.savefig(HERE/'figures/ensemble_ablation.pdf');fig.savefig(HERE/'figures/ensemble_ablation.png',dpi=200);plt.close(fig)


def main():
    source=HERE/'data/ensemble_ablation_predictions.csv'
    f=pd.read_csv(source)
    f=pd.concat([apply_labels(f.loc[f.dataset.eq(name)],name) for name in DATASETS],ignore_index=True)
    stats={}
    for name in DATASETS:
        g=f[f.dataset.eq(name)].reset_index(drop=True)
        stats[name]={key:score(g.y,g[key],g[key+'_reverse']) for key in ['E_original','E_anti','E_phi']}
        if name=='PTmul-NR':
            stats[name]['E_whole_only']=score(g.y,g.E_whole_only,-g.E_whole_only)
            stats[name]['paired_bootstrap_whole_control']=bootstrap(g.assign(E_original=g.E_whole_only))
        stats[name]['paired_bootstrap']=bootstrap(g)
        stats[name]['paired_bootstrap_phi']={k.replace('anti','phi'):v for k,v in bootstrap(g.assign(E_anti=g.E_phi)).items()}
    physical=paths()
    result=dict(protocol=dict(checkpoints='Same published frozen Janus and MSA checkpoints',
        label_policy=POLICY,
        ensemble='Fixed arithmetic mean, no fitting',single_ablation='Compare original head, published-bias antisymmetric head, and zero-hidden-bias EvoPhi-DDG',
        multiple_primary='Directional WT raw-star versus component-and-whole projection: also changes reverse background',
        multiple_control='Raw-head whole projection versus component/whole projection: same star backgrounds and whole projection',
        physical='Identical cached S669 states and ten paths per N; Janus cache held fixed',
        bootstrap='5000 paired WT-group draws, seed42; delta-r anti minus original, error gain original minus anti',
        classification='Experimental |y|>0.5, score>0, same retained variants; no threshold fitting'),datasets=stats,
        physical_nine=physical[physical.n_intermediates.eq(9)].to_dict(orient='records'))
    (HERE/'generated/ensemble_ablation_statistics.json').write_text(json.dumps(result,indent=2))
    table(stats);figure(stats,physical)
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
