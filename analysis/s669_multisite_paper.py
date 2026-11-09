"""S669 multisite paths: figures, tables and paired protein-bootstrap summaries."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from figure_formatting import format_axes,correlation
from table_formatting import performance_rows
HERE=Path(__file__).resolve().parent.parent
MODELS=['MSAesm_ddG original','MSAesm_ddG antisymmetric','EvoPhi-DDG','JanusDDG','Ensemble EvoPhi-DDG']
NAMES={'MSAesm_ddG original':'MSAesm_ddG (orig.)','MSAesm_ddG antisymmetric':'MSAesm_ddG (AS)',
       'EvoPhi-DDG':'EvoPhi-DDG','JanusDDG':'JanusDDG','Ensemble EvoPhi-DDG':'Ensemble'}
COLORS=dict(zip(MODELS,['#777d88','#c98b68','#b86829','#315da8','#168274']))


def load():
    frames=[pd.read_csv(HERE/f'data/s669_multisite_{name}_paths.csv.gz') for name in ['msa','janus']]
    f=pd.concat(frames,ignore_index=True);f=f[f.model.isin(MODELS)].copy()
    groups=pd.read_csv(HERE/'data/s669.csv').set_index('code').group
    assert len(f)==5*33450
    f['protein_group']=groups.loc[f.code].to_numpy()
    assert f.protein_group.nunique()==87
    assert f.groupby(['model','n_intermediates','repeat']).size().eq(669).all()
    np.testing.assert_allclose(f.path_sum-f.direct,f.closure,rtol=0,atol=1e-12)
    return f


def calculate():
    f=load();records=[];repeat=[];paired={}
    for (model,n,rep),g in f.groupby(['model','n_intermediates','repeat'],sort=False):
        repeat.append(dict(model=model,n_intermediates=int(n),repeat=int(rep),
            pearson=float(pearsonr(g.direct,g.path_sum).statistic),closure_mae=float(abs(g.closure).mean()),
            closure_rmse=float(np.sqrt(np.mean(g.closure**2))),max_absolute_closure=float(abs(g.closure).max())))
    repeat=pd.DataFrame(repeat);repeat.to_csv(HERE/'data/s669_multisite_repeat_metrics.csv',index=False)
    for (model,n),g in repeat.groupby(['model','n_intermediates'],sort=False):
        records.append(dict(model=model,n_intermediates=int(n),pearson_mean=float(g.pearson.mean()),pearson_std=float(g.pearson.std()),
            closure_mae_mean=float(g.closure_mae.mean()),closure_rmse_mean=float(g.closure_rmse.mean()),closure_rmse_std=float(g.closure_rmse.std()),max_absolute_closure=float(g.max_absolute_closure.max())))
    summary=pd.DataFrame(records);summary.to_csv(HERE/'data/s669_multisite_summary.csv',index=False)
    # Bootstrap complete WT-protein groups; retain all ten repeats in each draw.
    for n in [2,4,6,8,10]:
        g=f[(f.n_intermediates==n)&f.model.isin(['JanusDDG','EvoPhi-DDG','Ensemble EvoPhi-DDG'])]
        proteins=sorted(g.protein_group.unique());models=['JanusDDG','EvoPhi-DDG','Ensemble EvoPhi-DDG']
        counts=np.zeros((87,10));sq=np.zeros((87,10,3))
        for i,protein in enumerate(proteins):
            for rep in range(10):
                a=g[(g.protein_group==protein)&(g['repeat']==rep)]
                counts[i,rep]=sum(a.model.eq(models[0]))
                for j,name in enumerate(models):sq[i,rep,j]=np.sum(a.loc[a.model.eq(name),'closure'].to_numpy()**2)
        assert np.all(counts.sum(0)==669)
        rng=np.random.default_rng(42);weights=rng.multinomial(87,np.ones(87)/87,size=5000)
        boot=np.sqrt(np.einsum('bg,grm->brm',weights,sq)/(weights@counts)[:,:,None]).mean(1)
        point=np.sqrt(sq.sum(0)/counts.sum(0)[:,None]).mean(0)
        pooled_protein=np.sqrt(sq.sum(1)/counts.sum(1)[:,None])
        detail={}
        for j,name in enumerate(models[:2]):
            gain=boot[:,j]-boot[:,2]
            detail[name]=dict(rmse_gain=float(point[j]-point[2]),ci95=np.quantile(gain,[.025,.975]).tolist(),
                percent_reduction=float(100*(point[j]-point[2])/point[j]),n_proteins_improved=int(sum(pooled_protein[:,2]<pooled_protein[:,j])))
        j=g[g.model.eq('JanusDDG')].set_index(['code','repeat']).sort_index()
        m=g[g.model.eq('EvoPhi-DDG')].set_index(['code','repeat']).sort_index()
        e=g[g.model.eq('Ensemble EvoPhi-DDG')].set_index(['code','repeat']).sort_index()
        for col in ['direct','path_sum','closure']:np.testing.assert_allclose(e[col],(j[col]+m[col])/2,rtol=0,atol=1e-12)
        detail['closure_error_correlation']=float(pearsonr(j.closure,m.closure).statistic)
        detail['ensemble_rmse_ci95']=np.quantile(boot[:,2],[.025,.975]).tolist()
        paired[str(n)]=detail
    result=dict(protocol=dict(dataset='S669',n_variants=669,n_protein_groups=87,n_contexts=92,n_states=196196,n_directed_edges=229558,n_paths=33450,
        n_intermediates=[2,4,6,8,10],repeats=10,seed=42,no_experimental_labels=True,
        uncertainty='5000 paired whole-WT-protein bootstrap resamples, seed42; mean of ten replicate RMSEs; all repeats retained together'),
        summary=records,paired=paired)
    (HERE/'generated/s669_multisite_statistics.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(paired['10'],indent=2))
    return result


def table(result):
    summary=pd.DataFrame(result['summary'])
    rows=[r'\begin{tabular}{rlrrrr}\toprule $N$ & Predictor & Pearson & Closure MAE & Closure RMSE & Max. defect\\\midrule']
    labels=[r'MSAesm\_ddG original',r'MSAesm\_ddG antisymmetrized','EvoPhi-DDG','JanusDDG','Equal-weight ensemble']
    for n in [2,4,6,8,10]:
        g=summary[summary.n_intermediates.eq(n)].set_index('model').loc[MODELS]
        values=g[['pearson_mean','closure_mae_mean','closure_rmse_mean','max_absolute_closure']].to_numpy()
        lines=performance_rows(labels,values,['corr','.2f','.2f','.2f'],['max','min','min','min'])
        rows.extend(str(n)+' & '+line for line in lines);rows.append(r'\addlinespace')
    rows.append(r'\bottomrule\end{tabular}')
    (HERE/'generated/multisite_rows.tex').write_text('\n'.join(rows)+'\n')


def figures():
    f=load();stats=json.loads((HERE/'generated/s669_multisite_statistics.json').read_text());summary=pd.DataFrame(stats['summary'])
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.spines.top':False,'axes.spines.right':False,'axes.titlesize':9,'pdf.fonttype':42})
    fig,axes=plt.subplots(2,3,figsize=(8.0,5.8),layout='constrained')
    for ax,name,letter in zip(axes[0],['EvoPhi-DDG','JanusDDG','Ensemble EvoPhi-DDG'],'abc'):
        g=f[f.model.eq(name)&f.n_intermediates.eq(10)];s=summary[summary.model.eq(name)&summary.n_intermediates.eq(10)].iloc[0]
        lo=min(g.direct.min(),g.path_sum.min())-.15;hi=max(g.direct.max(),g.path_sum.max())+.15
        ax.scatter(g.direct,g.path_sum,s=4,alpha=.17,color=COLORS[name],rasterized=True,edgecolors='none');ax.plot([lo,hi],[lo,hi],'k--',lw=.7)
        ax.set(xlabel='Direct score (kcal/mol)',ylabel='Path sum (kcal/mol)' if letter=='a' else '',title=f'{letter}  {NAMES[name]}: 10 intermediates')
        ax.text(.04,.96,f'Mean r = {correlation(s.pearson_mean)}\nMean closure RMSE = {s.closure_rmse_mean:.2f}',transform=ax.transAxes,va='top',fontsize=7)
    for name in MODELS:
        g=summary[summary.model.eq(name)]
        axes[1,0].errorbar(g.n_intermediates,1-g.pearson_mean,yerr=g.pearson_std,marker='o',ms=3,capsize=2,color=COLORS[name],label=NAMES[name],lw=1)
        axes[1,1].errorbar(g.n_intermediates,g.closure_rmse_mean,yerr=g.closure_rmse_std,marker='o',ms=3,capsize=2,color=COLORS[name],label=NAMES[name],lw=1)
    axes[1,0].set(xlabel='Intermediate sequence states',ylabel='Pearson deficit (1 − r)',yscale='log',title='d  Multisite correlation deficit')
    axes[1,1].set(xlabel='Intermediate sequence states',ylabel='Closure RMSE (kcal/mol)',yscale='log',title='e  Multisite closure error')
    handles,labels=axes[1,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,-.01),ncol=3,frameon=False,fontsize=7)
    axes[1,1].set_ylim(bottom=.02)
    for ax in axes[1,:2]:ax.set_xticks([2,4,6,8,10]);ax.grid(alpha=.15)
    a=f[f.model.eq('JanusDDG')&f.n_intermediates.eq(10)].set_index(['code','repeat']).sort_index()
    b=f[f.model.eq('EvoPhi-DDG')&f.n_intermediates.eq(10)].set_index(['code','repeat']).sort_index()
    ax=axes[1,2];ax.scatter(a.closure,b.closure,s=4,alpha=.17,color='#168274',edgecolors='none',rasterized=True);ax.axhline(0,c='grey',ls='--',lw=.7);ax.axvline(0,c='grey',ls='--',lw=.7)
    ax.set(xlabel='JanusDDG closure error (kcal/mol)',ylabel='EvoPhi-DDG closure error (kcal/mol)',title=f'f  Closure-error r = {stats["paired"]["10"]["closure_error_correlation"]:.2f}')
    format_axes(fig,skip=[(ax,'x') for ax in axes[1,:2]])
    fig.savefig(HERE/'figures/multisite_transitivity.pdf',bbox_inches='tight');fig.savefig(HERE/'figures/multisite_transitivity.png',dpi=220,bbox_inches='tight');plt.close(fig)


def main():
    result=calculate();table(result);figures()

if __name__=='__main__':main()
