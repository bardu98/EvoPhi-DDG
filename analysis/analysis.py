"""Create manuscript tables, figures and protein-bootstrap statistics."""
from pathlib import Path
import json

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from table_formatting import best_cells, bold, performance_rows, formatted
from figure_formatting import label as display_label, format_axes, scientific_two_decimals
from matplotlib.ticker import FuncFormatter
from label_policy import apply_labels
from evophi_policy import apply_model

HERE=Path(__file__).resolve().parent.parent
ROOT=HERE.parent
REP=HERE/'external/PLM4Muts/reproduction'
for name in ['data','generated','figures']:(HERE/name).mkdir(exist_ok=True)
COLORS={'JanusDDG':'#315da8','EvoPhi-DDG':'#b86829','Ensemble':'#168274','MSAesm original':'#777d88','MSAesm anti':'#c98b68'}
DATASET_COLORS={'S669':'#315da8','S461':'#87549a','S96':'#bd8525','PTmul-NR':'#168274'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.spines.top':False,
 'axes.spines.right':False,'axes.labelsize':8,'axes.titlesize':8.5,'legend.fontsize':7,
 'pdf.fonttype':42,'savefig.bbox':'tight'})


def score(y,p):
    y=np.asarray(y);p=np.asarray(p)
    return dict(n=len(y),pearson=float(pearsonr(y,p).statistic),
        spearman=float(spearmanr(y,p).statistic),mae=float(abs(y-p).mean()),rmse=float(np.sqrt(np.mean((y-p)**2))))


def corr(y,p):
    a=y-y.mean();b=p-p.mean()
    return np.dot(a,b)/np.sqrt(np.dot(a,a)*np.dot(b,b))


def bootstrap(frame):
    rng=np.random.default_rng(42);groups=list(frame.groupby('group',sort=False).indices.values())
    y=frame.y.to_numpy();j=frame.J.to_numpy();m=frame.M.to_numpy();e=(j+m)/2
    samples=[]
    for _ in range(5000):
        idx=np.concatenate([groups[i] for i in rng.integers(len(groups),size=len(groups))])
        if np.std(y[idx])==0:continue
        samples.append([corr(y[idx],e[idx])-corr(y[idx],j[idx]),corr(y[idx],e[idx])-corr(y[idx],m[idx])])
    samples=np.asarray(samples);np.save(HERE/'generated'/f'bootstrap_{frame.attrs["name"]}.npy',samples)
    return dict(n_groups=len(groups),n_resamples=len(samples),
        ensemble_minus_J=dict(delta=float(corr(y,e)-corr(y,j)),ci95=np.quantile(samples[:,0],[.025,.975]).tolist()),
        ensemble_minus_M=dict(delta=float(corr(y,e)-corr(y,m)),ci95=np.quantile(samples[:,1],[.025,.975]).tolist()))


def load():
    s=apply_model(apply_labels(pd.read_csv(HERE/'data/s669.csv'),'S669'),'S669')
    t=apply_model(apply_labels(pd.read_csv(HERE/'data/ptmul_nr.csv'),'PTmul-NR'),'PTmul-NR')
    physics=pd.read_csv(HERE/'data/transitivity.csv')
    s.attrs['name']='S669';t.attrs['name']='PTmulNR'
    return s,t,physics


def dataset_frames(s,t):
    s=apply_model(apply_labels(s,'S669'),'S669');s.attrs['name']='S669'
    t=apply_model(apply_labels(t,'PTmul-NR'),'PTmul-NR');t.attrs['name']='PTmulNR'
    curated=apply_model(apply_labels(pd.read_csv(HERE/'data/s461.csv'),'S461'),'S461');curated.attrs['name']='S461'
    extra=apply_model(apply_labels(pd.read_csv(HERE/'data/s96.csv'),'S96'),'S96');extra.attrs['name']='S96'
    return [('S669',s),('S461',curated),('S96',extra),('PTmul-NR',t)]


def tables(s,t,physics,stats):
    for name,frame in dataset_frames(s,t)[:-1]:
        labels=[];values=[]
        for label,col in [('JanusDDG','J'),('MSAesm original','M_original'),('MSAesm antisymmetrized','M_anti'),('EvoPhi-DDG','M'),('Equal-weight ensemble','E')]:
            p=(frame.J+frame.M)/2 if col=='E' else frame[col]
            r=score(frame.y,p);labels.append(label)
            values.append([r['pearson'],r['spearman'],r['mae'],r['rmse']])
        rows=performance_rows(labels,values,['.4f','.4f','.3f','.3f'],['max','max','min','min'])
        (HERE/f'generated/{name.lower()}_rows.tex').write_text('\n'.join(rows)+'\n')
    labels=[];values=[]
    for label,col in [('JanusDDG','J'),('MSAesm: whole-variant antisymmetry','M_anti'),('EvoPhi-DDG','M'),('Equal-weight ensemble','E'),
                      ('EvoPhi-DDG: sequential, both directions','M_sequence'),('EvoPhi-DDG: all mutation orders','M_orders')]:
        p=(t.J+t.M)/2 if col=='E' else t[col];r=score(t.y,p)
        filtered=score(t.loc[t.retained,'y'],p[t.retained])
        labels.append(label);values.append([r['pearson'],r['mae'],r['rmse'],filtered['pearson']])
    rows=performance_rows(labels,values,['.4f','.3f','.3f','.4f'],['max','min','min','max'])
    (HERE/'generated/ptmul_rows.tex').write_text('\n'.join(rows)+'\n')
    rows=[]
    for n in [1,3,5,7,9]:
        vals={r.model:r for r in physics[physics.n_intermediates.eq(n)].itertuples()}
        keys=['MSAesm original','MSAesm_anti','EvoPhi-DDG','JanusDDG','Mean']
        values=[[vals[k].pearson_mean,vals[k].mae_mean,vals[k].rmse_mean] for k in keys]
        labels=[r'MSAesm\_ddG original',r'MSAesm\_ddG antisymmetrized','EvoPhi-DDG','JanusDDG','Equal-weight ensemble']
        for row in performance_rows(labels,values,['corr','.2f','.2f'],['max','min','min']):
            rows.append(str(n)+' & '+row)
    (HERE/'generated/transitivity_rows.tex').write_text('\n'.join(rows)+'\n')
    rows=[]
    for name,data in dataset_frames(s,t):
        labels=[];values=[];intervals=[]
        for base in ['J','M']:
            r=stats[name]['bootstrap']['ensemble_minus_'+base]
            labels.append('JanusDDG' if base=='J' else 'EvoPhi-DDG');values.append([r['delta']])
            intervals.append('['+formatted(r['ci95'][0],'+.2f')+', '+formatted(r['ci95'][1],'+.2f')+']')
        cells,winning=best_cells(values,['+.4f'],['max'])
        for label,cell,win,interval in zip(labels,cells,winning,intervals):
            rows.append(' & '.join([name,bold(label) if win else label,cell[0],interval])+r' \\')
    (HERE/'generated/bootstrap_rows.tex').write_text('\n'.join(rows)+'\n')
    headers={
        's669_rows':r'\begin{tabular}{lrrrr}\toprule Predictor & Pearson & Spearman & MAE & RMSE\\\midrule',
        's461_rows':r'\begin{tabular}{lrrrr}\toprule Predictor & Pearson & Spearman & MAE & RMSE\\\midrule',
        's96_rows':r'\begin{tabular}{lrrrr}\toprule Predictor & Pearson & Spearman & MAE & RMSE\\\midrule',
        'ptmul_rows':r'\begin{tabular}{lrrrr}\toprule Predictor & Pearson$_{82}$ & MAE$_{82}$ & RMSE$_{82}$ & Pearson$_{67}$\\\midrule',
        'transitivity_rows':r'\begin{tabular}{rlrrr}\toprule $N$ & Predictor & Pearson & Closure MAE & Closure RMSE\\\midrule',
        'bootstrap_rows':r'\begin{tabular}{llrr}\toprule Dataset & Comparator & $\Delta r$ & 95\% interval\\\midrule'}
    for name,header in headers.items():
        p=HERE/f'generated/{name}.tex';p.write_text(header+'\n'+p.read_text()+r'\bottomrule\end{tabular}'+'\n')


def accuracy_figure(s,t):
    datasets=dataset_frames(s,t)
    fig,axes=plt.subplots(len(datasets),3,figsize=(7.2,9.7),layout='constrained')
    for row,(name,f) in enumerate(datasets):
        if name=='S96':
            values=f[['y','J','M']].to_numpy()
            lim=(np.floor(values.min()-0.5),np.ceil(values.max()+0.5))
        else:
            lim=(-9,8) if name=='PTmul-NR' else (-10,7)
        for col,(label,pred) in enumerate([('JanusDDG',f.J),('EvoPhi-DDG',f.M),('Ensemble',(f.J+f.M)/2)]):
            a=axes[row,col];a.scatter(f.y,pred,s=25 if name=='PTmul-NR' else 13,alpha=.45,color=COLORS[label],edgecolors='none',rasterized=True)
            assert min(f.y.min(),pred.min())>=lim[0] and max(f.y.max(),pred.max())<=lim[1]
            a.plot(lim,lim,'k--',lw=.8);a.set(xlim=lim,ylim=lim)
            r=score(f.y,pred);a.text(.04,.94,f'r = {r["pearson"]:.2f}\nMAE = {r["mae"]:.2f}',transform=a.transAxes,va='top',fontsize=8,
                bbox=dict(facecolor='white',alpha=.8,edgecolor='none'))
            a.set_title(f'{chr(97+row*3+col)}  {name}: {display_label(label)}')
            if col==0:a.set_ylabel('Predicted $\\Delta\\Delta G$ (kcal/mol)')
            a.set_xlabel('Experimental $\\Delta\\Delta G$ (kcal/mol)')
    format_axes(fig)
    fig.savefig(HERE/'figures/accuracy.pdf');fig.savefig(HERE/'figures/accuracy.png',dpi=200);plt.close(fig)


def physics_figure(s,physics):
    fig,ax=plt.subplots(1,3,figsize=(7.2,3.3),layout='constrained')
    ax[0].scatter(s.M_original,s.M_reverse,s=9,alpha=.35,color=COLORS['MSAesm original'],label=display_label('MSAesm original'),rasterized=True)
    ax[0].scatter(s.M,-s.M,s=6,alpha=.4,color=COLORS['EvoPhi-DDG'],label=display_label('EvoPhi-DDG'),rasterized=True)
    ax[0].plot([-6,6],[6,-6],'k--',lw=.8);ax[0].set(xlim=(-6,6),ylim=(-6,6),xlabel='Forward prediction (kcal/mol)',ylabel='Reverse prediction (kcal/mol)',title='a  Reversal consistency')
    ax[0].legend(loc='upper right')
    for key,name in [('MSAesm original','MSAesm original'),('MSAesm_anti','MSAesm anti'),('EvoPhi-DDG','EvoPhi-DDG'),('JanusDDG','JanusDDG'),('Mean','Ensemble')]:
        g=physics[physics.model.eq(key)]
        ax[1].errorbar(g.n_intermediates,1-g.pearson_mean,yerr=g.pearson_std,marker='o',ms=4,capsize=2,label=display_label(name),color=COLORS[name],lw=1.5)
        ax[2].plot(g.n_intermediates,g.rmse_mean,marker='o',ms=4,color=COLORS[name],label=display_label(name))
    ax[1].set(xlabel='Intermediate residues',ylabel='Pearson deficit (1 − r)',yscale='symlog',ylim=(0,.3),title='b  Same-site correlation deficit')
    ax[2].set(xlabel='Intermediate residues',ylabel='Closure RMSE (kcal/mol)',yscale='symlog',ylim=(0,1),title='c  Absolute pathway error')
    for a in ax[1:]:a.set_xticks([1,3,5,7,9]);a.grid(alpha=.15)
    ax[1].set_yscale('symlog',linthresh=1e-7)
    ax[2].set_yscale('symlog',linthresh=1e-7)
    ax[2].set_yticks([0,1e-6,1e-4,1e-2,1])
    ax[1].set_yticks([0,1e-6,1e-5,1e-4,1e-3,1e-2,.1])
    handles,labels=ax[1].get_legend_handles_labels()
    fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.5,1.02),ncol=3,frameon=False,fontsize=7)
    fig.text(.5,-.03,'EvoPhi-DDG has exact same-site transitivity; nonzero closure values reflect floating-point roundoff.',ha='center',fontsize=6)
    format_axes(fig,skip=[(ax[1],'x'),(ax[2],'x')])
    for axis in ax[1:]:axis.yaxis.set_major_formatter(FuncFormatter(lambda value,_:'0' if value==0 else scientific_two_decimals(value) if value<.01 else f'{value:.2f}'))
    fig.savefig(HERE/'figures/physics.pdf');fig.savefig(HERE/'figures/physics.png',dpi=200);plt.close(fig)


def figures(s,t,physics,stats):
    physics_figure(s,physics)
    accuracy_figure(s,t)
    f=dict(dataset_frames(s,t))['S96']
    fig,axes=plt.subplots(1,3,figsize=(7.2,2.6),layout='constrained')
    lower=np.floor(min(f.y.min(),f.J.min(),f.M.min())-0.5)
    upper=np.ceil(max(f.y.max(),f.J.max(),f.M.max())+0.5)
    for col,(label,pred) in enumerate([('JanusDDG',f.J),('EvoPhi-DDG',f.M),('Ensemble',(f.J+f.M)/2)]):
        a=axes[col];a.scatter(f.y,pred,s=20,alpha=.55,color=COLORS[label],edgecolors='none',rasterized=True)
        a.plot([lower,upper],[lower,upper],'k--',lw=.8)
        a.set(xlim=(lower,upper),ylim=(lower,upper),xlabel='Experimental $\\Delta\\Delta G$ (kcal/mol)',title=f'{chr(97+col)}  S96: {display_label(label)}')
        if col==0:a.set_ylabel('Predicted $\\Delta\\Delta G$ (kcal/mol)')
        r=score(f.y,pred);a.text(.04,.95,f'r = {r["pearson"]:.2f}\nMAE = {r["mae"]:.2f}',transform=a.transAxes,va='top',
            bbox=dict(facecolor='white',alpha=.8,edgecolor='none'))
    format_axes(fig)
    fig.savefig(HERE/'figures/accuracy_s96.pdf');fig.savefig(HERE/'figures/accuracy_s96.png',dpi=200);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(7.2,3.8),layout='constrained')
    ej=s.J-s.y;em=s.M-s.y
    axes[0].hexbin(ej,em,gridsize=35,mincnt=1,cmap='viridis',linewidths=0,rasterized=True)
    axes[0].axhline(0,color='white',lw=.8);axes[0].axvline(0,color='white',lw=.8)
    axes[0].set(xlabel='JanusDDG residual (kcal/mol)',ylabel='EvoPhi-DDG residual (kcal/mol)',title='a  Residual complementarity')
    axes[0].text(.04,.95,f'Residual r = {stats["S669"]["residual_correlation"]:.2f}',transform=axes[0].transAxes,va='top',bbox=dict(facecolor='white',alpha=.9,edgecolor='none'))
    labels=['vs JanusDDG','vs EvoPhi-DDG','vs both'];locations=np.arange(3)
    for offset,(name,data) in zip([-.30,-.10,.10,.30],dataset_frames(s,t)):
        e=abs((data.J+data.M)/2-data.y);j=abs(data.J-data.y);m=abs(data.M-data.y)
        values=[(e<j).mean(),(e<m).mean(),(e<np.minimum(j,m)).mean()]
        axes[1].bar(locations+offset,100*np.array(values),width=.18,label=name,color=DATASET_COLORS[name])
    axes[1].set(xticks=locations,xticklabels=labels,ylim=(0,75),ylabel='Variants with lower absolute error (%)',title='b  Per-variant ensemble benefit');axes[1].legend()
    labels=[];ys=[]
    for name in ['S669','S461','S96','PTmul-NR']:
        for base in ['J','M']:
            b=stats[name]['bootstrap']['ensemble_minus_'+base];y=len(labels)
            axes[2].errorbar(b['delta'],y,xerr=[[b['delta']-b['ci95'][0]],[b['ci95'][1]-b['delta']]],fmt='o',capsize=3,color=DATASET_COLORS[name])
            labels.append(f'{name}:\nvs '+('JanusDDG' if base=='J' else 'EvoPhi-DDG'))
    axes[2].axvline(0,color='#777777',ls='--',lw=.8);axes[2].set(yticks=range(len(labels)),yticklabels=labels,xlabel='Ensemble change\nin Pearson $r$',title='c  Paired protein bootstrap');axes[2].invert_yaxis()
    axes[1].tick_params(axis='x',labelsize=6.5)
    plt.setp(axes[1].get_xticklabels(),rotation=35,ha='right',rotation_mode='anchor')
    axes[2].tick_params(axis='y',labelsize=7)
    format_axes(fig,skip=[(axes[1],'x'),(axes[2],'y')])
    fig.savefig(HERE/'figures/complementarity.pdf');fig.savefig(HERE/'figures/complementarity.png',dpi=200);plt.close(fig)


def main():
    s,t,p=load();stats={}
    for name,f in dataset_frames(s,t):
        e=(f.J+f.M)/2;ej=f.J-f.y;em=f.M-f.y
        stats[name]=dict(J=score(f.y,f.J),M=score(f.y,f.M),E=score(f.y,e),bootstrap=bootstrap(f),
            residual_correlation=float(pearsonr(ej,em).statistic),
            opposite_residual_fraction=float((ej*em<0).mean()),
            improves_both_fraction=float((abs(e-f.y)<np.minimum(abs(ej),abs(em))).mean()),
            mse_decomposition=dict(J=float(np.mean(ej**2)),M=float(np.mean(em**2)),
                cross=float(np.mean(ej*em)),E=float(np.mean((e-f.y)**2))),
            filtered={k:score(f.loc[f.retained,'y'],v[f.retained]) for k,v in [('J',f.J),('M',f.M),('E',e)]})
        if 'M_original' in f:stats[name]['original']=score(f.y,f.M_original)
        if name in ['S461','S96']:stats[name].pop('filtered')
    (HERE/'generated/statistics.json').write_text(json.dumps(stats,indent=2))
    tables(s,t,p,stats);figures(s,t,p,stats)
    print(json.dumps(stats,indent=2))


if __name__=='__main__':main()
