"""Janus-style classification: exclude experimental |ddG| <= 0.5, predict sign.

Uses bundled manuscript scores; no inference, refitting or threshold search.
ROC area is checked independently with a rank-based calculation.
"""
from pathlib import Path
import json

import numpy as np
import pandas as pd
from scipy.stats import rankdata
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from table_formatting import performance_rows
from figure_formatting import label as display_label, format_axes
from label_policy import apply_labels, POLICY
from evophi_policy import apply_model

HERE=Path(__file__).resolve().parent.parent
COLORS={'JanusDDG':'#315da8','MSAesm original':'#777d88','MSAesm anti':'#c98b68','EvoPhi-DDG':'#b86829','Ensemble':'#168274'}
METRICS=['precision','recall','specificity','balanced_accuracy','roc_auc','mcc','f1','average_precision']


def roc(y,s):
    y=np.asarray(y,dtype=bool);s=np.asarray(s,dtype=float)
    order=np.argsort(-s,kind='stable');ys=y[order];ss=s[order]
    end=np.r_[np.flatnonzero(np.diff(ss)!=0),len(s)-1]
    tp=np.r_[0,np.cumsum(ys)[end]];fp=np.r_[0,np.cumsum(~ys)[end]]
    return fp/(~y).sum(),tp/y.sum(),tp,fp


def score(y,s):
    y=np.asarray(y,dtype=bool);s=np.asarray(s,dtype=float);p=s>0
    tp=int((p&y).sum());fp=int((p&~y).sum());tn=int((~p&~y).sum());fn=int((~p&y).sum())
    precision=tp/(tp+fp) if tp+fp else 0.;recall=tp/(tp+fn);specificity=tn/(tn+fp)
    fpr,tpr,tps,fps=roc(y,s)
    auc=float(np.trapz(tpr,fpr))
    ranks=rankdata(s,method='average');np_=int(y.sum());nn=int((~y).sum())
    rank_auc=float((ranks[y].sum()-np_*(np_+1)/2)/(np_*nn))
    assert np.isclose(auc,rank_auc,atol=1e-12)
    pred_precision=np.divide(tps[1:],tps[1:]+fps[1:])
    ap=float(np.sum(np.diff(tpr)*pred_precision))
    den=np.sqrt(float((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    return dict(n=len(y),positive=np_,negative=nn,predicted_positive=int(p.sum()),
        tp=tp,fp=fp,tn=tn,fn=fn,precision=precision,recall=recall,specificity=specificity,
        balanced_accuracy=(recall+specificity)/2,roc_auc=auc,
        mcc=(tp*tn-fp*fn)/den if den else 0.,f1=2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else 0.,
        average_precision=ap,accuracy=float((p==y).mean()))


def verify():
    y=np.array([False,False,True,True])
    a=score(y,[-2,-1,1,2]);assert a['roc_auc']==a['balanced_accuracy']==a['precision']==a['recall']==1.
    b=score(y,[2,1,-1,-2]);assert b['roc_auc']==b['balanced_accuracy']==0.
    c=score(y,[0,0,0,0]);assert c['roc_auc']==c['balanced_accuracy']==.5 and c['predicted_positive']==0


def evaluate(frame,label='y'):
    selected=frame[label].abs()>.5
    data=frame.loc[selected].copy().reset_index(drop=True);truth=data[label].to_numpy()>.5
    models={'JanusDDG':data.J.to_numpy()}
    if 'M_original' in data:models['MSAesm original']=data.M_original.to_numpy()
    models['MSAesm anti']=data.M_anti.to_numpy();models['EvoPhi-DDG']=data.M.to_numpy();models['Ensemble']=(data.J.to_numpy()+data.M.to_numpy())/2
    scores={name:score(truth,p) for name,p in models.items()}
    groups=list(data.groupby('group',sort=False).indices.values());rng=np.random.default_rng(42)
    boot={name:[] for name in models};valid=0
    for _ in range(1000):
        idx=np.concatenate([groups[i] for i in rng.integers(len(groups),size=len(groups))])
        if len(np.unique(truth[idx]))<2:continue
        valid+=1
        for name,p in models.items():
            r=score(truth[idx],p[idx]);boot[name].append([r[k] for k in METRICS])
    intervals={name:{key:np.quantile(np.array(values)[:,i],[.025,.975]).tolist() for i,key in enumerate(METRICS)} for name,values in boot.items()}
    difference={}
    for name in models:
        if name=='Ensemble':continue
        delta=np.array(boot['Ensemble'])-np.array(boot[name])
        difference[name]={key:dict(value=scores['Ensemble'][key]-scores[name][key],
            ci95=np.quantile(delta[:,i],[.025,.975]).tolist()) for i,key in enumerate(METRICS)}
    data['class_stabilizing']=truth
    for name,p in models.items():
        short={'JanusDDG':'janus','MSAesm original':'msa_original','MSAesm anti':'msa_anti','EvoPhi-DDG':'evophi','Ensemble':'ensemble'}[name]
        data[short+'_score']=p;data[short+'_predicted_stabilizing']=p>0
    result=dict(original_n=len(frame),retained_n=len(data),excluded_n=int((~selected).sum()),
        boundary_n=int((frame[label].abs()==.5).sum()),stabilizing=int(truth.sum()),destabilizing=int((~truth).sum()),
        n_protein_groups=len(groups),positive_prevalence=float(truth.mean()),
        bootstrap_draws=1000,bootstrap_valid=valid,metrics=scores,ci95=intervals,
        ensemble_minus_component=difference)
    return result,data,truth,models


def make_figure(bundles,stats,filename='classification'):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':9,
        'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    fig,axes=plt.subplots(len(bundles),2,figsize=(7.2,2.25*len(bundles)),squeeze=False,layout='constrained')
    for row,(dataset,(_,data,y,models)) in enumerate(bundles.items()):
        models={k:v for k,v in models.items() if k!='MSAesm anti'}
        ax=axes[row,0]
        for name,p in models.items():
            fpr,tpr,_,_=roc(y,p);auc=stats[dataset]['metrics'][name]['roc_auc']
            ax.plot(fpr,tpr,color=COLORS[name],lw=1.4,label=f'{display_label(name)} ({auc:.2f})')
        ax.plot([0,1],[0,1],'k--',lw=.7)
        ax.set(xlim=(0,1),ylim=(0,1.02),xlabel='False-positive rate',ylabel='True-positive rate',
            title=f'{chr(97+row*2)}  {dataset}: ROC')
        ax.legend(loc='lower right',fontsize=6.5,framealpha=.9)
        ax=axes[row,1];keys=['precision','recall','balanced_accuracy'];x=np.arange(3);width=.78/len(models)
        for i,name in enumerate(models):
            values=np.array([stats[dataset]['metrics'][name][k] for k in keys])
            bounds=np.array([stats[dataset]['ci95'][name][k] for k in keys])
            # Percentile bootstrap intervals need not contain the sample estimate.
            # Draw intervals independently instead of forcing nonnegative yerr.
            xx=x+(i-(len(models)-1)/2)*width
            ax.bar(xx,values,width=width,color=COLORS[name],alpha=.9)
            ax.vlines(xx,bounds[:,0],bounds[:,1],color='#222222',lw=.65)
            ax.hlines(bounds[:,0],xx-width*.22,xx+width*.22,color='#222222',lw=.65)
            ax.hlines(bounds[:,1],xx-width*.22,xx+width*.22,color='#222222',lw=.65)
        ax.set(xticks=x,xticklabels=['Precision','Recall','Balanced\naccuracy'],ylim=(0,1.05),ylabel='Classification metric',
            title=f'{chr(98+row*2)}  {dataset}: sign classification')
        ax.grid(axis='y',alpha=.15)
    format_axes(fig,skip=[(ax,'x') for ax in axes[:,1]])
    fig.savefig(HERE/f'figures/{filename}.pdf');fig.savefig(HERE/f'figures/{filename}.png',dpi=200);plt.close(fig)


def classification_table(stats):
    rows=[r'\begin{tabular}{lrrrrr}\toprule Predictor & Pred.$+$ & Precision & Recall & Bal. acc. & ROC-AUC\\\midrule']
    for name,r in stats.items():
        rows.append(r'\multicolumn{6}{l}{\textit{'+f'{name}: {r["retained_n"]} variants, {r["stabilizing"]} stabilizing / {r["destabilizing"]} destabilizing'+r'}}\\')
        keys=['predicted_positive','precision','recall','balanced_accuracy','roc_auc']
        labels=list(r['metrics']);values=[[r['metrics'][model][key] for key in keys] for model in labels]
        display_labels=[{'MSAesm original':r'MSAesm\_ddG original','MSAesm anti':r'MSAesm\_ddG antisymmetrized'}.get(label,label) for label in labels]
        rows.extend(performance_rows(display_labels,values,['d','.3f','.3f','.3f','.3f'],[None,'max','max','max','max']))
        rows.append(r'\addlinespace')
    rows.append(r'\bottomrule\end{tabular}')
    (HERE/'generated/classification_rows.tex').write_text('\n'.join(rows)+'\n')


def main():
    verify();stats={};bundles={}
    for name,filename in [('S669','s669.csv'),('S461','s461.csv'),('S96','s96.csv'),('PTmul-NR','ptmul_nr.csv')]:
        frame=apply_model(apply_labels(pd.read_csv(HERE/'data'/filename),name),name);bundle=evaluate(frame)
        stats[name]=bundle[0];bundles[name]=bundle
        bundle[1].to_csv(HERE/'data'/f'classification_{filename}',index=False)
    result=dict(protocol=dict(experimental_filter='Strict y > +0.5 or y < -0.5; exclude both boundaries',
        label_policy=POLICY,
        positive_class='stabilizing, experimental ddG > +0.5 kcal/mol',
        predicted_positive='continuous model prediction > 0; zero maps to destabilizing',
        prediction_threshold_fitted=False,bootstrap='1000 paired protein-group draws, seed42',
        zero_division='Precision is set to 0 for resamples with no positive predictions',
        dataset_note='Janus Fig6 criterion on S669, curated S461, S96 and PTmul-NR; S461 and S96 overlap S669 and are not fully independent',
        verification='ROC trapezoidal area independently checked against Mann-Whitney rank AUC; perfect/reversed/tied fixtures pass'),
        datasets=stats)
    (HERE/'generated/classification_statistics.json').write_text(json.dumps(result,indent=2))
    classification_table(stats)
    flat=[]
    for name,r in stats.items():
        for model,m in r['metrics'].items():flat.append(dict(dataset=name,model=model,**m))
    pd.DataFrame(flat).to_csv(HERE/'generated/classification_metrics.csv',index=False)
    make_figure(bundles,stats)
    print(pd.DataFrame(flat)[['dataset','model','n','positive','precision','recall','balanced_accuracy','roc_auc','mcc']].to_string(index=False))
    print('Intervals:',json.dumps({n:s['ensemble_minus_component']['JanusDDG'] for n,s in stats.items()},indent=2))


if __name__=='__main__':main()
