"""Regenerate the EvoPhi-DDG head ablation from cached S2450 predictions."""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from scipy.stats import pearsonr,spearmanr
from analysis import HERE,ROOT
from figure_formatting import correlation,format_axes


def main():
    frame=pd.read_csv(HERE/'data/model_controls/S2450_predictions.csv')
    x=frame.original_AS.to_numpy();y=frame.zero_hidden_bias_AS.to_numpy();delta=y-x
    stats=dict(n=len(x),pearson=float(pearsonr(x,y).statistic),spearman=float(spearmanr(x,y).statistic),
               mae=float(abs(delta).mean()),max_abs=float(abs(delta).max()))
    fig,axes=plt.subplots(1,2,figsize=(8,3.7),layout='constrained')
    lo=min(x.min(),y.min())-.3;hi=max(x.max(),y.max())+.3
    ax=axes[0];ax.scatter(x,y,s=8,alpha=.4,color='#168274',rasterized=True)
    ax.plot([lo,hi],[lo,hi],'k--',lw=.8)
    ax.set(xlim=(lo,hi),ylim=(lo,hi),xlabel='MSAesm_ddG (AS) prediction (kcal/mol)',ylabel='EvoPhi-DDG prediction (kcal/mol)',title='a  Prediction agreement')
    ax.text(.04,.96,f'n = {len(x)}\nPearson r = {correlation(stats["pearson"])}\nSpearman ρ = {correlation(stats["spearman"])}',transform=ax.transAxes,va='top',fontsize=8)
    ax=axes[1];ax.scatter(x,delta,s=8,alpha=.4,color='#b86829',rasterized=True);ax.axhline(0,color='k',ls='--',lw=.8)
    ax.set(xlabel='MSAesm_ddG (AS) prediction (kcal/mol)',ylabel='EvoPhi-DDG − MSAesm_ddG (AS) (kcal/mol)',title='b  Effect of hidden-bias removal')
    ax.text(.04,.96,f'Mean |change| = {stats["mae"]:.2e}\nMax |change| = {stats["max_abs"]:.2e} kcal/mol',transform=ax.transAxes,va='top',fontsize=8)
    format_axes(fig,skip=[(axes[1],'y')]);axes[1].yaxis.set_major_formatter(FuncFormatter(lambda v,_:f'{v:.2e}'))
    for suffix in ['pdf','png']:fig.savefig(HERE/f'figures/S2450_bias_comparison.{suffix}',dpi=250,bbox_inches='tight')
    plt.close(fig)
    (HERE/'generated/bias_comparison_statistics.json').write_text(json.dumps(stats,indent=2)+'\n')

if __name__=='__main__':main()
