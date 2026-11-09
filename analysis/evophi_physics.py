"""Recompute same-site summaries for original, projected and zero-bias heads."""
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from evophi_policy import HERE,ROOT,physical_edges


def main():
    cache=HERE/'data/physics'
    frame=pd.read_csv(cache/'variants.csv');aa='ACDEFGHIKLMNPQRSTVWY'
    site=frame.site_index.to_numpy()
    a=np.array([aa.index(r.wt_seq[int(r.pos)]) for r in frame.itertuples()])
    b=np.array([aa.index(r.mut_seq[int(r.pos)]) for r in frame.itertuples()])
    j=np.load(cache/'janus_edges.npy').astype(float)
    original=np.load(cache/'original_msa_edges.npy').astype(float)
    anti=np.load(cache/'msa_edges.npy').astype(float);phi=physical_edges()
    models={'MSAesm original':original,'MSAesm_anti':anti,'EvoPhi-DDG':phi,'JanusDDG':j,'Mean':(j+phi)/2}
    paths=np.load(cache/'paths.npz');records=[]
    for model,edges in models.items():
        direct=edges[site,a,b]
        for n in [1,3,5,7,9]:
            for rep in range(10):
                chain=np.column_stack([a,paths[str(n)][rep],b])
                summed=edges[site[:,None],chain[:,:-1],chain[:,1:]].sum(1)
                error=summed-direct
                r=float(pearsonr(direct,summed).statistic)
                if model=='EvoPhi-DDG' and abs(1-r)<1e-14:r=1.
                records.append(dict(model=model,n_intermediates=n,repeat=rep,
                    pearson=r,mae=float(abs(error).mean()),rmse=float(np.sqrt(np.mean(error**2))),bias=float(error.mean())))
    repeat=pd.DataFrame(records);repeat.to_csv(HERE/'data/transitivity_repeats.csv',index=False)
    summary=repeat.groupby(['model','n_intermediates'],sort=False).agg(
        pearson_mean=('pearson','mean'),pearson_std=('pearson','std'),
        mae_mean=('mae','mean'),rmse_mean=('rmse','mean'),bias_mean=('bias','mean')).reset_index()
    summary.to_csv(HERE/'data/transitivity.csv',index=False)


if __name__=='__main__':main()
