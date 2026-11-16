"""Independent standard-library/numpy/pandas OLS-product reference, not CAD-SPOT."""
import numpy as np
import pandas as pd


def fit_residuals(table,covariates):
    numeric=table[covariates].to_numpy(float)
    sd=numeric.std(axis=0);sd[sd==0]=1
    numeric=(numeric-numeric.mean(axis=0))/sd
    sections=pd.get_dummies(table.sample_id,drop_first=True,dtype=float).to_numpy()
    x=np.column_stack([np.ones(len(table)),numeric,sections])
    y=table[['target_score','reference_score']].to_numpy(float)
    return y-x@np.linalg.lstsq(x,y,rcond=None)[0]


def standardize_sections(table,values):
    result=np.empty_like(values,dtype=float)
    for indices in table.groupby('sample_id',sort=True).indices.values():
        subset=values[indices];sd=subset.std(axis=0);sd=np.where(sd>0,sd,1.)
        result[indices]=(subset-subset.mean(axis=0))/sd
    return result


def rank_reference(table,covariates):
    residuals=fit_residuals(table,covariates)
    z=standardize_sections(table,residuals);priority=np.maximum(z[:,0],0)*np.maximum(z[:,1],0)
    chosen=np.zeros(len(table),dtype=bool)
    valid=np.isfinite(table[['x','y']].to_numpy(float)).all(axis=1)
    for indices in table.groupby('sample_id',sort=True).indices.values():
        indices=indices[valid[indices]]
        if len(indices)<7:continue
        p=priority[indices];center=np.median(p);scale=1.4826*np.median(np.abs(p-center))
        if scale<=0:scale=p.std() or 1.
        chosen[indices]=(z[indices,0]>0)&(z[indices,1]>0)&((p-center)/scale>=3)
    return dict(target_resid=residuals[:,0],reference_resid=residuals[:,1],
        target_resid_z=z[:,0],reference_resid_z=z[:,1],priority_score=priority,cadspot_candidate=chosen)
