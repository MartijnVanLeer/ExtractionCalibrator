import pandas as pd
import os
import xarray as xr 
import numpy as np
from tqdm import tqdm
modelname = 'Vlen'
results = pd.read_csv(os.path.join('..', 'Results', 'Vlen', 'RMSE_all.csv'), names = ['sim', 'RMSE', 'KGE', 'alpha', 'beta', 'r','xcorlen', 'zcorlen','frac', 'cc'])
cal_results = pd.read_csv(os.path.join('..', 'Results', modelname, f'Calibration_Performance_{modelname}.csv'))
RMSE = cal_results['RMSE'].values[-1]
KGE = cal_results['KGE'].values[-7]
Best = results.loc[(results.RMSE < RMSE) | (results.KGE > KGE)]
Best.reset_index(inplace = True, drop = True)
Best.index.name = 'index'
realizations = xr.Dataset.from_dataframe(Best)
row = Best.iloc[0]
df = pd.read_hdf(os.path.join('..', 'Results', modelname, 'KfieldsQC',f'xcorlens~{int(row.xcorlen)}', f'zcorlens~{row.zcorlen}', f'fracs~{row.frac}', 'k.h5'), key = 'c')
realizations = realizations.expand_dims({'x' : df.x.unique(), 'y' : df.y.unique()})
realizations['k'] = xr.DataArray(coords = (Best.index.values, realizations.x, realizations.y))
def harmonic_mean_func(values, dim):
    return values.count(dim=dim) / (1 / values).sum(dim=dim)

def harmonic_mean(group):
    return len(group) / (1 / group).sum()

for index, row in tqdm(Best.iterrows(), total=Best.shape[0]):
    df = pd.read_hdf(os.path.join('..', 'Results', modelname, 'KfieldsQC',f'xcorlens~{int(row.xcorlen)}', f'zcorlens~{row.zcorlen}', f'fracs~{row.frac}', 'k.h5'), key = 'c')
    Vals = df[['x', 'y', 'z']]
    Vals.index.name = 'index'
    Vals.loc[:,'k'] =  10**df.loc[:,f'K_{int(row.sim+1)}_{row.cc}']
    Vals = Vals.set_index(['x', 'y','z'])
    harmonic_mean_df = Vals.groupby(['x', 'y'])['k'].apply(lambda group: harmonic_mean(group)).reset_index()
    harmonic_mean_df = harmonic_mean_df.set_index(['x', 'y'], drop = True)
    harmonic_mean_xr = xr.Dataset.from_dataframe(harmonic_mean_df)
    realizations['k'].loc[index] = harmonic_mean_xr.k


realizations.to_netcdf(os.path.join('..', 'Results', modelname, 'OriginalBestRealizations.nc'))