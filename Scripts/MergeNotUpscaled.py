import pandas as pd
import os
import xarray as xr 
import numpy as np
modelname = 'Vlen'
results = pd.read_csv(os.path.join('..', 'Results', 'Vlen', 'RMSE_all.csv'), names = ['sim', 'RMSE', 'KGE', 'alpha', 'beta', 'r','xcorlen', 'zcorlen','frac', 'cc'])
cal_results = pd.read_csv(os.path.join('..', 'Results', modelname, f'Calibration_Performance_{modelname}.csv'))
RMSE = cal_results['RMSE'].values[-1]
KGE = cal_results['KGE'].values[-7]
Best = results.loc[(results.RMSE < RMSE) | (results.KGE > KGE)]
Best.reset_index(inplace = True)

realizations = xr.Dataset.from_dataframe(Best)
row = Best[0]
df = pd.read_hdf(os.path.join('..', 'Results', modelname, 'KfieldsQC',f'xcorlens~{int(row.xcorlen)}', f'zcorlens~{row.zcorlen}', f'fracs~{row.frac}', 'k.h5'), key = 'c')
realizations = realizations.expand_dims({'x' : df.x.unique(), 'y' : df.y.unique()})
# realizations['k'] = xr.DataArray(coords = (realizations.index, realizations.x, realizations.y, realizations.z))
def harmonic_mean_func(values, dim):
    return values.count(dim=dim) / (1 / values).sum(dim=dim)

for index, row in Best.iterrows():
    df = pd.read_hdf(os.path.join('..', 'Results', modelname, 'KfieldsQC',f'xcorlens~{int(row.xcorlen)}', f'zcorlens~{row.zcorlen}', f'fracs~{row.frac}', 'k.h5'), key = 'c')
    Vals = df[['x', 'y', 'z']]
    Vals['k'] =  10**df[f'K_{row.sim+1}_{row.cc}']
    ds = Vals.set_index(['index', 'x', 'y', 'z']).to_xarray()
    harmonic_mean = xr.apply_ufunc(
        harmonic_mean_func, ds['k'],
        input_core_dims=[['z']],
        kwargs={'dim': 'z'})
    realizations['k'] = harmonic_mean


realizations.to_netcdf(os.path.join('..', 'Results', modelname, 'OriginalBestRealizations.nc'))