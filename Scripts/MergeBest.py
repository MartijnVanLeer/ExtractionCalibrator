#Merge well fitting realizations and results into netCDF
import pandas as pd
import os
import xarray as xr 
import numpy as np

modelname = snakemake.params.modelname

#read csv
results = pd.read_csv(os.path.join('..', 'Results', modelname, 'RMSE_all.csv'), names = ['sim', 'RMSE', 'KGE', 'alpha', 'beta', 'r','xcorlen', 'zcorlen','frac', 'cc'])
ResidualsBest = pd.read_csv(os.path.join('..', 'Results', modelname,f'Residuals_{modelname}.csv'), index_col = "Time")
cal_results = pd.read_csv(os.path.join('..', 'Results', modelname, f'Calibration_Performance_{modelname}.csv'))
RMSE = cal_results['RMSE'].values[-1]
KGE = cal_results['KGE'].values[-7]

#select best and make xr ds
Best = results.loc[(results.RMSE < RMSE) | (results.KGE > KGE)]
if len(Best) == 0:
    print('No improved realizations')
    print(f'min RMSE = {results.RMSE.min()}')
    print(f'max KGE = {results.KGE.max()}')
    Best = results.nsmallest(500, 'RMSE')

realizations = xr.Dataset.from_dataframe(Best.reset_index())

#add cellid as dim
for index, row in Best.iterrows():
    TempDS = xr.open_dataset(os.path.join('..', 'Results', modelname, 'KfieldsQC',f'xcorlens~{int(Best.xcorlen.values[0])}', f'zcorlens~{Best.zcorlen.values[0]}', f'fracs~{Best.frac.values[0]}', 'UpscaledK.nc'))
realizations = realizations.expand_dims({'icell2d' : TempDS.icell2d.values})

#add empty dataarray k
realizations['k'] = xr.DataArray(coords = (realizations.index, realizations.icell2d))

#read k values from other nc files and move to realizations ds
for index, row in Best.iterrows():
    TempDS = xr.open_dataset(os.path.join('..', 'Results', modelname, 'KfieldsQC',f'xcorlens~{int(row.xcorlen)}', f'zcorlens~{row.zcorlen}', f'fracs~{row.frac}', 'UpscaledK.nc'))
    Vals = TempDS.sel(sim = row.sim, cc = row.cc)
    realizations['k'].loc[index] = Vals.k.values

realizations.to_netcdf(os.path.join('..', 'Results', modelname, 'BestRealizations.nc'))

