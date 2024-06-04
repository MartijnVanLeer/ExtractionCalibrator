import pandas as pd
import os
import xarray as xr 
import numpy as np

modelname = snakemake.params.modelname

#read csv
results = pd.read_csv(os.path.join('..', 'Results', modelname, 'RMSE_all.csv'), names = ['sim', 'RMSE', 'xcorlen', 'zcorlen','frac', 'cc'])
ResidualsBest = pd.read_csv(os.path.join('..', 'Results', modelname,f'Residuals_{modelname}.csv'), index_col = "Time")
#calc RMSE
residuals =ResidualsBest.to_numpy()
residuals = residuals[~np.isnan(residuals)]
RMSE = np.sqrt(np.mean(residuals**2))
#select best and make xr ds
Best = results[results.RMSE < RMSE]
realizations = xr.Dataset.from_dataframe(Best.reset_index())

#add cellid as dim
for index, row in Best.iterrows():
    TempDS = xr.open_dataset(os.path.join('..', 'Results', modelname, 'KfieldsQC',f'xcorlens~{int(Best.xcorlen.values[0])}', f'zcorlens~{int(Best.zcorlen.values[0])}', f'fracs~{Best.frac.values[0]}', 'UpscaledK.nc'))
realizations = realizations.expand_dims({'icell2d' : TempDS.icell2d.values})

#add empty dataarray k
realizations['k'] = xr.DataArray(coords = (realizations.index, realizations.icell2d))

#read k values from other nc files and move to realizations ds
for index, row in Best.iterrows():
    TempDS = xr.open_dataset(os.path.join('..', 'Results', modelname, 'KfieldsQC',f'xcorlens~{int(row.xcorlen)}', f'zcorlens~{int(row.zcorlen)}', f'fracs~{row.frac}', 'UpscaledK.nc'))
    Vals = TempDS.sel(sim = row.sim, cc = row.cc)
    realizations['k'].loc[index] = Vals.k.values

realizations.to_netcdf(os.path.join('..', 'Results', modelname, 'BestRealizations.nc'))

