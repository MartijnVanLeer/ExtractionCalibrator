import pandas as pd
import os
import xarray as xr 
import numpy as np

modelname = snakemake.params.modelname

results = pd.read_csv(os.path.join('..', 'Results', modelname, 'RMSE_all.csv'), names = ['sim', 'RMSE', 'xcorlen', 'zcorlen','frac', 'cc'])
ResidualsBest = pd.read_csv(os.path.join('..', 'Results', modelname,f'Residuals_{modelname}.csv'), index_col = "Time")
residuals =ResidualsBest.to_numpy()
residuals = residuals[~np.isnan(residuals)]
RMSE = np.sqrt(np.mean(residuals**2))
Best = results[results.RMSE < RMSE]


Best['k'] = np.nan

realizations = xr.Dataset.from_dataframe(Best.reset_index())
print(realizations)
for index, row in Best.iterrows():
    TempDS = xr.open_dataset(os.path.join('..', 'Results', modelname, 'KfieldsQC',f'xcorlens~{int(row.xcorlen)}', f'zcorlens~{int(row.zcorlen)}', f'fracs~{row.frac}', 'UpscaledK.nc'))
    Vals = TempDS.sel(sim = row.sim, cc = row.cc)
    if 'icell2d' not in realizations.dims.values():
        realizations = realizations.expand_dims({'icell2d' : TempDS.icell2d.values})
    realizations.loc[index, :]['k'] = Vals.k.values

realizations.to_netcdf(os.path.join('..', 'Results', modelname, 'BestRealizations.nc'))

