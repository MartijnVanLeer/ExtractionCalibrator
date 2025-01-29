import xarray as xr
import os

ds = xr.open_dataset(os.path.join('..', 'Results', 'Vlen', 'OriginalBestRealizations.nc'))
ds.load()
print('read')
ds = ds.assign(k = ds['k'].rename({'dim_0' : 'index'}))
ds = ds.drop_vars('dim_0')
# ds.to_netcdf(os.path.join('..', 'Results', 'Budel Output','Vlen', 'SelectedOriginal.nc'))
print('dropped dim_0')
dsopt = ds.where(ds.cc == 1.5, drop = True)
print('cc')
dsopt = dsopt.where(dsopt.xcorlen.isin([700,800,1000,1200]), drop = True)
print('x')
dsopt = dsopt.where(dsopt.zcorlen == 7.5, drop = True)
print('z')
dsopt = dsopt.where(dsopt.frac.isin([-0.05, 0]), drop = True)

dsopt.to_netcdf(os.path.join('..', 'Results', 'Budel Output','Vlen', 'OptimalOriginal.nc'))