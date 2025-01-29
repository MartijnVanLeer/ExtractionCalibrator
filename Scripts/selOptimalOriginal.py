import xarray as xr
import os

ds = xr.open_dataset(os.path.join('..', 'Results', 'Vlen', 'OriginalBestrealizations.nc'))
print('read')
dsopt = ds.where(ds.cc == 1.5, drop = True)
print('cc')
dsopt = dsopt.where(dsopt.xcorlen.isin([700,800,1000,1200]), drop = True)
print('x')
dsopt = dsopt.where(dsopt.zcorlen == 7.5, drop = True)
print('z')
dsopt = dsopt.where(dsopt.frac.isin([-0.05, 0]), drop = True)


dsopt.to_netcdf(os.path.join('..', 'Results', 'Budel Output','Vlen', 'OptimalOriginal.nc'))