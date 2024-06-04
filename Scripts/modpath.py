import os 
import nlmod
import flopy 
import OptimisationFuncs
import ModpathFuncs as mf
import xarray as xr
import pandas as pd

modelname = snakemake.params.modelname
orgFolder = os.path.join('..','Results',f'{modelname}', f'{modelname}_ss', 'Fitter')
destFolder = os.path.join('..','Results',f'{modelname}', f'{modelname}_ss', 'ModpathRuns')
OptimisationFuncs.copyOriginalFolder(modelname + '_ss', orgFolder ,destFolder , 'Fitter\\' )

#Load model and ds
ds = xr.open_dataset(os.path.join('..','Results',f'{modelname}', f'{modelname}_ss', f'{modelname}_ss.nc'))
ds.attrs['model_ws'] = destFolder

npfk, npfk33 = mf.load_calibrated_npf(modelname)
sim, npf = mf.load_ss(destFolder,ds, npfk, npfk33)

#run ref model
flowfrac_ref = mf.run_modpath_ref_bw(modelname, sim, ds, npf, layer)
dist_ref = mf.run_modpath_ref_fw(modelname, sim, ds, npf, layer)

#load realizations
rds = xr.open_dataset(os.path.join('..', 'Results', modelname, 'BestRealizations.nc'))

#run modpath for realizations
flowfrac, dist = mf.run_modpath_realizations(modelname,sim,ds,npf, rds, layer)
flowfrac  = pd.DataFrame({'Flowfrac' : flowfrac, 'Realization' : 'Realizations'})
flowfrac.append({'Flowfrac' : flowfrac_ref, 'Realization' : 'Reference'}, ignore_index=True)
flowfrac.to_csv(os.path.join('..','Results',f'{modelname}','flowfrac.csv'))

dist  = pd.DataFrame({'dist' : dist, 'Realization' : 'Realizations'})
dist.append({'Flowfrac' : dist_ref, 'Realization' : 'Reference'}, ignore_index=True)
dist.to_csv(os.path.join('..','Results',f'{modelname}','TT_dist.csv'))

