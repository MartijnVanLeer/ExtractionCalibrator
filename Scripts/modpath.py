import os 
from functions import ModpathFuncs as mf
import xarray as xr
import pandas as pd
import nlmod
from shutil import copytree
modelname = snakemake.params.modelname
layer = snakemake.params.layer
Name = snakemake.params.name

orgFolder = os.path.join('..','Results',f'{modelname}', f'{modelname}_ss', 'Fitter')
destFolder = os.path.join('..','Results',f'{modelname}', f'{modelname}_ss', 'ModpathRuns')
copytree(orgFolder ,destFolder, dirs_exist_ok = True)

#Load model and ds
ds = xr.open_dataset(os.path.join('..','Results',f'{modelname}', f'{modelname}_ss', f'{modelname}_ss.nc'))
ds.attrs['model_ws'] = destFolder
npfk, npfk33 = mf.load_calibrated_npf(modelname)
sim, npf = mf.load_ss(destFolder,ds, npfk, npfk33)

#get well nodes
ExWells = pd.read_csv(os.path.join('..','Data','dawaco',f'winputten_WG_{Name}.csv'))
wellxy = list(zip(ExWells['x_coordinaat'], ExWells['y_coordinaat']))

#run ref model
flowfrac_ref = mf.run_modpath_ref_bw(modelname, sim, ds, npf, layer, wellxy)
dist_ref = mf.run_modpath_ref_fw(modelname, sim, ds, npf, layer, wellxy)

#load realizations
rds = xr.open_dataset(os.path.join('..', 'Results', modelname, 'BestRealizations.nc'))

#run modpath for realizations
flowfrac, dist = mf.run_modpath_realizations(modelname,sim,ds,npf, rds, layer, wellxy)
flowfrac  = pd.DataFrame({'Flowfrac' : flowfrac, 'Realization' : 'Realizations'})
flowfrac = pd.concat([flowfrac, pd.DataFrame({'Flowfrac' : flowfrac_ref, 'Realization' : 'Reference'}, index = [0])], ignore_index=True)
flowfrac.to_csv(os.path.join('..','Results',f'{modelname}','flowfrac.csv'))

dist  = pd.DataFrame({'dist' : dist, 'Realization' : 'Realizations'})
dist =pd.concat([dist, pd.DataFrame({'dist' : dist_ref, 'Realization' : 'Reference'}, index = list(range(len(dist_ref))))], ignore_index=True)
dist.to_csv(os.path.join('..','Results',f'{modelname}','TT_dist.csv'))

