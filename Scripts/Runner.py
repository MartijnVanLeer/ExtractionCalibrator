import xarray as xr 
import flopy 
import OptimisationFuncs
import pandas as pd
import os
import nlmod
import numpy as np

model_name = snakemake.params.modelname
ds = xr.open_dataset(snakemake.input[1])
Layer = snakemake.params.simlayer

idx = pd.read_csv(os.path.join('..','Results',f'{model_name}',f'idx_SS_{model_name}.csv'))
ObsHeads =pd.read_csv(os.path.join('..','Results',f'{model_name}',f'ObsHead_{model_name}.csv'))
ObsWells = pd.read_csv(os.path.join('..','Results',f'{model_name}',f'ObsForCalibration_{model_name}_SS.csv'))
layno = idx[idx.SensLayers == Layer].idx.values

mds = xr.open_dataset(os.path.join('..','Results',f'{model_name}', f'{model_name}_t',f'{model_name}_t.nc'))

orgFolder = os.path.join('..','Results',f'{model_name}', f'{model_name}_t','Fitter','')
destFolder = os.path.join('..','Results',f'{model_name}', f'{model_name}_t', 'Runner','')
OptimisationFuncs.copyOriginalFolder(model_name + '_t', orgFolder ,destFolder , 'Runner\\' )
ds.attrs['model_ws'] = destFolder


sim = flopy.mf6.mfsimulation.MFSimulation.load('mfsim', sim_ws = destFolder, exe_name = mds.exe_name)
gwf = sim.get_model()
npf = gwf.get_package('NPF')


RMSE = []
for simno in ds.sim.values:
    npf.k[layno] = ds.sel(sim = simno).k.values
    npf.write()
    sim.run_simulation(silent = True)
    head = nlmod.gwf.get_heads_da(ds)
    df = pd.DataFrame(index = pd.DatetimeIndex(ObsHeads.index))
    for index, well in ObsWells.iterrows():
        modheads = head.isel(layer = int(well.Layno)).sel(icell2d = int(well.CellID)).sel(time = slice(pd.to_datetime(ObsHeads.index[0]),pd.to_datetime(ObsHeads.index[-1]) ))
        df[f'{well["putcode"]}'] = modheads.values

    residuals = df - ObsHeads
    residuals = residuals.to_numpy().flatten()
    residuals = residuals[~np.isnan(residuals)]
    residuals = sum(residuals**2)
    RMSE.append(np.sqrt(residuals))

RMSEdf = pd.DataFrame({'sim' : ds.sim.values, 'RMSE' : RMSE})
RMSEdf.to_csv(snakemake.output[0])







