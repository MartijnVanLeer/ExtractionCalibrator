import xarray as xr
import flopy 
import OptimisationFuncs
import pandas as pd
import os
import nlmod
import numpy as np
from tqdm import tqdm

model_name = snakemake.params.modelname
ds = xr.open_dataset(snakemake.input[1])
Layer = snakemake.params.simlayer
destFolder = snakemake.params.ws

idx = pd.read_csv(os.path.join('..','Results',f'{model_name}',f'idx_SS_{model_name}.csv'))
ObsHeads =pd.read_csv(os.path.join('..','Results',f'{model_name}',f'ObsHead_{model_name}.csv'), index_col = 'Time')
ObsHeads.index = pd.to_datetime(ObsHeads.index)
ObsWells = pd.read_csv(os.path.join('..','Results',f'{model_name}',f'ObsForCalibration_{model_name}_SS.csv'))
layno = idx[idx.SensLayers == Layer].idx.values[0]

mds = xr.open_dataset(os.path.join('..','Results',f'{model_name}', f'{model_name}_t',f'{model_name}_t.nc'))

orgFolder = os.path.join('..','Results',f'{model_name}', f'{model_name}_t','Fitter','')
OptimisationFuncs.copyOriginalFolder(model_name + '_t', orgFolder ,destFolder , 'Runner\\' )
mds.attrs['model_ws'] = destFolder


sim = flopy.mf6.mfsimulation.MFSimulation.load('mfsim', sim_ws = destFolder, exe_name = mds.exe_name)

gwf = sim.get_model()
npf = gwf.get_package('NPF')


RMSE = []
for simno in tqdm(ds.sim.values):
    data33 = npf.k33.data
    data33[layno] = ds.sel(sim = simno).k.values
    npf.k33.set_data(data33)
    data = npf.k.data
    data[layno] = ds.sel(sim = simno).k.values
    npf.k.set_data(data)
    npf.write()
    sim.run_simulation(silent = True)
    head = nlmod.gwf.get_heads_da(fname = os.path.join(destFolder,'Upscaler', 'Upscaler.hds'))
    df = pd.DataFrame(index = pd.DatetimeIndex(ObsHeads.index))
    for index, well in ObsWells.iterrows():
        modheads = head.isel(layer = int(well.Layno)).sel(icell2d = int(well.CellID)).sel(time = slice(pd.to_datetime(ObsHeads.index[0]),pd.to_datetime(ObsHeads.index[-1]) ))
        df[f'{well["putcode"]}'] = modheads.values

    residuals = df - ObsHeads

    residuals = residuals.to_numpy().flatten()
    residuals = residuals[~np.isnan(residuals)]
    residuals = sum(residuals**2)
    print(residuals)
    RMSE.append(np.sqrt(residuals))

RMSEdf = pd.DataFrame({'sim' : ds.sim.values, 'RMSE' : RMSE})
RMSEdf.to_csv(snakemake.output[0])







