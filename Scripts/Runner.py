import xarray as xr
import flopy 
import OptimisationFuncs
import pandas as pd
import os
import nlmod
import numpy as np
from tqdm import tqdm
import shutil 
modelname = snakemake.params.modelname
ds = xr.open_dataset(snakemake.input[1])
Layer = snakemake.params.simlayer
sim = snakemake.params.simulation
cc = snakemake.params.cc
xcorlens = sim['xcorlens']
zcorlens = sim['zcorlens']
fracs = sim['fracs']

#destFolder = snakemake.params.ws
destFolder = os.path.join(os.path.dirname(snakemake.input[1]), 'Runner')

idx = pd.read_csv(os.path.join('..','Results',f'{modelname}',f'idx_SS_{modelname}.csv'))
ObsHeads =pd.read_csv(os.path.join('..','Results',f'{modelname}',f'ObsHead_{modelname}.csv'), index_col = 'Time')
ObsHeads.index = pd.to_datetime(ObsHeads.index)
ObsWells = pd.read_csv(os.path.join('..','Results',f'{modelname}',f'ObsForCalibration_{modelname}_SS.csv'))
layno = idx[idx.SensLayers == Layer].idx.values[0]

#shutil.copyfile(os.path.join('..','Results',f'{modelname}', f'{modelname}_t',f'{modelname}_t.nc'), os.path.join(destFolder, 'ds.nc'))
#mds = xr.open_dataset(os.path.join(destFolder, 'ds.nc'))
mds = xr.open_dataset(os.path.join('..','Results',f'{modelname}', f'{modelname}_t',f'{modelname}_t.nc'))

orgFolder = os.path.join('..','Results',f'{modelname}', f'{modelname}_t','Fitter','')
OptimisationFuncs.copyOriginalFolder(modelname + '_t', orgFolder ,destFolder , 'Runner\\' )
mds.attrs['model_ws'] = destFolder


sim = flopy.mf6.mfsimulation.MFSimulation.load('mfsim', sim_ws = destFolder, exe_name = mds.exe_name)
gwf = sim.get_model()
npf = gwf.get_package('NPF')

thickness = mds.isel(layer = layno-1).botm-mds.sel(layer = Layer).botm

RMSE = []
cc_ls = []
simno_ls = []
for simno in tqdm(ds.sim.values):
    for corfac in ds.cc.values: 
        data33 = npf.k33.array
        data33[layno] = thickness.values / ds.sel(sim = simno, cc = corfac).k.values 
        npf.k33.set_data(data33)
        data = npf.k.array
        data[layno] = thickness.values / ds.sel(sim = simno, cc = corfac).k.values
        npf.k.set_data(data)
        npf.write()
        success, buff = sim.run_simulation(silent = True)
        df = pd.DataFrame(index = pd.DatetimeIndex(ObsHeads.index))
        if success:
            head = nlmod.gwf.get_heads_da(mds)
            for index, well in ObsWells.iterrows():
                modheads = head.isel(layer = int(well.Layno)).sel(icell2d = int(well.CellID)).sel(time = slice(pd.to_datetime(ObsHeads.index[0]),pd.to_datetime(ObsHeads.index[-1]) ))
                df[f'{well["putcode"]}'] = modheads.values

            residuals = df - ObsHeads
            residuals = residuals.to_numpy().flatten()
            residuals = residuals[~np.isnan(residuals)]
            residuals = np.mean(residuals**2)
            RMSE.append(np.sqrt(residuals))
        else:
            RMSE.append(np.nan)
            print('Modflow crashed')
        cc_ls.append(corfac)
        simno_ls.append(simno)


RMSEdf = pd.DataFrame({'sim' : simno_ls, 'RMSE' : RMSE, 'xcorlen' : xcorlens, 'zcorlen' : zcorlens, 'frac' : fracs, 'cc' : cc_ls})
RMSEdf.to_csv(snakemake.output[0], header = False)







