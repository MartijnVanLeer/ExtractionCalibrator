# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 14:08:46 2023

@author: leermdv
"""

import xarray as xr
import glob
import os
import pandas as pd
import numpy as np
import Upscale_funcs as uf
import nlmod
from nlmod.dims.grid import xyz_to_cid
from tqdm import tqdm
import flopy 

filename = snakemake.input[1]
model_name = snakemake.params.modelname
layer = snakemake.params.simlayer
ens_no = snakemake.params.ens_no
real_dx = snakemake.params.dx
ws = snakemake.params.ws

#load k realizations and move to ds
df = pd.read_csv(filename)

# xdf = xr.Dataset.from_dataframe(df)

#load model ds
mds = xr.open_dataset(os.path.join('..','Results',f'{model_name}', f'{model_name}_t',f'{model_name}_t.nc'))

orgFolder = os.path.join('..','Results',f'{model_name}', f'{model_name}_t','Fitter','')
sim = flopy.mf6.mfsimulation.MFSimulation.load('mfsim', sim_ws = orgFolder, exe_name = mds.exe_name)
gwf = sim.get_model()

#init result xarray
ids = mds.icell2d.values
result = xr.Dataset(data_vars=dict( k = (['sim', 'icell2d'], np.zeros((ens_no, len(ids))))), coords =  dict(sim = range(ens_no), icell2d = ids))
def add_cellid(Kfields,ds, layer):
    cellids = [] 
    for index, row in tqdm(Kfields.iterrows()):
        layer, cellid = gwf.modelgrid.intersect(row.x,row.y)
        cellids.append(cellid)
    Kfields['cellid'] = cellids
    return Kfields
print (df.iloc[35950:35953][['x','y','z']])
print (mds.extent)

df = add_cellid(df, mds, layer)

df.set_index(['x', 'y', 'z'], inplace = True)
test = df[df.cellid == 12]
print(np.count_nonzero(np.isnan(ds)))
testk= xr.Dataset.from_dataframe(test)
testk.to_netcdf('test.nc')

df.set_index(['x', 'y', 'z'], inplace = True)
#run modflow for modelcell for all realizations
for cellid in tqdm(ids):
    cell = df[df.cellid == cellid]
    cellk = xr.Dataset.from_dataframe(cell)
    clean  = cellk #.dropna('x', how = 'any').dropna('y', how = 'any').dropna('z', how = 'any')
    for sim in range(ens_no):
        k = clean[f"K_{sim+1}"].values
        if np.isnan(k).any():
            print(f'Nans spotted, cellid = {cellid}')
            raise Exception(f'Nan gevonden in cellid {cellid}')
        if (k.shape[0]) == 1 or (k.shape[1] == 1):
            result.loc[dict(icell2d = cellid, sim = sim)] = 10**(np.mean(k))
        else:
            fieldK = uf.Run_MF_WholeField(10**(k),
                            Lx = k.shape[0] *real_dx,
                            Ly = k.shape[1] *real_dx,
                            Lz = k.shape[2],
                            dx = real_dx,dy = real_dx,dz = 1, mds = mds, ws = ws)
            result.loc[dict(icell2d = cellid, sim = sim)] = fieldK
            


result.to_netcdf(snakemake.output[0])





