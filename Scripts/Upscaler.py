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

filename = snakemake.input[1]
model_name = snakemake.params.modelname
layer = snakemake.params.simlayer
ens_no = snakemake.params.ens_no
real_dx = snakemake.params.dx
ws = snakemake.params.ws
#load k realizations and move to ds
df = pd.read_csv(filename, index_col=['x','y','z'])
df.drop('Unnamed: 0', axis = 1)


#load model ds
mds = xr.open_dataset(os.path.join('..','Results',f'{model_name}', f'{model_name}_t',f'{model_name}_t.nc')).sel(layer = layer)
kds = xr.Dataset.from_dataframe(df)

ids = mds.icell2d.values
result = xr.Dataset(data_vars=dict( k = (['sim', 'icell2d'], np.zeros((ens_no, len(ids))))), coords =  dict(sim = range(ens_no), icell2d = ids))
for cellid in ids:
    cell = mds.sel(icell2d = cellid)
    dx = np.sqrt(cell.area.values)/2
    #cellk = kds.sel(x = slice(cell.x.values - dx, cell.x.values + dx),y = slice(cell.y.values - dx, cell.y.values + dx))
    cellk = kds.where(kds.cellid == cellid, drop = True,)
    cellk = cellk.dropna('z')
    print(cellk['K_sim1'].shape)
    for sim in range(ens_no):
        k = cellk[f"K_{sim+1}"].values
        fieldK = uf.Run_MF_WholeField(10**(k),
                                      Lx = k.shape[0] *real_dx,
                                      Ly = k.shape[1] *real_dx,
                                      Lz = k.shape[2],
                                      dx = real_dx,dy = real_dx,dz = 1, mds = mds, ws = ws)
        result.loc[dict(icell2d = cellid, sim = sim)] = fieldK

result.to_netcdf(snakemake.output[0])





