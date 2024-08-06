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
modelname = snakemake.params.modelname
layer = snakemake.params.simlayer
ens_no = snakemake.params.ens_no
real_dx = snakemake.params.dx
ws = snakemake.params.ws
cc = snakemake.params.cc

#load k realizations and cellids
df = pd.read_hdf(filename, key = 'c')
cellids = pd.read_csv(os.path.join('..','Results',f'{modelname}','cellids.csv'))


#load model ds
mds = xr.open_dataset(os.path.join('..','Results',f'{modelname}', f'{modelname}_t',f'{modelname}_t.nc'))

#init result xarray
ids = mds.icell2d.values
result = xr.Dataset(data_vars=dict( k = (['sim', 'icell2d','cc'], np.zeros((ens_no, len(ids), len(cc))))), coords =  dict(sim = range(ens_no), icell2d = ids, cc = cc))

#add cellids to df
df.set_index(['x', 'y', 'z'], inplace = True)
cellids.set_index(['x','y'], inplace = True)
df = df.merge(cellids, left_index=True, right_index=True, how = 'outer')
print (df.head(5))

#run modflow for modelcell for all realizations
for cellid in tqdm(ids):
    cell = df[df.cellid == cellid]
    cellk = xr.Dataset.from_dataframe(cell)
    for col in df.columns[1:]:
            sim = float(col.split('_')[1]) -1
            corfac = float(col.split('_')[2])
            k = cellk[col].values
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
            result.loc[dict(icell2d = cellid, sim = sim, cc = corfac)] = k.shape[2] / fieldK
            result['k'] = result.k.interpolate_na(dim = 'sim', method = 'nearest')


result.to_netcdf(snakemake.output[0])





