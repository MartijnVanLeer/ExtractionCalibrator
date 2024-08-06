
import numpy as np
import flopy 
import xarray as xr
import os
import pandas as pd
from tqdm import tqdm

modelname = snakemake.params.modelname
dx = dy = snakemake.params.dx

ds = xr.open_dataset(os.path.join('..','Results',f'{modelname}', f'{modelname}_t',f'{modelname}_t.nc'))
orgFolder = os.path.join('..','Results',f'{modelname}', f'{modelname}_t')
sim = flopy.mf6.mfsimulation.MFSimulation.load('mfsim', sim_ws = orgFolder, exe_name = ds.exe_name)
gwf = sim.get_model()

xmin = ds.extent[0] + 0.5 * dx
xmax = ds.extent[1] + 0.5 * dx
ymin = ds.extent[2] + 0.5 * dx
ymax = ds.extent[3] + 0.5 * dx

x = np.arange(xmin, xmax, dx)
y = np.arange(ymin, ymax, dx)

xnew = []
ynew = []
for xval in x: 
    for yval in y: 
        xnew.append(xval)
        ynew.append(yval)

def add_cellid(Kfields,gwf):
    cellids = [] 
    for index, row in tqdm(Kfields.iterrows(),'Intersecting grid..'):
        cellid = gwf.modelgrid.intersect(row.x,row.y)
        cellids.append(cellid)
    Kfields['cellid'] = cellids
    return Kfields

df = pd.DataFrame({'x' :xnew, 'y' : ynew})
cellids = [] 
for index, row in tqdm(df.iterrows(),'Intersecting grid..', miniters = 100, mininterval = 20):
    cellid = gwf.modelgrid.intersect(row.x,row.y)
    cellids.append(cellid)
df['cellid'] = cellids
df.to_csv(os.path.join('..','Results',f'{modelname}','cellids.csv'))





