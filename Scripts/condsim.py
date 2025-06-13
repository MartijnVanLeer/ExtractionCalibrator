# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 09:53:38 2023

@author: leermdv
"""


#%%

from functions import SISIM_R
import os
import pickle
import shutil
import xarray as xr
from functions import Heterogeniteit
import numpy as np
import flopy
import pandas as pd

def read_snakemake_rule(path, rule: str) -> "snakemake.rules.Rule": # type: ignore
    import snakemake as sm
    workflow = sm.Workflow(snakefile="snakefile")
    workflow.include(path)
    return workflow.get_rule(rule)

if "snakemake" not in globals():
    snakemake = read_snakemake_rule('snakefile','condsim')
    
Location = snakemake.params.Name
modelname = snakemake.params.modelname
Layer = snakemake.params.simlayer
sim = snakemake.params.simulation
cc = np.array(snakemake.params.cc)

xcorlens = sim['xcorlens']
zcorlens = sim['zcorlens']
fracedit = sim['fracs']
ens_no = snakemake.params.ens_no
dx = dy = snakemake.params.dx

with open(os.path.join('..','Results',f'{modelname}','boreholeindicators.pkl'), 'rb') as f:
    boringen = pickle.load(f)
frac = boringen.list.i[boringen.list.i > 0.5].count()/len(boringen.list)
print(f'Frac = {frac}')
fracs =  frac + fracedit
ds = xr.open_dataset(os.path.join('..','Results',f'{modelname}', f'{modelname}_t',f'{modelname}_t.nc'))

xmin = ds.extent[0]
ymin = ds.extent[2]
Lx = ds.extent[1] - ds.extent[0]
Ly = ds.extent[3] - ds.extent[2]


Lz = boringen.list.z.max() - boringen.list.z.min()
zmin = boringen.list.z.min()

res = SISIM_R.Cond_SISIM(boringen.list[['x','y','z','i']],
            xmin = xmin,ymin = ymin,zmin = zmin,
            Lx=Lx,  Ly=Ly,  Lz =Lz,
            dx =dx,dy =dy,dz =1,
            xcorlen =xcorlens, zcorlen = zcorlens,
            ens_no = ens_no, frac =fracs, 
            nmax = 50, seed = xcorlens*zcorlens*frac)

orgFolder = os.path.join('..','Results',f'{modelname}', f'{modelname}_t')
sim = flopy.mf6.mfsimulation.MFSimulation.load('mfsim', sim_ws = orgFolder, exe_name = ds.exe_name)
gwf = sim.get_model()

Kfields = boringen.add_k(res, ens_no,cc)
Kfields.to_hdf(snakemake.output[0], key = 'c', complevel = 9, mode = 'w')
