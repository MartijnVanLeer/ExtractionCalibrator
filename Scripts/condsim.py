# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 09:53:38 2023

@author: leermdv
"""


#%%

import SISIM_R
import os
import pickle
import shutil
import xarray as xr
import Heterogeniteit
import numpy as np

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

xcorlens = sim['xcorlens']
zcorlens = sim['zcorlens']
fracedit = sim['fracs']
ens_no = snakemake.params.ens_no
dx = dy = snakemake.params.dx

with open(os.path.join('..','Results',f'{modelname}','boreholeindicators.pkl'), 'rb') as f:
    boringen = pickle.load(f)
frac = boringen.list.i[boringen.list.i > 0.5].count()/len(boringen.list)
fracs =  frac + fracedit


ds = xr.open_dataset(os.path.join('..','Results',f'{modelname}', f'{modelname}_t',f'{modelname}_t.nc'))

xmin = ds.extent[0] - 0.5*dx
ymin = ds.extent[2] - 0.5*dx
Lx = ds.extent[1] - ds.extent[0] - dx 
Ly = ds.extent[3] - ds.extent[2] - dx 


Lz = boringen.list.z.max() - boringen.list.z.min()
zmin = boringen.list.z.min()


a,res = SISIM_R.Cond_SISIM(boringen.list[['x','y','z','i']],
            xmin = xmin,ymin = ymin,zmin = zmin,
            Lx=Lx,  Ly=Ly,  Lz =Lz,
            dx =dx,dy =dy,dz =1,
            xcorlen =xcorlens, zcorlen = zcorlens,
            ens_no = ens_no, frac =frac, 
            nmax = 100, seed = xcorlens*zcorlens*frac)




res = Heterogeniteit.add_cellid(res,ds, Layer)
Kfields = boringen.add_k(res)

Kfields.to_csv(snakemake.output[0])
