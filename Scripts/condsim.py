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

def read_snakemake_rule(path, rule: str) -> "snakemake.rules.Rule":
    import snakemake as sm
    workflow = sm.Workflow(snakefile="snakefile")
    workflow.include(path)
    return workflow.get_rule(rule)

if "snakemake" not in globals():
    snakemake = read_snakemake_rule('snakefile','condsim')
    
Location = snakemake.params.Name
modelname = snakemake.params.modelname
xcorlens = snakemake.params.simulation['xcorlens']
zcorlens = snakemake.params.simulation['zcorlens']
fracedit = snakemake.params.simulation['fracs']
ens_no = snakemake.params.ens_no
dx = dy = snakemake.params.dx

with open(os.path.join('..','Results',f'{modelname}','boreholeindicators.pkl'), 'rb') as f:
    boringen = pickle.load(os.path.join('..','Results',f'{modelname}','boreholeindicators.pkl'))
frac = boringen.list.i[boringen.list.i > 0.5].count()/len(boringen.list)
fracs =  frac + fracedit

xmin = boringen.layermodel.extent[0]
ymin = boringen.layermodel.extent[2]
Lx = boringen.layermodel.extent[1] - boringen.layermodel.extent[0]
Ly = boringen.layermodel.extent[3] - boringen.layermodel.extent[2]
zmin = boringen.list.z.min()



a,res = SISIM_R.Cond_SISIM(boringen.list[['x','y','z','i']],
            xmin = boringen.layermodel.extent[0],ymin = boringen.layermodel.extent[2],zmin = boringen.list.z.min(),
            Lx=boringen.layermodel.extent[1] - boringen.layermodel.extent[0],
            Ly=boringen.layermodel.extent[3] - boringen.layermodel.extent[2],
            Lz =abs(boringen.list.z.min()-boringen.list.z.max()),
            dx =dx,dy =dy,dz =1,
            xcorlen =xcorlens, zcorlen = zcorlens,
            ens_no = ens_no, frac =frac, nmax = 100, seed = 1337)

Kfields = boringen.add_k(res)
if not os.path.isdir(os.path.join('..' ,'Results',f'{modelname}','KfieldsQC')):
    os.mkdir(os.path.join('..' ,'Results',f'{modelname}','KfieldsQC'))
else:
    shutil.rmtree(os.path.join('..' ,'Results',f'{modelname}','KfieldsQC'))
    os.mkdir(os.path.join('..' ,'Results',f'{modelname}','KfieldsQC'))
res.to_csv(os.path.join('..' ,'Results',f'{modelname}','KfieldsQC',f'Ind_x{xcorlens}_z{zcorlens}_f{frac}_n{ens_no}.csv'))
Kfields.to_csv(os.path.join('..' ,'Results',f'{modelname}','KfieldsQC',f'K_x{xcorlens}_z{zcorlens}_f{frac}_n{ens_no}.csv'))
