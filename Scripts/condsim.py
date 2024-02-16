# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 09:53:38 2023

@author: leermdv
"""
import Heterogeniteit as h
import xarray as xr
import SISIM_R
from timeit import default_timer as timer
from os.path import join
import os
from tqdm import tqdm

def read_snakemake_rule(path, rule: str) -> "snakemake.rules.Rule":
    import snakemake as sm
    workflow = sm.Workflow(snakefile="snakefile")
    workflow.include(path)
    return workflow.get_rule(rule)

if "snakemake" not in globals():
    snakemake = read_snakemake_rule('snakefile','condsim')

Location = snakemake.params.Name
modelname = snakemake.params.modelname
Layer = snakemake.params.simlayer
TI = True if Location == 'Vlijmen' else False
borefolder = join("..",'Data',f'Boringen {Location}',r'Boormonsterprofiel_Geologisch booronderzoek')
ncpath =join('..','Results',f'{modelname}',f'{modelname}_t', 'layer_model.nc')
TIpath = join("..",'Data','AllDB.csv')
RegCat = join('..', 'Data','REGIS_catalogus_export18122015.csv' )



boringen = h.boringen(borefolder)
md = boringen.metadata
#voeg REGIS laag toe aan boringen
boringen.add_layers(xr.open_dataset(ncpath))
#selecteer trajecten binnen een laag en explode naar meterschaal
boringen.select_layer(Layer)
#Maak lijst met xyz en laag
boringen.listify()
cond = boringen.indicators(g1 = ['v', 'k', 'kz'])
#plot hoofdgrondsoort 
boringen.plot_lith_dist()

if TI:
    Kcore = h.KDist(TIpath, Layer)
    Kcore.plot_dist()
else:
    Kcore = h.Kreg(RegCat, Layer)

rn = boringen.plot_K_weighted(Kcore, TI)
#Trek K uit Kcore en voeg toe aan boringlijst
# cond.add_k(Kcore)

boringen.list.to_csv(join('..','Results',f'{modelname}_t','boreholeindicators.csv'))

#%%


xmin = boringen.layermodel.extent[0]
ymin = boringen.layermodel.extent[2]
Lx = boringen.layermodel.extent[1] - boringen.layermodel.extent[0]
Ly = boringen.layermodel.extent[3] - boringen.layermodel.extent[2]
zmin = cond.z.min()

xcorlens = [200,400,600,800]
zcorlens = [0.5,4,8]
frac = cond.i[cond.i > 0.5].count()/len(cond)
fracs = [frac,frac+0.1,frac +0.2]
ens_no = 1
dx = dy = 400
for xcorlen in tqdm(xcorlens, position = 1):
    for zcorlen in tqdm(zcorlens, position = 2):
        for frac in tqdm(fracs, position = 3):

            
            
            start = timer() 
            a,res = SISIM_R.Cond_SISIM(cond[['x','y','z','i']],
                        xmin = boringen.layermodel.extent[0],ymin = boringen.layermodel.extent[2],zmin = cond.z.min(),
                        Lx=boringen.layermodel.extent[1] - boringen.layermodel.extent[0],
                        Ly=boringen.layermodel.extent[3] - boringen.layermodel.extent[2],
                        Lz =abs(cond.z.min()-cond.z.max()),
                        dx =dx,dy =dy,dz =1,
                        xcorlen =xcorlen, zcorlen = zcorlen,
                        ens_no = ens_no, frac =frac, nmax = 100, seed = 1337)
            print(f'{round(timer() - start,1)} s for conditional simulation')
            
            Kfields = boringen.add_k(res)
            if not os.path.isdir(os.path.join('..' ,'Results',f'{modelname}_t','KfieldsQC')):
                os.mkdir(os.path.join('..' ,'Results',f'{modelname}_t','KfieldsQC'))
            res.to_csv(os.path.join('..' ,'Results',f'{modelname}_t','KfieldsQC',f'Ind_x{xcorlen}_z{zcorlen}_f{frac}_n{ens_no}'))
            Kfields.to_csv(os.path.join('..' ,'Results',f'{modelname}_t','KfieldsQC',f'K_x{xcorlen}_z{zcorlen}_f{frac}_n{ens_no}'))
