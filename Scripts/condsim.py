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
from tqdm import tqdm

Location = 'Vlijmen'
modelname = 'Budel400'
TI = False
borefolder = join("..",'Data',f'Boringen {Location}',r'Boormonsterprofiel_Geologisch booronderzoek')
ncpath =join('..','Results',modelname,'cache', 'layer_ds.nc')
TIpath = join("..",'Data','AllDB.csv')
RegCat = join('..', 'Data','REGIS_catalogus_export18122015.xlsx' )

Layer = 'WAk2'

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

rn = boringen.plot_K_weighted(Kcore)
#Trek K uit Kcore en voeg toe aan boringlijst
# cond.add_k(Kcore)

boringen.list.to_csv(f'../Results/{Location}/boreholeindicators.csv')

#%%


xmin = boringen.layermodel.extent[0]
ymin = boringen.layermodel.extent[2]
Lx = boringen.layermodel.extent[1] - boringen.layermodel.extent[0]
Ly = boringen.layermodel.extent[3] - boringen.layermodel.extent[2]
zmin = cond.z.min()

xcorlens = [200,400,600,800]
zcorlens = [0.5,4,8]
fracs = [0.5,0.55,0.60,0.65]#[0.5,0.6,0.7,0.8]

for xcorlen in tqdm(xcorlens, position = 1):
    for zcorlen in tqdm(zcorlens, position = 2):
        for frac in tqdm(fracs, position = 3):
            ens_no = 50
            dx = dy = 50
            
            
            start = timer() 
            a,res = SISIM_R.Cond_SISIM(cond[['x','y','z','i']],
                        xmin = boringen.layermodel.extent[0],ymin = boringen.layermodel.extent[2],zmin = cond.z.min(),
                        Lx=boringen.layermodel.extent[1] - boringen.layermodel.extent[0],
                        Ly=boringen.layermodel.extent[3] - boringen.layermodel.extent[2],
                        Lz =abs(cond.z.min()-cond.z.max()),
                        dx =dx,dy =dy,dz =1,
                        xcorlen =xcorlen, zcorlen = zcorlen,
                        ens_no = 50, frac =frac, nmax = 100, seed = 1337)
            print(f'{round(timer() - start,1)} s for conditional simulation')
            
            Kfields = boringen.add_k(res)
            if not os.path.isdir(f'../Results/{Location}/KfieldsQC'):
                os.mkdir(f'../Results/{Location}/KfieldsQC')
            res.to_csv(f'../Results/{Location}/KfieldsQC/Ind_x{xcorlen}_z{zcorlen}_f{frac}_n{ens_no}')
            Kfields.to_csv(f'../Results/{Location}/KfieldsQC/K_x{xcorlen}_z{zcorlen}_f{frac}_n{ens_no}')
