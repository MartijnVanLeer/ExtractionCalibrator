# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 09:53:38 2023

@author: leermdv
"""
import Heterogeniteit as h
import xarray as xr
from os.path import join
import pickle

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
ds = xr.open_dataset(join('..','Results',f'{modelname}', f'{modelname}_t',f'{modelname}_t.nc'))
TIpath = join("..",'Data','AllDB.csv')
RegCat = join('..', 'Data','REGIS_catalogus_export18122015.csv' )



boringen = h.boringen(borefolder)
md = boringen.metadata
#voeg REGIS laag toe aan boringen
#boringen.add_layers(xr.open_dataset(ncpath))
#selecteer trajecten binnen een laag en explode naar meterschaal
# boringen.select_layer(Layer)

boringen.select_range(Layer, ds)

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

boringen.list.to_csv(join('..','Results',f'{modelname}','boreholeindicators.csv'))
with open(join('..','Results',f'{modelname}','boreholeindicators.pkl'), 'wb') as f:
    pickle.dump(boringen, f)