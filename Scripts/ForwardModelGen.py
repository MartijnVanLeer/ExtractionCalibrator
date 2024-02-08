# -*- coding: utf-8 -*-
"""
Created on Mon May 15 16:02:46 2023

Forward model setup script with pure nlmod

@author: leermdv
"""

import nlmod
import flopy
import pandas as pd
import Helper
import os
import numpy as np

#%%Settings (model)
Name = 'Budel' #puttenveld Budel, Schijf, Vlijmen
modelname ='biervoor2018'

Range = 4000 #m Distance from grid border to centre of wells 3200
GHBrange = 3200 
delr = delc = 800
refineranges = {200:4}

steady_state = True

warmup = False if steady_state else True
use_geotop = False
NLzuid = False if Name == 'Vlijmen' else True
use_ahn = False #True if Name == 'Budel' else False
use_knmi = True


modelname = modelname +  '_ss' if steady_state else modelname + '_t'
    
startdate = '2014-12-31 00:00:00' #eerste is '2014-12-31 00:00:00', '2018-06-01 00:00:00' is grensdatum droogte budel
enddate = '2018-06-01 00:00:00' #laatste '2023-12-30 22:00:00'

#%% Loading data


print('Loading files..')
#Wells
ExWells = pd.read_csv(f'..\\Data\\dawaco\\winputten_WG_{Name}.csv')
ObsWells = pd.read_csv(f'..\\Data\\dawaco\\waarnemingsputten_WG_{Name}.csv')
Discharge_original = pd.read_csv(f'..\\Data\\pidata\\pidata_WG_{Name}.csv')
ObsHeads = pd.read_csv(f'..\\Data\\dawaco\\stijghoogtereeksen_WG_{Name}.csv')
lhmpath = '..\\Data\\lhm.nc'


Discharge = Helper.FixTS_Discharge(Discharge_original, startdate, enddate)
ObsHeads = Helper.FixTS_Obs(ObsHeads, startdate, enddate,ObsWells, Name)


if warmup:
    ExWells_warmup = pd.read_csv(f'..\\Data\\Langjarig\\gwo_data_{Name}_wells.csv', sep = ';')
    ExWellsAll = Helper.Add_warmup_wells(ExWells,ExWells_warmup)
    Discharge_warmup = pd.read_csv(f'..\\Data\\Langjarig\\gwo_data_{Name}_debieten.csv', sep = ';')
    Discharge_warmup_fixed = Helper.FixTS_Discharge_warmup(Discharge_warmup, startdate)
    
    Discharge = pd.concat([Discharge, Discharge_warmup_fixed])
    Discharge = Discharge.fillna(0)
    Discharge = Discharge.sort_index()
    ExWellsAll = ExWellsAll[ExWellsAll.putcode.isin(list(set(ExWells.putcode).intersection(Discharge.columns)))]
else:
    ExWellsAll = ExWells
    
ExWellsAll, Discharge = Helper.add_refresco(ExWellsAll,Discharge)
ExWellsAll, Discharge = Helper.add_bier(ExWellsAll,Discharge)
# ExWellsAll, _ = Helper.add_refresco(ExWellsAll, Discharge)
WellGdf = Helper.make_gdf(ExWells,ObsWells)    
    
    
    
model_ws = f"..\\Results\\{modelname}"

if not os.path.isdir(model_ws):
    os.mkdir(model_ws)



if not nlmod.util.check_presence_mfbinaries():
    nlmod.util.download_mfbinaries()
# figdir, cachedir = nlmod.util.get_model_dirs(model_ws)
cachedir=None


#%% Read source data
#To redownload uncomment this

# nlmod.cache.clear_cache(cachedir)
# for f in os.listdir(cachedir):
#     os.remove(os.path.join(cachedir, f))

#Regis 
print('Loading REGIS..')
extent = Helper.GetExtent(ExWells, Range)
layer_model = Helper.layermodel(extent, NLzuid,cachedir, "layer_ds")
# print(f'spreidingslengte:{Helper.Spreidingslengte(extent, "PZWAz4","WAk1", layer_model)}') 
#%% Grid functions
print('Making grid..')
ds = nlmod.to_model_ds(layer_model, modelname, model_ws, delr=delr)
refinements = Helper.refiner(ds, refineranges, WellGdf)
ds = nlmod.grid.refine(ds,refinement_features=refinements)
ds = Helper.resample(ds, layer_model, NLzuid)



if use_ahn:
    ahn_ds = nlmod.read.ahn.get_ahn(ds, cachedir=ds.cachedir, cachename="ahn")
    ds.update(ahn_ds)

print('Loading KNMI..')

if use_knmi and not steady_state:
        ds = nlmod.time.set_ds_time(ds, start =1, time= Discharge.index, steady_start=True)
        knmi_ds = nlmod.read.knmi.get_recharge(ds, cachedir=ds.cachedir, cachename="recharge", most_common_station=True)
        knmi_ds['recharge'].loc['time' == 'startdate'] = knmi_ds.recharge.mean(dim = 'time')
        ds.update(knmi_ds)
elif use_knmi and steady_state:        
        ds = nlmod.time.set_ds_time(ds, start =Discharge.index[0], time=  Discharge.index[-1], steady=steady_state,steady_start=True)
        knmi_ds = nlmod.read.knmi.get_recharge(ds, cachedir=cachedir, cachename='recharge')
        ds.update(knmi_ds)
elif not use_knmi and not steady_state:
        ds = nlmod.time.set_ds_time(ds, start =1, time= Discharge.index, steady_start=True)
elif not use_knmi and steady_state:        
        ds = nlmod.time.set_ds_time(ds, start =Discharge.index[0], time=  Discharge.index[-1], steady=steady_state,steady_start=True)
        




#%% Create sim 
print('Making sim..')
# create simulation
sim = nlmod.sim.sim(ds)

# create time discretisation
if not steady_state:
    tdis = nlmod.sim.tdis(ds, sim)
    ds.attrs['startdate'] = startdate
    ds.attrs['enddate'] = enddate
else:
    tdis = flopy.mf6.modflow.mftdis.ModflowTdis(
        sim,
        time_units=ds.time.time_units,
        nper=1,
    )

# create ims
ims = nlmod.sim.ims(sim, complexity = 'COMPLEX')

# create groundwater flow model
gwf = nlmod.gwf.gwf(ds, sim)

# Create discretization
dis = nlmod.gwf.dis(ds, gwf)

# create node property flow
npf = nlmod.gwf.npf(ds, gwf)

# Create the initial conditions package
ic = nlmod.gwf.ic(ds, gwf, starting_head=0)

if not steady_state:
    sto = nlmod.gwf.sto(ds, gwf)

print('WEL package..')
ds['idomain'] = nlmod.layers.get_idomain(ds)
#create wel package

wel_spd = Helper.get_wel_spd(ds, gwf,Discharge, ExWellsAll,steady_state)
wel = flopy.mf6.ModflowGwfwel(
        gwf,
        stress_period_data=wel_spd)

print('GHB package..')
# add constant head cells at model boundaries
ds.update(nlmod.grid.mask_model_edge(ds))
ghb = Helper.ghb(ds, gwf,cachedir,NLzuid, GHBrange, lhmpath = lhmpath, delr = delr)
Helper.plot_map(ds, gwf, 'ghb_head', 'KIz3')
#Create drain packakge from AHN 
drn = nlmod.gwf.surface_drain_from_ds(ds, gwf, resistance=1, elev = 'top')

riv = Helper.riv(ds,gwf)

print('WEL package..')
#Create recharge packagefrom KNMI
if use_knmi:
    if not steady_state:
        rch = nlmod.gwf.rch(ds, gwf)
    else:
        rch = Helper.steady_rch(ds,gwf)



# Create the output control package
oc = nlmod.gwf.oc(ds, gwf)


#%%run
nlmod.sim.write_and_run(sim, ds, write_ds=True) 
#%% post-processing
ModHeads = Helper.GetHeadsAtObs(ds, ObsWells, gwf)

# ModHeads.to_csv(f"{cachedir}\\Modelled_obs.csv")

#%%
putcodes = ObsWells[ObsWells['filter_bovenkant_ref'] >25].sort_values('filter_bovenkant_ref')['putcode']
putcodes = putcodes[putcodes.isin(ModHeads.columns)]
if not steady_state:
    Helper.ComparePlot(ObsHeads[startdate:], ModHeads[startdate:], putcodes, ObsWells)

#%%
# plot head in extraction layer

if not steady_state:   
    Helper.PlotPumpedHeads(ds,gwf,400, tstep = '2019-05-25 00:00:00')
else:
    Helper.PlotPumpedHeads(ds,gwf,0, tstep = '2015-01-01')
#%%
import Helper
Helper.plot_mean_res(ObsHeads,ModHeads, ObsWells, steady_state, depth = 350, tail = 365*1)

#%%
import Helper
Helper.plot_map(ds, gwf, 'ghb_head', 'KIz5')




