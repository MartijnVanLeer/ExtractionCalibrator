# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 14:32:52 2023

@author: leermdv
"""

import xarray as xr
import flopy 
import os
import OptimisationFuncs
import numpy as np
import pandas as pd
from tqdm import tqdm
import lmfit
import nlmod
import seaborn as sns
#Setting

def read_snakemake_rule(path, rule: str) -> "snakemake.rules.Rule":
    import snakemake as sm
    workflow = sm.Workflow(snakefile="snakefile")
    workflow.include(path)
    return workflow.get_rule(rule)

if "snakemake" not in globals():
    snakemake = read_snakemake_rule('snakefile','calibration_ss')

Location = snakemake.params.Name
modelname = snakemake.params.modelname
SensLayers = snakemake.params.SensLayers
CorLayers = snakemake.params.CorLayers
ghbCal = snakemake.params.ghbCal # 'obs', 'Single', None
KCal = snakemake.params.KCal
Weighted = snakemake.params.Weighted
BadWells =snakemake.params.BadWells # alleen begin'B57E0082_7','B57E0147_4'
Lambda = snakemake.params.Lambda
method = snakemake.params.method


#%% Load model and dataset
# Make new folder where files are edited 
orgFolder = os.path.join('..','Results',f'{modelname}', f'{modelname}_ss')
destFolder = os.path.join(orgFolder, 'Fitter', '')
OptimisationFuncs.copyOriginalFolder(modelname + '_ss', orgFolder ,destFolder , 'Fitter\\' )

#Load model and obs
ds = xr.open_dataset(os.path.join(orgFolder, f'{modelname}_ss.nc'))
ds.attrs['model_ws'] = destFolder

#Load packages
sim = flopy.mf6.mfsimulation.MFSimulation.load('mfsim', sim_ws = destFolder, exe_name =ds.exe_name)
gwf = sim.get_model()
npf = gwf.get_package('NPF')
npfk = npf.k.data
npfk33 = npf.k33.data
ghb = gwf.get_package('GHB')
ghb_spd = ghb.stress_period_data.data[0]

#Get some info from model 
idx = OptimisationFuncs.layerid(SensLayers,ds)
ObsWells = OptimisationFuncs.readObs(Location,gwf,ds)


#%% Sensitivity analysis

problem = OptimisationFuncs.problem(SensLayers, ObsWells)
Sens = np.ones((len(problem['outputs'])))
if Weighted:
    problem.sample_latin(100,seed = 3942)
    
    X = problem.samples
    Y = np.zeros((X.shape[0],len(problem['outputs'])))
    
    for i, x in tqdm(enumerate(X)):
        results = OptimisationFuncs.run_model(x, sim , idx , npf, npfk,npfk33, ObsWells,ds)
        Y[i,:] = results
    
    #correlation between observations
    cordf = pd.DataFrame(np.corrcoef(Y, rowvar = False), columns = ObsWells.putcode, index = ObsWells.putcode)
    # ObsWells = ObsWells.merge(cordf, left_on = 'putcode', right_index = True)

    for i, x in enumerate(problem['outputs']):
        Sens[i] = np.std(Y[:,i])#S1['ST']#*np.var(Y[:,i])*10000

ObsWells['Sensitivity'] = Sens


Best = ObsWells[ObsWells['Sensitivity'] > 0.01]
Best = Best[~Best.putcode.isin(BadWells)]
if Weighted: 
    Best = Best.drop_duplicates(subset = ['Sensitivity'])
Best.to_csv(os.path.join('..','Results',f'{modelname}',f'ObsForCalibration_{Location}.csv'))
ObsWells.to_csv(os.path.join(destFolder,'ObsWellsSens.csv'))

#%%Calibrate ss

ObsWells, ObsHeads = OptimisationFuncs.GetObs(modelname, Location, idx,ds)
ObsHeads = OptimisationFuncs.fix_heads(ObsHeads, minhead = 10, maxhead = None)

'''
Layer K values powers of 2 (-1 = *0.5, 1 = *2)
Initital simplex contains fac, -fac, 0 
Those are determined by the min/max
'''

params = OptimisationFuncs.init_params(idx,CorLayers, ghbCal, KCal, method)   

if method == 'NM':
    initsimplex = OptimisationFuncs.initsimplex(params, fac = 0.3)
    NMoptions = {'adaptive': True,
                'maxfev' :1000,
                'initial_simplex' : initsimplex,
                'xatol' : 0.1, #both xatol and fatol needed for termination
                'fatol' : 0.01
                }
    options = {'options': NMoptions,} 
    fitter = lmfit.Minimizer(OptimisationFuncs.run_calibration_ss, params,
                            fcn_args = (sim,gwf, idx ,npf,npfk, npfk33, ghb,ghb_spd, ObsWells, ObsHeads,ds,CorLayers,ghbCal, KCal,Lambda, method), iter_cb=OptimisationFuncs.per_iteration)
    result = fitter.minimize('Nelder-Mead',**options)
elif method == 'LM':
    LMoptions = {'ftol' : 1e-8,
                'xtol' : 1e-8,
                'max_nfev' : 1000,
                'factor' : 50}
    fitter = lmfit.Minimizer(OptimisationFuncs.run_calibration_ss, params,
                            fcn_args = (sim,gwf, idx ,npf,npfk, npfk33, ghb,ghb_spd, ObsWells, ObsHeads,ds,CorLayers,ghbCal, KCal,Lambda, method), iter_cb=OptimisationFuncs.per_iteration)
    result = fitter.minimize('leastsq',**LMoptions)

print(lmfit.fit_report(result))
#%% Run with best fit for residuals
best_params = {}
for x in range(len(result.x)):
    best_params[result.var_names[x]] = result.x[x]
    # best_params[result.var_names[x]] = 0
ObsWells = OptimisationFuncs.run_calibration_ss_result(best_params, sim ,gwf, idx ,npf,npfk, npfk33, ghb,ghb_spd,ObsWells, ObsHeads,ds,CorLayers, ghbCal, KCal )
print(f'MAE : {ObsWells.Residual.abs().mean()}')

# #%% plot
import matplotlib.pyplot as plt
# head = nlmod.gwf.get_heads_da(ds)
# pmv = flopy.plot.PlotMapView(modelgrid=gwf.modelgrid)
# array = head.isel(layer=22)
# ar = pmv.plot_array(array)
# cb = plt.colorbar(ar)

# #%%
# layno = 20
# fig, ax = plt.subplots()
# sns.scatterplot(data = ObsWells[ObsWells.Layno == layno], x = 'x_coordinaat', y = 'y_coordinaat', hue = 'Residual')
# ax.set_xlim((ds.extent[0], ds.extent[1]))
# ax.set_ylim((ds.extent[2], ds.extent[3]))
# # ObsHeads.loc[:,ObsWells[ObsWells.Layno == layno].putcode].plot()
fig, ax = plt.subplots()
sns.scatterplot(data = ObsWells, x = 'ObsHeadsSS', y = 'ModHead', hue = 'Layno', ax =ax)
ax.axline([min(ObsWells.ObsHeadsSS), min(ObsWells.ObsHeadsSS)],[max(ObsWells.ObsHeadsSS), max(ObsWells.ObsHeadsSS)])
fig.savefig(os.path.join('..','Results',f'{modelname}', 'Calibration_ss.png'))

#%%
best_paramdf = pd.DataFrame.from_dict(best_params, orient = 'index', columns = ['Value'])
ObsWells.to_csv(os.path.join('..','Results',f'{modelname}',f'ObsForCalibration_{modelname}_SS.csv'))
best_paramdf.to_csv(os.path.join('..','Results',f'{modelname}',f'BestParams_SS_{modelname}.csv'))
idx.to_csv(os.path.join('..','Results',f'{modelname}',f'idx_SS_{modelname}.csv'))
ObsHeads.to_csv(os.path.join('..','Results',f'{modelname}',f'ObsHeads_SS_{modelname}.csv'))

