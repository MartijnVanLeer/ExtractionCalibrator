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
import matplotlib.pyplot as plt
#Setting

def read_snakemake_rule(path, rule: str) -> "snakemake.rules.Rule":
    import snakemake as sm
    workflow = sm.Workflow(snakefile="snakefile")
    workflow.include(path)
    return workflow.get_rule(rule)

if "snakemake" not in globals():
    snakemake = read_snakemake_rule('snakefile','calibration_t')

Location = snakemake.params.Name
model_name = snakemake.params.modelname

BestParams = pd.read_csv(os.path.join('..','Results',f'{model_name}',f'BestParams_SS_{model_name}.csv'))
BestK = BestParams[BestParams['Unnamed: 0'].str[-1] != 'b']
BestGhb = BestParams[BestParams['Unnamed: 0'].str[-1] == 'b']
idx = pd.read_csv(os.path.join('..','Results',f'{model_name}',f'idx_SS_{model_name}.csv'))
ObsWells = pd.read_csv(os.path.join('..','Results',f'{model_name}',f'ObsForCalibration_{model_name}_SS.csv'))
ObsHeads = pd.read_csv(os.path.join('..','Results',f'{model_name}',f'ObsHeads_SS_{model_name}.csv'), index_col = 'Time')
ObsHeads.index =  pd.DatetimeIndex(ObsHeads.index)

CorLayers = snakemake.params.CorLayers
SensLayers = idx.SensLayers.values





#%% Load model and dataset
# Make new folder where files are edited 
orgFolder = os.path.join('..','Results',f'{model_name}', f'{model_name}_t','')
destFolder = os.path.join(orgFolder, 'Fitter','')
OptimisationFuncs.copyOriginalFolder(model_name + '_t', orgFolder ,destFolder , 'Fitter\\' )

#Load model and obs
ds = xr.open_dataset(os.path.join(orgFolder, f'{model_name}_t.nc'))
ds.attrs['model_ws'] = destFolder
ObsHeads = ObsHeads.loc[ds.startdate:ds.enddate]


#Load packages
sim = flopy.mf6.mfsimulation.MFSimulation.load('mfsim', sim_ws = destFolder, exe_name =ds.exe_name)
gwf = sim.get_model()
npf = gwf.get_package('NPF')
npfk = npf.k.data
npfk33 = npf.k33.data
sto = gwf.get_package('STO')
stoss = sto.ss.data
ghb = gwf.get_package('GHB')

npfk, npfk33 = OptimisationFuncs.npf_best_values(npfk, npfk33, idx, BestK, CorLayers)
ghb = OptimisationFuncs.ghb_best_values(ghb,gwf, BestGhb, idx, CorLayers)




#%%
params = OptimisationFuncs.init_params(idx,CorLayers, ghbCal = None, KCal = True, Transient= True)

initsimplex = OptimisationFuncs.initsimplex(params, fac = 0.25)


NMoptions = {'adaptive': True,
              'maxfev' :10,
              'initial_simplex' : initsimplex,
             'xatol' : 0.1, #both xatol and fatol needed for termination
             'fatol' : 10
              }
options = {'options': NMoptions,} 




fitter = lmfit.Minimizer(OptimisationFuncs.run_model_calibration_transient, params, fcn_args = (sim, idx,ObsWells,ObsHeads,ds, CorLayers, npfk, npfk33, stoss, npf, sto), iter_cb=OptimisationFuncs.per_iteration)
result = fitter.minimize('nelder-mead', **options)

print(lmfit.fit_report(result))

#%%
best_params_t = {}
for x in range(len(result.x)):
    best_params_t[result.var_names[x]] = result.x[x]
    
residuals, df, ObsHeads = OptimisationFuncs.run_best_result_transient(best_params_t,sim, idx,ObsWells,ObsHeads,ds, CorLayers, npfk, npfk33, stoss, npf, sto)


df.to_csv(os.path.join('..\\Results',f'{model_name}',f'ModHead_{model_name}.csv'))
ObsHeads.to_csv(os.path.join('..\\Results',f'{model_name}',f'ObsHead_{model_name}.csv'))
best_paramdf = pd.DataFrame.from_dict(best_params_t, orient = 'index', columns = ['Value'])
best_paramdf.to_csv(os.path.join('..\\Results',f'{model_name}',f'BestParams_t_{model_name}.csv'))

#%%plot 
for lay in idx[idx.laytype =='z'].idx.values:
    dfsel = df[ObsWells[ObsWells['Layno'] == lay].putcode]
    obssel = ObsHeads[ObsWells[ObsWells['Layno'] == lay].putcode]

    fig, axes = plt.subplots(nrows=len(dfsel.columns), ncols=1, figsize=(10, 2 * (len(dfsel.columns))))
    fig.set_dpi(600)
    for i, column in enumerate(dfsel.columns):  # Exclude the timestamp column
        ax = axes[i]
        ax.plot(dfsel.index, dfsel[column], label='Mod')
        ax.plot(obssel.index, obssel[column], label='Obs')
        ax.set_title(column)
        ax.legend()
        
    fig.suptitle(lay)
    
    plt.tight_layout()
    fig.save_fig(os.path.join('..\\Results',f'{model_name}',f'{lay}_Obs_Mod_heads.png'))
#%%

fig = plt.figure(figsize=(30,10))    
ax = sns.scatterplot(data = ObsWells, x = 'x_coordinaat', y = 'y_coordinaat', hue = 'Residual',legend = None)
def plotlabel(xvar, yvar, label):
    ax.text(xvar+0.002, yvar, label, fontsize = 20)

ObsWells.apply(lambda x: plotlabel(x['x_coordinaat'], x['y_coordinaat'], x['putcode'][:8]), axis = 1)
