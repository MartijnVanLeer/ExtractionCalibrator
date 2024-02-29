# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 15:08:09 2023

@author: leermdv
"""
import os
from shutil import copyfile
import numpy as np
import pandas as pd
import SALib
import nlmod
import lmfit
import flopy 
from time import sleep

def zscore(s, window, thresh=3, return_all=False):
    roll = s.rolling(window=window, min_periods=1, center=True)
    avg = roll.mean()
    std = roll.std(ddof=0)
    z = s.sub(avg).div(std)   
    m = z.between(-thresh, thresh)

    return s.where(m, np.nan)

def fix_heads(ObsHeads, minhead = None, maxhead = None):
    for series_name, series in ObsHeads.items():
        if minhead != None:
            ObsHeads[series_name] = series.where(series>minhead)
        if maxhead != None:
            ObsHeads[series_name] = series.where(series<maxhead)
        ObsHeads[series_name] = zscore(ObsHeads[series_name], window = 50, thresh = 4)
    return ObsHeads
        
    
def copyOriginalFolder(modelName, orgFolder ,destFolder , mfsimPrefix):
    if not os.path.isdir(destFolder):
        os.makedirs(destFolder, exist_ok=True)
    namLines = open(os.path.join(orgFolder, modelName + '.nam')).readlines()
    mfsimLines = open(os.path.join(orgFolder, 'mfsim.nam')).readlines()
    mfsimDestFile = open(os.path.join(destFolder,'mfsim.nam'), 'w')
    
    #Writing mfsim
    for line in mfsimLines:
        mfsimDestFile.write(line)
    mfsimDestFile.close()
    
    #packages
    packageFiles = []
    for line in namLines:
        if modelName + '.' in line and line.endswith('lst\n') == False:
            line = line.replace('\n','')
            line = line.split(' ')
            packageFiles.append(line[4])
    for line in mfsimLines:
        if modelName+'.' in line:
            line = line.replace('\'', '').replace('MODFLOW','').replace('\n', '')
            packageFiles.append(line.split(' ')[4])
    if os.path.exists(os.path.join(orgFolder,f'{modelName}.rch_ts')):
        packageFiles.append(f'{modelName}.rch_ts')
    for package in packageFiles:
        copyfile(os.path.join(orgFolder,package),os.path.join(destFolder, package))
        
def readObs(Location,gwf,ds):
    ObsWells = pd.read_csv(os.path.join('..','Data','dawaco',f'waarnemingsputten_WG_{Location}.csv'))
    ObsWells['putcode'] = ObsWells['putcode'].astype(str) + '_' +ObsWells['filter_nummer'].astype(str)
    ObsWells['putcode'] = ObsWells.putcode.str.rstrip('0').str.rstrip('.')   
    ObsWells = add_cellid(ObsWells,gwf,ds)
    return ObsWells

def GetObs(model_name, Location, idx,ds):
    ObsWells = pd.read_csv(os.path.join('..','Results',f'{model_name}',f'ObsForCalibration_{Location}.csv'))
    ObsHeads = pd.read_csv(os.path.join('..','Data','Preprocessed',f'stijghoogtereeksen_{Location}.csv'), index_col= 'Time')
    ObsWells = ObsWells[ObsWells['Layno'].isin(idx.idx)]
    columns = list(set(ObsWells.putcode).intersection(ObsHeads.columns))
    ObsHeads = ObsHeads[columns]
    ObsHeads = ObsHeads.set_index(pd.to_datetime(ObsHeads.index))
    ObsWells = ObsWells[ObsWells.putcode.isin(columns)]
    ObsHeads = ObsHeads[ObsWells.putcode]
    ObsHeads = ObsHeads.loc[ds.time.start:ds.time.values[0]]
    return ObsWells, ObsHeads

def layerid(SensLayers,ds):
    df = pd.DataFrame({'SensLayers' : SensLayers})
    idx = []
    for x in range(len(ds.layer)):
        if ds.layer[x].values in SensLayers:
            idx.append(x)
                
    df['idx'] = idx 
    df['laytype'] = df.SensLayers.str[-2]
    return df 

def problem(SensLayers, ObsWells):

    problem = SALib.ProblemSpec({'num_vars' : len(SensLayers),
            'names' : SensLayers,
            'bounds' : [[0.5,2]] * len(SensLayers),
            'outputs' : ObsWells.putcode.values})
    return problem

def add_cellid(ObsWells,gwf, ds):
    CellID = []
    Layno = []
    laytypes = ds.layer.str[-2]
    for index,well in ObsWells.iterrows():
            try:
                idx = gwf.modelgrid.intersect(well['x_coordinaat'] ,well['y_coordinaat'], -(well['filter_onderkant_ref'] - well['filter_referentiepunt_NAP']))
                CellID.append(idx[1])
                if laytypes[idx[0]] == 'k':
                    if well.putcode == 'B57E0064_3':
                        Layno.append(idx[0]-1)
                    else:
                        Layno.append(idx[0]+1)
                else:
                Layno.append(idx[0])
            except Exception:
                print(f'{well["putcode"]} outside model area')
                CellID.append(np.nan)
                Layno.append(np.nan)
                    
    ObsWells['CellID'] = CellID
    ObsWells['Layno'] = Layno

    ObsWells.dropna(inplace = True)
    return ObsWells

def run_model(X, sim , idx ,npf, npfk,npfk33, ObsWells,ds):
    newk = npfk.copy()
    newk33 = npfk33.copy()
    for i,layer in idx.iterrows():

        newk[layer.idx] = npfk[layer.idx] * X[i]
        newk33[layer.idx] = npfk33[layer.idx] * X[i]
            
    npf.k = newk
    npf.k33 = newk33
    npf.write()
    sim.run_simulation(silent = True)
    result = []
    head = nlmod.gwf.get_heads_da(ds)
    for index, well in ObsWells.iterrows():
        obsheads = head[:,int(well.Layno), int(well.CellID)]
        result.append(obsheads.values[0])
    return result

def init_params(idx,CorLayers, ghbCal, KCal, Transient = False):
    params = lmfit.Parameters()
    for index, lay in idx.iterrows():
        if lay['SensLayers'] not in CorLayers.values():
            if KCal:
                params.add(name = lay['SensLayers'], value = 0)# min  = -1, max = 1)
            if ghbCal == 'SensLayers':
                if 'z' in lay['SensLayers'][-2:]:
                    params.add(name = lay['SensLayers'] + "_ghb", value = 0)
    if ghbCal == 'Single':
        params.add(name = "ghb", value = 0)   
    if Transient:
        params.add(name = 'SS', value = 0)
    return params

def run_calibration_ss(p, sim ,gwf, idx ,npf, npfk,npfk33, ghb,ghb_spd,ObsWells, ObsHeads,ds,CorLayers,ghbCal, KCal):
    #K
    if KCal:
        newk = npfk.copy()
        newk33 = npfk33.copy()
        for i,layer in idx.iterrows():
            if layer.SensLayers not in CorLayers.values():
                newk[layer.idx] =npfk[layer.idx] * 2**p[layer.SensLayers]
                newk33[layer.idx] = npfk33[layer.idx] * 2**p[layer.SensLayers]
            if layer.SensLayers in CorLayers.keys():
                layno = idx[idx['SensLayers'] == CorLayers[layer.SensLayers]].idx.values[0]
                newk[layno] = npfk[layno] * 2**p[layer.SensLayers]
                newk33[layno] = npfk33[layno] * 2**p[layer.SensLayers]
                
        npf.k = newk
        npf.k33 = newk33
        while True:
            if os.access(os.path.join(sim.sim_path,npf.path[0] + '.' + npf.path[1]), os.W_OK):
                npf.write()
                break
    #ghb

    if ghbCal != None:
        new_spd = []
        for row in ghb_spd:
            # if row[0][0] in list(range(min(ObsWells.Layno.astype('int')),max(ObsWells.Layno.astype('int')+1))):
                if row[0][0] in idx.idx.values:
                        if idx[idx.idx == row[0][0]].laytype.values =='z':
                            if ghbCal == 'SensLayers':
                                if idx[idx.idx == row[0][0]].SensLayers.values not in CorLayers.values():
                                    new_spd.append([row[0],row[1] + p[f'{idx[idx.idx == row[0][0]].SensLayers.values[0]}_ghb'], row[2]])
                                else:                                
                                    new_spd.append([row[0],row[1] + p[f'{invdict(CorLayers, idx[idx.idx == row[0][0]].SensLayers.values)}_ghb'], row[2]])
                            elif ghbCal =='Single':
                                new_spd.append([row[0],row[1] + p['ghb'], row[2]])
                else: 
                    new_spd.append([row[0], row[1], row[2]])
        gwf.remove_package('GHB')
        ghb = flopy.mf6.ModflowGwfghb(gwf, stress_period_data=new_spd,save_flows=True,maxbound = len(new_spd))
        while True:
            if os.access(os.path.join(sim.sim_path,ghb.path[0] + '.' + ghb.path[1]), os.W_OK):
                ghb.write()
                break
    
    #run        
    sim.run_simulation(silent = True)
    
    #results
    head = nlmod.gwf.get_heads_da(ds)
    df = pd.DataFrame()
    ObsHeadsSS = ObsHeads.mean()
    
    Weights = ObsHeadsSS.copy()
    for index, well in ObsWells.iterrows():
        modheads = head[:,int(well.Layno), int(well.CellID)]
        df[f'{well["putcode"]}'] = modheads.values
        Weights[f'{well["putcode"]}'] = well.Sensitivity

        
    residuals = df - ObsHeadsSS
    residuals = residuals * Weights
    
    residuals =residuals.to_numpy()[0]
    residuals = residuals[~np.isnan(residuals)]
    residuals = sum(residuals**2)
    sleep(0.5)
    return residuals 

def run_calibration_ss_result(p, sim ,gwf, idx ,npf, npfk,npfk33, ghb,ghb_spd,ObsWells, ObsHeads,ds,CorLayers,ghbCal,KCal):
    
    #K
    
    if KCal:
        newk = npfk.copy()
        newk33 = npfk33.copy()
        for i,layer in idx.iterrows():
            if layer.SensLayers not in CorLayers.values():
                newk[layer.idx] =npfk[layer.idx] * 2**p[layer.SensLayers]
                newk33[layer.idx] = npfk33[layer.idx] * 2**p[layer.SensLayers]
            if layer.SensLayers in CorLayers.keys():
                layno = idx[idx['SensLayers'] == CorLayers[layer.SensLayers]].idx.values[0]
                newk[layno] = npfk[layno] * 2**p[layer.SensLayers]
                newk33[layno] = npfk33[layno] * 2**p[layer.SensLayers]
                
        npf.k = newk
        npf.k33 = newk33
        npf.write()
    #ghb

    if ghbCal != None:
        new_spd = []
        for row in ghb_spd:
            # if row[0][0] in list(range(min(ObsWells.Layno.astype('int')),max(ObsWells.Layno.astype('int')+1))):
                if row[0][0] in idx.idx.values:
                        if idx[idx.idx == row[0][0]].laytype.values =='z':
                            if ghbCal == 'SensLayers':
                                if idx[idx.idx == row[0][0]].SensLayers.values not in CorLayers.values():
                                    new_spd.append([row[0],row[1] + p[f'{idx[idx.idx == row[0][0]].SensLayers.values[0]}_ghb'], row[2]])
                                else:                                
                                    new_spd.append([row[0],row[1] + p[f'{invdict(CorLayers, idx[idx.idx == row[0][0]].SensLayers.values)}_ghb'], row[2]])
                            elif ghbCal =='Single':
                                new_spd.append([row[0],row[1] + p['ghb'], row[2]])
                else: 
                    new_spd.append([row[0], row[1], row[2]])
        gwf.remove_package('GHB')
        ghb = flopy.mf6.ModflowGwfghb(gwf, stress_period_data=new_spd,save_flows=True,maxbound = len(new_spd))
        ghb.write()
    
            
    sim.run_simulation(silent = True)
    head = nlmod.gwf.get_heads_da(ds)
    df = pd.DataFrame()
    ObsHeadsSS = ObsHeads.mean()
    
    for index, well in ObsWells.iterrows():
        modheads = head[:,int(well.Layno), int(well.CellID)]
        df[f'{well["putcode"]}'] = modheads.values
        
    ObsWells['Residual'] = (df - ObsHeadsSS).values[0]
    ObsWells['ObsHeadsSS'] = ObsHeadsSS.values
    
    ObsWells['ModHead'] = df.values[0]
    return ObsWells
def per_iteration(pars, iteration, resid, *args, **kws):
    print(" Iteration:  ", iteration, [f"{p.name} = {p.value:.2f}" for p in pars.values()], round(resid,5))
    
def invdict(dic, search_item):
    for key, item in dic.items():
        if item == search_item:
            return key
    
def generate_initial_simplex(midparams,fac = 1):
    n_dim = len(midparams)
    initial_simplex = np.zeros((n_dim + 1, n_dim))
    
    # Choose the first vertex as the midpoint between min and max for each dimension
    initial_simplex[0] = np.zeros((len(midparams)))
    # For the remaining vertices, vary one dimension at a time
    for i in range(n_dim-1):
        vertex = initial_simplex[0].copy()
        
        vertex[i] = fac
        vertex[i+1] = -fac
        initial_simplex[i + 1] = vertex
    initial_simplex[-1][-1] = fac
    initial_simplex[0][0] = -fac
    
    return initial_simplex

def initsimplex(params, fac, best_params = None):
    midparams = []
    for p in params:
        midparams.append(params[p].value)
    initsimplex = generate_initial_simplex( midparams, fac)
    return initsimplex

def npf_best_values(npfk, npfk33, idx, BestK, CorLayers):
    for i,layer in idx.iterrows():
        if layer.SensLayers not in CorLayers.values():
            npfk[layer.idx] =npfk[layer.idx] * 2**BestK.iloc[i,1]
            npfk33[layer.idx] = npfk33[layer.idx] * 2**BestK.iloc[i,1]
        if layer.SensLayers in CorLayers.keys():
            layno = idx[idx['SensLayers'] == CorLayers[layer.SensLayers]].idx.values[0]
            npfk[layno] = npfk[layno] * 2**BestK.iloc[i,1]
            npfk33[layno] = npfk33[layno] * 2**BestK.iloc[i,1]
    return npfk, npfk33

def ghb_best_values(ghb, gwf, BestGhb, idx, CorLayers):
    
    ghb_spd = ghb.stress_period_data.data[0]
    new_spd = []
    for row in ghb_spd:
        # if row[0][0] in list(range(min(ObsWells.Layno.astype('int')),max(ObsWells.Layno.astype('int')+1))):
            if row[0][0] in idx.idx.values:
                    if idx[idx.idx == row[0][0]].laytype.values =='z':
                            if idx[idx.idx == row[0][0]].SensLayers.values not in CorLayers.values():
                                new_spd.append([row[0],row[1] + BestGhb[BestGhb['Unnamed: 0'] ==  f'{idx[idx.idx == row[0][0]].SensLayers.values[0]}_ghb'].Value.values[0], row[2]])
                            else:                                
                                new_spd.append([row[0],row[1] + BestGhb[BestGhb['Unnamed: 0'] ==  f'{invdict(CorLayers, idx[idx.idx == row[0][0]].SensLayers.values[0])}_ghb'].Value.values[0], row[2]])

            else: 
                new_spd.append([row[0], row[1], row[2]])
    gwf.remove_package('GHB')
    ghb = flopy.mf6.ModflowGwfghb(gwf, stress_period_data=new_spd,save_flows=True,maxbound = len(new_spd))
    ghb.write()
    
def run_model_calibration_transient(p, sim, idx,ObsWells,ObsHeads,ds, CorLayers, npfk, npfk33, stoss, npf, sto):
    newk = npfk.copy()
    newk33 = npfk33.copy()
    
    for i,layer in idx.iterrows():
        if layer.SensLayers not in CorLayers.values():
            newk[layer.idx] =npfk[layer.idx] * 2**p[layer.SensLayers]
            newk33[layer.idx] = npfk33[layer.idx] * 2**p[layer.SensLayers]
        if layer.SensLayers in CorLayers.keys():
            layno = idx[idx['SensLayers'] == CorLayers[layer.SensLayers]].idx.values[0]
            newk[layno] = npfk[layno] * 2**p[layer.SensLayers]
            newk33[layno] = npfk33[layno] * 2**p[layer.SensLayers]
            
    npf.k = newk
    npf.k33 = newk33
    npf.write()
    sto.ss = stoss * 10**p['SS']
    sto.write()
    sim.run_simulation(silent = True)
    head = nlmod.gwf.get_heads_da(ds)
    df = pd.DataFrame(index = pd.DatetimeIndex(ObsHeads.index))
    for index, well in ObsWells.iterrows():
        modheads = head.isel(layer = int(well.Layno)).sel(icell2d = int(well.CellID)).sel(time = slice(pd.to_datetime(ObsHeads.index[0]),pd.to_datetime(ObsHeads.index[-1]) ))
        df[f'{well["putcode"]}'] = modheads.values
        
    # df = df - df.mean()
    # ObsHeads = ObsHeads - ObsHeads.mean()
    
    residuals = df - ObsHeads

    residuals = residuals.to_numpy().flatten()
    residuals = residuals[~np.isnan(residuals)]
    residuals = sum(residuals**2)
    return residuals

def run_best_result_transient(p, sim, idx,ObsWells,ObsHeads,ds, CorLayers, npfk, npfk33, stoss, npf, sto):
    newk = npfk.copy()
    newk33 = npfk33.copy()
    
    for i,layer in idx.iterrows():
        if layer.SensLayers not in CorLayers.values():
            newk[layer.idx] =npfk[layer.idx] * 2**p[layer.SensLayers]
            newk33[layer.idx] = npfk33[layer.idx] * 2**p[layer.SensLayers]
        if layer.SensLayers in CorLayers.keys():
            layno = idx[idx['SensLayers'] == CorLayers[layer.SensLayers]].idx.values[0]
            newk[layno] = npfk[layno] * 2**p[layer.SensLayers]
            newk33[layno] = npfk33[layno] * 2**p[layer.SensLayers]
            
    npf.k = newk
    npf.k33 = newk33
    npf.write()
    sto.ss = stoss * 10**p['SS']
    sto.write()
    sim.run_simulation(silent = True)
    head = nlmod.gwf.get_heads_da(ds)
    df = pd.DataFrame(index = pd.DatetimeIndex(ObsHeads.index))
    for index, well in ObsWells.iterrows():
        modheads = head.isel(layer = int(well.Layno)).sel(icell2d = int(well.CellID)).sel(time = slice(pd.to_datetime(ObsHeads.index[0]),pd.to_datetime(ObsHeads.index[-1]) ))
        df[f'{well["putcode"]}'] = modheads.values
        
    # df = df - df.mean()
    # ObsHeads = ObsHeads - ObsHeads.mean()
    
    residuals = df - ObsHeads

    return residuals, df, ObsHeads