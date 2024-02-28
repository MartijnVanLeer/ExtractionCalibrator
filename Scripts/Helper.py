# -*- coding: utf-8 -*-
"""
Created on Tue May 16 11:24:22 2023
Helper functions to keep Model scripts organized
@author: leermdv
"""
import pandas as pd
import numpy as np
import nlmod
import matplotlib.pyplot as plt
import geopandas as gpd
import shapely
import flopy
import xarray as xr 
from tqdm import tqdm
import os


def GetExtent(df,Range):
    x = df.x_coordinaat
    y = df.y_coordinaat
    xmin = round( x.mean() - Range)
    xmax = round(x.mean() + Range)
    ymin = round(y.mean() - Range)
    ymax = round(y.mean() + Range)
    
    return [xmin, xmax, ymin, ymax]

def get_ss_discharge(Discharge):
    Discharge_mean = Discharge.mean(axis = 'rows').to_frame()
    Discharge_mean = Discharge_mean.transpose()
    Discharge_mean['index'] = Discharge.index[0]
    Discharge_mean.set_index('index', inplace = True)
    return Discharge_mean

def Spreidingslengte(extent,layerkh,layerc,layer_model):
    x = sum(extent[:2])/2.
    y = sum(extent[2:])/2.
    Point = layer_model.sel(x = x, y= y, method = 'nearest')
    K = Point.kh.sel(layer = layerkh).values
    print(f'K:{K}')
    D = Point.top.sel(layer = layerkh).values - Point.botm.sel(layer = layerkh).values
    print(f'D:{D}')
    c = (Point.top.sel(layer = layerc).values - Point.botm.sel(layer = layerc).values)/Point.kv.sel(layer = layerc).values
    print(f'c:{c}')
    return np.sqrt(K*D*c)

def FixTS_Discharge(TS, startdate, enddate):
        TS['Time'] = pd.to_datetime(TS['datumtijd'], utc =True)
        TS.set_index(['Time'], inplace = True)
        TS.index = TS.index.tz_convert(None)
        TS = TS.pivot(columns = 'putcode', values = 'debietmeter')
        
        TS = TS.resample('D').sum()
        TS = TS.loc[startdate:enddate]
        
        return TS

def FixTS_Discharge_warmup(TS, startdate):
        TS = TS.dropna()
        TS.loc[:,'n'] = list(range(len(TS)))
        TS['Time'] = pd.to_datetime(TS['DateTime'], utc =True)
        TS.set_index('Time', inplace = True)
        TS.index = TS.index.tz_convert(None)
        TS = TS.rename(columns = {'Name' : 'putcode', 'Volume' : 'debietmeter'})
        TS.replace('NA', 0, inplace = True)
        TS = TS.pivot(columns = 'putcode', values = 'debietmeter')
        TS = TS.loc[:startdate]
        return TS
    
def Add_warmup_wells(ExWells,ExWells_warmup):
    ExWells_warmup = ExWells_warmup.rename(columns = {'Name' : 'putcode', 
                                   'XCoordinate' : 'x_coordinaat',
                                   'YCoordinate' : 'y_coordinaat', 
                                   'WellTopLevel' : 'filter_referentiepunt_NAP', 
                                   'FilterTopLevel':'filter_bovenkant_ref', 
                                   'FilterBottomLevel' : 'filter_onderkant_ref'})
    ExWells_warmup.loc[:,'filter_bovenkant_ref'] = -  ExWells_warmup.loc[:,'filter_bovenkant_ref']
    ExWells_warmup.loc[:,'filter_onderkant_ref'] = -  ExWells_warmup.loc[:,'filter_onderkant_ref']                 
    for x in range(len(ExWells_warmup)):
        if ExWells_warmup['putcode'][x] not in ExWells.putcode.values:
            ExWells = pd.concat([ExWells,ExWells_warmup[ExWells_warmup['putcode'] == ExWells_warmup['putcode'][x]]])
    return ExWells
    
def FixTS_Obs(TS, startdate, enddate, ObsWells, Name):
        TS['filter_nummer'] = TS['filter_nummer'].astype(int)
        ObsWells['filter_nummer'] = ObsWells['filter_nummer'].astype(int)
        TS['putcode'] = TS['putcode'].astype(str) + '_' + TS['filter_nummer'].astype(str)
        ObsWells['putcode'] = ObsWells['putcode'].astype(str) + '_' + ObsWells['filter_nummer'].astype(str)
        TS = TS.merge(ObsWells[['putcode', 'filter_referentiepunt_NAP']], on = ['putcode'])
        # TS['meting_ref'] = TS['meting_NAP'] - TS['filter_referentiepunt_NAP']
        TS['Time'] = pd.to_datetime(TS['datumtijd'], utc =True)
        TS.set_index(['Time'], inplace = True)
        TS.index = TS.index.tz_convert(None)
        TS = TS.pivot(columns = 'putcode', values = 'meting_NAP')
        
        TS = TS.resample('D').mean()
        TS = TS.loc[startdate:enddate]
        if not os.path.isfile(os.path.join('..','Data','Preprocessed',f'stijghoogtereeksen_{Name}.csv')):
            TS.to_csv(os.path.join('..','Data','Preprocessed',f'stijghoogtereeksen_{Name}.csv'), index = True)
        return TS

def add_refresco(ExWells,Discharge):
    wells = pd.read_csv(os.path.join('..','Data','dawaco','winputten_refresco.csv'), delimiter= ';', index_col = False)
    ExWells = pd.concat([ExWells,wells], ignore_index = True)
    for well in wells.putcode.values:
        Discharge[well] = np.nan
        if well in ['PP1', 'PP3', 'PP4']:
            Discharge.loc['2014-12-31 00:00:00':'2018-12-31 00:00:00', well] = 383000/365/3
            Discharge.loc['2019', well] = 300000/365/3
            Discharge.loc['2019-12-31 00:00:00':, well] = 125000/365/3
        else:
            Discharge.loc['2014-12-31 00:00:00':'2016-12-31 00:00:00', well] = 117000/365/5
            Discharge.loc['2016-12-31 00:00:00':'2018-12-31 00:00:00', well] = 250000/365/5
            Discharge.loc['2019', well] = 450000/365/5
            Discharge.loc['2019-12-31 00:00:00':, well] = 875000/365/5
    Discharge.fillna(0, inplace = True)
    return ExWells, Discharge

def add_bier(ExWells,Discharge):
    wells = pd.read_csv(os.path.join('..','Data','dawaco','winputten_budel.csv'), delimiter= ';', index_col = False)
    ExWells = pd.concat([ExWells,wells], ignore_index=True)
    for well in wells.putcode.values:
        Discharge[well] = np.nan
        Discharge.loc[:, well] = 150000/365

    Discharge.fillna(0, inplace = True)
    return ExWells, Discharge
            
        
        
    
    
def make_gdf(ExWells,ObsWells):
    ExWellGdf = gpd.GeoDataFrame(ExWells, geometry = gpd.points_from_xy(ExWells['x_coordinaat'], ExWells['y_coordinaat']))
    ExWellGdf['type'] = 'Ex'
    ObsWellGdf = gpd.GeoDataFrame(ObsWells, geometry = gpd.points_from_xy(ObsWells['x_coordinaat'], ObsWells['y_coordinaat']))
    ObsWellGdf['type'] = 'Obs'
    WellGdf = gpd.GeoDataFrame(pd.concat([ExWellGdf,ObsWellGdf]))
    return WellGdf

def layermodel(extent, NLzuid,nlzuidpad = os.path.join('..','Data', 'NLZuidmodel.nc')):
    if NLzuid: 
        layer_model_full = xr.open_dataset(nlzuidpad)
        layer_model_sel = layer_model_full.sel(x = slice(extent[0], extent[1]), y = slice(extent[3], extent[2]))
        # layer_model = layer_model.load()
        layer_model = layer_model_sel.dropna(dim = 'layer', how = 'all')
        layer_model['meantop'] = layer_model.botm.mean(dim = ['x', 'y'], skipna = True)
        layer_model = layer_model.sortby('meantop', ascending = False)
        layer_model.attrs['extent'] = extent
        # layer_model.transpose('layer', 'y', 'x')
    else:
        layer_model = nlmod.read.regis.get_combined_layer_models(extent,use_geotop=False)
    return layer_model

def refiner(ds, refineranges, WellGdf):
    refinements = []
    for r,level in refineranges.items():
        s= gpd.GeoSeries(WellGdf[WellGdf['type'] == 'Ex'].buffer(r).unary_union)
        refinements.append((s.buffer(1,cap_style = 3), level))
    return refinements

def resample(ds, layer_model, NLzuid):
    # if not NLzuid:
    kv = nlmod.resample.structured_da_to_ds(layer_model.kv.sel(layer = ds.layer), ds, method='average')
    kh = nlmod.resample.structured_da_to_ds(layer_model.kh.sel(layer = ds.layer), ds, method='average')
    ds['kh'], ds['kv'] = nlmod.layers.get_kh_kv(kh, kv, anisotropy = 10)
    ds['top'] = nlmod.resample.structured_da_to_ds(layer_model.top.sel(layer = ds.layer), ds, method='average')
    ds['botm'] = nlmod.resample.structured_da_to_ds(layer_model.botm.sel(layer = ds.layer), ds, method='average')
    return ds
    

def steady_rch(ds,gwf):
    spd = []
    avg = ds.recharge.mean(dim = 'time')
    for cell in ds.recharge.icell2d.values:
        spd.append([cell,avg.sel(icell2d = cell).values])
        
    
def get_wel_spd(ds, gwf, Discharge, ExWells, steady_state):
    if steady_state:
        Discharge = get_ss_discharge(Discharge)

    Discharge.loc[Discharge.index[0]] = Discharge.mean() #mean steady state start
    CellIDtop = []
    CellIDbot = []
    welspd_dict = {}
    
    for t in tqdm(range(len(Discharge))):
        welspd = []
        for index,well in ExWells.iterrows():
                if t == 0:
                    CellIDtop.append(gwf.modelgrid.intersect(well['x_coordinaat'] ,well['y_coordinaat'], -(well['filter_bovenkant_ref'] - well['filter_referentiepunt_NAP'])))
                    CellIDbot.append(gwf.modelgrid.intersect(well['x_coordinaat'] ,well['y_coordinaat'], -(well['filter_onderkant_ref'] - well['filter_referentiepunt_NAP'])))
                kt, icell2d = CellIDtop[index]
                kb =  CellIDbot[index][0]
                if kt != kb:
                    potwlayers = np.arange(kt, kb)
                else:
                    potwlayers = [kt]
                wlayers = []
                for k in potwlayers:
                    if  (ds.idomain.sel(icell2d = icell2d, layer = ds.layer[k]) == 1)&(ds.kh.sel(icell2d = icell2d, layer = ds.layer[k]) > 1):
                        wlayers.append(k)
                        
                        
                for k in wlayers:
                    if (ds.idomain.sel(icell2d = icell2d, layer = ds.layer[k]) == 1):
                            wdata = [(k, icell2d), -Discharge.loc[Discharge.index[t],well['putcode']] / len(wlayers)]
                            welspd.append(wdata) #naar m³/d
                        
            
        welspd_dict[t] = welspd
    return welspd_dict


def get_LHM_heads(lhmpath, GHBextent, NLzuid):
    LHMclip = xr.open_dataset(lhmpath).sel(x = slice(GHBextent[0], GHBextent[1])).sel(y = slice(GHBextent[3], GHBextent[2]))
    if NLzuid:
        for x in range(8):
            LHMclip[f'head_mean_l{x+1}'].rio.write_nodata(np.nan, inplace = True)
            LHMclip[f'head_mean_l{x+1}'] = LHMclip[f'head_mean_l{x+1}'].rio.write_crs("epsg:28992")
            LHMclip[f'head_mean_l{x+1}'] = LHMclip[f'head_mean_l{x+1}'].rio.interpolate_na()
    return LHMclip

def get_LHM_heads2(lhmpath, ds):
    LHMclip = xr.open_dataset(lhmpath).sel(x = slice(ds.extent[0], ds.extent[1])).sel(y = slice(ds.extent[3], ds.extent[2]))

    return LHMclip

def ghb(ds, gwf,cachedir, NLzuid, GHBrange =1000, lhmpath = '..\\Data\\lhm.nc', delr = 100):
    #Get bigger REGIS extent
    GHBextent = [ds.extent[0] - GHBrange-200,ds.extent[1] + GHBrange+200,  ds.extent[2] - GHBrange-200,ds.extent[3] + GHBrange +200 ]
    CellExt = [ds.x.min().values,ds.x.max().values, ds.y.min().values, ds.y.max().values]
    layer_model = layermodel(GHBextent,NLzuid)
    #load LHM data
    
    LHM = get_LHM_heads2(lhmpath, ds)
    
    
    stage = np.empty((len(ds.layer),len(ds.icell2d)))
    cond = np.empty((len(ds.layer),len(ds.icell2d)))
    cellids = np.where((ds['edge_mask'] != 0) & (ds['idomain'] == 1) & (ds['kh'] > 1.))
    ghb_spd = []
    for i in tqdm(range(len(cellids[0]))):
            cell = ds.sel(layer = ds.layer[cellids[0][i]],icell2d = cellids[1][i])
            stage_cell = gethead2(cell, LHM, CellExt, GHBrange)
            cond_cell = getcond(cell, layer_model, CellExt,  GHBrange, delr)
            if np.isnan(stage_cell):
                stage_cell = stage[cellids[0][i],cellids[1][i-1]]
            if np.isnan(cond_cell):
                cond_cell = cond[cellids[0][i],cellids[1][i-1]]

            stage[cellids[0][i],cellids[1][i]] = stage_cell
            cond[cellids[0][i],cellids[1][i]] = cond_cell
            ghb_spd.append([(cellids[0][i], cellids[1][i]), stage_cell, cond_cell])
    
    ds['ghb_head'] = xr.DataArray(stage, coords = [ds.layer.values, ds.icell2d.values], dims = ['layer','icell2d'])
    ds['ghb_head'] = ds['ghb_head'].where(ds['edge_mask'] != 0)
    ds['ghb_cond'] = xr.DataArray(cond, coords = [ds.layer.values, ds.icell2d.values], dims = ['layer','icell2d'])
    ds['ghb_cond'] = ds['ghb_cond'].where(ds['edge_mask'] != 0)
    ghb = flopy.mf6.ModflowGwfghb(
        gwf,
        stress_period_data=ghb_spd,
        save_flows=True,)
    return ghb 
            
def gethead(cell, LHM, ext, GHBrange,delr):
    x,y = ghbXY(cell, ext, GHBrange,delr)
    Headcell = LHM.sel(x = x, y = y, method = 'nearest')
    for x in range(8):
        if cell.botm.values >= Headcell[f'wvp{x+1}-bot'].values:
            head = Headcell[f'head_mean_l{x+1}'].values
            break
        else:
            head = Headcell['head_mean_l8'].values
    return head

def gethead2(cell, LHM, ext, GHBrange):
    xmid = (ext[1] + ext[0])/2
    ymid = (ext[3] + ext[2])/2
    Headcell = LHM.sel(x = xmid, y = ymid, method = 'nearest')
    
    for x in range(8):
        if cell.botm.values >= Headcell[f'wvp{x+1}-bot'].values:
            lay = x+1 
            break
        else:
            lay = 8
    dy = LHM[f'head_mean_l{lay}'].diff('y').mean().values
    dx = LHM[f'head_mean_l{lay}'].diff('x').mean().values
    mean = LHM[f'head_mean_l{lay}'].mean().values
    if cell.x.values == ext[0]:
        head  = dx * ((ext[1] - ext[0])/2+GHBrange)/250 + mean - dy * (cell.y.values - ymid)/250 
    elif cell.x == ext[1]:
        head  = dx * -((ext[1] - ext[0])/2+GHBrange)/250+ mean - dy * (cell.y.values - ymid)/250 
    if cell.y == ext[2]:
        head  = dy * ((ext[3] - ext[2])/2+GHBrange)/250 + mean - dx * (cell.x.values - xmid)/250 
    elif cell.y == ext[3]:
        head  = dy * -((ext[3] - ext[2])/2+GHBrange)/250+ mean - dx * (cell.x.values - xmid)/250 
    return head


def ghbXY(cell, ext, GHBrange, delr):
    if cell.x.values == ext[0]:
        x  = cell.x.values-GHBrange
        y = cell.y.values
    elif cell.x == ext[1]:
        x = cell.x.values+GHBrange
        y = cell.y.values
    else:
        x = cell.x.values
        
    if cell.y == ext[2]:
        y  = cell.y.values-GHBrange
    elif cell.y == ext[3]:
        y = cell.y.values+GHBrange
    else:
        y = cell.y.values
    return x, y

def getcond(cell, layer_model, ext, GHBrange, delr):
    x,y = ghbXY(cell, ext, GHBrange,delr)
    if cell.x.values > x :
        line = layer_model.sel(y = y, method = 'nearest').sel(layer = cell.layer.values).sel(x=slice(x,cell.x.values))
    elif cell.x.values < x:
        line = layer_model.sel(y = y, method = 'nearest').sel(layer = cell.layer.values).sel(x=slice(cell.x.values,x))
    elif cell.y.values < y:
        line = layer_model.sel(x=x, method = 'nearest').sel(layer = cell.layer.values).sel(y = slice(y,cell.y.values))
    elif cell.y > y:
        line = layer_model.sel(x=x, method = 'nearest').sel(layer = cell.layer.values).sel(y = slice(cell.y.values,y))
    else:
        line = None
        print('error')
    logk = np.log(line['kh'].values)
    Kav  = np.exp(logk.mean())
    Dav = abs(line['top'].mean().values - line['botm'].mean().values)
    cond = Kav*Dav*delr/GHBrange
    return float(cond)

def riv(ds,gwf):
    bgt = nlmod.read.bgt.get_bgt(ds.extent)
    la = nlmod.gwf.surface_water.download_level_areas(bgt, extent=ds.extent, raise_exceptions=False)
    sw = nlmod.gwf.surface_water.add_stages_from_waterboards(bgt, la=la)
    # set stage for model to mean of summer and winter levels
    sw["stage"] = sw[["winter_stage", "summer_stage"]].mean(axis=1)
    
    # use a water depth of 0.5 meter
    sw["rbot"] = sw["stage"] - 0.5
    
    # we need to mask out the NaNs
    sw.drop(sw.index[sw["rbot"].isna()], inplace=True)
    # intersect surface water with grid
    if len(sw != 0):
        sfw_grid = nlmod.grid.gdf_to_grid(sw, gwf)
        # add bed resistance to calculate conductance
        bed_resistance = 1.0  # days
        sfw_grid["cond"] = sfw_grid.area / bed_resistance
        sfw_grid.set_index("cellid", inplace=True)
        
        # build stress period data for RIV package
        riv_spd = nlmod.gwf.surface_water.build_spd(sfw_grid, "RIV", ds)
        riv = flopy.mf6.ModflowGwfriv(gwf, stress_period_data=riv_spd)
    else:
        riv = None
    return riv

def ComparePlot(df1, df2, putcode, ObsWells):
    fig, ax = plt.subplots(len(putcode), 1)
    fig.set_dpi(300)
    fig.set_size_inches(10,20)
    if len(putcode)== 1:
        df1[putcode].plot.line(ax = ax)
        df2[putcode].plot.line(ax = ax)
        ax.legend(["Observed", "Modelled"])
    else:
        for x in range(len(putcode)):
            try:
                df1[putcode[x]].plot.line(ax = ax[x])
                df2[putcode[x]].plot.line(ax = ax[x])
                plt.text(-.01, .5, ObsWells[ObsWells['putcode'] == putcode[x]]['filter_bovenkant_ref'].values[0], ha='right', va='top', transform=ax[x].transAxes)
                ax[x].set_yticks([])
                if ax[x] != ax[-1]:
                    ax[x].set_xticks([])
            except KeyError:
                pass

            
        ax[0].legend(["Observed", "Modelled"])

def DischargePlot(Discharge):
    fig, ax = plt.subplots(len(Discharge.columns), 1)
    fig.set_dpi(300)
    fig.set_size_inches(10,20)

    for x in range(len(Discharge.columns)):
        Discharge[Discharge.columns[x]].plot.line(ax = ax[x])
    
        if ax[x] != ax[-1]:
            ax[x].set_xticks([])
                            

                    
def PlotPumpedHeads(ds,gwf,zoomedrange, tstep = 120):
    fig, ax = plt.subplots()
    head = nlmod.gwf.output.get_heads_da(ds)
    zoom = [ds.extent[0]+zoomedrange,ds.extent[1] - zoomedrange, ds.extent[2]+zoomedrange,ds.extent[3] - zoomedrange]

    pmv = flopy.plot.PlotMapView(modelgrid=gwf.modelgrid, extent = zoom, ax=ax)
    layer = ds.layer.values[gwf.wel.stress_period_data.array[0][0][0][0]]
    array = head.sel(layer=layer, time = tstep)
    pmv.plot_array(array, alpha=1, vmin = array.min().values, vmax = array.max().values)  
    # pmv.plot_bc('WEL', package = gwf.wel, color ='red', zorder = 1)              

def plot_mean_res(ObsHeads,ModHeads, WellGdf, steady_state, depth, tail):
    fig, ax = plt.subplots()
    fig.set_dpi(1200)
    depthcolumns = WellGdf[WellGdf['filter_bovenkant_ref'] < depth].putcode.values
    columns = list(set(depthcolumns).intersection(ModHeads.columns))
    columns = list(set(columns).intersection(ObsHeads.columns))
    if not steady_state:
        res = ObsHeads[columns] - ModHeads[columns]
        res = res.tail(n=tail)
    else:
        res = ObsHeads[columns].mean() - ModHeads[columns]

    res.mean().plot(kind = 'bar', ax = ax)
    labels = []
    for c in columns:
        labels.append(WellGdf[WellGdf['putcode'] == c].filter_bovenkant_ref.values[0])
    ax.set_xticklabels(labels)
    ax.set_title('mean residual per observation')
    ax.set_ylabel ('mean residual (meters)')
    
def plot_map(ds,gwf,array, layer):
    # ds_masked = ds.where(ds['var'] != -9999.) 
    # layer = ds.layer.values[gwf.wel.stress_period_data.array[0][0][0][0]]
    fig, ax = plt.subplots()
    ax = nlmod.plot.modelgrid(ds)
    layerarray = ds[array].sel(layer = layer)
    mapview = nlmod.plot.data_array(layerarray, ds = ds, ax = ax)
    cb = fig.colorbar(mapview, ax = ax)
    ax.set_title(array)

def GetHeadsAtObs(ds, ObsWells, gwf):
        df = pd.DataFrame(index = ds.time.data)
        head = nlmod.gwf.get_heads_da(ds)
        for index,well in ObsWells.iterrows():
                try:
                    CellIDbot = gwf.modelgrid.intersect(well['x_coordinaat'] ,well['y_coordinaat'], -(well['filter_onderkant_ref'] - well['filter_referentiepunt_NAP']))
                except Exception:
                    print(f'{well["putcode"]} outside model area')
                    continue
                if len(CellIDbot) == 3:
                    obsheads = head[:,CellIDbot[0], CellIDbot[1], CellIDbot[2]].values
                elif len(CellIDbot) == 2:
                    obsheads = head[:,CellIDbot[0], CellIDbot[1]]
                df[f'{well["putcode"]}'] = obsheads

        return df                   
                   
                   
                   
                   