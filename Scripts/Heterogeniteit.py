# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 16:14:09 2023

@author: leermdv
"""

import pandas as pd
import glob 
import statistics
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import scipy.stats
from tqdm import tqdm
from nlmod.dims.grid import xyz_to_cid
from nlmod.dims.grid import xy_to_icell2d
import io
import os

class boringen():
    def __init__(self, folder,ds):
        self.get_boringen(folder,ds)
        self.selection = None
    def get_boringen(self,folder, ds):
        csv_files = glob.glob(folder + "/*.txt")
        self.boreholes = {}
        metadata = []
        for file in tqdm(csv_files, 'Read borehole files..'):
            base = os.path.basename(file)
            name = os.path.splitext(base)[0]
            #metadata
            with open(file) as f:
                text = f.read()
            sections = text.split('LITHOLOGIE LAGEN')
            mdtext = io.StringIO(sections[0])
            md = pd.read_csv(mdtext,sep = ': \t', engine='python',)
            metadata.append(md)
            #boreholes
            bhtext = io.StringIO(sections[1].split('LITHOLOGIE SUBLAGEN')[0])
            bh = pd.read_csv(bhtext, sep = '\t', skiprows = [0], index_col = False, engine='python')
            bh['Bovenkant laag (m beneden maaiveld)'] = bh['Bovenkant laag (m beneden maaiveld)'].round()
            bh['Onderkant laag (m beneden maaiveld)'] = bh['Onderkant laag (m beneden maaiveld)'].round()
            result = []
            for ix, row in bh.iterrows():
                top = row['Bovenkant laag (m beneden maaiveld)']
                bot = row['Onderkant laag (m beneden maaiveld)']
                for range_ in range(int(top), int(bot)):
                    row['z'] = -range_ + round(float(md.loc['Maaiveldhoogte (meter t.o.v. NAP)'].values[0]))
                    row['x'] = int(md.loc['X-coordinaat (m)'].values[0])
                    row['y'] = int(md.loc['Y-coordinaat (m)'].values[0])
                    row['quality'] = md.loc['Kwaliteitcode beschrijving lithologie', 'ALGEMENE GEGEVENS BORING']
                    row['boringnr'] = ix
                    result.append(row.copy())
            self.boreholes[name] = pd.DataFrame(result)
        self.metadata = pd.concat(metadata, axis = 1)
        self.metadata.columns = self.metadata.iloc[0]

    def add_layers(self, layermodel):
        self.layermodel = layermodel
        for boringnr, df in tqdm(self.boreholes.items(), 'adding layers..'):
            cell = self.layermodel.sel(x = df['x'].values[0], method = 'nearest').sel(y = df['y'].values[0], method = 'nearest')
            dummy = []
            for idx, row in df.iterrows():
                layer = cell.layer.where((cell.top >= row['z']) & (cell.botm <= row['z']), drop = True,).layer.values
                if len(layer) == 0:
                    layer = None
                else:
                    layer = layer[0]
                row['layer'] = layer
                dummy.append(row.copy())
            self.boreholes[boringnr] = pd.DataFrame(dummy)

    def select_layer(self, layer):
        self.selection = self.boreholes
        for boringnr, df in self.boreholes.items():
            df = df[df['layer'] == layer]
            self.selection[boringnr] = df

    def select_range(self, layer, ds):
        self.selection = self.boreholes
        idx = list(ds.layer).index(layer)
        thickness = ds.isel(layer = idx-1).botm-ds.sel(layer = layer).botm
        condrange = np.ceil(thickness.max().values)
        ds['Layermid'] = (ds.isel(layer = idx-1).botm + ds.sel(layer = layer).botm)/2
        print(len(self.boreholes.items()))
        for boringnr, df in self.boreholes.items():
            tdf = df.copy()
            cellid = xy_to_icell2d
            mid = ds.sel(icell2d = xy_to_icell2d((tdf['x'].values[0],tdf['y'].values[0]), ds)).Layermid.values
            minrange = mid - 0.5 * condrange
            maxrange = mid + 0.5 * condrange
            tdf = tdf[tdf['z'] >= minrange]
            tdf = tdf[tdf['z'] <= maxrange]
            tdf['z'] = -np.arange(len(tdf))
            self.selection[boringnr] = tdf
    
    def listify(self):        
        if self.selection != None:
            self.list = pd.concat(self.selection.values())
        else:
            self.list = pd.concat(self.boreholes.values())
        self.list.reset_index(inplace = True)
        lith = []
        for index, row in self.list.iterrows():
            if row.Hoofdgrondsoort == 'klei':
                if row['Bijmenging silt'] in ['zwak siltig', 'matig siltig','---']:
                    if row['Bijmenging zand'] != '---':
                        lith.append('kz')
                    else:
                        lith.append('k')
                else:
                    lith.append('kz')
            elif row.Hoofdgrondsoort == 'zand':
                if row['Bijmenging klei'] != '---':
                    lith.append('kz')
                else:
                    if row.Zandmediaanklasse in ['uiterst grof (O)','grove  categorie (O)','zeer grof (O)','uiterst grof', 'zeer grof']:
                        lith.append('zg')
                    elif row.Zandmediaanklasse in ['matig grof (O)','matig fijn (O)','matig grof', 'zandmediaan onduidelijk', '---']:
                        lith.append('zm')
                    else:
                        lith.append('zf')
            elif row.Hoofdgrondsoort == 'leem':
                lith.append('kz')
            elif row.Hoofdgrondsoort in ['veen', 'bruinkool']:
                lith.append('v')
            elif row.Hoofdgrondsoort == 'grind':
                lith.append('zg')
        self.list['Lithoclass'] = lith
        return self.list
    
    def indicators(self, g1):
        self.group = g1
        self.list['i (no weight)'] = np.where(self.list.Lithoclass.isin(g1),1,0)
        indicator = []
        for index, row in self.list.iterrows():
            
            if row.quality == 'A':
                if row.Lithoclass in g1:
                    indicator.append(1)
                else:
                    indicator.append(0)
            if row.quality == 'B':
                if row.Lithoclass in g1:
                    indicator.append(0.9)
                else:
                    indicator.append(0.1)
            if row.quality == 'C':
                if row.Lithoclass in g1:
                    indicator.append(0.8)
                else:
                    indicator.append(0.2)
            if row.quality == 'D':
                if row.Lithoclass in g1:
                    indicator.append(0.7)
                else:
                    indicator.append(0.3)
            if row.quality == 'E':
                if row.Lithoclass in g1:
                    indicator.append(0.6)
                else:
                    indicator.append(0.4)
        self.list['i']  = indicator
        return self.list
    
    def plot_lith_dist(self):
        fig, ax = plt.subplots()
        fig.set_dpi = 300
        # sns.histplot(self.list, x = 'Hoofdgrondsoort',ax = ax)
        sns.histplot(self.list, x = 'Lithoclass',ax = ax)
        ax.set_title('Local lithoclass distribution')
        
    def plot_K_weighted(self, Kcore, TI):
        rn = []
        lth = []
        rng = np.random.default_rng()
        for col in self.list.Lithoclass.unique():
            count = self.list.Lithoclass[self.list.Lithoclass == col].count()
            if TI:
                sel = Kcore.dist[Kcore.dist['LTH_klass_toegekend'] == col]['logk']
                draw = rng.normal(sel.mean(), sel.std(), count)
            else:
                sel = Kcore.dist[Kcore.dist['LITHO_CLASS_CD'] == col]
                draw = pert(sel.MIN_KV.values, sel.MEAN_KV.values,sel.MAX_KV.values, count)
            for x in draw:
                rn.append(x)
                lth.append(col)
        self.Kweighted = pd.DataFrame({'Lithoclass' : lth, 'K':rn})
        fig,ax = plt.subplots()
        fig.set_dpi(300)
        ax.set_title('Local K distribution')
        
        sns.histplot(rn, ax =ax, stat = 'density')
        x0,x1 = ax.get_xlim()
        x_pdf = np.linspace(x0,x1,100)
        
        self.mu1, self.std1 = scipy.stats.norm.fit(self.Kweighted[self.Kweighted['Lithoclass'].isin(self.group)]['K'])
        y_pdf1 = scipy.stats.norm.pdf(x_pdf, self.mu1, self.std1)
        ax.plot(x_pdf, y_pdf1, c = 'r')
        self.mu2, self.std2 = scipy.stats.norm.fit(self.Kweighted[self.Kweighted['Lithoclass'].isin(self.group) == False]['K'])
        y_pdf2 = scipy.stats.norm.pdf(x_pdf, self.mu2, self.std2)
        ax.plot(x_pdf, y_pdf2, c = 'green')
        return rn
    
    def add_k(self, res, ens_no):
        rng = np.random.default_rng()
        K1 = rng.normal(self.mu1,self.std1, len(res))
        K2 = rng.normal(self.mu2,self.std2, len(res))
        Kfield = res[['x','y','z','cellid']]
        for x in range(ens_no):
            Kfield.loc[:,f"K_{x+1}"] = np.where(res[f'sim{x+1}'] == 1, K1,K2)
            np.random.shuffle(K1)
            np.random.shuffle(K2)
        return Kfield
    
    

    
class KDist():
    def __init__(self,folder, layer):
        self.type = 'TI'
        stratls = [char for char in layer if char.isupper()]
        self.strat = "".join(stratls)
        self.full = pd.read_csv(folder)
        self.layer = layer
        self.dist = self.full[self.full['Strat'] == self.strat] 
        self.dist['logk'] = np.log10(self.dist['K (m/d 10C)'])
    def gm(self):
        return self.dist.logk.mean()
    def hm(self):
        return statistics.harmonic_mean(self.dist['K (m/d 10C)'])
    def am(self):
        return self.dist['K (m/d 10C)'].mean()
    def plot_dist(self):

        fp = sns.FacetGrid(self.dist, col = 'LTH_klass_toegekend', col_wrap = 3, col_order = ['v','k','kz','zf','zm','zg'],
                           xlim = (-7,2))
        fp.map_dataframe(sns.histplot, 'logk',binwidth=0.5,  stat = 'density',)
        fp.map_dataframe(annotate)
        fp.set_titles(col_template = '{col_name}')

        fp.figure.set_dpi(300)
        fp.figure.suptitle('Core K distribution')
        for ax in fp.axes.ravel():
            var = ax.get_title()
            data = self.dist[self.dist.LTH_klass_toegekend.eq(var)]
            mu, std = scipy.stats.norm.fit(data['logk'])
            x0,x1 = ax.get_xlim()
            x_pdf = np.linspace(x0,x1,100)
            y_pdf = scipy.stats.norm.pdf(x_pdf, mu, std)
            ax.plot(x_pdf, y_pdf, c = 'r')
        fp.savefig()

class Kreg():
    def __init__(self,folder, layer):
        self.type = 'REG'
        stratls = [char for char in layer if char.isupper()]
        self.strat = "".join(stratls)
        self.full = pd.read_csv(folder, delimiter = ';')
        self.layer = layer
        dist = self.full[self.full['UNIT_CD'] == self.strat]
        self.dist = dist[dist['HYD_VERSION'] == 'v02r2']
        self.dist = self.dist[['LITHO_CLASS_CD','MIN_KV','MEAN_KV','MAX_KV']]

        for col in ['MIN_KV','MEAN_KV','MAX_KV']:
            self.dist[col] = np.log10(self.dist[col])
        self.dist['minfac'] = self.dist.MEAN_KV/self.dist.MIN_KV
        self.dist['maxfac'] = self.dist.MAX_KV/self.dist.MEAN_KV
        self.dist.replace('kk', 'k', inplace = True)
        self.dist.replace('b', 'v', inplace = True)

def pert(mini, mean, maxi, size=1, lamb=4):
    r = maxi - mini
    alpha = 1 + lamb * (mean - mini) / r
    beta = 1 + lamb * (maxi - mean) / r
    return mini + np.random.beta(alpha, beta, size=size) * r
    
def annotate(data, **kws):
    n = len(data)
    ax = plt.gca()
    ax.text(.1, .9, f"N = {n}", transform=ax.transAxes, bbox = dict(facecolor = 'white', edgecolor ='black'))     
    
def trim(Kfields,ds, Layer):
    keep = []
    cellids = [] 
    for index, row in Kfields.iterrows():
        celllayer, cellid = xyz_to_cid((row.x,row.y, row.z), ds)
        cellids.append(cellid)
        if ds.layer[celllayer] == Layer:
            keep.append(True)
        else: 
            keep.append(False)
    Kfields['keep'] = keep
    Kfields['cellid'] = cellids
    return Kfields[Kfields.keep==True] 

def add_cellid(Kfields,ds):
    cellids = [] 
    for index, row in Kfields.iterrows():
        cellid = xy_to_icell2d((row.x,row.y), ds)
        cellids.append(cellid)
    Kfields['cellid'] = cellids
    return Kfields
    
    
    
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
