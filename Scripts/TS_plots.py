# -*- coding: utf-8 -*-
"""
Created on Wed May  3 11:26:22 2023

@author: leermdv
"""

import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt

Location = 'Vlijmen'

ExWells = pd.read_csv(f'../Data/dawaco/winputten_WG_{Location}.csv')
ObsWells = pd.read_csv(f'../Data/dawaco/waarnemingsputten_WG_{Location}.csv')
Obs = pd.read_csv(f'../Data/dawaco/stijghoogtereeksen_WG_{Location}.csv')

def plotmap(ExWells,ObsWells, Location):
    Obsgdf = gpd.GeoDataFrame(ObsWells, geometry=gpd.points_from_xy(ObsWells.x_coordinaat, ObsWells.y_coordinaat))
    Exgdf = gpd.GeoDataFrame(ExWells, geometry=gpd.points_from_xy(ExWells.x_coordinaat, ExWells.y_coordinaat))
    fig, ax = plt.subplots()
    fig.set_dpi(500)
    Exgdf.plot(ax = ax, color = 'red', edgecolor=  'black')
    Obsgdf.plot(ax = ax, color = 'blue', marker = '+')
    ax.set_title(Location)

plotmap(ExWells, ObsWells, Location)

#%%
def plot_obs(ObsWells):
    fig,ax = plt.subplots()
    fig.set_dpi(500)
    ax.scatter(ObsWells.putcode, -ObsWells.filter_bovenkant_ref)
        
plot_obs(ExWells)
