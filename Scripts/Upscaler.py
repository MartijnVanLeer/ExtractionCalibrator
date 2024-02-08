# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 14:08:46 2023

@author: leermdv
"""

import xarray as xr
import glob
import os
import pandas as pd
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import matplotlib.cm as cmx

# path = r'C:\Users\leermdv\OneDrive - TNO\Documents\Python Scripts\Brabant water\Results\Vlijmen\Kfields'
path = r'C:\Users\leermdv\OneDrive - TNO\Documents\Python Scripts\Brabant water\Results\Vlijmen\KfieldsQC'
all_files = os.listdir(path)
kfiles = [i for i in all_files if i.startswith('K_')]

#%%
# for filename in kfiles:
filename = kfiles[-2]

df = pd.read_csv(os.path.join(path,filename), index_col=['x','y','z'])
df.drop('Unnamed: 0', axis = 1, )
ds = xr.Dataset.from_dataframe(df)
ds.K_1.isel(y = 30).T.plot(norm = colors.TwoSlopeNorm(vmin =  -7, vcenter = -1, vmax = 2))
plt.title(filename)

