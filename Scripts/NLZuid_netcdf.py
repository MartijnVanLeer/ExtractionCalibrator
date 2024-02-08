# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 15:54:49 2023

@author: leermdv
"""

import xarray as xr
import rioxarray
import os 
from tqdm import tqdm
readfolder = r"\\Tsn.tno.nl\data\projects\060\5\52146\Werkdocumenten\03_grondwaterdynamiek\03_testlocatieZeeland\PhDopschaling\Brabantwater\ZNL_tiff_v6"

kh = []
kv = []
kd = []
t = [] 
b = []
cw = []
d = []
khl = []
kvl= []
kdl = []
tl = []
bl = []
cwl = []
dl = []
for filename in tqdm(os.listdir(readfolder)):
    if filename.endswith('tif'):
        filename = filename.replace('.tif', '')
        if filename.endswith('kh'):
            kh.append(rioxarray.open_rasterio(os.path.join(readfolder, filename + '.tif')))
            khl.append(filename.split('_')[1])
        if filename.endswith('kv'):
            kv.append(rioxarray.open_rasterio(os.path.join(readfolder, filename + '.tif')))
            kvl.append(filename.split('_')[1])
        if filename.endswith('t'):
            t.append(rioxarray.open_rasterio(os.path.join(readfolder, filename + '.tif')))
            tl.append(filename.split('_')[1])
        if filename.endswith('b'):
            b.append(rioxarray.open_rasterio(os.path.join(readfolder, filename + '.tif')))
            bl.append(filename.split('_')[1])
        # if filename.endswith('cw'):
        #     cw.append(rioxarray.open_rasterio(os.path.join(readfolder, filename + '.tif')))
        #     cwl.append(filename.split('_')[1])
        # if filename.endswith('d'):
        #     if filename.endswith('kd'):
        #         kd.append(rioxarray.open_rasterio(os.path.join(readfolder, filename + '.tif')))
        #         kdl.append(filename.split('_')[1])
        #     else:
        #         d.append(rioxarray.open_rasterio(os.path.join(readfolder, filename + '.tif')))   
        #         dl.append(filename.split('_')[1])

                
folder = r'\\Tsn.tno.nl\data\projects\060\5\52146\Werkdocumenten\03_grondwaterdynamiek\03_testlocatieZeeland\PhDopschaling\Brabantwater\netcdfs'
names = ['kh', 'kv', 't', 'b',
         # 'cw', 'd', 'kd'
         ]
layernames= [khl, kvl, tl, bl]
for idx, x in tqdm(enumerate([kh,kv,t,b,cw,d, kd])):
    da = xr.concat(x, dim = 'layer')
    da = da.squeeze()
    da = da.drop_vars(['band','spatial_ref'])
    da = da.assign_coords({'layer' : layernames[idx]})
    da = da.where(da != da.max())
    da.to_netcdf(f'{folder}\\{names[idx]}.nc')
