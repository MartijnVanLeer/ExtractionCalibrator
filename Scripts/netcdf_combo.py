# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 12:03:51 2023

@author: leermdv
"""

import xarray as xr
import os

folder = r"\\Tsn.tno.nl\data\projects\060\5\52146\Werkdocumenten\03_grondwaterdynamiek\03_testlocatieZeeland\PhDopschaling\Brabantwater\netcdfs"

# ds = xr.open_dataset(f'{folder}\\NLZuidmodel.nc')

kh_da = xr.open_dataset(os.path.join(folder,'kh.nc'))
kh_da = kh_da.rename({'__xarray_dataarray_variable__' : 'kh'})

kv_da= xr.open_dataset(os.path.join(folder,'kv.nc'))
kv_da = kv_da.rename({'__xarray_dataarray_variable__' : 'kv'})

t_da= xr.open_dataset(os.path.join(folder,'t.nc'))
t_da = t_da.rename({'__xarray_dataarray_variable__' : 'top'})

b_da= xr.open_dataset(os.path.join(folder,'b.nc'))
b_da = b_da.rename({'__xarray_dataarray_variable__' : 'botm'})

# cw_da= xr.open_dataset(os.path.join(folder,'cw.nc')).astype('f8')
# cw_da = cw_da.rename({'__xarray_dataarray_variable__' : 'c'})

# d_da= xr.open_dataset(os.path.join(folder,'d.nc')).astype('f4')
# d_da = d_da.rename({'__xarray_dataarray_variable__' : 'd'})

# kd_da= xr.open_dataset(os.path.join(folder,'kd.nc')).astype('f8')
# kd_da = kd_da.rename({'__xarray_dataarray_variable__' : 'kd'})


das = [ t_da,b_da,kh_da, kv_da,
        # # cw_da,
        # d_da,#kd_da
        ] 


ds = xr.merge(das)

ds.to_netcdf(r"\\Tsn.tno.nl\data\projects\060\5\52146\Werkdocumenten\03_grondwaterdynamiek\03_testlocatieZeeland\PhDopschaling\Brabantwater\netcdfs\\NLZuidmodel.nc")
