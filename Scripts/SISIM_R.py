# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 12:41:44 2022

@author: Martijn van Leer
"""

import os
import platform

if platform.platform().startswith('Win'):
    os.environ['R_HOME'] = r'C:\ProgramData\anaconda3\pkgs\r-base-4.1.3-hdca333a_12\lib\R'
    os.environ["PATH"]   =  r'C:\ProgramData\anaconda3\pkgs\r-base-4.1.3-hdca333a_12\lib\R\bin\x64' + ";" + os.environ["PATH"]
else: 
   os.environ['R_HOME'] = '/home/4120973/.conda/pkgs/r-base-4.3.2-hb8ee39d_2/lib/R'
   os.environ["PATH"]   =  '/home/4120973/.conda/pkgs/r-base-4.3.2-hb8ee39d_2/lib/R/bin/Rscript' + ";" + os.environ["PATH"] 
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
gstat = importr('gstat')
gstat = importr('sp')
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None
from rpy2.robjects import pandas2ri
pandas2ri.activate()


def SISIM(Lx=100,Ly=100,Lz =0.5,dx =2,dy =2,dz =0.05, xcorlen =40, zcorlen = 0.1, ens_no = 50, frac =0.2, nmax = 20, seed = 1337):
    """
    Parameters
    ----------
    Lx : float, optional
        x-length of field. The default is 100.
    Ly : float,, optional
        y-length of field. The default is 100.
    Lz : float,, optional
        z-length of field. The default is 0.5.
    dx : float, optional
        cell size in x direction. The default is 2.
    dy : float,, optional
        cell size in y direction. The default is 2.
    dz : float,, optional
        cell size in z direction.. The default is 0.05.
    xcorlen : float, optional
        Horizontal correlation length (range). The default is 40.
    zcorlen : float, optional
        Vertical correlation length (range) The default is 0.1.
    ens_no : int, optional
        Number of realizations. The default is 50.
    frac : int, optional
        Fraction of binary distribution. The default is 0.2.
    nmax : int, optional
        Number of close points gstat takes into account. The default is 20.

    Returns
    -------
    binfieldList : list
        List containing all generated field as numpy arrays.

    """
    robjects.r('''
               # r code for SISIM function
               SIS <- function (Lx,Ly,Lz,dx,dy,dz, xcorlen, zcorlen, ens_no, frac, nmax=20, seed){
               set.seed(seed)
               xyz <- expand.grid(seq(dx,Lx,dx), seq(dy,Ly,dy), seq(dz,Lz,dz))
               names(xyz) <- c("x","y", "z")
               g <- gstat(NULL, 'var1', formula = z~1, locations = ~x+y+z, dummy = TRUE, beta = frac,model = vgm(1,"Sph",xcorlen, anis = c(0,0,0,1,zcorlen/xcorlen)),nmin = 1,nmax = nmax, maxdist = as.numeric(xcorlen), force = TRUE)
               predict(g, newdata = xyz, nsim = ens_no, indicators = TRUE)
               }
               ''')
    SIS = robjects.r['SIS'] 
    res = SIS(Lx,Ly,Lz,dx,dy,dz, xcorlen, zcorlen, ens_no, frac, nmax, seed = (xcorlen*100+zcorlen)*1337) #run R function
    binfieldList = []
    for idx in range(ens_no):
        binfield = np.array(res[idx+3])#(res[f'sim{idx+1}']) # convert pd.dataframe to np.array, 
        binfieldList.append(binfield.reshape(int(Lz/dz),int(Lx/dx),int(Ly/dy)).T) #reshape list to 3D box
    return binfieldList

def SGSIM(Lx=100,Ly=100,Lz =0.5,dx =2,dy =2,dz =0.05, xcorlen =40, zcorlen = 0.1, ens_no = 50, LogMean = 0, LogVar = 1, nmax = 20):
    """
    Parameters
    ----------
    Lx : float, optional
        x-length of field. The default is 100.
    Ly : float,, optional
        y-length of field. The default is 100.
    Lz : float,, optional
        z-length of field. The default is 0.5.
    dx : float, optional
        cell size in x direction. The default is 2.
    dy : float,, optional
        cell size in y direction. The default is 2.
    dz : float,, optional
        cell size in z direction.. The default is 0.05.
    xcorlen : float, optional
        Horizontal correlation length (range). The default is 40.
    zcorlen : float, optional
        Vertical correlation length (range) The default is 0.1.
    ens_no : int, optional
        Number of realizations. The default is 50.
    LogMean : float, optional
        mean of the logtransformed K distribution. The default is 0.
    LogVar : float, optional
        Variance of the log distribution of K. The default is 1.
    nmax : TYPE, optional
        DESCRIPTION. The default is 20.

    Returns
    -------
    GausfieldList : List
        List containing all generated field as numpy arrays.

    """
    xcorlen = 0.5*xcorlen
    zcorlen = 0.5*zcorlen
    robjects.r('''
               # r code for SGSIM function
               SGS <- function (Lx,Ly,Lz,dx,dy,dz, xcorlen, zcorlen, ens_no, nmax,LogMean, LogVar){
               xyz <- expand.grid(seq(dx,Lx,dx), seq(dy,Ly,dy), seq(dz,Lz,dz))
               names(xyz) <- c("x","y", "z")
               g <- gstat(NULL, 'var1', formula = z~1, locations = ~x+y+z, dummy = TRUE, beta = LogMean,model = vgm(LogVar,"Sph",xcorlen, anis = c(0,0,0,1,zcorlen/xcorlen)), nmax = nmax)
               predict(g, newdata = xyz, nsim = ens_no)
               }
               ''')
    SIS = robjects.r['SGS'] 
    res = SIS(Lx,Ly,Lz,dx,dy,dz, xcorlen, zcorlen, ens_no, nmax, LogMean, LogVar) #run R function
    GausfieldList = []
    for idx in range(ens_no):
        binfield = np.array(res[idx+3])#(res[f'sim{idx+1}']) # convert pd.dataframe to np.array, 
        GausfieldList.append(np.e**(binfield.reshape(int(Lz/dz),int(Lx/dx),int(Ly/dy)).T)) #reshape list to 3D box
    return GausfieldList

def Cond_SISIM(borelogs_grid_df_r,xmin,ymin,zmin,Lx=100,Ly=100,Lz =0.5,dx =2,dy =2,dz =0.05, xcorlen =40, zcorlen = 0.1, ens_no = 50, frac =0.2, nmax = 20, seed = 1337):
    """
  

    """
    robjects.r('''
               # r code for SGSIM function
                    SIS <- function (borelogs_grid_df_r,xmin,ymin,zmin,Lx,Ly,Lz,dx,dy,dz, xcorlen, zcorlen, ens_no, frac, nmax,seed){
                    set.seed(seed)
                    xyz_sel <- borelogs_grid_df_r
                    coordinates(xyz_sel) = ~x+y+z
                                   
                    xyz <- expand.grid(seq(dx/2,Lx,dx),seq(dy/2,Ly,dy),seq(dz/2,Lz,dz))
                    names(xyz) <- c("x","y","z")
                    xyz[, 1] = xyz[, 1] + xmin
                    xyz[, 2] = xyz[, 2] + ymin
                    xyz[, 3] = xyz[, 3] + zmin
                    coordinates(xyz) = ~x+y+z
                

                    g <- gstat(formula = i~1,id = "i", data = xyz_sel, locations = xyz_sel,dummy = FALSE,beta = frac, model = vgm(frac,"Sph",xcorlen, anis = c(0,0,0,1,zcorlen/xcorlen)), maxdist = as.numeric(xcorlen), nmax = nmax,vdist = TRUE)
                    sdf <- predict(g, newdata = xyz, nsim = ens_no, indicators = TRUE, debug.level = -1)
                    as.data.frame(sdf)
                    }          
'''
               )
    SIS = robjects.r['SIS'] 
    res_r = SIS(borelogs_grid_df_r,xmin,ymin,zmin,Lx,Ly,Lz,dx,dy,dz, xcorlen, zcorlen, ens_no, frac, nmax,seed)
    # with robjects.default_converter + pandas2ri.converter:
    res = robjects.conversion.get_conversion().rpy2py(res_r)
    return res

def Cond_SGSIM(borelogs_grid_df_r,LogMean, LogVar, xmin,ymin,zmin,Lx=100,Ly=100,Lz =0.5,dx =2,dy =2,dz =0.05, xcorlen =40, zcorlen = 0.1, ens_no = 50, nmax = 20, seed = 1337):
    """
  

    """
    xcorlen = 0.5*xcorlen
    zcorlen = 0.5*zcorlen
    robjects.r('''
               # r code for SGSIM function
                    SGS <- function (borelogs_grid_df_r,xmin,ymin,zmin,Lx,Ly,Lz,dx,dy,dz, xcorlen, zcorlen, ens_no, nmax,seed,LogMean, LogVar){
                    set.seed(seed)
                    xyz_sel <- borelogs_grid_df_r
                    coordinates(xyz_sel) = ~x+y+z
                                   
                    xyz <- expand.grid(seq(dx/2,Lx,dx),seq(dy/2,Ly,dy),seq(dz/2,Lz,dz))
                    names(xyz) <- c("x","y","z")
                    xyz[, 1] = xyz[, 1] + xmin
                    xyz[, 2] = xyz[, 2] + ymin
                    xyz[, 3] = xyz[, 3] + zmin
                    coordinates(xyz) = ~x+y+z
                                   
                    g <- gstat(NULL, 'i',formula = i~1, data = xyz_sel,dummy = FALSE,model = vgm(LogVar,"Sph",xcorlen, anis = c(0,0,0,1,zcorlen/xcorlen)), nmax = nmax)
                    sdf <- predict(g, newdata = xyz, nsim = ens_no)
                    as.data.frame(sdf)
                    }          
'''
               )
    SGS = robjects.r['SGS'] 
    res_r = SGS(borelogs_grid_df_r,xmin,ymin,zmin,Lx,Ly,Lz,dx,dy,dz, xcorlen, zcorlen, ens_no, nmax,seed,LogMean, LogVar)
    # with robjects.default_converter + pandas2ri.converter:
    res = robjects.conversion.get_conversion().rpy2py(res_r)
    GausfieldList = []
    for idx in range(ens_no):
        binfield = np.array(res[f'sim{idx+1}']) # convert pd.dataframe to np.array, 
        GausfieldList.append((binfield.reshape(int(Lz/dz),int(Lx/dx),int(Ly/dy)).T)) #reshape list to 3D box
    return GausfieldList,res