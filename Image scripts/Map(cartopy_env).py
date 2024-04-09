# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 11:37:39 2023

@author: leermdv
"""
#%%

import geopandas as gpd
import matplotlib.pyplot as plt
import cartopy
import cartopy.mpl.geoaxes
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
from matplotlib_scalebar.scalebar import ScaleBar
from shapely.geometry.point import Point
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import pandas as pd

#default settings for text and lines
textsize = 8
linewidth = 0.3
plt.rcParams.update({'font.size': textsize,
                     'axes.linewidth':linewidth,
                     'xtick.major.width':linewidth,
                     'ytick.major.width':linewidth,
                     'xtick.minor.size':linewidth,
                     'xtick.minor.size':linewidth,
                     'patch.linewidth' :linewidth
                     })

#init fig
fig = plt.figure()

cm = 1/2.54 #unit converter


#load files
gdf = gpd.read_file(r"C:\Users\leermdv\OneDrive - TNO\Documents\Python Scripts\Schouwen Scripts\Image scripts\Landsgrenzen\CountryBorders.shp")
Peilshp =  r"C:\Users\leermdv\OneDrive - TNO\Documents\Python Scripts\ExtractionCalibrator\ObsForCalibration_Budel.csv" #location shapefile with RD coordinates
Peildf = pd.read_csv(Peilshp )
Peildf = gpd.GeoDataFrame(Peildf, geometry = gpd.points_from_xy(Peildf.x_coordinaat, Peildf.y_coordinaat))
Loc = Peildf.dissolve().centroid #select single row from file

#load terrain background img
class StadiaStamen(cimgt.Stamen):
    def _image_url(self, tile):
         x,y,z = tile
         url = f"https://tiles.stadiamaps.com/tiles/stamen_terrain_background/{z}/{x}/{y}.png"
         return url
     
stamen_terrain = StadiaStamen('terrain-background')
#convert borders to same crs as terrain map 
crs = stamen_terrain.crs
gdf = gdf.to_crs(crs.proj4_init)

ax = fig.add_subplot(1,1,1, projection = crs)
fig.set_size_inches(8.6*cm, 9*cm)
fig.set_dpi(1200)

def plot_loc(ax, crs, Loc):
    Location = Loc.set_crs('epsg:28992')
    Location = Location.to_crs(crs.proj4_init)
    Location.plot(ax=ax, color = 'red', edgecolor = 'black', label = 'Extraction site', zorder = 10)
    return ax, Location

def plot_scalebar(ax):
    points = gpd.GeoSeries([Point(3, 51), Point(4, 51)], crs=crs)  # Geographic WGS 84 - degrees
    points = points.to_crs(28992) # Projected RD in meters
    distance_meters = points[0].distance(points[1])
    ax.add_artist(ScaleBar(distance_meters, location = 'lower left', frameon = False ))
    x, y, arrow_length = 0.96, 0.98, 0.05
    ax.annotate('N', xy=(x, y), xytext=(x, y-arrow_length-0.01),
                arrowprops=dict(facecolor='black', width=0, headwidth=7, headlength = 7),
                ha='center', va='center', fontsize=8,
                xycoords=ax.transAxes)
    return ax

def plot_scalebar_s(ax):
    points = gpd.GeoSeries([Point(3, 51), Point(4, 51)], crs=crs)  # Geographic WGS 84 - degrees
    points = points.to_crs(28992) # Projected RD in meters
    distance_meters = points[0].distance(points[1])
    ax.add_artist(ScaleBar(distance_meters, location = 'lower left', frameon = False ))

    return ax
#terrain background map
ax.set_extent([2.0, 7.5, 50.5, 54], crs=ccrs.Geodetic())
ax.add_feature(cartopy.feature.LAND, edgecolor='black')
# ax.add_image(stamen_terrain, 6)
ax.add_feature(cartopy.feature.LAKES, edgecolor='black')
ax.add_feature(cartopy.feature.RIVERS)
ax.add_feature(cartopy.feature.OCEAN)

#plot location from gdf
ax, Location = plot_loc(ax, crs, Loc)

#plot scalebar
ax = plot_scalebar(ax)

#plot borders
gdf.plot(ax = ax, zorder = 3, color = 'black', linestyle = '--', linewidth = 1, alpha = 0.6)

#init inset with settings
# axins = inset_axes(ax, width="40%", height="35%", loc='upper left',axes_class=cartopy.mpl.geoaxes.GeoAxes, 
#                   axes_kwargs=dict(projection=crs))
#inset background map (openstreetmap)
# axins.add_image(cimgt.OSM(),10)
# axins.set_extent([3.76,4.0,51.66,51.79], crs=ccrs.Geodetic()) #extent of inset map (latitude/longitude coordinates)

#Plot location in inset
# axins, Location = plot_loc(axins, crs, Loc)
# axins = plot_scalebar_s(axins)
#add country names (locations in percentage of axis)
ax.annotate('The Netherlands', xy = (0.45,0.45), xycoords = 'axes fraction', zorder =10,  fontweight = 'light')
ax.annotate('Germany', xy = (0.8,0.3), xycoords = 'axes fraction', zorder = 11, fontweight = 'light')
ax.annotate('Belgium', xy = (0.4,0.1), xycoords = 'axes fraction', zorder = 12, fontweight = 'light')

#add legend
ax.legend(loc = 'lower right', fancybox = True, framealpha = 0.8, title = '$\\bf{Legend}$')

#settings for marking the inset location
# box, c1, c2 = mark_inset(ax, axins, loc1=3, loc2=4, fc="none", ec="0.0") 

fig.tight_layout()
# fig.savefig(r'..\Images\Figure1.pdf')
# fig.savefig(r'..\Images\Figure1.png')
