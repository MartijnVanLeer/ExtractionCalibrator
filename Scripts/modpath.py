# -*- coding: utf-8 -*-
"""
Created on Tue May  7 09:29:54 2024

@author: leermdv
"""

import nlmod
import xarray as xr 
import flopy 
import os
import sys 
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


ws = r"C:\Users\leermdv\OneDrive - TNO\Documents\Python Scripts\ExtractionCalibrator\Data\fullrun_ss"
ds = xr.open_dataset(os.path.join(ws, 'fullrun_ss.nc'))
ds.attrs['model_ws'] = r"C:\Users\leermdv\OneDrive - TNO\Documents\Python Scripts\ExtractionCalibrator\Data\fullrun_ss"
exe_name = os.path.join(os.path.dirname(nlmod.__file__), "bin", "mf6")

if sys.platform.startswith("win"):
    exe_name += ".exe"
ds.attrs['exe_name'] = exe_name
sim = flopy.mf6.mfsimulation.MFSimulation.load('mfsim', sim_ws =ws, exe_name =ds.exe_name)


gwf = sim.get_model()
npf = gwf.get_package('NPF')
npf.save_flows = True
npf.write()

nlmod.sim.write_and_run(sim, ds, write_ds = False)
#%%
layno = list(ds.layer).index('KIz2')
mpf = nlmod.modpath.mpf(gwf)
_mpfbas = nlmod.modpath.bas(mpf,)
nodes = nlmod.modpath.layer_to_nodes(mpf, layno)


# nodes = np.array(nodes)[ds.icell2d.where(ds.y == ds.y.values[369], drop = True).values.astype(int)]
nodes = nodes
# pg = nlmod.modpath.pg_from_fdt(nodes, divisions = 3)
pg = nlmod.modpath.pg_from_pd(nodes, localx=0.5, localy=0.5, localz=0.0)
mpsim = nlmod.modpath.sim(mpf, pg, "forward")
nlmod.modpath.write_and_run(mpf)

#%%
pdata = nlmod.modpath.load_pathline_data(mpf)


fpth = os.path.join(ws,'modpath', "mp7_fullrun_ss.mpend")
e = flopy.utils.EndpointFile(fpth)

welnodes = nlmod.modpath.package_to_nodes(gwf, 'WEL', mpf)
well_epd = e.get_destination_endpoint_data(dest_cells=welnodes)


#%%


f = plt.figure(dpi = 1000)
ax = f.add_subplot(111, projection='3d')
ax.set_xlim((ds.extent[0], ds.extent[1]))
ax.set_ylim((ds.extent[2], ds.extent[3]))

for pid in well_epd["particleid"]:
    pf = pdata[pdata["particleid"] == pid]
    ax.plot(pf["x"],pf["y"], pf["z"], linewidth = 0.2, linestyle = '-')
    # ax.scatter(pf["x"],pf["y"], pf["z"], s= 0.1, c  = pf['time'])
# ax.view_init(10, 30, 0)
# ax.plot(pf["x"], pf["z"], color="k", linewidth=0.5, label="pathline")
#%%
f = plt.figure(dpi = 300)
ax = f.add_subplot()
mm = flopy.plot.PlotMapView(modelgrid=gwf.modelgrid, ax=ax)


mm.plot_endpoint(well_epd, direction="starting", colorbar=True, s =0.5, cmap = 'hot', zorder =10)

mm.plot_grid(alpha = 0.2)
zr = 2000


#%%
f = plt.figure(dpi = 300)
ax = f.add_subplot()
sns.histplot(well_epd.time/365, ax = ax, bins = 40)
ax.set_title('Particles released in centre of aquifer')
ax.set_xlim(0, 350000/365)
ax.set_xlabel('Travel time (years)')

#%%
f, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 4))

for i, pid in enumerate(well_epd["particleid"]):
    pf = pdata[pdata["particleid"] == pid]
    x0, y0, z0 = pf[["x", "y", "z"]][0]
    distance = np.sqrt((pf["x"] - x0) ** 2 + (pf["y"] - y0) ** 2 + (pf["z"] - z0) ** 2)
    ax.plot(pf["time"] / 365.25, distance, label=pid)

ax.set_ylabel("distance [m]")
ax.set_xlabel("time [year]")
ax.set_title("distance travelled per particle")
ax.grid()
#%%
import numpy as np
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
f = plt.figure(dpi = 300)
ax = f.add_subplot()
mm = flopy.plot.PlotMapView(modelgrid=gwf.modelgrid, ax=ax)
mm.plot_endpoint(well_epd, direction="starting", colorbar=True, s =0.5, cmap = 'hot',zorder = 10)
mm.plot_grid(alpha = 0.2)
points = np.array(list(zip(well_epd['x0'],well_epd['y0'])))   # 30 random points in 2-D
hull = ConvexHull(points)

cent = np.mean(points, axis=0)
pts = points[hull.vertices]

k = 1.1

poly = Polygon(k*(pts - cent) + cent, closed=True,
               capstyle='round',facecolor = 'white', edgecolor = 'black', alpha=0.5, zorder = 1)
plt.gca().add_patch(poly)


