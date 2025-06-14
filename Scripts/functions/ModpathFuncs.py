import os 
import nlmod
import flopy 
import numpy as np
from tqdm import tqdm

#Load packages from transient model for calibrated k values
def load_calibrated_npf(modelname):
    tfolder = os.path.join('..','Results',f'{modelname}', f'{modelname}_t', 'Fitter')
    tmod = flopy.mf6.mfsimulation.MFSimulation.load('mfsim', sim_ws = tfolder)
    gwft = tmod.get_model()
    npft = gwft.get_package('NPF')
    npfk = npft.k.data
    npfk33 = npft.k33.data
    return npfk, npfk33


#get initial steady state calibrated model and update K values
def load_ss(destFolder,ds, npfk, npfk33):
    sim = flopy.mf6.mfsimulation.MFSimulation.load('mfsim', sim_ws = destFolder, exe_name =ds.exe_name)
    gwf = sim.get_model()
    npf = gwf.get_package('NPF')
    npf.k = npfk
    npf.k33 = npfk33
    npf.save_flows = True
    return sim, npf

#run 'homogeneous' model with updated k values
def run_modpath_ref_bw(modelname, sim, ds, npf, layer,wellxy):
    npf.write()
    nlmod.sim.write_and_run(sim, ds, write_ds = False, silent = True)
    flowfrac = run_bw(modelname, sim, ds, layer,wellxy)
    return flowfrac

def run_modpath_ref_fw(modelname, sim, ds, npf, layer,wellxy):
    npf.write()
    nlmod.sim.write_and_run(sim, ds, write_ds = False, silent = True)
    dist = run_fw(modelname, sim, ds, layer,wellxy)
    return dist

def run_modpath_realizations(modelname,sim,ds,npf, rds, layer,wellxy):
    flowfrac = []
    dist = []
    layno = list(ds.layer).index(layer)
    for i in tqdm(range(len(rds.index))):
        data33 = npf.k33.array
        data33[layno] = rds.isel(index = i).k.values
        npf.k33.set_data(data33)
        data = npf.k.array
        data[layno] = rds.isel(index = i).k.values
        npf.k.set_data(data)
        npf.write()
        nlmod.sim.write_and_run(sim, ds, write_ds = False, silent = True)
        dist.extend(run_fw(modelname, sim, ds, layer, wellxy))
        flowfrac.append(run_bw(modelname,sim,ds,layer,wellxy))
    return flowfrac, dist



def run_bw(modelname, sim, ds, layer, wellxy):
    gwf = sim.get_model()
    #BW tracking

    mpf_bw = nlmod.modpath.mpf(gwf)
    mpfbas = nlmod.modpath.bas(mpf_bw,)
    layno = list(ds.layer).index(layer)
    layernodes = []
    for lay in range(layno-1):
        layernodes += nlmod.modpath.layer_to_nodes(mpf_bw, lay)
    welnodes = nlmod.modpath.xy_to_nodes(wellxy,mpf_bw,ds,layno+1)
    pg_bw = nlmod.modpath.pg_from_fdt(welnodes, divisions = 10)
    mpsim = nlmod.modpath.sim(mpf_bw, pg_bw, "backward", gwf = gwf, weaksinkoption = 'stop_at', weaksourceoption= 'stop_at')
    nlmod.modpath.write_and_run(mpf_bw, silent = True)
    #Get fraction of endpoints through aquitard
    fpth = os.path.join(ds.model_ws, 'modpath', f"mp7_{modelname.lower()}_ss.mpend")
    e = flopy.utils.EndpointFile(fpth)
    aquitard_epd = e.get_destination_endpoint_data(dest_cells = layernodes)
    aquitard_epd = aquitard_epd[np.isin(aquitard_epd.initialcellface, [1,2,3,4])]
    all_epd = e.get_alldata()
    all_epd = all_epd[np.isin(all_epd.initialcellface, [1,2,3,4])]
    flowfrac =  sum(1/aquitard_epd.time)/sum(1/all_epd.time)
    return flowfrac

def run_fw(modelname, sim, ds, layer, wellxy):
    layno = list(ds.layer).index(layer)
    gwf = sim.get_model()
    mpf = nlmod.modpath.mpf(gwf)
    mpfbas = nlmod.modpath.bas(mpf,)
    layernodes = nlmod.modpath.layer_to_nodes(mpf, layno-1)
    # pg = nlmod.modpath.pg_from_fdt(nodes, divisions = 3)
    pg = nlmod.modpath.pg_from_pd(layernodes, localx=0.5, localy=0.5, localz=0.1)
    mpsim = nlmod.modpath.sim(mpf, pg, "forward")
    nlmod.modpath.write_and_run(mpf, silent = True)

    fpth = os.path.join(ds.model_ws, 'modpath', f"mp7_{modelname.lower()}_ss.mpend")
    e = flopy.utils.EndpointFile(fpth)
    welnodes = nlmod.modpath.xy_to_nodes(wellxy,mpf,ds,layno+1)
    all_epd = e.get_destination_endpoint_data(dest_cells = welnodes)
    dist = all_epd.time/365
    return dist
