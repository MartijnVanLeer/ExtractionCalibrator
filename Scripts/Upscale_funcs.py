import os 
import flopy 
from tqdm import tqdm
import numpy as np

def setup_mf(Lx,Ly,Lz,dx,dy,dz,z, mds,ws):
    mf_exe = mds.exe_name
    workspace = ws
    if not os.path.isdir(workspace):
        os.makedirs(workspace, exist_ok=True)
    
    sim = flopy.mf6.MFSimulation(sim_name='Upscaler', version="mf6", exe_name=mf_exe, sim_ws=workspace)
    tdis = flopy.mf6.ModflowTdis(sim)
    gwf = flopy.mf6.ModflowGwf(sim, modelname='Upscaler', model_nam_file="{}.nam".format('Upscaler'),save_flows=True)
    botm = np.zeros((int(Lz/dz),int(Lx/dx), int(Ly/dy)))
    for lay in range((int(Lz/dz))):
        botm[lay,:,:] = -z[lay]
    dis = flopy.mf6.ModflowGwfdis(gwf, nlay = int(Lz/dz), nrow = int(Lx/dx), ncol = int(Ly/dy), delr = dx, delc = dy, top = 0, botm = botm )
    ims = flopy.mf6.ModflowIms(sim)
    ic = flopy.mf6.ModflowGwfic(gwf, strt =0)
    chd_spd = [] 
    for row in range(int(Lx/dx)):
        for col in range(int(Ly/dy)):
             chd_spd.append([(0,row,col),0 ])
             chd_spd.append( [(int(Lz/dz)-1,row,col),Lz-dz])
    chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd, save_flows=True)
    oc = flopy.mf6.ModflowGwfoc(gwf,head_filerecord='Upscaler' + ".hds", budget_filerecord="{}.cbc".format('Upscaler'), saverecord=[('BUDGET', 'ALL'), ('HEAD', 'ALL')])
    npf = flopy.mf6.ModflowGwfnpf(gwf, k = 10, k33 = 1,  k33overk = False, save_specific_discharge=True)
    sim.write_simulation(silent = True)
    return sim

def Run_MF_WholeField(Kfields,Lx,Ly,Lz,dx,dy,dz, mds, ws):
    z = np.linspace(dz,Lz, int(Lz/dz))
    sim = setup_mf(Lx,Ly,Lz,dx,dy,dz,z, mds,ws)
    
    if Kfields.ndim == 4:
        keffWholeField = []
        for field in tqdm(Kfields, 'Running MODFLOW on total generated fields..'):
            keffWholeField.append(run_mf(sim, field))
    elif Kfields.ndim ==3:
        keffWholeField = run_mf(sim, Kfields, mds,ws)
    return keffWholeField

def run_mf(sim, Kfield,mds, ws):
    ''' Change npf package, run MF, extract K value'''
    gwf = sim.get_model()
    npf = gwf.get_package('NPF')
    npf.k33.set_data(Kfield.transpose(2,0,1))
    npf.k.set_data(Kfield.transpose(2,0,1))
    npf.write()
    success, buff = sim.run_simulation(silent = True)
    if not success:
        print(Kfield)
        raise Exception('Modflow crashed')
    cbb = flopy.utils.CellBudgetFile(os.path.join(ws, f"{gwf.name}.cbc"))
    qs = cbb.get_data(text='DATA-SPDIS')[0]
    qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(qs, gwf)
    K = abs(qz[:, :,:].mean())
    return K
    