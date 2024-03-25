import os
import pandas as pd 

sim = snakemake.params.simulation
path = snakemake.params.resultpath

if not os.path.isfile(path):
    df = pd.DataFrame({'sim' : [],'RMSE' : [], 'xcorlen' : [] , 'zcorlen' : [], 'frac' : []})
else:
    df = pd.read_csv(path)



x = sim['xcorlens']
df['xcorlen'] = x
z = sim['zcorlens']
df['zcorlen'] = z
f = sim['fracs']
df['frac'] = f


Alldf = pd.read_csv(path)
Alldf = pd.concat(Alldf, df)
Alldf = pd.to_csv(path)