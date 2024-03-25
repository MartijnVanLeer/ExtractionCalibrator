import os
import pandas as pd 

if not os.path.isfile(snakemake.output[0]):
    df = pd.DataFrame({'sim' : [],'RMSE' : [], 'xcorlen' : [] , 'zcorlen' : [], 'frac' : []})
else:
    df = pd.read_csv(snakemake.input[0])


sim = snakemake.params.simulation


df['xcorlen'] = sim['xcorlens']
df['zcorlen'] = sim['zcorlens']
df['frac'] = sim['fracs']


Alldf = pd.read_csv(snakemake.output[0])
Alldf = pd.concat(Alldf, df)
Alldf = pd.to_csv(snakemake.output[0])