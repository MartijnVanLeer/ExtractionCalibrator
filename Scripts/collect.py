import os
import pandas as pd 

if not os.path.isfile(snakemake.output[0]):
    df = pd.DataFrame({'sim' : [],'RMSE' : [], 'xcorlen' : [] , 'zcorlen' : [], 'frac' : []})
    df.to_csv[snakemake.output[0]]

sim = snakemake.params.simulation

df = pd.read_csv(snakemake.input[0])
df['xcorlen'] = sim['xcorlen']
df['zcorlen'] = sim['zcorlen']
df['frac'] = sim['fracs']


Alldf = pd.read_csv(snakemake.output[0])
Alldf = pd.concat(Alldf, df)
Alldf = pd.to_csv(snakemake.output[0])