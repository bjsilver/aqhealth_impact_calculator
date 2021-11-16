#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 10:31:42 2021

@author: eebjs
"""


import pandas as pd
import matplotlib.pyplot as plt

nofire = pd.read_csv('/nfs/a340/eebjs/acrobear/hia_results/CAMS_weighted_hia.csv', index_col=[0,1])

fire = pd.read_csv('/nfs/a340/eebjs/acrobear/hia_results/CAMS_hia.csv', index_col=[0,1])

countries_lookup = pd.read_csv('/nfs/a340/eebjs/hiadata/countries_isocodes.csv',
                               index_col='alpha-3')

def reformat_hiadf(hiadf):
    # convert to thousands
    hiadf = hiadf / 1000
    # sum age groups
    hiadf = hiadf.groupby(hiadf.index.get_level_values(0)).sum()
    # replace country isocodes with names
    # hiadf.index = countries_lookup.loc[hiadf.index, 'name']
    
    return hiadf
    
nofire = reformat_hiadf(nofire)
fire = reformat_hiadf(fire)
    


var = 'Non-accidental function (Non-Communicable + LRI deaths)'

# calculate diff
diff  = fire - nofire
diff = diff.where(fire[var]>5).dropna()

# calculate percent difference
pdiff = ((fire / nofire) * 100) - 100
pdiff = pdiff.where(fire[var]>5).dropna()

#%% plot bar comparison %

fig, ax = plt.subplots(figsize=[10,6])
pdiff[var].plot.bar()
ax.set_ylabel('change in premature mortalities attributed to fires (%)')
ax.set_xlabel('country')

#%% plot bar absolute deaths

fig, ax = plt.subplots(figsize=[10,6])
diff[var].plot.bar()
ax.set_ylabel('change in premature mortalities attributed to fires (thousands)')
ax.set_xlabel('country')