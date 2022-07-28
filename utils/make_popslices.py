#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 20 16:12:33 2022

@author: eebjs
"""

import pandas as pd 
from regrid_population_count import load_popds
from hia import config
from tqdm import tqdm
import pickle
import xarray as xr

countries_lookup = pd.read_csv('./lookups/country_lookup.csv',
                               index_col='ISOCODE').squeeze('columns')

popds = load_popds()

countries = popds['National Identifier Grid, v4.11 (2010): National Identifier Grid']

#%%
das = []
for year in list(range(2000,2025,5)):
    
    da = popds[f'Population Count, v4.11 ({year})']
    da.name = 'count'
        
    das.append(da)



#%%
popds = xr.concat(das, dim='year')
popds = popds.assign_coords({'year':list(range(2000,2025,5))})

#%%

das = []
for year in range(2000, 2021):
    print('interpolating', year)
    if year in popds.year:
        da = popds.loc[{'year':year}]
    else:
        da = popds.interp({'year':year})
    das.append(da)
            
popds = xr.concat(das, dim='year')


#%%
    
for year in range(2000, 2021):

    print('making popslices', year)

    popcount = popds.loc[{'year':year}]
        
    das = {}
    for country_isocode in tqdm(countries_lookup.index):
        
        cvalue = countries_lookup.loc[country_isocode]
        
        country_popda = popcount.where(countries == cvalue)
        country_popda.name = country_isocode
        
        # drop empty lon and lat bands
        country_popda = country_popda.dropna(dim='latitude', how='all').dropna(dim='longitude', how='all')
        
        # set nan to 0
        country_popda = country_popda.where(country_popda.notnull(), 0)
        country_popda = country_popda.drop('year')
        
        das[country_isocode] = country_popda
        
    with open(f'./grids/pop_count_slices_{year}.P', 'wb') as handle:
        pickle.dump(das, handle, protocol=pickle.HIGHEST_PROTOCOL)

