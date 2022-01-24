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

countries_lookup = pd.read_csv('./lookups/country_lookup.csv',
                               index_col='ISOCODE', squeeze=True)

popds = load_popds()
popcount = popds['Population Count, v4.11 ('+\
                          str(config['population_year'])+')']
countries = popds['National Identifier Grid, v4.11 (2010): National Identifier Grid']

das = {}
for country_isocode in tqdm(countries_lookup.index):
    
    cvalue = countries_lookup.loc[country_isocode]
    
    country_popda = popcount.where(countries == cvalue)
    country_popda.name = country_isocode
    
    # drop empty lon and lat bands
    country_popda = country_popda.dropna(dim='latitude', how='all').dropna(dim='longitude', how='all')
    
    # set nan to 0
    country_popda = country_popda.where(country_popda.notnull(), 0)
    
    das[country_isocode] = country_popda
    
with open('./grids/pop_count_slices.P', 'wb') as handle:
    pickle.dump(das, handle, protocol=pickle.HIGHEST_PROTOCOL)

