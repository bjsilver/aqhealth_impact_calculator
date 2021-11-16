#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 17:22:01 2021

@author: eebjs
"""

import pandas as pd
from fuzzywuzzy import process

dpath = '/nfs/a340/eebjs/hiadata/'

countries_lookup = pd.read_csv(dpath+'countries_isocodes.csv',
                               index_col='name')

bhdf = pd.read_csv(dpath+'baseline_health_data/raw/IHME-GBD_2019_DATA-31c9bf16-1.csv')

who_countries = bhdf['location_name'].unique()
who_countries_ids = bhdf['location_id'].unique()
    
#%%
results = pd.Series(index=who_countries, dtype=object)
for cstr in who_countries:
    print(cstr)

    if cstr in countries_lookup.index:
        results.loc[cstr] = countries_lookup.loc[cstr]['alpha-3']
    else:
        match = process.extractOne(cstr, countries_lookup.index)
        
        results.loc[cstr] = countries_lookup.loc[match[0]]['alpha-3']
    
# replace results index with location_id
results.index = who_countries_ids
    
results.to_csv(dpath+'WHO_country_isocode_mapper.csv')
