#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 20 16:12:33 2022

@author: eebjs
"""

import pandas as pd 
import pickle

dpath = '/nfs/a340/eebjs/hiadata/population_structure/'


mapper = pd.read_csv('../lookups/WHO_country_isocode_mapper.csv',
                     index_col=0).squeeze('columns')

for year in range(2001, 2020):
    
    print(year)
    
    with open(f'../grids/pop_count_slices_{year}.P', 'rb') as handle:
        popslices = pickle.load(handle)
    
    df = pd.read_csv(dpath+f'IHME_GBD_2019_POP_{year}_Y2020M10D15.CSV')
    
    df = df.where(df['sex_name']=='both').dropna()
    
    # reformat age group names
    df.loc[:, 'age_group_name'] = df['age_group_name'].str.replace(' to ','-')
    df.loc[:, 'age_group_name'] = df['age_group_name'].str.replace(' plus','+')
    
    
    age_groups = ['25-29', '30-34', '35-39', '40-44', '45-49', '50-54', '55-59', '60-64',
           '65-69', '70-74', '75-79', '80+']
    
    uncertainties = ['lower', 'mid', 'upper']
    
    mindex = pd.MultiIndex.from_product([uncertainties, age_groups, mapper.values],
                                        names=['uncertainty', 'age_group_name', 'country_isocode'])
  
    
    srs = []
    for location_id in mapper.index:
        
        country_isocode = mapper.loc[location_id]
        
        total_population = int(popslices[country_isocode].sum())
        
        cdf = df.loc[df['location_id'] == location_id]
        cdf = cdf[cdf['age_group_name'].isin(age_groups)]
        cdf = cdf[['age_group_name', 'val','upper','lower']].set_index('age_group_name').stack()
        cdf.index = cdf.index.set_names('uncertainty', level=1)
        
        # divide by total population
        cdf = cdf / total_population
        
        sr = pd.concat([cdf], keys=[country_isocode], names=['country_isocode'])
        srs.append(sr)
        
        
    age_structures = pd.concat(srs)
    age_structures.to_pickle(f'../lookups/age_structure_{year}.P')