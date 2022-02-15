#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 1 10:02:28 2021

@author: eebjs
"""

import pandas as pd
from hia import config
import xarray as xr


def read_hiadf(path):
    # open df
    hiadf = pd.read_csv(path, index_col=[0,1,2])
    # convert to thousands
    hiadf = hiadf / 1000
    # sum age groups
    hiadf = hiadf.groupby(hiadf.index.get_level_values(0)).sum()
    # replace country isocodes with names
    # hiadf.index = countries_lookup.loc[hiadf.index, 'name']
    
    return hiadf

def load_popds():
    
    # open the netcdf contents description csv
    contents = pd.read_csv(config['popdata_dpath']+\
                           config['popdata_contents_fname'])
    # keep the important rows
    contents = contents.where(contents.file_name=='gpw_v4_population_count_rev11').dropna()
    # increment index
    contents.index = contents.index + 1
    
    # open GPW population grid
    popds = xr.open_dataset(config['popdata_dpath']+\
                            config['popdata_fname'])
    
    ds = xr.Dataset(
                      coords=popds.drop_dims('raster').coords,
                      attrs=popds.drop_dims('raster').attrs)
    
    for rast in popds.raster.values:
        
        da = popds['Population Count, v4.11 (2000, 2005, 2010, 2015, 2020): 2.5 arc-minutes'][rast-1]
        da_name = contents.loc[rast, 'raster_name']
        ds[da_name] = da.drop('raster')
        
    return ds