#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 20 16:12:33 2022

@author: eebjs
"""

import pandas as pd 
from tqdm import tqdm
import pickle
import xarray as xr

# get mapper for isocode to country name
countries_lookup = pd.read_csv('../lookups/country_lookup.csv',
                               index_col='ISOCODE').squeeze('columns')


def load_popds():

    # open the netcdf contents description csv
    contents = pd.read_csv('/nfs/a340/eebjs/hiadata/population_count/raw/gpw_v4_netcdf_contents_rev11.csv')
    # keep the important rows
    contents = contents.where(contents.file_name=='gpw_v4_population_count_rev11').dropna()
    # increment index
    contents.index = contents.index + 1

    # open GPW population grid
    popds = xr.open_dataset('/nfs/a340/eebjs/hiadata/population_count/raw/gpw_v4_population_count_rev11_2pt5_min.nc')

    ds = xr.Dataset(
                      coords=popds.drop_dims('raster').coords,
                      attrs=popds.drop_dims('raster').attrs)

    for rast in popds.raster.values:

        da = popds['Population Count, v4.11 (2000, 2005, 2010, 2015, 2020): 2.5 arc-minutes'][rast-1]
        da_name = contents.loc[rast, 'raster_name']
        ds[da_name] = da.drop('raster')

    return ds

popds = load_popds()

countries = popds['National Identifier Grid, v4.11 (2010): National Identifier Grid']

#%% get the 5-yearly population data for 2000 to 2020
das = []
for year in list(range(2000,2025,5)):
    
    da = popds[f'Population Count, v4.11 ({year})']
    da.name = 'count'
        
    das.append(da)

#%% create a dataset to store the results
popda = xr.concat(das, dim='year')
popda = popda.assign_coords({'year':pd.date_range('2000-01-01', '2021-01-01', freq='5YS')})

#%% interpolate for every year

popda = popda.resample({'year':'YS'}).interpolate()
# save a copy
popda.to_netcdf('/nfs/a340/eebjs/hiadata/population_count/interpolated/popcount.nc',
                encoding={'count':{'zlib':True,'complevel':7}})

#%% make slices for each year, country
    
for year in popda.year:

    print('making popslices', year)

    popcount = popda.loc[{'year':year}]
        
    das = {} #store slices in this dictonary
    for country_isocode in tqdm(countries_lookup.index):
        
        cvalue = countries_lookup.loc[country_isocode]
        
        # mask using the national identifier grid from GPW
        country_popda = popcount.where(countries == cvalue)
        country_popda.name = country_isocode
        
        # drop empty lon and lat bands
        country_popda = country_popda.dropna(dim='latitude', how='all').dropna(dim='longitude', how='all')
        
        # set nan to 0
        country_popda = country_popda.where(country_popda.notnull(), 0)
        country_popda = country_popda.drop('year')
        
        das[country_isocode] = country_popda
    
    # pickle the dictionary
    with open(f'./grids/pop_count_slices_{pd.to_datetime(year.values).year}.P', 'wb') as handle:
        pickle.dump(das, handle, protocol=pickle.HIGHEST_PROTOCOL)

