#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 15:50:10 2021

script to reformat and regrid data downloaded from:
https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-count-rev11/data-download

this script only sums population count onto a coarser grid and cannot do more
complex regridding

@author: eebjs
"""

import xarray as xr
import pandas as pd
import numpy as np
from tqdm import tqdm
from scipy import stats

# load global target grid:
gridto = xr.open_dataset('/nfs/see-fs-02_users/eebjs/acrobear/scripts/hia/grids/common_grid.nc')

# path to directory containing gpw netcdf and mapping csv
pop_path = '/nfs/a340/eebjs/hiadata/population_count/raw/'
fname = 'gpw_v4_population_count_rev11_2pt5_min.nc'

# the netcdf is badly formatted
popds = xr.open_dataset(pop_path+fname)

# load and format a lookup table that maps to the raster coordinate in the netcd
lookup = pd.read_csv(pop_path+'gpw_v4_netcdf_contents_rev11.csv')
lookup = lookup.where(lookup.file_name=='gpw_v4_population_count_rev11').dropna()
lookup.index = lookup.index + 1

# load countries lookup:
countries = pd.read_csv(pop_path+'gpw_v4_national_identifier_grid_rev11_lookup.txt', error_bad_lines=False, sep='\t')
countries = countries[['Value', 'ISOCODE']]


#%% reformat
ds = xr.Dataset(
                  coords=popds.drop_dims('raster').coords,
                  attrs=popds.drop_dims('raster').attrs)

for rast in popds.raster.values:
    
    da = popds['Population Count, v4.11 (2000, 2005, 2010, 2015, 2020): 2.5 arc-minutes'][rast-1]
    da_name = lookup.loc[rast, 'raster_name']
    ds[da_name] = da.drop('raster')
    
cpc = pd.DataFrame([countries['Value'], countries['ISOCODE']]).T
cpc = cpc.set_index('Value')
cpc.to_csv('/nfs/a340/eebjs/hiadata/population_count/regridded/country_lookup.csv')


#%% coarsen to target grid function

def coarsen_to_gridto(var, gridto=gridto):

    da = ds[var]
    # replace nan with zero
    # da = da.fillna(0)
    
    # create a dataarray to fill
    popda = xr.DataArray(coords=gridto.drop(['lat_b', 'lon_b']).coords)
    for y in tqdm(range(gridto.dims['latitude'])):
        
        # get lat coordinate (centre, lower and upper)
        clat = float(gridto.coords['latitude'][y])
        llat = float(gridto.coords['lat_b'][y])
        ulat = float(gridto.coords['lat_b'][y+1])
        
        for x in range(gridto.dims['longitude']):
            
            # get lon coordinate (centre, lower and upper)
            clon = float(gridto.coords['longitude'][x])
            llon = float(gridto.coords['lon_b'][x])
            ulon = float(gridto.coords['lon_b'][x+1])
            
            # sum population counts in box
            popsum = float(da.loc[{'longitude':slice(llon, ulon), 'latitude':slice(llat, ulat)}].sum())
            
            # assign summed value to coarsened dataarray
            popda.loc[{'longitude':clon, 'latitude':clat}] = popsum
            
    return popda


        
#%% coarsen the population count dataset

counts = {}
# for year in ['2000', '2005', '2010', '2015', '2020']:
for year in ['2020']:
    da = coarsen_to_gridto(var='Population Count, v4.11 ('+year+')')
    counts[year] = da
    
master_ds = xr.Dataset(counts)
master_ds.to_netcdf('/nfs/a340/eebjs/hiadata/population_count/regridded/population_count.nc')

#%% coarsen the country mask
var = 'National Identifier Grid, v4.11 (2010): National Identifier Grid'
da = ds[var]
# replace nan with zero
# da = da.fillna(0)

# create a dataarray to fill
countries = xr.DataArray(coords=gridto.drop(['lat_b', 'lon_b']).coords)
for y in tqdm(range(gridto.dims['latitude'])):

    
    # get lat coordinate (centre, lower and upper)
    clat = float(gridto.coords['latitude'][y])
    llat = float(gridto.coords['lat_b'][y])
    ulat = float(gridto.coords['lat_b'][y+1])
    
    
    
    for x in range(gridto.dims['longitude']):
        
        # get lon coordinate (centre, lower and upper)
        clon = float(gridto.coords['longitude'][x])
        llon = float(gridto.coords['lon_b'][x])
        ulon = float(gridto.coords['lon_b'][x+1])
        
        # slice the country grid for box
        cda = da.loc[{'longitude':slice(llon, ulon), 'latitude':slice(llat, ulat)}]
        # find the mode
        countries.loc[{'longitude':clon, 'latitude':clat}] = float(stats.mode(cda.values.flatten())[0])
        
        
        
countries = countries.to_dataset(name='mask')
# countries.attrs['lookup'] = cpc.to_dict()['ISOCODE']
        
countries.to_netcdf('/nfs/a340/eebjs/hiadata/population_count/regridded/country_mask.nc')
