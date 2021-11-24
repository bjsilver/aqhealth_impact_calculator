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
from tqdm import tqdm
from scipy import stats
import yaml






#%% reformat
def reformat_GPW_ds(popds, lookup):
    
    ds = xr.Dataset(
                      coords=popds.drop_dims('raster').coords,
                      attrs=popds.drop_dims('raster').attrs)
    
    for rast in popds.raster.values:
        
        da = popds['Population Count, v4.11 (2000, 2005, 2010, 2015, 2020): 2.5 arc-minutes'][rast-1]
        da_name = lookup.loc[rast, 'raster_name']
        ds[da_name] = da.drop('raster')
        
    return ds

#%% coarsen to target grid function

def coarsen_to_gridto(da, gridto):

    # replace nan with zero
    # da = da.fillna(0)
    
    # create a dataarray to fill
    popda = xr.DataArray(coords=gridto.drop(['lat_b', 'lon_b']).coords)
    for y in tqdm(range(gridto.dims['latitude']),
                  'regridding population data'):
        
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


        


#%% coarsen the country mask

def coarsen_country_grid(da, gridto):

    # create a dataarray to fill
    countries = xr.DataArray(coords=gridto.drop(['lat_b', 'lon_b']).coords)
    for y in tqdm(range(gridto.dims['latitude']),
                  'regridding countries mask'):
    
        
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
        
        
        
    return countries.to_dataset(name='mask')        



#%% main

def regrid_population_count(yamlfile):
    
    # load config file
    config = yaml.safe_load(open(yamlfile))
    
    ### LOAD DATA
    
    # load global target grid:
    gridto = xr.open_dataset('./grids/common_grid.nc')
    
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
    # reformat GPW ds
    ds = reformat_GPW_ds(popds, lookup=contents)

    
    # load countries lookup:
    countries = pd.read_csv(config['popdata_dpath']+\
                            config['popdata_lookup_fname'],
                            sep='\t')
    countries = countries[['Value', 'ISOCODE']]
    # reformat and save as csv
    cpc = pd.DataFrame([countries['Value'], countries['ISOCODE']]).T
    cpc = cpc.set_index('Value')
    # save lookup for GPW value -> ISO mapping
    cpc.to_csv('./lookups/country_lookup.csv')
    
    
    
    ### COARSEN POPULATION COUNT
    popda = coarsen_to_gridto(ds['Population Count, v4.11 ('+\
                              str(config['population_year'])+')'],
                              gridto=gridto)
    popda = popda.to_dataset(name=str(config['population_year']))
    popda.to_netcdf('./grids/population_count.nc')
    
    
    
    ### COARSEN COUNTRY GRID
    countryda = coarsen_country_grid(ds['National Identifier Grid, v4.11 (2010): National Identifier Grid'], 
                                     gridto=gridto)
    countryda.to_netcdf('./grids/country_mask.nc')
    
if __name__ == '__main__':
    regrid_population_count()