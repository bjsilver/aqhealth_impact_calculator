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
from hia import config
from utils.ufuncs import load_popds

def coarsen_to_gridto(da, gridto, lonname, latname):

    # replace nan with zero
    # da = da.fillna(0)
    
    # create a dataarray to fill
    popda = xr.DataArray(coords=gridto.drop(['lat_b', 'lon_b']).coords)
    for y in tqdm(range(gridto.dims[latname]),
                  'regridding population data'):
        
        # get lat coordinate (centre, lower and upper)
        clat = float(gridto.coords[latname][y])
        llat = float(gridto.coords['lat_b'][y])
        ulat = float(gridto.coords['lat_b'][y+1])
        
        for x in range(gridto.dims[lonname]):
            
            # get lon coordinate (centre, lower and upper)
            clon = float(gridto.coords[lonname][x])
            llon = float(gridto.coords['lon_b'][x])
            ulon = float(gridto.coords['lon_b'][x+1])
            
            # sum population counts in box
            popsum = float(da.loc[{'longitude':slice(llon, ulon), 'latitude':slice(llat, ulat)}].sum())
            
            # assign summed value to coarsened dataarray
            popda.loc[{lonname:clon, latname:clat}] = popsum
            
    popda.attrs = gridto.attrs
            
    return popda     




def coarsen_country_grid(da, gridto, lonname, latname):

    # create a dataarray to fill
    countries = xr.DataArray(coords=gridto.drop(['lat_b', 'lon_b']).coords)
    for y in tqdm(range(gridto.dims[latname]),
                  'regridding countries mask'):
    
        
        # get lat coordinate (centre, lower and upper)
        clat = float(gridto.coords[latname][y])
        llat = float(gridto.coords['lat_b'][y])
        ulat = float(gridto.coords['lat_b'][y+1])
        
        
        
        for x in range(gridto.dims[lonname]):
            
            # get lon coordinate (centre, lower and upper)
            clon = float(gridto.coords[lonname][x])
            llon = float(gridto.coords['lon_b'][x])
            ulon = float(gridto.coords['lon_b'][x+1])
            
            # slice the country grid for box
            cda = da.loc[{'longitude':slice(llon, ulon), 'latitude':slice(llat, ulat)}]
            # find the mode
            countries.loc[{lonname:clon, latname:clat}] = float(stats.mode(cda.values.flatten())[0])
        
    countries.attrs = gridto.attrs
        
    return countries.to_dataset(name='mask')        





def regrid_population_count():
    
    
    lonname, latname = config['lonname'], config['latname']
    
    ### LOAD DATA
    
    # load global target grid:
    gridto = xr.open_dataset('./grids/common_grid.nc')
    
    
    # reformat GPW ds
    ds = load_popds()

    
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
                              gridto=gridto,
                              lonname=lonname, latname=latname)
    popda = popda.to_dataset(name=str(config['population_year']))
    popda.to_netcdf('./grids/population_count.nc')
    
    
    
    ### COARSEN COUNTRY GRID
    countryda = coarsen_country_grid(ds['National Identifier Grid, v4.11 (2010): National Identifier Grid'], 
                                     gridto=gridto,
                                     lonname=lonname, latname=latname)
    countryda.to_netcdf('./grids/country_mask.nc')
