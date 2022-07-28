#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 15:07:44 2021

@author: eebjs
"""

import xarray as xr
import numpy as np
from hia import config
import os

def make_common_grid():
    
    # crop to the model domain
    crop_to = xr.load_dataset(config['model_path'])
       
    gridin = xr.open_dataset(config['popdata_dpath']+\
                                 config['popdata_fname'])
            
    gridto = make_grid_from_popdata(gridin, crop_to)
        
    os.remove('./grids/common_grid.nc')
    gridto.to_netcdf('./grids/common_grid.nc')
        
    
    print('common grid saved as ./grids/common_grid.nc')

def get_resolution(coords):
    res = np.unique(np.abs(np.diff(coords)))
    if len(res) > 1:
        raise ValueError('Dataset is not on a regular grid')
    else:
        res = res[0]
        return res
    
def make_grid_from_modeldata(gridin, crop_to):
    
    lats = crop_to.coords[config['latname']]
    lons = crop_to.coords[config['lonname']]
    
    # get lon and lat resolution
    latres = get_resolution(lats)
    lonres = get_resolution(lons)
    
    minlat, maxlat = lats.min(), lats.max()
    minlon, maxlon = lons.min(), lons.max()
    
    # lats = np.arange(minlat, maxlat+latres, latres)
    # lons = np.arange(minlon.min(), maxlon+lonres, lonres)
    
    # make dataset
    ds =  xr.Dataset(coords={key: gridin.coords[key] for key in [config['latname'], config['lonname']]})
    
    # add bound coordinates
    lat_b = np.arange(minlat-latres/2, maxlat+latres/2+latres, latres)[::-1]
    lon_b = np.arange(minlon-lonres/2, maxlon+lonres/2+lonres, lonres)
    ds.coords['lat_b'] = lat_b
    ds.coords['lon_b'] = lon_b

    ds.attrs['latitude resolution'] = latres
    ds.attrs['longitude resolution'] = lonres
    
    return ds 
