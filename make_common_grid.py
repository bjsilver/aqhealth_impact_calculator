#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 15:07:44 2021

@author: eebjs
"""

import xarray as xr
import numpy as np
import yaml

# load config file
config = yaml.safe_load(open("./yamls/acrobear_gemm.yml"))

def get_resolution(coords):
    res = np.unique(np.abs(np.diff(coords)))
    if len(res) > 1:
        print('inconsistent resolution along', coords.long_name)
    else:
        res = res[0]
        return res
    
def make_grid_from_modeldata(gridin):
    
    
    # get lon and lat resolution
    latres = get_resolution(gridin.latitude)
    lonres = get_resolution(gridin.longitude)
    
    minlat, maxlat = gridin.latitude.min(), gridin.latitude.max()
    minlon, maxlon = gridin.longitude.min(), gridin.longitude.max()
    
    # lats = np.arange(minlat, maxlat+latres, latres)
    # lons = np.arange(minlon.min(), maxlon+lonres, lonres)
    
    # make dataset
    ds =  xr.Dataset(coords=gridin.coords)
    
    # add bound coordinates
    lat_b = np.arange(minlat-latres/2, maxlat+latres/2+latres, latres)[::-1]
    lon_b = np.arange(minlon-lonres/2, maxlon+lonres/2+lonres, lonres)
    ds.coords['lat_b'] = lat_b
    ds.coords['lon_b'] = lon_b

    ds.attrs['latitude resolution'] = latres
    ds.attrs['longitude resolution'] = lonres
    
    return ds

def make_common_grid():
    
    gridin = xr.load_dataset(config['model_path'])
    
    gridto = make_grid_from_modeldata(gridin)

    gridto.to_netcdf('./grids/common_grid.nc')
    
    print('common grid saved as ./grids/common_grid.nc')