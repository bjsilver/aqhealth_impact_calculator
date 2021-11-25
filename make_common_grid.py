#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 15:07:44 2021

@author: eebjs
"""

import xarray as xr
import numpy as np
import yaml

def get_resolution(coords):
    res = np.unique(np.abs(np.diff(coords)))
    if len(res) > 1:
        raise ValueError('Dataset is not on a regular grid')
    else:
        res = res[0]
        return res
    
def make_grid_from_modeldata(gridin, latname, lonname):
    
    lats = gridin.coords[latname]
    lons = gridin.coords[lonname]
    
    # get lon and lat resolution
    latres = get_resolution(lats)
    lonres = get_resolution(lons)
    
    minlat, maxlat = lats.min(), lats.max()
    minlon, maxlon = lons.min(), lons.max()
    
    # lats = np.arange(minlat, maxlat+latres, latres)
    # lons = np.arange(minlon.min(), maxlon+lonres, lonres)
    
    # make dataset
    ds =  xr.Dataset(coords={key: gridin.coords[key] for key in [latname, lonname]})
    
    # add bound coordinates
    lat_b = np.arange(minlat-latres/2, maxlat+latres/2+latres, latres)[::-1]
    lon_b = np.arange(minlon-lonres/2, maxlon+lonres/2+lonres, lonres)
    ds.coords['lat_b'] = lat_b
    ds.coords['lon_b'] = lon_b

    ds.attrs['latitude resolution'] = latres
    ds.attrs['longitude resolution'] = lonres
    
    return ds

def make_common_grid(yamlfile):
    
    # load config file
    config = yaml.safe_load(open(yamlfile))
    
    gridin = xr.load_dataset(config['model_path'])
    
    gridto = make_grid_from_modeldata(gridin, 
                                      latname=config['latname'],
                                      lonname=config['lonname'])

    gridto.to_netcdf('./grids/common_grid.nc')
    
    print('common grid saved as ./grids/common_grid.nc')
    
