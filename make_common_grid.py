#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 15:07:44 2021

@author: eebjs
"""

import xarray as xr
import numpy as np

# infile = sys.argv[1]
infile = '/nfs/a340/eebjs/acrobear/cams/annual_means/cams_surface_pm25_2020.nc'

gridin = xr.load_dataset(infile)
# gridin = gridin.drop_dims('time')

#%%
def get_resolution(coords):
    res = np.unique(np.abs(np.diff(gridin.latitude)))
    if len(res) > 1:
        print('inconsistent resolution along', coords.long_name)
    else:
        res = res[0]
        return res
    
def make_common_grid(gridin):
    
    
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
    
#%%

# gridto = xeu.grid_global(d_lon=lonres, d_lat=latres)
gridto = make_common_grid(gridin)

gridto.to_netcdf('/nfs/see-fs-02_users/eebjs/acrobear/scripts/hia/grids/common_grid.nc')

