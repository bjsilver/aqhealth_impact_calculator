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


        
    
    

def get_resolution(coords):
    res = np.unique(np.abs(np.diff(coords)))
    if len(res) > 1:
        raise ValueError('Dataset is not on a regular grid')
    else:
        res = res[0]
        return res
    

def make_grid_from_popdata(gridin, crop_to):


    ds = xr.Dataset(coords={key: gridin.coords[key] for key in ['latitude', 'longitude']})

    lats = crop_to.coords[config['latname']]
    lons = crop_to.coords[config['lonname']]
    ds = ds.loc[{'latitude':slice(lats.max(), lats.min()),
                 'longitude':slice(lons.min(), lons.max())}]
    return ds



def make_common_grid():
    
    # crop to the model domain
    crop_to = xr.load_dataset(config['model_path'])
       
    gridin = xr.open_dataset(config['popdata_dpath']+\
                                 config['popdata_fname'])
            
    gridto = make_grid_from_popdata(gridin, crop_to)
        
    os.remove('./grids/common_grid.nc')
    gridto.to_netcdf('./grids/common_grid.nc')
    
    print('common grid saved as ./grids/common_grid.nc')