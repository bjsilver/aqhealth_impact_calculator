#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 17:44:47 2021

@author: eebjs
"""

import os
import xarray as xr
import numpy as np
from regrid_population_count import load_popds
from hia import config, SCENARIO_NAME
import xesmf as xe
import pickle
    
    

#%%
def regrid_model_to_popcount():
    
    if os.path.exists('./grids/model_regridded_'+SCENARIO_NAME+'.nc'):
        print('regridded model data found at \n'+\
             './grids/model_regridded_'+SCENARIO_NAME+'.nc' )
        return
    
    modelda = xr.load_dataset(config['model_path'])[config['pm25var_name']]
    countries = load_popds()
    
    lats = modelda.coords[config['latname']]
    lons = modelda.coords[config['lonname']]
    
    # get max and min lat and lon from regridded countries data
    maxlat, minlat = lats.max(), lats.min()
    maxlon, minlon = lons.max(), lons.min()

    countries = countries.loc[{'latitude':slice(maxlat, minlat),
                                      'longitude':slice(minlon, maxlon)}]

    popda = countries['Population Count, v4.11 ('+\
                              str(config['population_year'])+')']
    popda = popda.loc[{'latitude':slice(maxlat, minlat), 'longitude':slice(minlon, maxlon)}]
    popda.name = str(config['population_year'])
    popda.to_netcdf('./grids/population_count.nc')
    
    
    weights_fname = f'pop_{popda.shape[0]}x{popda.shape[1]}_to_mod_{modelda.shape[0]}x{modelda.shape[1]}_weights.nc'
    
    if os.path.exists(weights_fname):
        print('reusing regridder weights')
        regridder = xe.Regridder(modelda, popda, 'bilinear', weights=weights_fname)
    else:
        print('making regridder')
        regridder = xe.Regridder(modelda, popda, 'bilinear')
    
    da_regridded = regridder(modelda, keep_attrs=True)
    
    da_regridded.to_netcdf('./grids/model_regridded_'+SCENARIO_NAME+'.nc')
    

    
if __name__ == "__main__":
    regrid_model_to_popcount()