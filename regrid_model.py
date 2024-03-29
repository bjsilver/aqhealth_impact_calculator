#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 17:44:47 2021

@author: eebjs
"""

import os
import xarray as xr
from hia import config
import xesmf as xe
    
    

#%%
def regrid_model():

    SCENARIO_NAME = config['scenario_name']
    print(SCENARIO_NAME)
    
    if os.path.exists('./grids/model_regridded_'+SCENARIO_NAME+'.nc'):
        print('regridded model data found at \n'+\
              './grids/model_regridded_'+SCENARIO_NAME+'.nc' )
        return
    
    modelda = xr.load_dataset(config['model_path'])[config['pm25var_name']]
    common = xr.load_dataset(('./grids/common_grid.nc'))
    
    lats = modelda.coords[config['latname']]
    lons = modelda.coords[config['lonname']]
    
    # get max and min lat and lon from regridded countries data
    maxlat, minlat = lats.max(), lats.min()
    maxlon, minlon = lons.max(), lons.min()

    common = common.loc[{'latitude':slice(maxlat, minlat),
                                      'longitude':slice(minlon, maxlon)}]

    # popda = countries['Population Count, v4.11 (2000)']
    # popda = popda.loc[{'latitude':slice(maxlat, minlat), 'longitude':slice(minlon, maxlon)}]
    # popda.name = str(config['population_year'])
    # popda.to_netcdf('./grids/population_count.nc')
    
    
    weights_fname = f'pop_{len(common.longitude)}x{len(common.latitude)}_to_mod_{modelda.shape[0]}x{modelda.shape[1]}_weights.nc'
    
    if os.path.exists('./grids/'+weights_fname):
        print('reusing regridder weights')
        regridder = xe.Regridder(modelda, common, 'bilinear', weights='./grids/'+weights_fname)
    else:
        print('making regridder')
        regridder = xe.Regridder(modelda, common, 'bilinear')
        regridder.to_netcdf('./grids/'+weights_fname)
        
    da_regridded = regridder(modelda, keep_attrs=True)
    
    da_regridded.to_netcdf('./grids/model_regridded_'+SCENARIO_NAME+'.nc')
    

    
if __name__ == "__main__":
    regrid_model()