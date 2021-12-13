#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 17:44:47 2021

@author: eebjs
"""

import xarray as xr
import pandas as pd
from regrid_population_count import reformat_GPW_ds
from hia import config
import xesmf as xe

### LOAD POP DATA

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

popda = ds['Population Count, v4.11 (2020)']

modelda = xr.load_dataset(config['model_path'])['pm2p5']
lats = modelda.latitude
lons = modelda.longitude
minlat, maxlat = lats.min(), lats.max()
minlon, maxlon = lons.min(), lons.max()


popda = popda.loc[{'latitude':slice(maxlat, minlat), 'longitude':slice(minlon, maxlon)}]

#%%
regridder = xe.Regridder(modelda, popda, 'nearest_s2d')

da_regridded = regridder(modelda, keep_attrs=True)
