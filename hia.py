#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 13:34:36 2021

@author: eebjs
"""

import yaml

# load config file
config = yaml.safe_load(open("./yamls/acrobear_gemm.yml"))

print('Starting health impact assessment for:', config['project_name'])

#%% 1. make common grid based on input model data
print('1: making common grid')
from make_common_grid import make_common_grid
make_common_grid()
print('')


#%% 2. regrid population data to common grid
print('2: regridding population data')
from regrid_population_count import regrid_population_count
regrid_population_count()
print('')

#%% perform HIA
print('3. Health Impact assessment')
from gemm_hia import gemm_hia
gemm_hia()
print('')