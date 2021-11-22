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

#%% make common grid
print('making common grid')
from make_common_grid import make_common_grid
make_common_grid(infile=config['model_path'])

#%% regrid population data to common grid
print('regrid the population data')
from regrid_population_count import regrid_population_count
regrid_population_count(popdata_dpath=config['popdata_dpath'], 
                        popdata_fname=config['popdata_fname'], 
                        popdata_contents_fname=config['popdata_contents_fname'], 
                        popdata_lookup_fname=config['popdata_lookup_fname'], 
                        population_year=config['population_year'])