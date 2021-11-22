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
