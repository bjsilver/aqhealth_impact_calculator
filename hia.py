#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 13:34:36 2021

@author: eebjs
"""

import yaml
import sys


with open('./config.yaml', 'r') as stream:
    config = yaml.safe_load(stream)




def hia():
    
    with open('./config.yaml', 'r') as stream:
        config = yaml.safe_load(stream)
    
    print('Starting health impact assessment for:', config['scenario_name'])    
    
    #%% 1. make common grid based on input model data
    print('1: making common grid')
    from make_common_grid import make_common_grid
    make_common_grid()
    print('')
    
    #%% 2. regrid population or model data to common grid
    print('2: regridding population data')
    from regrid_model import regrid_model
    regrid_model()
    print('')
    
    #%% perform HIA
    print('3. Health Impact assessment')
    from hia_calculation_parallel import hia_calculation
    hia_calculation()
    print('')

if __name__ == "__main__":
    # config = yaml.safe_load()
    with open(sys.argv[1], 'r') as stream:
        config = yaml.safe_load(stream)
    with open('./config.yaml', 'w') as output:
        yaml.safe_dump(config, output)
    hia()
