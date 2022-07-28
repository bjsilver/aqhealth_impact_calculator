#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 13:34:36 2021

@author: eebjs
"""

import yaml
import os
import sys

config = yaml.safe_load('./config.yaml')


def hia():
    
    print('Starting health impact assessment for:', config['scenario_name'])    
    
    #%% 1. make common grid based on input model data
    print('1: making common grid')
    from make_common_grid import make_common_grid
    make_common_grid()
    print('')
    
    #%% 2. regrid population or model data to common grid
    print('2: regridding population data')
    if config['regrid_to'] == 'mod':
        from regrid_population_count import regrid_population_count
        regrid_population_count()
    elif config['regrid_to'] == 'pop':
        from regrid_model_to_popcount import regrid_model_to_popcount
        regrid_model_to_popcount()
    print('')
    
    #%% perform HIA
    print('3. Health Impact assessment')
    from hia_calculation import hia_calculation
    hia_calculation()
    print('')

if __name__ == "__main__":
    config = yaml.safe_load(open(sys.argv[1]))
    hia()
