#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 15:37:30 2022

@author: eebjs
"""

import yaml
import os
import importlib

config_path = '/nfs/see-fs-02_users/eebjs/acrobear/scripts/aqhealth_impact_calculator/config_files/'

for year in range(2001,2021):
    
    for scenario in ['fireon', 'fireoff']:
        
        if os.path.exists(f'./results/cesm_{year}_{scenario}/by_country_results.csv'):
            print('skipping', year, scenario, 'as results dir already exists')
    
            print(f'creating config file for ./results/cesm_{year}_{scenario}')
                

        config = {
         'scenario_name':f'cesm_{year}_{scenario}',
         'fire_scenario': 'fire scenario',
         'health_function': 'GEMM',
         'population_year': 2001,
         'regrid_to': 'pop',
         'model_path': f'/nfs/a340/eebjs/hiadata/pm25_annual_means/pm25_{year}_{scenario}.nc',
         'pm25var_name': 'pm25',
         'latname': 'lat',
         'lonname': 'lon',
         'dpath': '/nfs/a340/eebjs/hiadata/',
         'popdata_dpath': '/nfs/a340/eebjs/hiadata/population_count/raw/',
         'popdata_fname': 'gpw_v4_population_count_rev11_2pt5_min.nc',
         'popdata_contents_fname': 'gpw_v4_netcdf_contents_rev11.csv',
         'popdata_lookup_fname': 'gpw_v4_national_identifier_grid_rev11_lookup.txt',
         'bh_fpath': '/nfs/a340/eebjs/hiadata/baseline_health_data/raw/IHME-GBD_2019_DATA-31c9bf16-1.csv',
         'popstruct_fpath': '/nfs/a340/eebjs/hiadata/population_structure/IHME_GBD_2019_POP_2010_Y2020M10D15.CSV',
         'gemm_params_fpath': '/nfs/a340/eebjs/hiadata/gemm_function/gemm_parameters.csv',
         'hazard_ratio_function': 'gemm',
         'include_china': 'IN'
         }
        
        # make a folder for results if there isn't one
        results_fpath = './results/'+config['scenario_name']
        if not os.path.exists(results_fpath):
            os.mkdir(results_fpath)
            # add the config file as yaml
            fname = results_fpath + '/'+config['scenario_name']+'.yml'
            with open(fname, 'w') as outfile:
                yaml.dump(config, outfile, default_flow_style=False, sort_keys=False)
        
        # os.remove('./config.yaml')
        # with open('./config.yaml', 'w') as outfile:
        #     yaml.dump(config, outfile, default_flow_style=False)
 
        
        
