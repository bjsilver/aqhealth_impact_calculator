#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  9 09:22:34 2021

@author: eebjs
"""
#%% load packages and data

import xarray as xr
import pandas as pd
import numpy as np
from tqdm import tqdm
idx = pd.IndexSlice
import yaml

# dictionary that maps the gemm function cause names to GBD cause names
gemm_to_gbd = {'Non-accidental function (Non-Communicable + LRI deaths)':['Non-communicable diseases', 'Lower respiratory infections'],
               'Ischaemic Heart Disease':'Ischemic heart disease',
               'Strokes':'Stroke',
               'Chronic Obstructive Pulomonary Disease':'Chronic obstructive pulmonary disease',
               'Lung Cancer':'Tracheal, bronchus, and lung cancer',
               'Lower Respiratory Infections':'Lower respiratory infections'}

# dictionary to map the uncertainty names used in gemm to gbd columns
gbd_uncert = {'lower':'lower',
              'mid':'val',
              'upper':'upper'}


# slice gridded dataarray at country
def country_slice(da, country_isocode, countries, countries_lookup):
    
    cvalue = int(countries_lookup.loc[country_isocode])
    
    return da.where(countries == cvalue)

# population weight country mean
def get_popweight_mean(pm25, popcount, country_isocode, countries, countries_lookup):
    
    # get pm25 slice of country
    pm25slice = country_slice(pm25, country_isocode, countries, countries_lookup)
    
    # get population slice of country
    popslice = country_slice(popcount, country_isocode, countries, countries_lookup)
    
    # population weight
    popweighted = ((pm25slice * popslice) / popslice.sum()).sum()
    
    return float(popweighted)

# get baseline health metric
def get_bhm(bhdf, country_isocode, age_group, cause, isomap, measure, uncert):
    
    who_country_id = isomap.loc[country_isocode]
    gbd_cause = gemm_to_gbd[cause]
    if type(gbd_cause) is str:
        gbd_cause = [gbd_cause]
        
    age_group = age_group.replace('-', ' to ')
    age_group = age_group.replace('+', ' plus')
            
    result = float(bhdf.loc[idx[who_country_id, age_group,
                                gbd_cause, measure]][gbd_uncert[uncert]].sum())
    if result == 0:
        raise ValueError('Query returned zero')
            
    return result

# get mortality rate array
def gemm_hazard_ratio(pm25, alpha, mu, tau, theta):
    """Relative risk calculation in GEMM

    Keyword arguments:
    pm25     -- the population weighted pm2.5 concentration
    alpha, theta, tau and mu -- parameters from the gemm model
    they are specific to each cause and age cohort
    """
    
    # subtract counterfactual
    z = pm25 - 2.4
    # set values less than zero to zero
    z = z.where(z > 0, 0)
    
    gamma = np.log(1 + z / alpha) / (1 + np.exp((mu - z) / tau))
    
    hazard_ratio = np.exp(theta * gamma)
    
    return hazard_ratio
        

# function to get an array of the age group populaton by country by age group
def get_age_group_population(popstruct, popcount, country_isocode, age_group, isomap, countries, countries_lookup, uncert):
    
    who_country_id = isomap.loc[country_isocode]
    
    age_group = age_group.replace('-', ' to ')
    age_group = age_group.replace('+', ' plus')
    
    # get the data for the country
    cdata = popstruct.loc[who_country_id]
    total_population = float(cdata.loc[
        (cdata['age_group_name'] == 'All Ages') &\
        (cdata['sex_name'] == 'both')
        ][gbd_uncert[uncert]])
        
    age_group_population = cdata.loc[
        (cdata['age_group_name'] == age_group) &\
        (cdata['sex_name'] == 'both')
        ][gbd_uncert[uncert]]
    
    
    age_group_proportion = float(age_group_population / total_population)
    
    popslice = country_slice(popcount, country_isocode, countries, countries_lookup)
    
    age_group_population_sliced = popslice * age_group_proportion
    
    return age_group_population_sliced
    
                
def gemm_hia(yamlfile):
    
    # load config file
    config = yaml.safe_load(open(yamlfile))
    
    ### LOAD DATA
    
    # load regridded population count dataset
    popcount = xr.open_dataset('./grids/population_count.nc')
    popcount = popcount[str(config['population_year'])]
    # load the country mask
    countries = xr.open_dataset('./grids/country_mask.nc')['mask']

    countries_lookup = pd.read_csv('./lookups/country_lookup.csv',
                                   index_col='ISOCODE')
    # load baseline health data
    bhdf = pd.read_csv(config['bh_fpath'],
                       index_col=['location_id', 'age_name', 'cause_name', 'measure_name'])
    bhdf = bhdf.sort_index()
    
    # load PM data
    pm25 = xr.open_dataset(config['model_path'])
    pm25 = pm25[config['pm25var_name']]
    
    # load GEMM parameters
    gemm_params = pd.read_csv(config['gemm_params_fpath'], 
                              index_col=[0,1,2])
    # population age structure
    popstruct = pd.read_csv(config['popstruct_fpath'],
                            index_col='location_id')
    # exclude or include china (IN or EX)
    gemm_params = gemm_params.loc['All-regions, '+config['include_china']+'cluding China']
    # load ISO code mapper for WHO country names
    isomap = pd.read_csv('./lookups/WHO_country_isocode_mapper.csv',
                         index_col=1, squeeze=True)
    
    
    ### MAKE ITERABLES

    # get a list of countries in the country mask
    country_codes = np.unique(countries)
    country_codes = country_codes[country_codes < 1000] # remove NaN code
    countries_in = countries_lookup[countries_lookup['Value'].isin(country_codes)].index
    # get list of causes in GEMM
    causes = list(gemm_to_gbd.keys())
    # get list of age ranges
    age_groups = gemm_params.index.get_level_values(1).unique()
    # get an iterable of uncertainty
    uncertainties = ['lower', 'mid', 'upper']

    ### ITERATE THROUGH HIA

    # create a dataframe to store results
    mindex = pd.MultiIndex.from_product([countries_in, age_groups, uncertainties],
                                        names=['country', 'age_group', 'uncertainty'])
    results = pd.DataFrame(index=mindex, columns=causes)
    results = results.sort_index() # sort for faster indexing

    
            
    for cause in causes:

        print(cause+':')
        for age_group in age_groups:
            print(age_group, end='')
            
            # get GEMM parameters
            theta, theta_se, alpha, mu, tau = gemm_params.loc[cause, age_group]
            
            for uncert in uncertainties:
                    
                # modify theta if calculating an uncertainty bound
                if uncert == 'lower':
                    theta_uncert = theta - theta_se
                elif uncert == 'mid':
                    theta_uncert = theta
                elif uncert == 'upper':
                    theta_uncert = theta + theta_se
                    
                hazard_ratio = gemm_hazard_ratio(pm25, 
                            alpha=alpha, mu=mu, tau=tau, theta=theta_uncert)
                
                for country_isocode in countries_in:
                    if country_isocode not in isomap.index:
                        # print('skipping', country_isocode, 'due to no data')
                        continue
                    
                    # get baseline death rate
                    baseline_deaths = get_bhm(bhdf, country_isocode, age_group, cause,
                                              isomap, measure='Deaths',
                                              uncert=uncert)
                    
                    # get age group population
                    age_group_pop = get_age_group_population(popstruct, popcount, country_isocode, age_group, isomap, countries, countries_lookup, uncert=uncert)
                    
                    mortality_rate = (1 - 1 / hazard_ratio) * baseline_deaths
        
                    deaths = int((mortality_rate * age_group_pop / 100000).sum())
                    
                    results.loc[idx[country_isocode, age_group, uncert], cause] = deaths
                    
                    
            print(' ✓', end='  ')
        print('\n')
    
    results.to_csv('./results/'+config['project_name']+'_results.csv')
    print('results saved to '+'./results/'+config['project_name']+'_results.csv')
