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
idx = pd.IndexSlice

dpath = '/nfs/a340/eebjs/hiadata/'

# load regridded population count dataset
popcount = xr.open_dataset(dpath+'population_count/regridded/population_count.nc')['2020']

# load the country mask
countries = xr.open_dataset(dpath+'population_count/regridded/country_mask.nc')
countries_lookup = pd.read_csv(dpath+'population_count/regridded/country_lookup.csv',
                               index_col='ISOCODE')

# load baseline health data
bhdf = pd.read_csv(dpath+'baseline_health_data/raw/IHME-GBD_2019_DATA-31c9bf16-1.csv', index_col=['location_id', 'age_name', 'cause_name', 'measure_name']).sort_index()

# load PM data
# pm25 = xr.open_dataset('/nfs/a340/eebjs/acrobear/cams/annual_means/cams_surface_pm25_2019.nc')['pm2p5'] *1e9 #to convert kg-> ug

pm25 = xr.open_dataset('/nfs/a340/eebjs/acrobear/cams/weighted/background_GFAS_weighted.nc')['weighted t=.7'].mean('time')

# load GEMM parameters
gemm_params = pd.read_csv(dpath+'gemm_function/gemm_parameters.csv', index_col=[0,1,2])
# exclude or include china (IN or EX)
gemm_params = gemm_params.loc['All-regions, INcluding China']

# population age structure
popstruct = pd.read_csv(dpath+'population_structure/'+'IHME_GBD_2019_POP_2019_Y2020M10D15.CSV',
                        index_col='location_id')

# load ISO code mapper for WHO country names
isomap = pd.read_csv(dpath+'WHO_country_isocode_mapper.csv',
                     index_col=1, squeeze=True)

#%% mapping dictionaries

# dictionary that maps the gemm function cause names to GBD cause names
gemm_to_gbd = {'Non-accidental function (Non-Communicable + LRI deaths)':['Non-communicable diseases', 'Lower respiratory infections'],
               'Ischaemic Heart Disease':'Ischemic heart disease',
               'Strokes':'Stroke',
               'Chronic Obstructive Pulomonary Disease':'Chronic obstructive pulmonary disease',
               'Lung Cancer':'Tracheal, bronchus, and lung cancer',
               'Lower Respiratory Infections':'Lower respiratory infections'}

#%% functions

# get relative risk array
def gemm_deaths(pm, alpha, mu, tau, theta, baseline_deaths):
    """Relative risk calculation in GEMM

    Keyword arguments:
    pm     -- the population weighted pm2.5 concentration
    alpha, theta, tau and mu -- parameters from the gemm model
    they are specific to each cause and age cohort
    """
    
    z = max(0, (pm-2.4))
    
    gamma = np.log(1 + z / alpha) / (1 + np.exp((mu - z) / tau))
    
    hazard_ratio = np.exp(theta * gamma)
    
    deaths = (1 - 1 / hazard_ratio) * baseline_deaths

    return deaths

# slice gridded dataarray at country
def country_slice(da, isocode):
    
    cvalue = int(countries_lookup.loc[isocode])
    
    return da.where(countries['mask'] == cvalue)

# population weight country mean
def get_popweight_mean(pm25, isocode):
    
    # get pm25 slice of country
    pm25slice = country_slice(pm25, isocode)
    
    # get population slice of country
    popslice = country_slice(popcount, isocode)
    
    # population weight
    popweighted = float(((pm25slice * popslice) / popslice.sum()).sum())
    
    return popweighted

# get baseline health metric
def get_bhm(country_isocode, age_group, cause, measure):
    
    who_country_id = isomap.loc[country_isocode]
    gbd_cause = gemm_to_gbd[cause]
    if type(gbd_cause) is str:
        gbd_cause = [gbd_cause]
        
    age_group = age_group.replace('-', ' to ')
    age_group = age_group.replace('+', ' plus')
            
    result = float(bhdf.loc[idx[who_country_id, age_group,
                                gbd_cause, measure]]['val'].sum())
    if result == 0:
        raise ValueError('Query returned zero')
            
    return result
        

# function to get the population structure for a given country and age group
def get_age_group_population(country_isocode, age_group):
    
    who_country_id = isomap.loc[country_isocode]
    
    age_group = age_group.replace('-', ' to ')
    age_group = age_group.replace('+', ' plus')
    
    # get the data for the country
    cdata = popstruct.loc[who_country_id]
    total_population = float(cdata.loc[
        (cdata['age_group_name'] == 'All Ages') &\
        (cdata['sex_name'] == 'both')
        ]['val'])
        
    age_group_population = cdata.loc[
        (cdata['age_group_name'] == age_group) &\
        (cdata['sex_name'] == 'both')
        ]['val']
    
    
    age_group_proportion = float(age_group_population / total_population)
    
    popslice = country_slice(popcount, country_isocode)
    
    age_group_population_sliced = popslice * age_group_proportion
    
    return int(age_group_population_sliced.sum())
    
    # get population in country slice
    # using this rather than IHME total population in case only parts of country are in domain
        
    
    return float(age_group_population)
    
    
#%% get iterables

# get a list of countries in the country mask
country_codes = np.unique(countries['mask'])
country_codes = country_codes[country_codes < 1000] # remove NaN code
countries_in = countries_lookup[countries_lookup['Value'].isin(country_codes)].index

# get list of causes in GEMM
causes = list(gemm_to_gbd.keys())

# get list of age ranges
age_groups = gemm_params.index.get_level_values(1).unique()

# get an iterable of uncertainty
uncertainties = ['lower', 'mid', 'upper']

#%% HIA

# create a dataframe to store results
mindex = pd.MultiIndex.from_product([countries_in, age_groups, uncertainties])
results = pd.DataFrame(index=mindex, columns=causes)
results = results.sort_index() # sort for faster indexing


cause = causes[0]


for country_isocode in countries_in:
    print(country_isocode)
    if country_isocode not in isomap.index:
        print('skipping', country_isocode, 'due to no data')
        continue
    
    pm_popweighted = get_popweight_mean(pm25, country_isocode)
        
    for cause in causes:
    
        for age_group in age_groups:
            
            # get baseline death rate
            baseline_deaths = get_bhm(country_isocode=country_isocode,
                    age_group=age_group,
                    cause=cause,
                    measure='Deaths')
            
            age_group_pop = get_age_group_population(country_isocode=country_isocode,
                                                     age_group=age_group)
            
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
                
                mortality_rate = gemm_deaths(pm=pm_popweighted, 
                            alpha=alpha, mu=mu, tau=tau, theta=theta_uncert, 
                            baseline_deaths=baseline_deaths)
    
                # print(age_group, mortality_rate * age_group_pop / 100000)
                deaths = mortality_rate * age_group_pop / 100000
                
                results.loc[idx[country_isocode, age_group, uncert], cause] = deaths
            
# drop rows of results with no results
results = results.dropna(how='all')

results.to_csv('/nfs/a340/eebjs/acrobear/hia_results/CAMS_weighted_hia.csv')