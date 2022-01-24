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
from hia import config
from regrid_population_count import load_popds
import pickle


# dictionary to map the uncertainty names used in gemm to gbd columns
gbd_uncert = {'lower':'lower',
              'mid':'val',
              'upper':'upper'}

# load the common grid
gridto = xr.open_dataset('./grids/common_grid.nc')
# drop the bound coordinates if present
if 'lat_b' in gridto.coords:
    gridto = gridto.drop(['lat_b', 'lon_b'])
    
age_structure = pd.read_pickle('./lookups/age_structure.P')


# slice gridded dataarray at country
def country_slice(da, country_isocode, countries, countries_lookup):
    
    cvalue = countries_lookup.loc[country_isocode]
    
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
def get_bhm(bhdf, country_isocode, age_group, cause, isomap, measure, uncert,
            causemap):
    
    who_country_id = isomap.loc[country_isocode]
    gbd_cause = causemap[cause]
    if type(gbd_cause) is str:
        gbd_cause = [gbd_cause]
        

    result = float(bhdf.loc[idx[who_country_id, age_group,
                                gbd_cause, measure]][gbd_uncert[uncert]].sum())
    # if result == 0:
    #     raise ValueError('Query returned zero')
            
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
def get_age_group_population(country_isocode, age_group, uncert, 
                             popslices, age_structure):

    
    age_group_proportion = age_structure.loc[(uncert, age_group, country_isocode)]
    
    popslice = popslices[country_isocode].copy()
    
    age_group_population_sliced = popslice * age_group_proportion
    
    return age_group_population_sliced
    
                
def hia_calculation():
    
    ### LOAD DATA
    
    # depending on whether regridded to population of model data
    if config['regrid_to'] == 'mod':
        # then load regridded popcount/mask and raw model
        
        pm25 = xr.open_dataset(config['model_path'])
        pm25 = pm25[config['pm25var_name']]
        if 'time' in pm25.dims:
            pm25 = pm25.mean('time')
            
        # load the country mask
        countries = xr.open_dataset('./grids/country_mask.nc')['mask']
        
    elif config['regrid_to'] == 'pop':
        # then load raw popcount/mask and regridded model
        
        popds = load_popds()
        countries = popds['National Identifier Grid, v4.11 (2010): National Identifier Grid']
        
        pm25 = xr.open_dataset('./grids/model_regridded_'+\
                               config['scenario_name']+'.nc')
        pm25 = pm25[config['pm25var_name']]
        if 'time' in pm25.dims:
            pm25 = pm25.mean('time')
        
    with open('./grids/pop_count_slices.P', 'rb') as handle:
        popslices = pickle.load(handle)

    countries_lookup = pd.read_csv('./lookups/country_lookup.csv',
                                   index_col='ISOCODE', squeeze=True)
    # load and clean baseline health data
    bhdf = pd.read_csv(config['bh_fpath'])
    agecol = bhdf['age_name'].str.replace(' to ', '-')
    agecol = agecol.str.replace(' plus', '+')
    bhdf['age_name'] = agecol
    bhdf = bhdf.set_index(['location_id', 'age_name', 'cause_name', 'measure_name'])
    bhdf = bhdf.sort_index()
    


    
    # load ISO code mapper for WHO country names
    isomap = pd.read_csv('./lookups/WHO_country_isocode_mapper.csv',
                         index_col=1, squeeze=True)
    
    # use the correct hazard function
    if config['hazard_ratio_function'] == 'gemm':
        from health_functions.gemm import GEMM_Function as HealthFunc
        from health_functions.gemm import gemm_to_gbd as causemap
    # other health functions could be added later like...
    # elif config['hazard_ratio_function'] == 'gbd2019':
        # from health_functions.gemm import GBD2019_Function as HealthFunc
    
    ### MAKE ITERABLES

    # get a list of countries in the country mask
    # crop to pm25 array
    countries = countries.where(pm25)
    country_codes = np.unique(countries)
    country_codes = country_codes[country_codes < 1000] # remove NaN code
    countries_in = countries_lookup[countries_lookup.isin(country_codes)].index
    # get list of causes in GEMM
    causes = HealthFunc.causes
    # get list of age ranges
    age_groups = HealthFunc.age_groups
    # get an iterable of uncertainty
    uncertainties = ['lower', 'mid', 'upper']
    # uncertainties = ['mid'] # just one for testing
    

    ### ITERATE THROUGH HIA

    # create a dataframe to store results
    mindex = pd.MultiIndex.from_product([countries_in, age_groups, uncertainties],
                                        names=['country', 'age_group', 'uncertainty'])
    results = pd.DataFrame(index=mindex, columns=causes)
    results = results.sort_index() # sort for faster indexing

    # make ds to save gridded results in
    ds = gridto.copy()
            
    for cause in causes:
        print(cause+':')
        
        uncerts = []
        for uncert in uncertainties:
            print(uncert, end=' ')
        # country_deaths list
            age_das = []
            for age_group in age_groups:
                print(age_group, end=' ')
            
            
                    
                # get instance of health function for cause, age_group and uncert
                healthfunc = HealthFunc(cause=cause, age_group=age_group, uncert=uncert)
                hazard_ratio = healthfunc.calculate_hazard_ratio(pm25)
                
                age_da = xr.DataArray(name=age_group, coords=gridto.coords)
                age_da[:,:] = 0
                
                for country_isocode in countries_in:
                    if country_isocode not in isomap.index:
                        # print('skipping', country_isocode, 'due to no data')
                        continue
                    
                    # print(country_isocode)
                    # get baseline death rate
                    
                    baseline_deaths = get_bhm(bhdf, country_isocode, age_group, cause,isomap, measure='Deaths', uncert=uncert, causemap=causemap)
                    
          
                    # get age group population
                    age_group_pop = get_age_group_population(country_isocode, age_group, uncert, popslices, age_structure)
                 
                    hazard_ratio_slice = hazard_ratio.where(age_group_pop.notnull())
                    # calculate mortality rate from hazard ratio
                    mortality_rate = (1 - 1 / hazard_ratio_slice) * baseline_deaths
                    
                    # calculate premature mortalities
                    deaths = mortality_rate * age_group_pop / 100000
                    # deaths = deaths.expand_dims({'country':[country_isocode]})
                    
                    age_da.loc[deaths.coords] += deaths.values
                    
                    results.loc[idx[country_isocode, age_group, uncert], cause] = float(deaths.sum())
                    
                age_das.append(age_da)
                
            # sum across age groups
            deaths_summed = sum(age_das)
            
            # add uncertainty dimension
            deaths_summed = deaths_summed.expand_dims({'uncert':[uncert]})
            uncerts.append(deaths_summed)
            print(' ✓')
        
        da = xr.concat(uncerts, dim='uncert')
        ds[cause] = da
        print('\n')
        
    ds.to_netcdf('./results/'+config['scenario_name']+'/gridded_results.nc')
        
    
    results.to_csv('./results/'+config['scenario_name']+'/by_country_results.csv')
    print('./results/'+config['scenario_name']+'/by_country_results.csv')

#%%

if __name__ == "__main__":
    hia_calculation()