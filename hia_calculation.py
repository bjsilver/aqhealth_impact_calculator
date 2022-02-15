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
from hia import config, SCENARIO_NAME
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

class BaseLineHealthData:
    
    def __init__(self, config):
        self.config = config
        
        if self.config['hazard_ratio_function'] == 'gemm':
            from health_functions.gemm import gemm_to_gbd as causemap
        self.causemap = causemap

        # load and clean baseline health data
        bhdf = pd.read_csv(self.config['bh_fpath'])
        agecol = bhdf['age_name'].str.replace(' to ', '-')
        agecol = agecol.str.replace(' plus', '+')
        bhdf['age_name'] = agecol
        bhdf = bhdf.set_index(['location_id', 'age_name', 'cause_name', 'measure_name'])
        bhdf = bhdf.sort_index()
        self.bhdf = bhdf
        
        self.isomap = pd.read_csv('./lookups/WHO_country_isocode_mapper.csv',
                             index_col=1, squeeze=True)
     
    def lookup(self, country_isocode, cause, age_group,
               measure, uncert):
        
        who_country_id = self.isomap.loc[country_isocode]
        gbd_cause = self.causemap[cause]
        if type(gbd_cause) is str:
            gbd_cause = [gbd_cause]
            

        result = float(self.bhdf.loc[idx[who_country_id, age_group,
                                    gbd_cause, measure]][gbd_uncert[uncert]].sum())
        
        return result
    
class PopulationData():
    
    def __init__(self):
        with open('./grids/pop_count_slices.P', 'rb') as handle:
            self.popslices = pickle.load(handle)
        self.age_structure = pd.read_pickle('./lookups/age_structure.P')
        
    def lookup(self, country_isocode, age_group, uncert):
        age_group_proportion = self.age_structure.loc[(uncert, age_group, country_isocode)]
        popslice = self.popslices[country_isocode].copy()
        age_group_population_sliced = popslice * age_group_proportion
        return age_group_population_sliced
    
    # population weight country mean
    # def get_popweight_mean(pm25, popcount, country_isocode, countries, countries_lookup):
        
    #     # get pm25 slice of country
    #     pm25slice = country_slice(pm25, country_isocode, countries, countries_lookup)
        
    #     # get population slice of country
    #     popslice = country_slice(popcount, country_isocode, countries, countries_lookup)
        
    #     # population weight
    #     popweighted = ((pm25slice * popslice) / popslice.sum()).sum()
        
    #     return float(popweighted)
    
                
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
                               SCENARIO_NAME+'.nc')
        pm25 = pm25[config['pm25var_name']]
        if 'time' in pm25.dims:
            pm25 = pm25.mean('time')
        
    

    countries_lookup = pd.read_csv('./lookups/country_lookup.csv',
                                   index_col='ISOCODE', squeeze=True)

    


    
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
    
    # instatiate datasets
    bhd = BaseLineHealthData(config)
    popdata = PopulationData()
    

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
                    
                    baseline_deaths = bhd.lookup(country_isocode=country_isocode,
                                                 cause=cause, 
                                                 age_group=age_group, 
                                                 measure='Deaths', 
                                                 uncert=uncert)
                    
          
                    # get age group population
                    age_group_pop = popdata.lookup(country_isocode, age_group, uncert)
                 
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
        
    ds.to_netcdf('./results/'+SCENARIO_NAME+'/gridded_results.nc')
        
    
    results.to_csv('./results/'+SCENARIO_NAME+'/by_country_results.csv')
    print('./results/'+SCENARIO_NAME+'/by_country_results.csv')

#%%

if __name__ == "__main__":
    hia_calculation()