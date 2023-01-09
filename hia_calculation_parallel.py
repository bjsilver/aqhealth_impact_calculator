#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  9 09:22:34 2021

@author: eebjs
"""
#%% load packages and data

import xarray as xr
import pandas as pd
# import numpy as np
idx = pd.IndexSlice
from hia import config
import pickle
import health_functions
import multiprocessing as mp
from itertools import product
from tqdm import tqdm
from pyproj import Geod
from shapely.geometry import LineString, Point, Polygon


# use the correct hazard function
HealthFunc = getattr(health_functions, config['hazard_ratio_function'])

# dictionary to map the uncertainty names used in gemm to gbd columns
gbd_uncert = {'lower':'lower',
              'mid':'val',
              'upper':'upper'}

# load the common grid
gridto = xr.load_dataset('./grids/common_grid.nc')
# drop the bound coordinates if present
if 'lat_b' in gridto.coords:
    gridto = gridto.drop(['lat_b', 'lon_b'])

def get_latlon_areas(latcentres, loncentres, latres, lonres):

    geod = Geod(ellps="WGS84")
    
    areas = pd.Series(dtype=float, index=latcentres)
    for clat in latcentres:
        # get bounds from coords
        lat_south, lat_north = clat - latres/2, clat + latres/2
        lon_west, lon_east = 0 - lonres/2, 0 + lonres/2

        # create polygon
        poly = Polygon(

               LineString([
                   Point(lon_east, lat_south),
                   Point(lon_east, lat_north),
                   Point(lon_west, lat_north),
                   Point(lon_west, lat_south)
                   ])

            )

        area, _ = geod.geometry_area_perimeter(poly)
        areas.loc[clat] = area
        
    da = xr.DataArray(areas, dims={'latitude':latcentres})
    da = da.expand_dims({'longitude':loncentres})

    return da

area_weights = get_latlon_areas(gridto.latitude.values,
                                gridto.longitude.values,
                                latres=2.5/60, lonres=2.5/60)
    
# load ISO code mapper for WHO country names
isomap = pd.read_csv('./lookups/WHO_country_isocode_mapper.csv',
                     index_col=1).squeeze('columns')


class BaseLineHealthData:
    
    def __init__(self, config, HealthFunc=HealthFunc):
        self.config = config
        self.causemap = HealthFunc.gbd_mapper

        # load and clean baseline health data
        bhdf = pd.read_csv(self.config['bh_fpath'])
        agecol = bhdf['age_name'].str.replace(' to ', '-')
        agecol = agecol.str.replace(' plus', '+')
        bhdf['age_name'] = agecol
        bhdf = bhdf.set_index(['location_id', 'age_name', 'cause_name', 'measure_name'])
        bhdf = bhdf.sort_index()
        self.bhdf = bhdf
        
        self.isomap = pd.read_csv('./lookups/WHO_country_isocode_mapper.csv',
                             index_col=1).squeeze('columns')
     
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
    
    def __init__(self, year):
        with open(f'./grids/pop_count_slices_{year}.P', 'rb') as handle:
            self.popslices = pickle.load(handle)
        self.age_structure = pd.read_pickle(f'./lookups/age_structure_{year}.P')
        
    uncert_dict = {'lower':'lower',
                   'mid':'val',
                   'upper':'upper'}
        
    def lookup(self, country_isocode, age_group, uncert, uncert_dict=uncert_dict):
        age_group_proportion = self.age_structure.loc[(country_isocode, 
                                                       age_group, uncert_dict[uncert])]
        popslice = self.popslices[country_isocode].copy()
        age_group_population_sliced = popslice * age_group_proportion
        return age_group_population_sliced
    
    def population_array(self, country_isocode):
        return self.popslices[country_isocode]
    
    def get_countries_in(self, pm25):
        # check if pollution data contains countries
        countries_in = []
        for isocode, popslice in self.popslices.items():
            isin = any(popslice.longitude.isin(pm25.longitude)) & any(popslice.latitude.isin(pm25.latitude))
            if isin:
                countries_in.append(isocode)
        return countries_in
        

# instatiate datasets
bhd = BaseLineHealthData(config)
popdata = PopulationData(year=config['population_year'])


pm25 = xr.load_dataset('./grids/model_regridded_'+\
                       config['scenario_name']+'.nc')
pm25 = pm25[config['pm25var_name']]
if 'time' in pm25.dims:
    pm25 = pm25.mean('time')
# drop other coords 
if 'lev' in pm25.coords:
    pm25 = pm25.drop('lev')
    
# get a list of countries in the country mask
# crop to pm25 array
countries_in = popdata.get_countries_in(pm25)



def hiacalc(args):
    
    cause, uncert, age_group = args
    
    # print(cause, uncert, age_group)
    
    sr = pd.Series(index=pd.MultiIndex.from_product([countries_in, [age_group], [uncert]]),
                 name=cause, dtype=float)
                
    # get instance of health function for cause, age_group and uncert
    healthfunc = HealthFunc(cause=cause, age_group=age_group, uncert=uncert)
    hazard_ratio = healthfunc.calculate_hazard_ratio(pm25)
    
    da = xr.DataArray(.0, coords=gridto.coords, name=cause)
    
    for country_isocode in countries_in:
        if country_isocode not in isomap.index:
            continue
    
        
        baseline_deaths = bhd.lookup(country_isocode=country_isocode,
                                     cause=cause, 
                                     age_group=age_group, 
                                     measure='Deaths', 
                                     uncert=uncert)
        

        # get age group population
        age_group_pop = popdata.lookup(country_isocode, age_group, uncert)
     
        hazard_ratio_slice = hazard_ratio.where(age_group_pop)
        # calculate mortality rate from hazard ratio
        mortality_rate = (1 - 1 / hazard_ratio_slice) * baseline_deaths
        
        # calculate premature mortalities
        deaths = mortality_rate * age_group_pop / 100000
        deaths = deaths.where(deaths>0, 0)
        
        da.loc[deaths.coords] += deaths.values
        sr.loc[country_isocode, age_group, uncert] = float(deaths.sum())
        
    da = da.expand_dims({'uncert':[uncert], 'age_group':[age_group]})

    return da, sr
    

                
def hia_calculation():

    
    ### MAKE ITERABLES

    
    causes = HealthFunc.causes
    # get list of age ranges
    age_groups = HealthFunc.age_groups
    # get an iterable of uncertainty
    uncertainties = ['lower', 'mid', 'upper']
    
    
    with mp.Pool(mp.cpu_count()-2) as pool:
        fargs = list(product(causes, uncertainties, age_groups))
        outs = []
        for out in tqdm(pool.imap_unordered(hiacalc, fargs), total=len(fargs)):
            outs.append(out)
        
    # unpack results
    das, srs = map(list, zip(*outs)) 

    # combine gridded results
    mdas = []
    for cause in causes:
        cdas = [da for da in das if da.name == cause]
        das_summed = []
        for uncert in uncertainties:
            udas = [da.drop('uncert') for da in cdas if da.uncert == uncert]
            udas_sum = sum([da.values[0] for da in udas])
            da = xr.DataArray(udas_sum, coords={'uncert':[uncert], 'latitude':gridto.latitude, 'longitude':gridto.longitude})
            das_summed.append(da)
        uda = xr.concat(das_summed, dim='uncert')
        uda.name = cause
        mdas.append(uda)
        
    ds = xr.merge(mdas)
    comp = dict(zlib=True, complevel=5)
    encoding = {var: comp for var in ds.data_vars}
    ds.to_netcdf('./results/'+config['scenario_name']+'/gridded_results.nc', encoding=encoding)
            

    msrs = []
    for cause in causes:
        csrs = [sr for sr in srs if sr.name == cause]
        msrs.append(pd.concat(csrs).sort_index())
    results = pd.concat(msrs, axis=1)
    results.index.names = ['country', 'age_group', 'uncertainty']
    
    results.to_csv('./results/'+config['scenario_name']+'/by_country_age-group_uncertainty_results.csv')
    print('./results/'+config['scenario_name']+'/by_country_results.csv')

    # sum age groups
    results = results.groupby(['country', 'uncertainty']).sum()
    # replace country isocodes with names
    # hiadf.index = countries_lookup.loc[hiadf.index, 'name']
    # get mid uncertainty and drop level
    results = results.loc[:, 'mid', :]
    # hiadf.index = hiadf.index.droplevel(1)
    # convert to thousands
    results = results / 1000
        
    # calculate population weighted mean PM
    for country_isocode in countries_in:
        if country_isocode not in isomap.index:
            # print('skipping', country_isocode, 'due to no data')
            continue
        
        popslice = popdata.population_array(country_isocode)
        pm25slice = pm25.where(popslice)
        results.loc[country_isocode, 'mean PM2.5'] = pm25slice.weighted(area_weights.loc[pm25slice.coords]).mean()
        # slice popslice to pm25slice in case pm25 data doesn't cover entire country
        popslice = popslice.where(pm25slice)
        results.loc[country_isocode, 'population (within model area)'] = int(popslice.sum())
        popweighted = pm25slice.weighted(popslice).mean()
        results.loc[country_isocode, 'population-weighted PM2.5'] = float(popweighted)
        
    results.to_csv('./results/'+config['scenario_name']+'/by_country_results.csv')
    
#%%

if __name__ == "__main__":
    hia_calculation()
