#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 30 10:45:07 2021

@author: eebjs
"""

import yaml
import os
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import cartopy.feature as cfeature
import pandas as pd
from cartopy.io.shapereader import Reader
import cartopy.io.shapereader as shpreader
from matplotlib import cm
from matplotlib import colors
from matplotlib import colorbar
import matplotlib.ticker as ticker

def discretise_cmap(cmap, norm, interval=5):
    cols = []
    ranges = ticker.MaxNLocator(20).tick_values(norm.vmin, norm.vmax)
    for b in ranges:
        cols.append(cmap(norm(b)))
    cmap = colors.ListedColormap(cols)
    norm = colors.BoundaryNorm(ranges, cmap.N)
    return(cmap, norm)


shape_path = '/nfs/see-fs-02_users/eebjs/OneDrive/AIA_project/scripts/shapefiles/'

FIRE_SCENARIO_NAME = 'vandonk_2005'
NOFIRE_SCENARIO_NAME = 'vandonk_2005_cesmscaled'


fire_config = yaml.safe_load(open('../config_files/'+FIRE_SCENARIO_NAME+'.yml'))
nofire_config = yaml.safe_load(open('../config_files/'+NOFIRE_SCENARIO_NAME+'.yml'))

# make a folder for results if there isn't one
results_fpath = '../results/'+FIRE_SCENARIO_NAME +'__VS__'+ NOFIRE_SCENARIO_NAME +'/'
if not os.path.exists(results_fpath):
    os.mkdir(results_fpath)
    
def read_hiadf(path):
    # open df
    hiadf = pd.read_csv(path, index_col=[0,1,2])

    # sum age groups
    hiadf = hiadf.groupby(['country', 'uncertainty']).sum()
    # replace country isocodes with names
    # hiadf.index = countries_lookup.loc[hiadf.index, 'name']
    
    # get mid uncertainty and drop level
    hiadf = hiadf.loc[:, 'mid', :]
    hiadf.index = hiadf.index.droplevel(1)
    
    # convert to thousands
    hiadf = hiadf / 1000
    
    return hiadf


#%% plot means

def plot_annual_pm_mean():
    
    proj = ccrs.Orthographic(central_latitude=90,
                                                  central_longitude=140)
    
    fig, axes = plt.subplots(ncols=2, 
                                subplot_kw={'projection':proj})
    
    
    for ax, scenario in zip(axes,[NOFIRE_SCENARIO_NAME, FIRE_SCENARIO_NAME]):
        
        config = yaml.safe_load(open('../config_files/'+scenario+'.yml'))
    
        da = xr.open_dataset(config['model_path'])[config['pm25var_name']]
        if 'time' in da.dims:
            da = da.mean('time')
        da.name = 'PM2.5 (ug/m3)'
        
        # da = da[0:2]

        ax.coastlines()
        conts = da.plot.contourf(ax=ax, 
                                 transform=ccrs.PlateCarree(), cmap='magma_r',
                                 levels=np.arange(0,55, 5),
                                 add_colorbar=False)
        
        ax.set_extent([-180, 180, 45, 90],
                      ccrs.PlateCarree())
        ax.add_feature(cfeature.BORDERS, lw=.5, zorder=2, edgecolor='black')
        
        ax.set_title(config['fire_scenario'])
        
    fig.subplots_adjust(wspace=.25)
    cbar_ax=fig.add_axes([0.49, 0.25, 0.02, 0.5])
    cbar_ax.set_title('$\\mathrm{PM_{2.5}\ \mu g\ m^{-3}}$',
                      fontsize=10)
    fig.colorbar(conts, cax=cbar_ax)
    
    
    # get aspect and scale
    figheight = fig.get_figheight()
    figwidth = fig.get_figwidth()
    fig.set_figheight(figheight*1.7)
    fig.set_figwidth(figwidth*1.7)
    
    fig.savefig(results_fpath+'mean_PM_comparison.png',
                dpi=300, bbox_inches='tight')
    
    
plot_annual_pm_mean()


#%% plot annual mean pm diff

da_fire = xr.open_dataset(fire_config['model_path'])[fire_config['pm25var_name']]
da_nofire = xr.open_dataset(nofire_config['model_path'])[nofire_config['pm25var_name']]
da = da_fire - da_nofire
if 'time' in da.dims:
    da = da.mean('time')
bound = round(float(max(da.quantile(.01),da.quantile(.99))), 0)
interval = (bound//10)*2
da.name = 'PM2.5 (ug/m3)'
fig = plt.figure()
ax = plt.subplot(projection=ccrs.Orthographic(central_latitude=90,
                                              central_longitude=140))



ax.coastlines()
# da = da[0:2]
da.plot.contourf(ax=ax, transform=ccrs.PlateCarree(), cmap='RdBu_r',
                 levels=np.arange(-10, 12, 2))

ax.set_extent([-180, 180, 45, 90],
              ccrs.PlateCarree())
ax.add_feature(cfeature.BORDERS, lw=.5, zorder=2, edgecolor='black')

ax.set_title('Increase in annual mean $\\mathrm{PM_{2.5}}$ attributed to fire')


# get aspect and scale
figheight = fig.get_figheight()
figwidth = fig.get_figwidth()
fig.set_figheight(figheight*1.7)
fig.set_figwidth(figwidth*1.7)

fig.savefig(results_fpath+'_mean_PM_diff.png',
            dpi=300, bbox_inches='tight')

#%% plot population count
def plot_popcount(SCENARIO_NAME):
    config = yaml.safe_load(open('../config_files/'+SCENARIO_NAME+'.yml'))
    
    da = xr.open_dataset('../grids/population_count.nc')[str(config['population_year'])]
    da.name = 'population count (thousands)'
    da = da.where(da != 0, np.nan)
    
    fig = plt.figure()
    ax = plt.subplot(projection=ccrs.Orthographic(central_latitude=90,
                                                  central_longitude=140))

    
    
    ax.coastlines()
    da.plot(ax=ax, transform=ccrs.PlateCarree(), cmap='bone_r',
                     vmax=1000)
    
    ax.set_extent([-180, 180, 45, 90],
                  ccrs.PlateCarree())
    ax.add_feature(cfeature.BORDERS, lw=.5, zorder=2, edgecolor='black')
    ax.add_feature(cfeature.OCEAN)
    ax.add_feature(cfeature.LAND)
    ax.set_title('Population')
    
    # get aspect and scale
    figheight = fig.get_figheight()
    figwidth = fig.get_figwidth()
    fig.set_figheight(figheight*1.7)
    fig.set_figwidth(figwidth*1.7)
    
    fig.savefig(results_fpath+SCENARIO_NAME+'_population_count.png',
                dpi=300, bbox_inches='tight')
    
plot_popcount(FIRE_SCENARIO_NAME)

#%% plot premature mortalities change

# cause = 'Non-accidental function (Non-Communicable + LRI deaths)'
# da_fire = xr.open_dataset('../results/'+FIRE_SCENARIO_NAME+'/gridded_results.nc')[cause]
# da_nofire = xr.open_dataset('../results/'+NOFIRE_SCENARIO_NAME+'/gridded_results.nc')[cause]
# da = da_fire - da_nofire
# da = da[1].drop('uncert')
# da = da.where(da != 0, np.nan)


# fig = plt.figure()
# ax = plt.subplot(projection=ccrs.Orthographic(central_latitude=90,
#                                               central_longitude=140))

# bound = abs(float(max(da.quantile(.1), da.quantile(.9))))
# ax.coastlines()
# da.plot.contourf(ax=ax, transform=ccrs.PlateCarree(), cmap='RdBu_r',
#                  levels=np.arange(-10, 12, 2))

# ax.set_extent([-180, 180, 45, 90],
#               ccrs.PlateCarree())
# ax.add_feature(cfeature.BORDERS, lw=.5, zorder=2, edgecolor='black')
# ax.add_feature(cfeature.OCEAN)
# ax.add_feature(cfeature.LAND)

# # get aspect and scale
# figheight = fig.get_figheight()
# figwidth = fig.get_figwidth()
# fig.set_figheight(figheight*1.7)
# fig.set_figwidth(figwidth*1.7)

# fig.savefig(results_fpath+'mortality_change.png',
#             dpi=300, bbox_inches='tight')

#%% absolute deaths map

cause = 'Non-accidental function (Non-Communicable + LRI deaths)'
fire_df = read_hiadf('../results/'+FIRE_SCENARIO_NAME+'/by_country_results.csv')
nofire_df = read_hiadf('../results/'+NOFIRE_SCENARIO_NAME+'/by_country_results.csv')

diff = (fire_df[cause] - nofire_df[cause]) 

fig = plt.figure()
ax = plt.subplot(projection=ccrs.Orthographic(central_latitude=90,
                                              central_longitude=140))

# load global countries shapefile
reader = Reader(shape_path+'ne_50m_admin_0_countries.shp')
shapes = reader.records()

shpfilename = shpreader.natural_earth(resolution='110m',
                                          category='cultural',
                                          name='admin_0_countries')
reader = shpreader.Reader(shpfilename)
countries = reader.records()

# create normalise object
bound = round(abs((diff.quantile(.99))))
norm = colors.Normalize(vmin=0,vmax=bound)
# make colourmap
cmap = cm.get_cmap('OrRd')
cmap, norm = discretise_cmap(cmap, norm)

ax.set_extent([-180, 180, 45, 90],
              ccrs.PlateCarree())


for country in countries:

    country_isoname = country.attributes['ISO_A3']

    if country_isoname not in diff.index:
        continue
    
    mort = diff.loc[country_isoname]
    
    ax.add_geometries([country.geometry], ccrs.PlateCarree(),
                      facecolor=cmap(norm(mort)),
                      edgecolor='black',
                      linewidth=.25)
    
# colorbar
cax = fig.add_axes([.25, .05, .5, .03]) # position of colorbar
cbar = colorbar.ColorbarBase(cax, cmap=cmap, norm=norm,
                             orientation='horizontal')
cax.set_title('annual premature deaths attributed to fire (thousands)')
    
figheight = fig.get_figheight()
figwidth = fig.get_figwidth()
fig.set_figheight(figheight*1.7)
fig.set_figwidth(figwidth*1.7)

fig.savefig(results_fpath+'fire_deaths_abs_map.png',
            dpi=300, bbox_inches='tight')

#%% absolute deaths map global

cause = 'Non-accidental function (Non-Communicable + LRI deaths)'
fire_df = read_hiadf('../results/'+FIRE_SCENARIO_NAME+'/by_country_results.csv')
nofire_df = read_hiadf('../results/'+NOFIRE_SCENARIO_NAME+'/by_country_results.csv')

diff = (fire_df[cause] - nofire_df[cause]) 

fig = plt.figure()
ax = plt.subplot(projection=ccrs.PlateCarree())

# load global countries shapefile
reader = Reader(shape_path+'ne_50m_admin_0_countries.shp')
shapes = reader.records()

shpfilename = shpreader.natural_earth(resolution='110m',
                                          category='cultural',
                                          name='admin_0_countries')
reader = shpreader.Reader(shpfilename)
countries = reader.records()

# create normalise object
bound = round(abs((diff.quantile(.999))))
norm = colors.Normalize(vmin=0,vmax=bound)
# make colourmap
cmap = cm.get_cmap('OrRd')
cmap, norm = discretise_cmap(cmap, norm)

ax.set_extent([-180, 180, -58, 90],
              ccrs.PlateCarree())


for country in countries:

    country_isoname = country.attributes['ISO_A3']

    if country_isoname not in diff.index:
        continue
    
    mort = diff.loc[country_isoname]
    
    ax.add_geometries([country.geometry], ccrs.PlateCarree(),
                      facecolor=cmap(norm(mort)),
                      edgecolor='black',
                      linewidth=.25)
    
# colorbar
cax = fig.add_axes([.25, .17, .5, .03]) # position of colorbar
cbar = colorbar.ColorbarBase(cax, cmap=cmap, norm=norm,
                             orientation='horizontal')
cax.set_title('annual premature deaths attributed to fire (thousands)')
    
figheight = fig.get_figheight()
figwidth = fig.get_figwidth()
fig.set_figheight(figheight*1.7)
fig.set_figwidth(figwidth*1.7)

fig.savefig(results_fpath+'fire_deaths_abs_map_global.png',
            dpi=300, bbox_inches='tight')

#%% relative fire deaths map

cause = 'Non-accidental function (Non-Communicable + LRI deaths)'
fire_df = read_hiadf('../results/'+FIRE_SCENARIO_NAME+'/by_country_results.csv')
nofire_df = read_hiadf('../results/'+NOFIRE_SCENARIO_NAME+'/by_country_results.csv')

df_fire = read_hiadf('../results/'+FIRE_SCENARIO_NAME+'/by_country_results.csv')[cause]
df_fire.name = 'fire attributed'
df_nofire = read_hiadf('../results/'+NOFIRE_SCENARIO_NAME+'/by_country_results.csv')[cause]
df_nofire.name = 'no fire'

# get the % fire attr
perfire =  ((df_fire - df_nofire) / df_fire ) * 100

fig = plt.figure()
ax = plt.subplot(projection=ccrs.Orthographic(central_latitude=90,
                                              central_longitude=140))

# load global countries shapefile
reader = Reader(shape_path+'ne_50m_admin_0_countries.shp')
shapes = reader.records()

shpfilename = shpreader.natural_earth(resolution='110m',
                                          category='cultural',
                                          name='admin_0_countries')
reader = shpreader.Reader(shpfilename)
countries = reader.records()

# create normalise object
bound = round(abs((perfire.quantile(.9))))
norm = colors.Normalize(vmin=0,vmax=bound)
# make colourmap
cmap = cm.get_cmap('OrRd')
cmap, norm = discretise_cmap(cmap, norm)

ax.set_extent([-180, 180, 45, 90],
              ccrs.PlateCarree())


for country in countries:

    country_isoname = country.attributes['ISO_A3']

    if country_isoname not in perfire.index:
        continue
        
        
    mort = perfire.loc[country_isoname]
    
    ax.add_geometries([country.geometry], ccrs.PlateCarree(),
                      facecolor=cmap(norm(mort)),
                      edgecolor='black',
                      linewidth=.25)
    
# colorbar
cax = fig.add_axes([.25, .05, .5, .03]) # position of colorbar
cbar = colorbar.ColorbarBase(cax, cmap=cmap, norm=norm,
                             orientation='horizontal')
cax.set_title('% of total premature deaths attributed to fire')
    
figheight = fig.get_figheight()
figwidth = fig.get_figwidth()
fig.set_figheight(figheight*1.7)
fig.set_figwidth(figwidth*1.7)

fig.savefig(results_fpath+'fire_deaths_rel_map.png',
            dpi=300, bbox_inches='tight')

#%% relative fire deaths map global

cause = 'Non-accidental function (Non-Communicable + LRI deaths)'
fire_df = read_hiadf('../results/'+FIRE_SCENARIO_NAME+'/by_country_results.csv')
nofire_df = read_hiadf('../results/'+NOFIRE_SCENARIO_NAME+'/by_country_results.csv')

df_fire = read_hiadf('../results/'+FIRE_SCENARIO_NAME+'/by_country_results.csv')[cause]
df_fire.name = 'fire attributed'
df_nofire = read_hiadf('../results/'+NOFIRE_SCENARIO_NAME+'/by_country_results.csv')[cause]
df_nofire.name = 'no fire'

# get the % fire attr
perfire =  ((df_fire - df_nofire) / df_fire ) * 100

fig = plt.figure()
ax = plt.subplot(projection=ccrs.PlateCarree())

# load global countries shapefile
reader = Reader(shape_path+'ne_50m_admin_0_countries.shp')
shapes = reader.records()

shpfilename = shpreader.natural_earth(resolution='110m',
                                          category='cultural',
                                          name='admin_0_countries')
reader = shpreader.Reader(shpfilename)
countries = reader.records()

# create normalise object
bound = round(abs((perfire.quantile(.99))))
norm = colors.Normalize(vmin=0,vmax=bound)
# make colourmap
cmap = cm.get_cmap('OrRd')
cmap, norm = discretise_cmap(cmap, norm)

ax.set_extent([-180, 180, -58, 90],
              ccrs.PlateCarree())


for country in countries:

    country_isoname = country.attributes['ISO_A3']

    if country_isoname not in perfire.index:
        continue
        
        
    mort = perfire.loc[country_isoname]
    
    ax.add_geometries([country.geometry], ccrs.PlateCarree(),
                      facecolor=cmap(norm(mort)),
                      edgecolor='black',
                      linewidth=.25)
    
# colorbar
cax = fig.add_axes([.25, .17, .5, .03]) # position of colorbar
cbar = colorbar.ColorbarBase(cax, cmap=cmap, norm=norm,
                             orientation='horizontal')
cax.set_title('% of total premature deaths attributed to fire')
    
figheight = fig.get_figheight()
figwidth = fig.get_figwidth()
fig.set_figheight(figheight*1.7)
fig.set_figwidth(figwidth*1.7)

fig.savefig(results_fpath+'fire_deaths_rel_map_global.png',
            dpi=300, bbox_inches='tight')

#%% plot premature mortalities abs for each scenario
def death_map(SCENARIO_NAME):
    config = yaml.safe_load(open('../config_files/'+SCENARIO_NAME+'.yml'))
    
    cause = 'Non-accidental function (Non-Communicable + LRI deaths)'
    da = xr.open_dataset('../results/'+SCENARIO_NAME+'/gridded_results.nc')[cause]
    da = da[1].drop('uncert')
    da = da.where(da != 0, np.nan)
    
    # bound = abs(round(float(max(da.quantile(.01), da.quantile(.99)))))

    fig = plt.figure()
    ax = plt.subplot(projection=ccrs.Orthographic(central_latitude=90,
                                                  central_longitude=140))

    ax.coastlines()
    da.plot(ax=ax, transform=ccrs.PlateCarree(), cmap='Reds',
            vmin=0,vmax=10)

    ax.set_extent([-180, 180, 40, 90],
                  ccrs.PlateCarree())
    ax.add_feature(cfeature.BORDERS, lw=.5, zorder=2, edgecolor='black')
    ax.add_feature(cfeature.OCEAN)
    ax.add_feature(cfeature.LAND)

    # get aspect and scale
    figheight = fig.get_figheight()
    figwidth = fig.get_figwidth()
    fig.set_figheight(figheight*1.7)
    fig.set_figwidth(figwidth*1.7)
    
    fig.savefig(results_fpath+SCENARIO_NAME+'_mortality.png',
                dpi=300, bbox_inches='tight')
  
death_map(FIRE_SCENARIO_NAME)
death_map(NOFIRE_SCENARIO_NAME)

#%% death bar chart

arctic_countries = ['RUS', 'USA', 'CAN', 'GRL', 'ISL', 'NOR', 'SWE', 'FIN']

cause = 'Non-accidental function (Non-Communicable + LRI deaths)'
df_fire = read_hiadf('../results/'+FIRE_SCENARIO_NAME+'/by_country_results.csv')[cause]
df_nofire = read_hiadf('../results/'+NOFIRE_SCENARIO_NAME+'/by_country_results.csv')[cause]

def deaths_bar(df_fire, df_nofire):
    
    df_fire = df_fire.loc[arctic_countries]
    df_fire.name = 'fire'
    df_nofire = df_nofire.loc[arctic_countries]
    df_nofire.name = 'no fire'
    
    
    fig, ax = plt.subplots(figsize=[6,6])

    ax.bar(df_nofire.index, df_nofire, width=0.6, label='non-fire attributed',
           color='grey', edgecolor='white')
    ax.bar(df_fire.index, df_fire-df_nofire, width=0.6, bottom=df_nofire.values,
           label='fire-attributed', color='#FF662F', edgecolor='white')
    ax.set_ylabel('annual excess deaths  (thousands)')
    
    # convert to percent
    # df_fire = df_fire/df_fire*100
    # df_nofire = df_nofire/df_fire*100
    
    # ax2 = ax.twinx()
    # ax2.bar(df_nofire.index, df_nofire, width=0.2, label='without fire',
    #        color='grey', edgecolor='white', align='edge')
    # ax2.bar(df_fire.index, df_fire-df_nofire, width=0.2, bottom=df_nofire.values,
    #        label='with fire', color='#FF662F', edgecolor='white',
    #        align='edge')
    # ax2.margins(0)
    # ax2.set_ylabel('proportion excess death attribution (%)')
    ax.legend()
    ax.grid(alpha=.4)
    
    fig.savefig(results_fpath+'bar_mortality.png',
                dpi=300, bbox_inches='tight')

deaths_bar(df_fire, df_nofire)

#%% fire only death bar chart

arctic_countries = ['RUS', 'USA', 'CAN', 'GRL', 'ISL', 'NOR', 'SWE', 'FIN']
cause = 'Non-accidental function (Non-Communicable + LRI deaths)'
df_fire = read_hiadf('../results/'+FIRE_SCENARIO_NAME+'/by_country_results.csv')[cause]
df_nofire = read_hiadf('../results/'+NOFIRE_SCENARIO_NAME+'/by_country_results.csv')[cause]

def deaths_bar(df_fire, df_nofire):
    
    df_fire = df_fire.loc[arctic_countries]
    df_fire.name = 'fire attributed'
    df_nofire = df_nofire.loc[arctic_countries]
    df_nofire.name = 'no fire'
    
    # get the diff
    df_fire = df_fire - df_nofire
    
    
    fig, ax = plt.subplots(figsize=[6,6])

    ax.bar(df_fire.index, df_fire, width=0.6,
           label='fire-attributed', color='#FF662F', edgecolor='white')
    ax.set_ylabel('annual excess deaths  (thousands)')
    

    ax.legend()
    ax.grid(alpha=.4)
    
    fig.savefig(results_fpath+'bar_fire_mortality.png',
                dpi=300, bbox_inches='tight')

deaths_bar(df_fire, df_nofire)

#%% proportion of deaths due to fire

arctic_countries = ['RUS', 'USA', 'CAN', 'ISL', 'NOR', 'SWE', 'FIN']

cause = 'Non-accidental function (Non-Communicable + LRI deaths)'
df_fire = read_hiadf('../results/'+FIRE_SCENARIO_NAME+'/by_country_results.csv')[cause]
df_nofire = read_hiadf('../results/'+NOFIRE_SCENARIO_NAME+'/by_country_results.csv')[cause]

def deaths_bar(df_fire, df_nofire):
    
    df_fire = df_fire.loc[arctic_countries]
    df_fire.name = 'fire attributed'
    df_nofire = df_nofire.loc[arctic_countries]
    df_nofire.name = 'no fire'
    
    # get the % fire attr
    perfire =  ((df_fire - df_nofire) / df_fire ) * 100
    
    fig, ax = plt.subplots(figsize=[6,6])

    ax.bar(perfire.index, perfire, width=0.6,
           label='percent fire-attributed', color='#FF662F', edgecolor='white')
    ax.set_ylabel('%')
    

    ax.legend()
    ax.grid(alpha=.4)
    
    fig.savefig(results_fpath+'bar_fire_%.png',
                dpi=300, bbox_inches='tight')

deaths_bar(df_fire, df_nofire)

#%% scatter fire vs non-fire deaths

arctic_countries = ['RUS', 'USA', 'CAN', 'ISL', 'NOR', 'SWE', 'FIN']

cause = 'Non-accidental function (Non-Communicable + LRI deaths)'
df_fire = read_hiadf('../results/'+FIRE_SCENARIO_NAME+'/by_country_results.csv')[cause]
df_nofire = read_hiadf('../results/'+NOFIRE_SCENARIO_NAME+'/by_country_results.csv')[cause]
colours = df_fire.where(df_fire.index.isin(arctic_countries), 'blue')
colours = colours.where(colours== 'blue', 'orange')

fig, ax = plt.subplots(figsize=[6,6])

ax.scatter(x=df_fire.values, y=(df_fire-df_nofire).values,
           color=colours)

ax.set_ylim(0,15)
ax.set_xlim(0,220)
