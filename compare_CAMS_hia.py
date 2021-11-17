#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 10:31:42 2021

@author: eebjs
"""


import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
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


nofire = pd.read_csv('/nfs/a340/eebjs/acrobear/hia_results/CAMS_weighted_hia.csv', index_col=[0,1])

fire = pd.read_csv('/nfs/a340/eebjs/acrobear/hia_results/CAMS_hia.csv', index_col=[0,1])

countries_lookup = pd.read_csv('/nfs/a340/eebjs/hiadata/countries_isocodes.csv',
                               index_col='alpha-3')

def reformat_hiadf(hiadf):
    # convert to thousands
    hiadf = hiadf / 1000
    # sum age groups
    hiadf = hiadf.groupby(hiadf.index.get_level_values(0)).sum()
    # replace country isocodes with names
    # hiadf.index = countries_lookup.loc[hiadf.index, 'name']
    
    return hiadf
    
nofire = reformat_hiadf(nofire)
fire = reformat_hiadf(fire)
    


var = 'Non-accidental function (Non-Communicable + LRI deaths)'

# calculate diff
diff  = fire - nofire
# diff = diff.where(fire[var]>5).dropna()

# calculate percent difference
pdiff = ((fire / nofire) * 100) - 100
# pdiff = pdiff.where(fire[var]>5).dropna()

#%% plot bar comparison %

fig, ax = plt.subplots(figsize=[15,4])
pdiff[var].plot.bar()
ax.set_ylabel('change in premature mortalities\n attributed to fires (%)')
ax.set_xlabel('country')
ax.grid(axis='y')

fig.savefig('/nfs/b0122/eebjs/acrobear/figs/hia/fire_morts_%.png', dpi=300,
            bbox_inches='tight')

#%% plot bar absolute deaths

fig, ax = plt.subplots(figsize=[15,4])
diff[var].plot.bar()
ax.set_ylabel('change in premature mortalities \n attributed to fires (thousands)')
ax.set_xlabel('country')
ax.grid(axis='y')

fig.savefig('/nfs/b0122/eebjs/acrobear/figs/hia/fire_morts_abs.png', dpi=300,
            bbox_inches='tight')

#%% plot map comparison 



# calculate percent difference
pdiff = ((fire / nofire) * 100) - 100

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
bound = round(max(abs(pdiff[var].min()), abs(pdiff[var].max())))
norm = colors.Normalize(vmin=-100,vmax=100)
# make colourmap
cmap = cm.get_cmap('RdBu')
cmap, norm = discretise_cmap(cmap, norm)


for country in countries:

    country_isoname = country.attributes['ISO_A3']

    if country_isoname not in pdiff.index:
        continue
    
    mort = pdiff.loc[country_isoname, var]
    
    ax.add_geometries([country.geometry], ccrs.PlateCarree(),
                      facecolor=cmap(norm(mort)),
                      label='poo',
                      edgecolor='black',
                      linewidth=.25)
    
# set extent
ax.set_extent([-180, 180, 40, 80], crs=ccrs.PlateCarree())
    
# colorbar
cax = fig.add_axes([.25, .38, .5, .03]) # position of colorbar
cbar = colorbar.ColorbarBase(cax, cmap=cmap, norm=norm,
                             orientation='horizontal')

# get aspect and scale
figheight = fig.get_figheight()
figwidth = fig.get_figwidth()
fig.set_figheight(figheight*2.5)
fig.set_figwidth(figwidth*2.5)

fig.savefig('/nfs/b0122/eebjs/acrobear/figs/hia/fire_morts_%map.png', dpi=300,
            bbox_inches='tight')


#%% plot map comparison abs

# calculate diff
diff  = fire - nofire

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
bound = round(max(abs(diff[var].min()), abs(diff[var].max())))
norm = colors.Normalize(vmin=-bound,vmax=bound)
# make colourmap
cmap = cm.get_cmap('RdBu')
# cmap, norm = discretise_cmap(cmap, norm)


for country in countries:

    country_isoname = country.attributes['ISO_A3']

    if country_isoname not in diff.index:
        continue
    
    mort = diff.loc[country_isoname, var]
    
    ax.add_geometries([country.geometry], ccrs.PlateCarree(),
                      facecolor=cmap(norm(mort)),
                      edgecolor='black',
                      linewidth=.25)
    
# set extent
ax.set_extent([-180, 180, 40, 80], crs=ccrs.PlateCarree())
    
# colorbar
cax = fig.add_axes([.25, .38, .5, .03]) # position of colorbar
cbar = colorbar.ColorbarBase(cax, cmap=cmap, norm=norm,
                             orientation='horizontal')

# get aspect and scale
figheight = fig.get_figheight()
figwidth = fig.get_figwidth()
fig.set_figheight(figheight*2.5)
fig.set_figwidth(figwidth*2.5)

fig.savefig('/nfs/b0122/eebjs/acrobear/figs/hia/fire_morts_absmap.png', dpi=300,
            bbox_inches='tight')
