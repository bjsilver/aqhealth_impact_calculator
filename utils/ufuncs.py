#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 10:02:28 2021

@author: eebjs
"""

import pandas as pd


def read_hiadf(path):
    # open df
    hiadf = pd.read_csv(path, index_col=[0,1,2])
    # convert to thousands
    hiadf = hiadf / 1000
    # sum age groups
    hiadf = hiadf.groupby(hiadf.index.get_level_values(0)).sum()
    # replace country isocodes with names
    # hiadf.index = countries_lookup.loc[hiadf.index, 'name']
    
    return hiadf