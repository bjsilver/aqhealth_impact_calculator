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


bhdpath = '/nfs/a340/eebjs/hiadata/baseline_health_data/raw/'

acm = pd.read_csv(bhdpath+'/nfs/a340/eebjs/hiadata/baseline_health_data/raw')