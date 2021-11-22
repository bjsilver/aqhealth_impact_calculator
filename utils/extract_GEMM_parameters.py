#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 10:37:13 2021

@author: eebjs
"""

import pandas as pd

df = pd.read_excel('/nfs/a340/eebjs/hiadata/gemm_function/raw/GEMM Calculator (PNAS)_ab.xlsx', sheet_name='GEMM fit parameters')

# the locations of the top right corners of the parameter tables in the sheet
ilocs = ((3, 1),  (3, 8),  (3, 15),  (3, 22),  (3, 29),  (3, 36),
         (24, 1), (24, 8), (24, 15), (24, 22), (24, 29), (24, 36))

#%% make a function to cut and format the table parts from the sheet
def extract_table(iloc):

    irow, icol = iloc
    
    includes_china = df.iloc[iloc]
    cause = df.iloc[irow+2, icol+1].replace('\n', ' ')
    
    # extract the table part of the sheet relative to iloc
    cdf = df.iloc[irow+7:irow+19, icol:icol+6]
    cdf.columns = ['age_group', 'theta', 'theta_se' , 'alpha', 'mu', 'tau']
    
    # correct mistake in xlsx file by changing 30-35 to 30-34
    cdf = cdf.where(cdf!='30-35', '30-34')
    
    cdf = cdf.set_index('age_group')
    
    # create mindex
    mindex = pd.MultiIndex.from_product([[includes_china], [cause], cdf.index],
                                        names=['includes_china?',
                                               'cause',
                                               'age_group'])
    cdf.index = mindex
    
    return cdf

#%%

cdfs = []
for iloc in ilocs:
    
    cdf = extract_table(iloc)
    cdfs.append(cdf)
    
gemm_parameters = pd.concat(cdfs)

gemm_parameters.to_csv('/nfs/a340/eebjs/hiadata/gemm_function/gemm_parameters.csv')
