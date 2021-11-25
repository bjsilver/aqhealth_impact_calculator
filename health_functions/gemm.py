#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 10:35:31 2021

@author: eebjs
"""


import pandas as pd
from hia import config
import numpy as np

# dictionary that maps the gemm function cause names to GBD cause names
gemm_to_gbd = {'Non-accidental function (Non-Communicable + LRI deaths)':['Non-communicable diseases', 'Lower respiratory infections'],
                'Ischaemic Heart Disease':'Ischemic heart disease',
                'Strokes':'Stroke',
                'Chronic Obstructive Pulomonary Disease':'Chronic obstructive pulmonary disease',
                'Lung Cancer':'Tracheal, bronchus, and lung cancer',
                'Lower Respiratory Infections':'Lower respiratory infections'}

class GEMM_Function():
    
    # load GEMM params and apply options
    gemm_params = pd.read_csv(config['gemm_params_fpath'], 
                              index_col=[0,1,2])
    gemm_params = gemm_params.sort_index()
    
    # exclude or include china (IN or EX)
    gemm_params = gemm_params.loc['All-regions, '+config['include_china']+'cluding China']
    
    # create iterators for age groups and causes
    causes = gemm_params.index.get_level_values('cause').unique()
    age_groups = gemm_params.index.get_level_values('age_group').unique()
    
    def __init__(self, cause, age_group, uncert, 
                 gemm_params=gemm_params):

        
        self.theta, self.theta_se, self.alpha, self.mu, self.tau =\
            gemm_params.loc[cause, age_group]
            
        # modify theta if calculating an uncertainty bound
        if uncert == 'lower':
            self.theta_uncert = self.theta - self.theta_se
        elif uncert == 'mid':
            self.theta_uncert = self.theta
        elif uncert == 'upper':
            self.theta_uncert = self.theta + self.theta_se
            
        
        
    def calculate_hazard_ratio(self, pm25):
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
        
        gamma = np.log(1 + z / self.alpha) / (1 + np.exp((self.mu - z) / self.tau))
        
        hazard_ratio = np.exp(self.theta * gamma)
        
        return hazard_ratio
