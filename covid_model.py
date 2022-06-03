#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 31 06:48:21 2022

@author: diego
"""
import numpy as np
import EpiKit as epi
import dictfit

S         =            ['S',  'E',  'Ip',  'Ic',  'Isc', 'Rp', 'R']
dict_of_rel_inf   = {
    'adult': np.array([0, 0, 0.55, 1.00, 0.55, 0, 0, ]),
    'children': np.array( [ 0,    0,0.55*0.63, 0.63,0.55*0.63,  0,   0, ]),
    }


def beta_from_R0(R0, data, net, r_period, weekdays, generation=1):
    if  data == 'office' and weekdays==7:
        DictFit = dictfit.get_dict['office7']
    elif data == 'office' and r_period!=1:
        Dkey = {0.25:'office_r025',0.5:'office_r05',1.5:'office_r15',}
        DictFit = dictfit.get_dict[Dkey[args_read.r_period]]
    else:
        DictFit = dictfit.get_dict[ data]
    
    if net in DictFit[generation].keys():
        A,f = DictFit[generation][net]
    else:
        print('---> net {} not in the list'.format(args_read.net))
        A,f = DictFit[1]['DYN']
    return -np.log( 1 - R0/A)/ (f  /(60*24) )

def get_disease_model(beta,ParamSim, data):
    if ParamSim['vax'] == '75p':
        ParamVax = {'f_sigma': 1-0.73,
                    'f_c': 1-0.70,
                    'f_r': 1-0.24,
                    }
    
    elif ParamSim['vax'] == '93p':
        ParamVax = {'f_sigma': 1-0.85,
                    'f_c': 1-0.93,
                    'f_r': 1-0.5,
                    }
    if data in ['office','hospital']:
        relative_infectiousness = np.array(
              [dict_of_rel_inf['adult'],
               dict_of_rel_inf['adult'] * ParamVax['f_r']])
    elif data == 'school':
        relative_infectiousness = np.array(
              [dict_of_rel_inf['children'],
               dict_of_rel_inf['adult'] * ParamVax['f_r']])        
  
    S_inf = relative_infectiousness * beta 
    r_period = ParamSim['r_period'] 
    Param = {'epsilon' : 1/4,
              'mu_p' : (1/ r_period ) * 1/1.8,
              'mu' :  (1 /r_period )* 1/5,
              'mu_R' : 1/17.7,
              'p_sc' : np.array([0.5, 1 - 0.5 * ParamVax['f_c'] ]),
              'sigma' : np.array([1.0, 1.0 * ParamVax['f_sigma']  ]),
              'weekden': ParamSim['weekden'] ,
              'beta': beta,
              'g_epsilon' : 3,
              'g_mu_p' : 1,
              'g_mu' :   6,
              'g_mu_R' : 2,
              'r_period':r_period ,
              }
    Param.update(ParamVax)
    React = [ epi.Reaction('S',  'E', beta ), #0
              epi.Reaction('E', 'Ip', Param['epsilon'], gamma = Param['g_epsilon'] ), #1
              epi.Reaction('Ip',                          #2
                               {'Isc' : Param['p_sc'] ,
                                 'Ic' : 1 - Param['p_sc']},
                                        Param['mu_p'], gamma = Param['g_mu_p'] ),
              epi.Reaction('Ic',  'Rp', Param['mu'], gamma = Param['g_mu']  ), #3
              epi.Reaction('Isc', 'Rp', Param['mu'], gamma = Param['g_mu']  ), #4
              epi.Reaction('Rp',  'R',  Param['mu_R'],  gamma = Param['g_mu_R'] ),] #5


    Dict_states = {'healthy':'S',
                   'after_infection':'E',
                   'recovered':'R'}


    return Param, React, S_inf,Dict_states


    # S         =            ['S',  'E',  'Ip',  'Ic',  'Isc', 'Rp', 'R']
    # S_inf0    = np.array([ [ 0,    0,0.55*0.63, 0.63,0.55*0.63,  0,   0, ],
    #                        [ 0,    0,   0.55,   1.00,   0.55,    0,   0, ],
    #                        ])
    # SIGMA_STUD = 0.5
    # PROD_DET = 0.3
# Param_out, React, Sinf = get_covid_model(beta,Param)



