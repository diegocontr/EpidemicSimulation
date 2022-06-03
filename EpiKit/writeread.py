#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 09:03:36 2021

@author: diego
"""

import numpy as np
import pandas as pd
import networkx as nx
import json
import gzip


def load_json_dataset(filename):
    with gzip.open(filename, 'r') as infile:
        jread = infile.read()
        jread= jread.decode('utf-8')
    js2 = json.loads( jread  )
    
    #js2 = load_json_dataset(filename)
    contact_meta = pd.read_json( js2['metadata'])    
    g0 = nx.Graph()
    g0.add_nodes_from( np.arange(214)  )
    
    if js2['type'] == 'daily':
        Gt = []
        for elist in js2['graphs']:
            g = g0.copy()
            g.add_weighted_edges_from(elist)
            Gt.append(g)
        return Gt,  contact_meta, None , js2['type']    
            
    elif js2['type'] == 'temporal':
        Gt = []
        for gday in js2['graphs']:
            graphs_day = {}
            for k in gday:
                elist = gday[k]
                g = g0.copy()
                g.add_weighted_edges_from(elist)
                graphs_day[int(k)] = g
            Gt.append(graphs_day)   
        dt = js2['dt'] 
        return Gt,  contact_meta, dt, js2['type'] 
            
    elif js2['type'] == 'matrix':
        Gt = js2['graphs']      
        return Gt,  contact_meta, None , js2['type'] 

# #%%
# {
#  'type':str,
#  'metadata': str( pandas.to_json) ,
#  'graphs':list,
#  'dt':float,
#  }
# #%%

def format_results_out(sim, Lres:list, Params_list:list, level=0, test = True):

    if level == 2:
        keys_saving = ['biogroup', 'initial_time', 'final_size',  'R0', 
                       'states', 'Gedges', 'Gnodes', 'external_infections', 
                       'size(t)']
        if test:
            keys_saving += ['#test', '#qdays', 'quarantined(t)', 'test(t)',  'dict_qtime']
        sim_properties =  {'time': Lres[0].agg_in_time['time'].tolist() ,
                            'states': Lres[0].S,
                              }
        
    elif level == 1:
        keys_saving = ['initial_time', 'final_size',   'R0',
                        'Gedges', 'Gnodes', 'dict_qtime', 'external_infections', 
                       'size(t)']
        if test:
            keys_saving += ['#test', '#qdays', 'dict_qtime']  
        sim_properties =  {'time': Lres[0].agg_in_time['time'].tolist() ,
                            'states': Lres[0].S,
                              }
    elif level == 0:
        keys_saving = [ 'final_size', 'R0']
        if test:
            keys_saving += ['#test', '#qdays']  
        sim_properties =  {   }

    joinParams = {}
    for p in Params_list:
        joinParams.update(p)
    Param_out = {k: ( v if isinstance(v, (int, float)) else list(v)  )
                                    for k,v in joinParams.items() }

    res_dict = []
    for r in Lres:
        
        G = r.infection_tree
        for u in G.nodes:
            if not 'time' in G.nodes[u]:
                G.nodes[u]['time'] = np.nan
        R = [G.out_degree(u)  for u in r.external_infections  if G.nodes[u]['time']==0 ]
        
            
            
        test_stat0 =  {'biogroup': r.biogroup.tolist() ,
                      'initial_time': r.OutResult['time0'],  
                       'final_size': r.OutResult['size'],      
                       '#test': int( r.number_of_tests ),
                       '#qdays': float(r.days_in_quarantine), 
                       'R0': np.mean(R),
                       'states': r.agg_state_in_time.tolist(),
                       'Gedges': [ (int(i), int(j),v['time']) for (i,j),v in dict( G.edges ).items() ],
                       'Gnodes': [ ( int(k), v['time'] )
                                   for k,v in dict( G.nodes ).items() ],
                       'dict_qtime': np.array(r.quarantines_per_biogroup).tolist() ,
                       'external_infections':  np.array(r.external_infections).tolist(),
                       'size(t)': r.agg_in_time['size'].tolist(),
                       'quarantined(t)': r.agg_in_time['#qdays'].tolist(),
                       'test(t)': r.agg_in_time['#test'].tolist(),
                      }
        test_stat = {k:test_stat0[k] for k in keys_saving}
        res_dict.append( test_stat )
        
    return {'Param': Param_out, 
            'sim': sim_properties, 
            'runs': res_dict,
            'level': level,
            }    
#%%
# def format_results_out2(level=0):
#     if level == 2:
#         keys_saving = ['biogroup', 'initial_time', 'final_size', '#test', '#qdays', 'R0', 
#                        'states', 'Gedges', 'Gnodes', 'dict_qtime', 'external_infections', 
#                        'size(t)', 'quarantined(t)', 'test(t)']
#         sim_properties =  {'time': 1,
#                             'states': 1,
#                               }
        
#     elif level == 1:
#         keys_saving = ['initial_time', 'final_size', '#test', '#qdays',  'R0',
#                         'Gedges', 'Gnodes', 'dict_qtime', 'external_infections', 
#                        'size(t)']
#         sim_properties =  {'time': 1 ,
#                             'states': 1,
#                               }
#     elif level == 0:
#         keys_saving = [ 'final_size', '#test', '#qdays', 'R0']
#         sim_properties =  {   }


#     res_dict = []
#     for i in range(1):
        
#         test_stat = {k:1 for k in keys_saving}
#         res_dict.append( test_stat )
        
#     return {'Param': 1, 
#             'sim': sim_properties, 
#             'runs': res_dict,
#             'level': level,
#             }    
# format_results_out2(level=0)

#%%


def struct_to_json_file(filename, struct):
    js = json.dumps( struct  , separators=(',', ':') )
    json_bytes = js.encode('utf-8')
    with gzip.open(filename, 'w') as fout:
        fout.write(json_bytes)


def results_from_json_file(filename):
    if 'gzip' in filename or 1:
        with gzip.open(filename, 'r') as infile:
            jread = infile.read()
            jread= jread.decode('utf-8')
    else:
        with open(filename, 'r') as infile:
            jread = infile.read()

    js = json.loads( jread  )
    return rebuild_results_from_json(js)



