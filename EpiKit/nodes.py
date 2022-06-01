#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 00:50:16 2021

@author: diego
"""

import numpy as np
import matplotlib.pyplot as plt
import numpy.random as nprm
import networkx as nx



#from statistics import mode



class Agents_lite:
    def __init__(self,Nodes,Nt,Nq):
        self.S = Nodes.S
        self.N = Nodes.N
        self.Ns = Nodes.Ns
        self.Nr = Nodes.Nr
        self.ind_nodes = Nodes.ind_nodes
        self.infection_tree = Nodes.infection_tree
        self.node_category = Nodes.node_category
        self.nclass_teacher = Nodes.nodes_category_teacher
        self.biogroup = Nodes.biogroup

        self.symp_tested = Nodes.symp_tested
        self.tested = Nodes.tested
        self.time_isolated = Nodes.time_isolated
        self.time_detected = Nodes.time_detected
        self.days_quarantine = Nodes.days_quarantine
        self.tmax = Nodes.tmax

        self.Qcount = Nodes.Qcount

        self.number_of_tests = Nt
        self.quarantines_per_biogroup = Nq
        self.ExtInf = Nodes.ExtInf
        #self.external_infections = nodes.external_infections


class Agents:
    def __init__(self,S, S_inf, N, Nr, tmax):
        self.S = S
        self.Sinf = S_inf
        self.N = N
        self.Ns = len(S)
        self.Nr = Nr

        self.Sdict = {self.S[i]:i for i in range(self.Ns) }
        self.Dict_inf = [ {self.S[i]:self.Sinf[j][i] for i in range(self.Ns) }
                          for j in range(self.Sinf.shape[0]) ]

        self.ind_nodes = np.arange(self.N)
        #initial = nprm.choice(self.ind_nodes, I0 )
        #self.state = np.array(['S' if not( i in initial) else 'E' for i in range(N)],dtype='<U4')

        self.state = np.array(['S' for i in range(N)],dtype='<U4')
        #self.state[I0] = 'E'

        self.time = np.zeros((N)) + tmax+1
        self.rnext = np.zeros((N),int) + self.Nr
        self.tmax = tmax

        self.infection_tree =nx.DiGraph()
        self.infection_tree.add_nodes_from( self.ind_nodes  )

    def count_states(self):
         return np.array([ np.count_nonzero(self.state==s) for s in self.S ])




class AgentsNet: #DYN
    def __init__(self,Graphs,S, S_inf, sigma, Nr,tmax, precal=True):

        Gt = Graphs['Gt']
        meta = Graphs['meta']
        TeachersDict = Graphs['TeachersDict']
        N = meta.index.size
        #super().__init__(S, S_inf, N, Nr, tmax)
        self.shuffle = False
        self.tmax = tmax
        self.sigma = sigma
        self.S = S
        self.Sinf = S_inf
        self.N = N
        self.Ns = len(S)
        self.Nr = Nr
        self.Sdict = {self.S[i]:i for i in range(self.Ns) }
        self.Dict_inf = [ {self.S[i]:self.Sinf[j][i] for i in range(self.Ns) }
                          for j in range(self.Sinf.shape[0]) ]

        self.ind_nodes = np.arange(self.N)
        self.node_category = np.array(meta['class'])
        self.categories_of_teachers = [1, 3]


        for k in TeachersDict:
            self.node_category[k] = TeachersDict[k]

        self.unique_nclass = np.unique(self.node_category)

        self.nodes_category_teacher = np.array(meta['class'])

        self.biogroup = np.array([1 if self.nodes_category_teacher[i] == 'Teachers' else 0
                                             for i in self.ind_nodes] )

        self.susceptibility = np.array([ self.sigma[i] for i in self.biogroup] )
        self.tested = np.zeros( self.N )
        self.symp_tested  = np.zeros( self.N )
        self.positive_detected = np.zeros( self.N )
        self.time_isolated = np.zeros( (self.N, 2) )
        self.time_detected =  np.zeros( self.N ) + self.tmax +1
        self.days_quarantine = np.zeros( self.N )
        self.type_infection = 0 * self.biogroup

        self.Graphs = Graphs
        self.meta = meta
        #self.T = max(list(Gt[1].keys()) )
        self.dt_sec = Graphs['dt_sec']
        self.Ttot = Graphs.secf
        self.dt = self.dt_sec/(3600*24)
        self.dt0 = 20/(3600*24)

        self.ndays = len(Gt)

        self.Gt = Gt
        if precal:
           self.set_prob_infection()

        self.activitydays = (0,1,3,4)

        self.reset('S')

    def set_shuffle(self, val=True):
        self.shuffle = True

    def reset(self,I0):
        self.Qcount = []
        self.ExtInf = []
        self.qnodes = np.array(['F' for i in range(self.N)],dtype='<U4')
        self.qnext =  np.zeros((self.N,2),int) + self.Nr

        #self.external_infections = [I0]
        if type(I0) == str:
            self.state = np.array([I0 for i in range(self.N)],dtype='<U4')
        else:
            self.state = I0

        self.susceptibility = np.array([ self.sigma[i] for i in self.biogroup] )

        self.time = np.zeros((self.N)) + self.tmax + 1
        self.rnext = np.zeros((self.N),int) + self.Nr + 1
        self.tested = np.zeros((self.N),int)
        self.symp_tested = np.zeros((self.N),int)
        self.time_isolated = np.zeros((self.N,2))
        self.type_infection = 0 * self.biogroup

        self.infection_tree = nx.DiGraph()
        self.infection_tree.add_nodes_from( self.ind_nodes  )

        if self.shuffle:
            nprm.shuffle(self.Gt)

    def reset_node(self,node):
        self.rnext[node] = self.Nr + 1
        self.time[node] = self.tmax + 1
        self.state[node] = 'S'
        self.tested[node] = 0
        self.symp_tested[node] = 0
        self.time_isolated[node] = 0

    def activitytimes(self,t):
        return self.Graphs.activitytimes(t)

    def neighbors(self,Ii,it,d2):
        #print(dayindex_dataset,it)
        return self.Gt[d2][it].neighbors(Ii)

    def prob_infection2(self, iI,iS,it,d2):
        f_inf = self.Dict_inf[ self.biogroup[iI] ][ self.state[iI] ]
        sigma = self.susceptibility[ iS ]
        w = self.Gt[d2][it][iS][iI]['weight']
        return sigma * f_inf * w * self.dt0

    def prob_infection(self, iI,iS,it,d2):
        return self.Gt[d2][it][iS][iI]['p_inf'][iS][self.state[iI]]

    def set_prob_infection(self):
        for d in range(self.ndays):
          for it in self.Gt[d].keys():
            for i,j in self.Gt[d][it].edges:
                w = self.Gt[d][it][i][j]['weight']
                f_i = self.Dict_inf[ self.biogroup[i] ]
                f_j = self.Dict_inf[ self.biogroup[j] ]
                sigma_i, sigma_j = self.susceptibility[ i ] ,  self.susceptibility[ j ]
                # probability from j_susceptible i to infect j
                pij = {s: 1 - np.exp( - sigma_j * f_i[s] * w * self.dt0) for s in f_i }
                pji = {s: 1 - np.exp( - sigma_i * f_j[s] * w * self.dt0) for s in f_j }
                #pij = {s: sigma_j * f_i[s] * w * self.dt0 for s in f_i }
                #pji = {s: sigma_i * f_j[s] * w * self.dt0 for s in f_j }
                self.Gt[d][it][i][j]['p_inf'] = {i: pji, j: pij}

class Agents_TemporalNetwork(AgentsNet):
    def __init__(self,Graphs,S, S_inf, sigma, Nr,tmax):
        super().__init__(Graphs,S, S_inf, sigma, Nr,tmax, precal=False)
    def prob_infection(self, iI,iS,it,d):
        w = self.Gt[d][it][iI][iS]['weight']
        state = self.state[iI]
        f_I = self.Dict_inf[ self.biogroup[iI] ][state]
        sigma_S =   self.susceptibility[ iS ]
        return  1 - np.exp( - sigma_S * f_I * w * self.dt0)
       # return   sigma_S * f_I * w * self.dt0
      # return  1 - np.power(1-sigma_S*f_I , w * self.dt0)
# #%%
# p = 0.1
# x = np.linspace(0,10,100)
# plt.plot(x,p*x)
# plt.plot(x, 1 - np.power(1-p, x))
# #%%

class Agents_DailyNetworks(AgentsNet):
    def __init__(self,Graphs,S, S_inf, sigma, Nr,tmax):
        super().__init__(Graphs,S, S_inf, sigma, Nr,tmax, precal=False)
        self.dtT = self.dt_sec/self.Ttot
    def set_shuffle(self, val=True):
        self.shuffle = True

    def prob_infection(self, iI,iS,it,d):
        w = self.Gt[d][iI][iS]['weight'] * self.dtT
        state = self.state[iI]
        f_I = self.Dict_inf[ self.biogroup[iI] ][state]
        sigma_S =   self.susceptibility[ iS ]
        return  1 - np.exp( - sigma_S * f_I * w * self.dt0)
        #return   sigma_S * f_I * w * self.dt0
        #return  1 - np.power(1-sigma_S*f_I , w * self.dt0)

    def neighbors(self,Ii,it,d2):
        return self.Gt[d2].neighbors(Ii)



class Agents_ContactMatrix(AgentsNet):
    def __init__(self, Graphs, S, S_inf, sigma, Nr,tmax):
        super().__init__(Graphs,S, S_inf, sigma, Nr,tmax, precal=False)
       # super().__init__(None,S, S_inf, sigma, Nr,tmax, precal=False)
        self.dtT = self.dt_sec/self.Ttot
        self.ClassProps = {}

        for c in self.unique_nclass:
            nc = self.ind_nodes[self.node_category == c]
            sigma = self.susceptibility[nc]
            self.ClassProps[c] = (nc, sigma, nc.size)

    def reset(self,I0):
        super().reset(I0)
        self.ClassProps = {}
        for c in self.unique_nclass:
            nc = self.ind_nodes[self.node_category == c]
            sigma = self.susceptibility[nc]
            self.ClassProps[c] = (nc, sigma, nc.size)

    def set_shuffle(self, val=True):
        self.shuffle = True

    def prob_infection(self, iI,iS,it,d):
        return  1 #- np.exp( - sigma_S * f_I * w * self.dt0)

    def neighbors(self,iI,it,d2):
        nei = []
        cI  = self.nodes_category_teacher[iI]
        state = self.state[iI]
        f_I = self.Dict_inf[ self.biogroup[iI] ][state]
        for c in self.unique_nclass:
            w = self.Gt[d2][cI][c] * self.dtT
            #nc = self.ind_nodes[self.node_category==c]
            #sigma_S = self.susceptibility[nc]
            nc, sigma_S,N = self.ClassProps[c]
            prob = 1 - np.exp( - sigma_S * f_I * w * self.dt0)
            r = nprm.random(nc.size)
            nei += list(nc[prob>r])
            #nei += list(nprm.choice(nc,p))
        return nei


