#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 00:50:16 2021

@author: diego
"""

import numpy as np
import matplotlib.pyplot as plt
import numpy.random as nprm

from .nodes import AgentsNet, Agents_TemporalNetwork, Agents_DailyNetworks, Agents_ContactMatrix



class Result:
    def __init__(self,sim,times,jump_hours, saving_level=2):
        self.saving_level = saving_level
        self.saving_test = True

        self.jump = int(jump_hours * (3600 // sim.Nodes.dt_sec))  # hours * ( )
        self.times_out = times[::self.jump]
        nodes_zeros = np.zeros((self.times_out.size, sim.N), dtype='<U4')
        agg_zeros = np.zeros(self.times_out.size, int)


        self.t_index = 0
        
        self.S = sim.S
        self.biogroup = sim.Nodes.biogroup
        #self.df_in_time = pd.DataFrame(columns=['time','size','#test','#qdays'])
        self.agg_in_time = {k:agg_zeros.copy() for k in ['time','size','#test','#qdays']}
        self.state_history  = nodes_zeros.copy()
        self.quarantined_history  =   nodes_zeros.copy()
    def evolution_all(self):
        Xt = []
        for nstates in self.state_history:
            n = np.array([np.count_nonzero(nstates == s)
                         for s in self.S])
            Xt.append(n)
        return np.array(Xt)
    
    def update(self, sim, t):
        size = sim.OutResult['size']
        nt, qt = sim.Test.get_test_quarantine() 
       # dqt = [sum([(self.quarantines_per_biogroup[a][Q]) for Q in ['Qw', 'Qpi', 'Qi']])
        for k,val in zip(['time','size','#test','#qdays'], (t,size,nt, sum(qt) )):
            self.agg_in_time[k][self.t_index] = val
        if self.saving_level >0:
            self.state_history[self.t_index] = sim.Nodes.state.copy()
            self.quarantined_history[self.t_index] = sim.Nodes.qnodes.copy()
        
        self.t_index += 1 
        
    def finish(self, sim):
        self.OutResult = sim.OutResult
        self.external_infections = sim.Nodes.ExtInf
        self.infection_tree = sim.Nodes.infection_tree
        nt, qt = sim.Test.get_test_quarantine() 
        self.number_of_tests = nt #self.agg_in_time['#test'][-1]
        self.days_in_quarantine = sum(qt) #self.agg_in_time['#qdays'][-1]        
        if self.saving_level >1:
            self.agg_state_in_time = self.evolution_all()
            self.quarantines_per_biogroup = sim.Test.quarantines_per_biogroup
        
    def plot_t(self):
        tg = self.times_out
        Xt = self.evolution_all()
        plt.plot(tg, Xt[:, 0], 'C2-', label='S')
        plt.plot(tg, Xt[:, 1], 'C1-', label='E')
        plt.plot(tg, Xt[:, 4] + Xt[:, 3] + Xt[:, 2], 'C3-', label='I')
        plt.plot(tg, Xt[:, -1] + Xt[:, -2], 'C0-', label='R')
        plt.legend()


class Reaction:
    def __init__(self,r, p, k, gamma=None):
        self.index = 0
        self.r = r
        if gamma is None:
            self.time_distribution = lambda: nprm.exponential(k)
        else:
            tau = 1 / (k * gamma)
            self.time_distribution = lambda: nprm.gamma( gamma, tau)
        if isinstance(p, str):
            self.complex = False
            self.p = p
        if isinstance(p, dict):
            self.complex = True
            self.p = np.array([k for k in p.keys()])
            self.prob = np.array([p[k] for k in p.keys()])
            if self.prob.shape[1] > 1:
                self.accprob = [np.cumsum(self.prob.T[i])
                                for i in range(self.prob.shape[1])]
            else:
                self.accprob = [np.cumsum(self.prob)]
        self.k = k

    def product(self, nclass=np.nan):
        if self.complex:
            r = nprm.random()
            # this is a mess
            accprob = self.accprob[nclass] if nclass == nclass else self.accprob
            #print( accprob )
            ip = accprob[accprob < r].size
            return self.p[ip]
        else:
            return self.p
    def __repr__(self):
        return '{} -> {}   tau={}'.format(self.r, self.p, 1 / self.k)


# a = [Model.gR.R[2].time_distribution() for i in range(20000)]; print(np.mean(a), np.std(a) )

class Set_of_Reactions:
    def __init__(self, parent, React, S, tmax):
        self.parent = parent
        self.R = React
        self.S = S
        self.Ns = parent.Ns
        self.Nr = parent.Nr
        #self.time_distribution = nprm.exponential

        self.Sdict = {self.S[i]: i for i in range(self.Ns)}
        self.rates = np.array([r.k for r in self.R])

        self.tmax = tmax

        for ir in range(self.Nr):
            self.R[ir].index = ir
        self.prod_rnext = {r.r: r.index for r in self.R}
        rnext_keys = list(self.prod_rnext.keys())
        for s in S:
            if s not in rnext_keys:
                self.prod_rnext[s] = False

    def get_prop(self, node, t):
        return self.rates

    def infect_r(self, node, rnow):
        reaction = self.R[rnow]
        p = reaction.product()
        rnext = self.prod_rnext[p]
        reaction_next = self.R[rnext]
        time = reaction.time_distribution(1 / reaction_next.k)
        return p, rnext, time

    def get_rnext(self, state):
        rnext = self.prod_rnext[state]
        if rnext:
            reaction = self.R[rnext]
            time = reaction.time_distribution()
            return rnext, time
        else:
            return False, self.tmax + 1

    def set_and_rnext(self, inode):
        state = self.parent.Nodes.state[inode]
        return self.get_rnext(state)

    def product_and_rnext(self, inode, rnow):
        reaction = self.R[rnow]
        nclass = self.parent.Nodes.biogroup[inode]
        p = reaction.product(nclass)
        rnext, time = self.get_rnext(p)
        return p, rnext, time

class NoTest:
    def __init__(self):
        pass
    def reset(self, sim=None):
        pass
    def testing(self, sim, nodes, t, activity):
        pass
    def get_test_quarantine(self):
        return 0,0
    def counting_qteacher(self):
        return 0, 0



class AgentModel:

    def __init__(self, Graphs, Param:dict, React:list,
                 S:list, S_inf:np.array, S_dict:dict, N:int, tmax: float):
        '''
        AgentModel integrate all components of simulations.

        Parameters
        ----------
        Graphs : GraphData
            GraphData instance with the network to use.
        Param : dict
            Dictionary with parameters name as key, and its values.
            It requieres the keys: 'beta' and 'sigma'.
        React : list
            List of Reaction object describing possible transitions.
        S : list
            List of possible states.
        S_inf : np.array
            Relative infectiousness for state and biogroup.
        S_dict : dict
            Dictionary of key states. It must have the keys 'healthy',
            'after_infection', and 'recovered'.
        N : int
            Population number of agents.
        tmax : float
            Maximum time the simulation can reach.

        '''
        self.S_inf = S_inf
        self.S = S
        self.healthy_state = S_dict['healthy']
        self.state_after_infection =  S_dict['after_infection']
        self.recovered_state =  S_dict['recovered']
        self.tmax = tmax
        self.stop = 'at_tmax' # []
        self.N = N
        self.Ns = len(S)
        self.Nr = len(React)
        self.gR = Set_of_Reactions(self, React, S, tmax)
        self.vax_selection = False
        self.Test = NoTest()
        self.beta = Param['beta']
       # self.OutResult = {'time0': 0}

        if Graphs.type == 'daily':
            self.Nodes = Agents_DailyNetworks(Graphs,
                                              self.S, self.S_inf,
                                              Param['sigma'],
                                              self.Nr, self.tmax)

        if Graphs.type == 'matrix':
            self.Nodes = Agents_ContactMatrix(Graphs,
                                              self.S, self.S_inf,
                                              Param['sigma'],
                                              self.Nr, self.tmax)

        elif Graphs.type == 'temporal':
            self.Nodes = Agents_TemporalNetwork(Graphs,
                                                self.S, self.S_inf,
                                                Param['sigma'],
                                                self.Nr, self.tmax)

        self.weekden = Param['weekden']
        
    def set_test_protocol(self, Test):
        self.Test = Test
        
    def change_beta(self,beta_new:float):
        '''
        Change the value of the infection rate beta to a new value.

        Parameters
        ----------
        beta_new : float
            New value of the infection rate beta

        '''
        beta_old = self.beta
        self.S_inf = self.S_inf * beta_new/ beta_old
        self.Nodes.Dict_inf = [ {self.S[i]:self.S_inf[j][i] for i in range(self.Ns) }
                          for j in range(self.S_inf.shape[0]) ]
        self.beta = beta_new
        
    def find_transmissions(self, i_infected: int, t:float, it:int, dayindex_dataset:int):
        '''
        Search for healthy neighbors of i_infected and evaluate if they should
        get infected.

        Parameters
        ----------
        i_infected : int
            Index of infected node.
        t : float
            Current time.
        it : int
            Index of the current time.
        dayindex_dataset : int
            Index of the day in the dataset.

        '''
        for iS in self.Nodes.neighbors(i_infected, it, dayindex_dataset):
            # print('>>',iS)
            if self.Nodes.state[iS] == self.healthy_state:
                tq0, tqf = self.Nodes.time_isolated[iS]
                if ( not(tq0 < t < tqf) and
                     nprm.random() < self.Nodes.prob_infection(i_infected, iS, it, dayindex_dataset)):
                        self.set_infection(i_infected, iS, t)  # ,it,dayindex_dataset)


    def set_infection(self, i_infected:int, j_susceptible:int, t:float):
        '''
        Change the state of the node i_susceptible to the product of the first 
        reaction in set_of_reactions object.
        Updates the infection tree accordingly.

        Parameters
        ----------
        i_infected : int
            Index of infected node.
        j_susceptible : int
            Inex of the node to be infected.
        t : float
            Time of the infection.

        '''
        
        if i_infected == i_infected:
            self.Nodes.infection_tree.add_edge(i_infected, j_susceptible, time=t)
            self.Nodes.infection_tree.nodes[j_susceptible]['time'] = t
        product, rnext, time = self.gR.product_and_rnext(j_susceptible, 0)

        self.Nodes.time[j_susceptible] = t + time + 15 / (60 * 24)
        self.Nodes.rnext[j_susceptible] = 1
        self.Nodes.state[j_susceptible] = product
        self.OutResult['size'] += 1

    def change(self, i_node: int, rnow:int, t:float):
        '''
        Change the state of a node to the product of a reaction.

        Parameters
        ----------
        i_node : int
            index of the node to change.
        rnow : int
            index of the reaction which produces the new state.
        t : float
            Current time.

        '''
        product_now, rnext, time_next = self.gR.product_and_rnext(i_node, rnow)
        self.Nodes.time[i_node] = t + time_next
        self.Nodes.rnext[i_node] = rnext
        self.Nodes.state[i_node] = product_now

    def set_initial_condition(self, I0):
        '''
        Function that takes a numpy array and assign the state prescribed to
        each node. This functions first reset the nodes, and then go through
        each node assigning the state, and time of transition. The infection 
        tree and the member Nodes.ExtInf is updated with the nodes in the 
        state_after_infection.

        Parameters
        ----------
        I0 : np.array
            Array with the state to assign to each node.

        Returns
        -------
        None.

        '''
        self.Nodes.reset(I0)
        for inode in self.Nodes.ind_nodes:
            state = self.Nodes.state[inode]
            if state != self.healthy_state:
                rnext, time = self.gR.set_and_rnext(inode)
                # time= 2/24
                self.Nodes.time[inode] = time
                self.Nodes.rnext[inode] = rnext
                if state == self.state_after_infection:
                    self.Nodes.infection_tree.add_node(inode)
                    self.Nodes.infection_tree.nodes[inode]['time'] = 0
                    self.Nodes.ExtInf.append(inode)
                #print  (state,rnext,time)

    def set_node_state(self, inode, state):
        if state == self.healthy_state:
            self.Nodes.reset_node(inode)
        else:
            self.Nodes.state[inode] = state
            rnext, time = self.gR.set_and_rnext(inode)
            self.Nodes.time[inode] = time
            self.Nodes.rnext[inode] = rnext

    def set_introductions(self, time_0: float):
        '''
        Creates the memeber self.wtime_noactivity to introduce new infections.

        Parameters
        ----------
        time_0 : float
            Initial time of the simulation.

        '''
        self.ext_infection_time = 0  # 7*self.weekden #old simulations
        self.new_ext_inf = self.tmax + 1

        wtimes = np.arange(0, 7, self.Nodes.dt)
        wact, wtindex, wday2 = self.Nodes.activity_times(wtimes + time_0)
        self.wtimes_noactivity = wtimes[~wact]

    def check_introductions(self, t: float):
        '''
        Check if at this time corresponds to introduce a new infection. 
        If so, choses a new healthy node to infect.

        Parameters
        ----------
        t : float
            Current time in the simulation.

        '''    
        if t >= self.ext_infection_time:
            self.ext_infection_time += self.weekden * 7
            itime = nprm.choice(self.wtimes_noactivity, 1)[0]
            week = nprm.randint(0, self.weekden)
            self.new_ext_inf = t + 7 * week + itime
        if t >= self.new_ext_inf:
            self.new_ext_inf = self.tmax + 1
            Snodes = self.Nodes.ind_nodes[(self.Nodes.state == self.healthy_state) * np.logical_or(
                (self.Nodes.biogroup == 0), (self.Nodes.biogroup == 1))]
            if Snodes.size > 0:
                if self.vax_selection:
                    prob = self.Nodes.susceptibility[Snodes]
                    prob = prob / prob.sum()
                    inode = nprm.choice(Snodes, 1, p=prob)[0]
                else:
                    inode = nprm.choice(Snodes, 1)[0]
                self.set_infection(np.nan, inode, t)  # , 0, 0)
                self.Nodes.infection_tree.nodes[inode]['time'] = t
                self.Nodes.ExtInf.append(inode)
            else:
                self.ext_infection_time = self.tmax + 1

    def check_if_stop(self):
        '''
        Using the member value self.stop, it evaluates if the simulation should
        stop.

        Returns
        -------
        bool
            True if the condition to end the simulation is satisfied.

        '''
        if self.stop == 'at_tmax':
            return False
        elif self.stop == 'no_trans':
            return self.Nodes.time.min() > self.tmax     
        elif self.stop == 'end_gen1':  
            inode = self.Nodes.ExtInf[0]
            return self.Nodes.state[inode] == self.recovered_state 
            
    def solve(self, set_condition, jump_hours:float=24.0, saving_level:int = 2 ) -> Result:
        '''
        Function perform the time discrete simulation.

        Parameters
        ----------
        set_condition : function or np.array 
            Prescribe the initial state  of the nodes.
        jump_hours : float, optional
            Frequency in hours when to save the state of the system. 
            The default is 24.0.
        saving_level : int, optional    
            Level of detail to save. 
            The default is 2.  
        Returns
        -------
        Out : Result
            Data of the simulation result.

        '''
        # pick a random time
        times = np.arange(0, self.tmax + self.Nodes.dt, self.Nodes.dt)
        if self.Nodes.Graphs.net_cycle > 0:
            time_0 = nprm.choice(times[times < self.Nodes.Graphs.net_cycle])
        else:
            time_0 = 0

        self.times_displaced = times + time_0

        # choose vaccinated and initial state
        # if set_condition is a function we evaluate it and then asign it
        if callable(set_condition):
            s0 = set_condition(self)
            self.set_initial_condition(s0)
        # if it is an array, we asign it directly    
        elif isinstance(set_condition, np.array) :
            self.set_initial_condition(set_condition)  
        else:
            raise TypeError('Not a valid type for <set_condition> parameter. It must be a function or an numpy array.')

        self.OutResult = {'size': self.Nodes.ind_nodes[self.Nodes.state == self.state_after_infection].size,
                          'time0': time_0}

        self.Test.reset(sim=self)

        activitytime, tindex, dayindex_in_graph = self.Nodes.activity_times(self.times_displaced)
        Out = Result(self, times, jump_hours, saving_level = saving_level)
        jump = Out.jump

        to_continue = True

        self.set_introductions(time_0)

        for t, activity, dayindex, it, ind in zip(
                times, activitytime, dayindex_in_graph, tindex, range(times.size)):
            self.evol(t, activity, dayindex, it, to_continue)
            if ind % jump == 0:
                Out.update(self,t)
                if self.Nodes.time.min() > self.tmax:
                    to_continue = False
                else:
                    to_continue = True
                    
                if self.check_if_stop():
                    break
                    
        Out.finish(self)            
        return Out

    def evol(self, t:float, activity:bool, dayindex:int, it:int, to_continue:bool):
        '''
        Actions to take place in each time step. Namely,
        - Check for new infections
        - Apply testing protocols
        - Check if it is time to move an infected node to the next 
          stage of their infections.
        - Check for introductions of new infections


        Parameters
        ----------
        t : float
            Current time
        activity : bool
            If contacts are active.
        dayindex : int
            Index of networks in the dataset.
        it : int
            Index of current time.
        to_continue : bool
            If True, only new introductions and testing are evaluated.

        '''

        if to_continue and activity:
            for inode in self.Nodes.ind_nodes[
                              (self.Nodes.state != self.healthy_state) &
                              (self.Nodes.state != self.recovered_state)]:
                state = self.Nodes.state[inode]
                category = self.Nodes.biogroup[inode]
                if self.Nodes.Dict_inf[category][state]:
                    tq0, tqf = self.Nodes.time_isolated[inode]
                    if not (tq0 < t < tqf):
                        self.find_transmissions(inode, t, it, dayindex)

        self.Test.testing(self, self.Nodes, t, activity)
        if to_continue:
            for inode in self.Nodes.ind_nodes[(self.Nodes.time < self.tmax)]:
                rnext = self.Nodes.rnext[inode]
                if rnext < self.Nr:
                    if self.Nodes.time[inode] <= t:
                        self.change(inode, rnext, t)

        if self.weekden != 0:
            self.check_introductions(t)
