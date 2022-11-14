#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This part of the code is particularly released under the CRAPL 
license https://matt.might.net/articles/crapl/

The implementation of the different tests represents low re-usability code,
so this part was written with no major consideration on its maintainability
and re-usability.

The author of this code adheres to the values of open science, and releases 
this code in the interest of transparency. However, if my code is of interest 
to you, this part may need to be re-implemented. 

Simulations of epidemics without testing protocols or quarantines do not 
require this part of the code.

"""

import numpy as np
import numpy.random as nprm

### TODO : remove Q states or time: use one or the other
### TODO : remove repetition of code
### TODO : make a more coherent output
class Test_None:
    def __init__(self, Param):
        if 'Npop' in Param.keys():
            self.Npop = Param['Npop']
        else:
            raise(ValueError, 'Population "Npop" not informed in the dictionary of parameters' )
             
        self.Param = Param
        self.infectious_states = ('Ic', 'Isc', 'Ip')
        self.testable_states = ('Ic', 'Ip', 'Isc', 'Rp')
        self.p_detect = Param['p_detect']  # np.array([0.3, 0.5])
        self.type = Param['type']
        self.Delta_Q = Param['Delta_Q']
        if self.type == 'PCR':
            self.turn_around = 1
            self.sensitivity = {'Ic': 0.8, 'Ip': 0.7,
                                'Isc': 0.8, 'Rp': 0.0, 'R': 0.0}

        if self.type == 'antigen':
            self.turn_around = 30 / (60 * 24)
            self.sensitivity = {'Ic': 0.8, 'Ip': 0.7,
                                'Isc': 0.7, 'Rp': 0.0, 'R': 0.0}


        if self.type == 'auto':
            self.turn_around = 15 / (60 * 24)
            self.sensitivity = {'Ic': 0.8, 'Ip': 0.5,
                                'Isc': 0.7, 'Rp': 0.0, 'R': 0.0}

        if self.type == 'perfect':
            self.turn_around = 15 / (60 * 24)
            self.sensitivity = {'Ic': 1.0, 'Ip': 1.0,
                                'Isc': 1.0, 'Rp': 0.0, 'R': 0.0}

        self.sensitivity_Qteacher = {
            'Ic': 0.8, 'Ip': 0.7, 'Isc': 0.8, 'Rp': 0.8, 'R': 0.0}

        self.qtime_waiting = np.array([0, self.turn_around])
        self.qtime = np.array([0, self.Delta_Q + self.turn_around])
        self.qtime_reg = np.array(
            [self.turn_around, self.Delta_Q + self.turn_around])
        self.reset()

    def working_hours(t):
        return (0.375 < (t % 1) < 0.833)

    def testing(self, sim, nodes, t, activity):
        pass

    def reset(self, sim=None):
        self.number_of_tests = 0
        self.Ndetected = {0: 0, 1: 0, 2: 0, 3: 0}
        self.sim = sim

        self.quarantines_per_biogroup = {0: {'Qi': 0, 'Qpi': 0, 'Qw': 0},
                                         1: {'Qi': 0, 'Qpi': 0, 'Qw': 0},
                                         2: {'Qi': 0, 'Qpi': 0, 'Qw': 0},
                                         3: {'Qi': 0, 'Qpi': 0, 'Qw': 0}}


    def get_test_quarantine(self):
        d = [sum([(self.quarantines_per_biogroup[a][Q]) for Q in ['Qw', 'Qpi', 'Qi']])
             for a in [0, 1, 2, 3]]
        return int(self.number_of_tests), np.array(d)

    def testing_symptomatic(self, sim, Nodes, t, activity):
        positive = []
        for inode in Nodes.ind_nodes[Nodes.state == 'Ic']:
            if not Nodes.symp_tested[inode]:
                Nodes.symp_tested[inode] = 1
                adult = Nodes.biogroup[inode]
                if nprm.random() < self.p_detect[adult]:
                    Nodes.tested[inode] += 1
                    self.number_of_tests += 1
                    if nprm.random() < self.sensitivity['Ic']:
                        positive.append(inode)
                    else:
                        Nodes.qnodes[inode] = 'Ineg'

        return positive

    def testing_symptomatic_posneg(self, sim, Nodes, t, activity):
        positive = []
        negative = []
        for inode in Nodes.ind_nodes[Nodes.state == 'Ic']:
            if not Nodes.symp_tested[inode]:
                Nodes.symp_tested[inode] = 1
                adult = Nodes.biogroup[inode]
                if nprm.random() < self.p_detect[adult]:
                    Nodes.tested[inode] += 1
                    self.number_of_tests += 1
                    if nprm.random() < self.sensitivity['Ic']:
                        positive.append(inode)
                    else:
                        negative.append(inode)
        return positive, negative

    def testing_group(self, Nodes, to_test):
        positive = []
        for inode in to_test:
            state = Nodes.state[inode]
            Nodes.tested[inode] += 1
            self.number_of_tests += 1
            if state in self.testable_states:
                adult = Nodes.biogroup[inode]
                if nprm.random() < self.sensitivity[state]:
                    # nodes.tested[inode] += 1
                    #nodes.symp_tested[inode] = 1
                    positive.append(inode)
        return positive

    def send_to_quarantine(self, Nodes, inode, q_time, Q='Qi'):
        tfin = min([self.sim.tmax, q_time[1]])
        dqt = tfin - q_time[0]
        # print(Q,dqt)

        Nodes.qnodes[inode] = Q
        Nodes.qnext[inode] = q_time[0] + self.Delta_Q

        Nodes.time_isolated[inode] = q_time
        Nodes.days_quarantine[inode] += dqt
        Nodes.time_detected[inode] = q_time[0]
        if Nodes.state[inode] in self.infectious_states:
            Nodes.positive_detected[inode] = 1
        self.quarantines_per_biogroup[Nodes.biogroup[inode]][Q] += dqt
        if Q == 'Qi':
            self.Ndetected[Nodes.biogroup[inode]] += 1
        #self.Qtime_w[nodes.biogroup[inode]  ]  += self.qtime[1]

    def testing_before_F(self, Nodes, t):
        nodes = Nodes.ind_nodes[(Nodes.qnodes == 'Qi') *
                                (Nodes.time_isolated[:, 1] <= t)]
        positive = []
        for inode in nodes:
            adult = Nodes.biogroup[inode]
            state = Nodes.Nodes[inode]
            self.number_of_tests += 1
            if nprm.random() < self.sensitivity[state]:
                positive.append(inode)

        positive, negative = self.testing_symptomatic_posneg(
            sim, Nodes, t, activity)

        self.isolation_q_teacher(sim, Nodes, t, positive,
                                 qinterval=self.qtime, qstate='Qi')

    def set_qstate_F(self, Nodes, t):
       # self.testing_before_F( nodes, t)
        Nodes.qnodes[(Nodes.qnodes != 'F') *
                     (Nodes.time_isolated[:, 1] < t)] = 'F'

    def counting_qteacher(self):
        return [0, 0]


class Test_ST_office(Test_None):
    def __init__(self, Param):
        super().__init__(Param)
        self.reset()

    def reset(self, sim=None):
        super().reset(sim)
        self.time_test = 0.0  # midday
        self.periodicity_test = 0.25 / 24  # 15 m
        self.testing_teacher = []
        self.waiting_teacher = []

    def set_node_state(self, sim, Nodes, inode, state):
        Nodes.state[inode] = state
        rnext, time = sim.gR.set_and_rnext(inode)
        Nodes.time[inode] = time
        Nodes.rnext[inode] = rnext

    def teacher_replacement(self, sim, Nodes, t):
        pass

    def isolation_waiting_teacher(self, sim, Nodes, t):
        pass

    def isolation_q_teacher(self, sim, Nodes, t,
                            positive, qinterval, qstate='Qi'):
        q_time = qinterval + t
        for inode in positive:
            Nodes.Qcount.append((str(qstate), int(inode), float(q_time[0])))
            self.send_to_quarantine(Nodes, inode, q_time, qstate)

    def testing(self, sim, nodes, t, activity):
        if t > self.time_test:
            self.time_test = t + self.periodicity_test
            self.set_qstate_F(nodes, t)
            if activity:
                self.isolation_waiting_teacher(sim, nodes, t)
            if (not activity):  # : #and (0.375< (t%1) < 0.833):
                self.time_test = t + self.periodicity_test
                positive, negative = self.testing_symptomatic_posneg(
                    sim, nodes, t, activity)

                self.isolation_q_teacher(sim, nodes, t, positive,
                                         qinterval=self.qtime, qstate='Qi')
                self.isolation_q_teacher(
                    sim,
                    nodes,
                    t,
                    negative,
                    qinterval=self.qtime_waiting,
                    qstate='Qw')
                #self.isolation_symptomatics(sim, nodes, t, positive, negative)
                # self.isolation_positive_teacher(sim, nodes, t, positive)
                # self.isolation_negative(sim, nodes, t, negative)
                self.teacher_replacement(sim, nodes, t)


class Test_STRS_office(Test_ST_office):
    def __init__(self, Param):

        self.Delta_R = Param['Delta_R']
        self.adhesion = Param['adhesion']
        self.rdecay = Param['adhesion_decay']
        self.Delta_R2 = 4
        super().__init__(Param)

    # def reset(self,sim=None):
    #     super().reset(sim)
    #     self.reactive_tests = []
    #     self.random_adhesion = nprm.random(self.Npop)
    #     self.adhesion_t = self.adhesion

    #     self.number_of_tests = round(self.Npop * self.adhesion)
    #     self.ind_to_test = nprm.choice( np.arange(self.Npop), self.number_of_tests , replace=False )
    #     self.to_test = np.zeros(self.Npop,int)
    #     self.to_test[self.ind_to_test] +=1
    def N_to_test(self):
        return round(self.Nnovax * self.adhesion)

    def reset(self, sim=None):
        super().reset(sim)
        self.reactive_tests = []
        self.reactive_tests2 = []

        self.time_wteacher = 0
        self.random_adhesion = nprm.random(self.Npop)
        self.adhesion_t = self.adhesion
        if sim is not None:
            ages = sim.Nodes.biogroup
            self.nonvax = sim.Nodes.ind_nodes[ages == 0]
            self.Nnovax = self.nonvax.size
            self.Ntest = round(self.nonvax.size * self.adhesion)
            self.ind_to_test = nprm.choice(
                self.nonvax, self.Ntest, replace=False)
            self.to_test = np.zeros(self.Npop, int)
            self.to_test[self.ind_to_test] += 1
        else:
            self.Nnovax = self.Npop
            self.Ntest = self.N_to_test()  # round(self.Npop * self.adhesion)
            self.ind_to_test = nprm.choice(np.arange(self.Npop),
                                           self.Ntest, replace=False)
            self.to_test = np.zeros(self.Npop, int)
            self.to_test[self.ind_to_test] += 1

    def add_reactive_screaning(self, Nodes, t, positive):
        for inode in positive:
            c = Nodes.node_category[inode]
            if not(c in [c0 for t, c0 in self.reactive_tests2]):
                # treac =  t+self.Delta_R
                # while not( int(treac%7) in nodes.activitydays ):
                #     treac += 1
                self.reactive_tests.append((t + self.Delta_R, c))
                self.reactive_tests2.append((t + self.Delta_R2, c))

    def reactive_screaning(self, sim, Nodes, t):
        for reactive_tests in [self.reactive_tests,
                               self.reactive_tests2]:
            for (treac, nclass) in reactive_tests:
                if t > treac:
                    reactive_tests.remove((treac, nclass))
                    cnodes = Nodes.ind_nodes[(Nodes.node_category == nclass) *
                                             (Nodes.qnodes != 'Qi')]
                    r = self.to_test[cnodes]

                    to_test = cnodes[r == 1]
                    positive = self.testing_group(Nodes, to_test)
                    self.isolation_q_teacher(sim, Nodes, t, positive,
                                             qinterval=self.qtime, qstate='Qi')
                    for inode in positive:
                        Nodes.symp_tested[inode] = 1

    def testing(self, sim, nodes, t, activity):
        if t > self.time_test:
            self.set_qstate_F(nodes, t)
            self.time_test = t + self.periodicity_test

            if (not activity):  # and (0.375< (t%1) < 0.833):
                positive, negative = self.testing_symptomatic_posneg(
                    sim, nodes, t, activity)
                #self.isolation_symptomatics(sim, nodes, t, positive, negative)
                self.isolation_q_teacher(sim, nodes, t, positive,
                                         qinterval=self.qtime, qstate='Qi')
                self.isolation_q_teacher(
                    sim,
                    nodes,
                    t,
                    negative,
                    qinterval=self.qtime_waiting,
                    qstate='Qw')

                #self.teacher_replacement(sim, nodes, t)
                self.add_reactive_screaning(nodes, t, positive)
            if activity:
                self.reactive_screaning(sim, nodes, t)


class Test_Regular_office(Test_ST_office):
    def __init__(self, Param):
        if Param['W'] in [14.0, 7.0]:
            self.W = [Param['W']]
            self.nround = 1
            self.time_regtest0 = 1 + 8.5 / 24
        elif Param['W'] == 3.5:
            self.nround = 2
            self.W = [2, 5]
            self.time_regtest0 = 1 + 8.5 / 24
        elif Param['W'] == 1.75:
            self.nround = 4
            self.W = [1, 2, 1, 3]
            self.time_regtest0 = 8.5 / 24
        else:
            raise ValueError(
                'W={} is not a valid periodicity'.format(Param['W']))
        self.adhesion = Param['adhesion']
        self.rdecay = Param['adhesion_decay']
        super().__init__(Param)

    def N_to_test(self):
        return round(self.Nnovax * self.adhesion_t)

    def reset(self, sim=None):
        super().reset(sim)
        self.round = 0
        if sim is not None:
            ti = sim.OutResult['time0']
            self.round = 0
            self.time_regtest = -ti
            while True:
                self.time_regtest += self.W[self.round % self.nround]
                self.round = self.round + 1
                if self.time_regtest > 0:
                    break
        else:
            t0 = np.nan
            self.time_regtest = 'waiting'

        self.time_wteacher = 0
        self.random_adhesion = nprm.random(self.Npop)
        self.adhesion_t = self.adhesion
        if sim is not None:
            ages = sim.Nodes.biogroup
            #self.nonvax = sim.nodes.ind_nodes[np.logical_or(ages==0,ages==1 )]
            self.nonvax = sim.Nodes.ind_nodes[ages == 0]
            self.Nnovax = self.nonvax.size
            self.Ntest = round(self.nonvax.size * self.adhesion)
            self.ind_to_test = nprm.choice(
                self.nonvax, self.Ntest, replace=False)
            self.to_test = np.zeros(self.Npop, int)
            self.to_test[self.ind_to_test] += 1
        else:
            self.Nnovax = self.Npop
            self.Ntest = self.N_to_test()  # round(self.Npop * self.adhesion)
            self.ind_to_test = nprm.choice(np.arange(self.Npop),
                                           self.Ntest, replace=False)
            self.to_test = np.zeros(self.Npop, int)
            self.to_test[self.ind_to_test] += 1

    def testing(self, sim, nodes, t, activity):
        super().testing(sim, nodes, t, activity)
        if activity and t >= self.time_regtest:
            # self.time_regtest += self.W
            self.time_regtest += self.W[self.round % self.nround]
            self.round = self.round + 1
            #print(t%7,self.time_regtest, self.round,self.round % self.nround)
            #self.set_qstate_F(nodes, t)
            to_test = nodes.ind_nodes[(self.to_test == 1) *
                                      (nodes.qnodes != 'Qi') *
                                      (nodes.qnodes != 'Qpi')]

            positive = self.testing_group(nodes, to_test)
            self.isolation_q_teacher(sim, nodes, t, positive,
                                     qinterval=self.qtime_reg, qstate='Qi')
            # if self.rdecay != 1.0:
            #     self.adhesion_t *= self.rdecay
            #     self.N_to_test()

            #     self.number_of_tests = self.N_to_test() #round(self.Nnovax * self.adhesion_t)
            #     self.ind_to_test = nprm.choice( self.ind_to_test, self.number_of_tests , replace=False )
            #     self.to_test = np.zeros(self.Npop,int)
            #     self.to_test[self.ind_to_test] +=1



class Test_TW_office(Test_ST_office):
    def __init__(self, Param, n=5):
        self.split = n
        self.saveTW = False
        super().__init__(Param)

    def reset(self, sim=None):
        super().reset(sim)
        groups = np.arange(self.Npop)
        nprm.shuffle(groups)
        d = self.Npop // self.split
        #self.groups = groups.reshape((self.split, self.Npop//self.split ))
        self.groups = []
        for i in range(self.split):
            self.groups.append(groups[d * i: d * (i + 1)])
        self.q_telework = np.array([0, 1])
        self.time_setq = 0 / 24
        self.ind_group = 0

        if sim is not None:
            ti = sim.OutResult['time0']
            self.time0 = ti

        else:
            print('Warning: no sim object found in Test_TW_office')
            self.time0 = 'waiting'

    def send_to_telework(self, Nodes, t):
        q_telework = self.q_telework + t
        self.ind_group = (self.ind_group + 1) % self.split
        group = self.groups[self.ind_group]
        for inode in group:
            tq0, tqf = Nodes.time_isolated[inode]
            if not(tq0 < t < tqf):
                Nodes.time_isolated[inode] = q_telework
                if self.saveTW:
                    Nodes.qnodes[inode] = 'TW'

    def testing(self, sim, nodes, t, activity):
        super().testing(sim, nodes, t, activity)
        if t >= self.time_setq:
            if int((t + self.time0) % 7) in [0, 1, 2, 3, 4]:
                # if activity:
                self.time_setq = t + 1
                self.send_to_telework(nodes, t)


class Test_TW_office_2tw(Test_ST_office):
    def __init__(self, Param, n=5):
        self.split = n
        self.saveTW = False
        super().__init__(Param)

    def reset(self, sim=None):
        super().reset(sim)
        groups = np.arange(self.Npop)
        nprm.shuffle(groups)
        d = self.Npop // self.split
        #self.groups = groups.reshape((self.split, self.Npop//self.split ))
        self.groups = []
        for i in range(self.split):
            self.groups.append(groups[d * i: d * (i + 1)])
        self.q_telework = np.array([0, 1])
        self.time_setq = 0 / 24
        self.ind_group = 0

        if sim is not None:
            ti = sim.OutResult['time0']
            self.time0 = ti

        else:
            print('Warning: no sim object found in Test_TW_office')
            self.time0 = 'waiting'

    def send_to_telework(self, Nodes, t):
        q_telework = self.q_telework + t
        self.ind_group = (self.ind_group + 1) % self.split
        self.ind_group2 = (self.ind_group + 4) % self.split
        for group in [self.groups[self.ind_group],
                      self.groups[self.ind_group2]]:
            for inode in group:
                tq0, tqf = Nodes.time_isolated[inode]
                if not(tq0 < t < tqf):
                    Nodes.time_isolated[inode] = q_telework
                    if self.saveTW:
                        Nodes.qnodes[inode] = 'TW'

    def testing(self, sim, nodes, t, activity):
        super().testing(sim, nodes, t, activity)
        if t >= self.time_setq:
            if int((t + self.time0) % 7) in [0, 1, 2, 3, 4]:
                # if activity:
                self.time_setq = t + 1
                self.send_to_telework(nodes, t)


class Test_RegTW_office(Test_Regular_office):
    def __init__(self, Param, n=5):
        self.split = n
        self.saveTW = False
        super().__init__(Param)

    def reset(self, sim=None):
        super().reset(sim)
        groups = np.arange(self.Npop)
        nprm.shuffle(groups)
        d = self.Npop // self.split
        self.groups = []
        for i in range(self.split):
            self.groups.append(groups[d * i: d * (i + 1)])
        self.q_telework = np.array([0, 1])
        self.time_setq = 8 / 24
        self.ind_group = 0

        if sim is not None:
            ti = sim.OutResult['time0']
            self.time0 = ti
        else:
            self.time0 = 'waiting'

    def send_to_telework(self, Nodes, t):
        q_telework = self.q_telework + t
        self.ind_group = (self.ind_group + 1) % self.split
        group = self.groups[self.ind_group]
        for inode in group:
            tq0, tqf = Nodes.time_isolated[inode]
            if not(tq0 < t < tqf):
                Nodes.time_isolated[inode] = q_telework
                # if self.saveTW:
                # nodes.qnodes[inode] = 'TW'

    def testing(self, sim, nodes, t, activity):
        super().testing(sim, nodes, t, activity)
        if t >= self.time_setq:
            if int((t + self.time0) % 7) in [0, 1, 2, 3, 4]:
                # if activity:
                self.time_setq = t + 1
                self.send_to_telework(nodes, t)

##########################


class Test_ST(Test_ST_office):
    def __init__(self, Param):
        super().__init__(Param)
        self.reset()

    def teacher_replacement(self, sim, Nodes, t):
        to_remove = []
        for test_now in range(len(self.testing_teacher)):

            iprof, tprof, nq = self.testing_teacher[test_now]
            if t > tprof:
                #self.testing_teacher.remove(  test_now  )
                to_remove.append(self.testing_teacher[test_now])
                if nq == 1:
                    self.number_of_tests += 1
                    Nodes.tested[iprof] += 1
                    if nprm.random() < self.sensitivity_Qteacher['Rp']:
                        self.testing_teacher.append(
                            (iprof, self.qtime[1] + t, 2))
                        self.quarantines_per_biogroup[1]['Qi'] += 1
                    else:
                        self.set_node_state(sim, Nodes, iprof, 'Rp')
                else:
                    self.set_node_state(sim, Nodes, iprof, 'Rp')
        if to_remove:
            #print( len(self.testing_teacher) ,  len(to_remove))
            for elem in to_remove:
                self.testing_teacher.remove(elem)
            #print( len(self.testing_teacher) )

    def isolation_q_teacher(self, sim, Nodes, t,
                            positive, qinterval, qstate='Qi'):
        q_time = qinterval + t
        for inode in positive:
            Nodes.Qcount.append((str(qstate), int(inode), float(q_time[0])))
            if not Nodes.biogroup[inode] in Nodes.categories_of_teachers:
                self.send_to_quarantine(Nodes, inode, q_time, qstate)
            else:
                self.waiting_teacher.append((inode, *q_time))

    def isolation_waiting_teacher(self, sim, Nodes, t):
        for wt in self.waiting_teacher:
            inode, t0, tf = wt
            if t > t0:
                self.testing_teacher.append((inode, tf, 1))
                self.waiting_teacher.remove(wt)
                self.set_node_state(sim, Nodes, inode, 'S')

    def counting_qteacher(self):
        return [len(self.waiting_teacher), len(self.testing_teacher)]


class Test_STPI1(Test_ST):
    def __init__(self, Param):
        #self.cases_reac = n
        super().__init__(Param)

    def reset(self, sim=None):
        super().reset(sim)
        # end of quarantine for class
        self.qclass = {c: 0 for c in self.Lclass}

    def preventive_isolation(self, Nodes, t, positive):
        q_time = self.qtime + t
        for inode in positive:
            c = Nodes.node_category[inode]
            if (self.qclass[c] < t):
                self.qclass[c] = q_time[1]
                for jnode in Nodes.ind_nodes[(Nodes.node_category == c)]:
                    Q = 'Qi' if Nodes.qnodes[jnode] == 'Qi' else 'Qpi'
                    self.send_to_quarantine(Nodes, jnode, q_time, Q)

    def testing(self, sim, nodes, t, activity):
        if t > self.time_test:
            if (not activity):  # and (0.375< (t%1) < 0.833):
                self.set_qstate_F(nodes, t)
                self.time_test = t + self.periodicity_test
                positive, negative = self.testing_symptomatic_posneg(
                    sim, nodes, t, activity)
                #self.isolation_symptomatics(sim, nodes, t, positive, negative)

                self.isolation_q_teacher(sim, nodes, t, positive,
                                         qinterval=self.qtime, qstate='Qi')
                self.isolation_q_teacher(
                    sim,
                    nodes,
                    t,
                    negative,
                    qinterval=self.qtime_waiting,
                    qstate='Qw')

                self.preventive_isolation(nodes, t, positive)


class Test_STPIn(Test_ST):
    def __init__(self, Param, n=1):
        self.cases_reac = n
        super().__init__(Param)

    def reset(self, sim=None):
        super().reset(sim)
        # end of quarantine for class
        self.qclass = {c: 0 for c in self.Lclass}

    def preventive_isolation(self, Nodes, t, positive):
        q_time = self.qtime + t
        for inode in positive:
            c = Nodes.node_category[inode]
            infected_same_class = Nodes.ind_nodes[(Nodes.node_category == c) *
                                                  # (nodes.time_isolated[:,1]>  t)   ].size
                                                  (Nodes.qnodes == 'Qi')].size

            if (infected_same_class >= self.cases_reac) and (
                    self.qclass[c] < t):
                self.qclass[c] = q_time[1]
                for jnode in Nodes.ind_nodes[(Nodes.node_category == c)]:
                    Q = 'Qi' if Nodes.qnodes[jnode] == 'Qi' else 'Qpi'
                    self.send_to_quarantine(Nodes, jnode, q_time, Q)

    def testing(self, sim, nodes, t, activity):
        if t > self.time_test:
            if (not activity):  # and (0.375< (t%1) < 0.833):
                self.set_qstate_F(nodes, t)
                self.time_test = t + self.periodicity_test
                positive, negative = self.testing_symptomatic_posneg(
                    sim, nodes, t, activity)
                #self.isolation_symptomatics(sim, nodes, t, positive, negative)

                self.isolation_q_teacher(sim, nodes, t, positive,
                                         qinterval=self.qtime, qstate='Qi')
                self.isolation_q_teacher(
                    sim,
                    nodes,
                    t,
                    negative,
                    qinterval=self.qtime_waiting,
                    qstate='Qw')

                self.preventive_isolation(nodes, t, positive)


class Test_STRS(Test_ST):
    def __init__(self, Param):

        self.Delta_R = Param['Delta_R']
        self.adhesion = Param['adhesion']
        self.rdecay = Param['adhesion_decay']
        self.Delta_R2 = 4
        super().__init__(Param)

    def N_to_test(self):
        return round(self.Nnovax * self.adhesion)

    def reset(self, sim=None):
        super().reset(sim)
        self.reactive_tests = []
        self.reactive_tests2 = []

        self.time_wteacher = 0
        self.random_adhesion = nprm.random(self.Npop)
        self.adhesion_t = self.adhesion
        if sim is not None:
            ages = sim.Nodes.biogroup
            self.nonvax = sim.Nodes.ind_nodes[np.logical_or(
                ages == 0, ages == 1)]
            self.Nnovax = self.nonvax.size
            self.Ntest = round(self.nonvax.size * self.adhesion)
            self.ind_to_test = nprm.choice(
                self.nonvax, self.Ntest, replace=False)
            self.to_test = np.zeros(self.Npop, int)
            self.to_test[self.ind_to_test] += 1
        else:
            self.Nnovax = self.Npop
            self.Ntest = self.N_to_test()  # round(self.Npop * self.adhesion)
            self.ind_to_test = nprm.choice(np.arange(self.Npop),
                                           self.Ntest, replace=False)
            self.to_test = np.zeros(self.Npop, int)
            self.to_test[self.ind_to_test] += 1

    def add_reactive_screaning(self, Nodes, t, positive):
        for inode in positive:
            c = Nodes.node_category[inode]
            if not(c in [c0 for t, c0 in self.reactive_tests]):
                # treac =  t+self.Delta_R
                # while not( int(treac%7) in nodes.activitydays ):
                #     treac += 1
                self.reactive_tests.append((t + self.Delta_R, c))
                self.reactive_tests.append((t + self.Delta_R2, c))

    def reactive_screaning(self, sim, Nodes, t):

        for (treac, nclass) in self.reactive_tests:
            if t > treac:
                self.reactive_tests.remove((treac, nclass))
                cnodes = Nodes.ind_nodes[(Nodes.node_category == nclass) *
                                         (Nodes.qnodes != 'Qi')]
                r = self.to_test[cnodes]

                to_test = cnodes[r == 1]
                positive = self.testing_group(Nodes, to_test)
                self.isolation_q_teacher(sim, Nodes, t, positive,
                                         qinterval=self.qtime, qstate='Qi')
                for inode in positive:
                    Nodes.symp_tested[inode] = 1

    def testing(self, sim, nodes, t, activity):
        if t > self.time_test:
            self.set_qstate_F(nodes, t)
            self.time_test = t + self.periodicity_test

            if (not activity):  # and (0.375< (t%1) < 0.833):
                positive, negative = self.testing_symptomatic_posneg(
                    sim, nodes, t, activity)
                #self.isolation_symptomatics(sim, nodes, t, positive, negative)
                self.isolation_q_teacher(sim, nodes, t, positive,
                                         qinterval=self.qtime, qstate='Qi')
                self.isolation_q_teacher(
                    sim,
                    nodes,
                    t,
                    negative,
                    qinterval=self.qtime_waiting,
                    qstate='Qw')

                self.teacher_replacement(sim, nodes, t)
                self.add_reactive_screaning(nodes, t, positive)
            # if activity:
            self.reactive_screaning(sim, nodes, t)


class Test_Regular(Test_ST):
    def __init__(self, Param):
        if Param['W'] in [14.0, 7.0]:
            self.W = [Param['W']]
            self.nround = 1
            self.time_regtest0 = 1 + 8.5 / 24
        elif Param['W'] == 3.5:
            self.nround = 2
            self.W = [2, 5]
            self.time_regtest0 = 1 + 8.5 / 24
        elif Param['W'] == 1.75:
            self.nround = 4
            self.W = [1, 2, 1, 3]
            self.time_regtest0 = 8.5 / 24
        else:
            raise ValueError(
                'W={} is not a valid periodicity'.format(Param['W']))
        self.adhesion = Param['adhesion']
        self.rdecay = Param['adhesion_decay']
        super().__init__(Param)

    def N_to_test(self):
        return round(self.Nnovax * self.adhesion_t)

    def reset(self, sim=None):
        super().reset(sim)
        self.round = 0
        if sim is not None:
            ti = sim.OutResult['time0']
            self.round = 0
            self.time_regtest = -ti
            while True:
                self.time_regtest += self.W[self.round % self.nround]
                self.round = self.round + 1
                if self.time_regtest > 0:
                    break
        else:
            t0 = np.nan
            self.time_regtest = 'waiting'

        self.time_wteacher = 0
        self.random_adhesion = nprm.random(self.Npop)
        self.adhesion_t = self.adhesion
        if sim is not None:
            ages = sim.Nodes.biogroup
            self.nonvax = sim.Nodes.ind_nodes[np.logical_or(
                ages == 0, ages == 1)]
            self.Nnovax = self.nonvax.size
            self.Ntest = round(self.nonvax.size * self.adhesion)
            self.ind_to_test = nprm.choice(
                self.nonvax, self.Ntest, replace=False)
            self.to_test = np.zeros(self.Npop, int)
            self.to_test[self.ind_to_test] += 1
        else:
            self.Nnovax = self.Npop
            self.Ntest = self.N_to_test()  # round(self.Npop * self.adhesion)
            self.ind_to_test = nprm.choice(np.arange(self.Npop),
                                           self.Ntest, replace=False)
            self.to_test = np.zeros(self.Npop, int)
            self.to_test[self.ind_to_test] += 1

    def testing(self, sim, nodes, t, activity):
        super().testing(sim, nodes, t, activity)
        if activity and t >= self.time_regtest:
            # self.time_regtest += self.W
            self.time_regtest += self.W[self.round % self.nround]
            self.round = self.round + 1
            #print(t%7,self.time_regtest, self.round,self.round % self.nround)
            #self.set_qstate_F(nodes, t)
            to_test = nodes.ind_nodes[(self.to_test == 1) *
                                      (nodes.qnodes != 'Qi') *
                                      (nodes.qnodes != 'Qpi')]

            positive = self.testing_group(nodes, to_test)
            self.isolation_q_teacher(sim, nodes, t, positive,
                                     qinterval=self.qtime_reg, qstate='Qi')
            # if self.rdecay != 1.0:
            #     self.adhesion_t *= self.rdecay
            #     self.N_to_test()

            #     self.number_of_tests = self.N_to_test() #round(self.Nnovax * self.adhesion_t)
            #     self.ind_to_test = nprm.choice( self.ind_to_test, self.number_of_tests , replace=False )
            #     self.to_test = np.zeros(self.Npop,int)
            #     self.to_test[self.ind_to_test] +=1
