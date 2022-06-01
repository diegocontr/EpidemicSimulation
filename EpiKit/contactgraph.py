#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 00:50:16 2021

@author: diego
"""

import numpy as np
import pandas as pd
import networkx as nx


def get_contact_graph(file, contact_meta, dt_min, Gall0=False,
                      duration_day=31080,
                      tmin=31220):

    contact = pd.read_csv(file, sep='\t',
                          names=['t', 'i', 'j', 'Ci', 'Cj'])
    return get_contact_graph_c(contact, contact_meta, dt_min, Gall0,
                               duration_day,  tmin)


def get_contact_graph_c(contact, contact_meta, dt_min, Gall0=False,
                        duration_day=31080,
                        tmin=31220):


    dt_G = int(dt_min*60)
    tfinal = tmin + duration_day

    utimes = np.arange(0, duration_day, dt_G)
    G0 = nx.Graph()
    G0.add_nodes_from(np.arange(contact_meta.shape[0]))
    Gdict = contact_meta.copy().set_index('new index').to_dict()
    nx.set_node_attributes(G0, Gdict['class'], "class")
    #nx.set_node_attributes(G0, Gdict['id'], "id")

    Gt = {t: G0.copy() for t in utimes}
    Gall = Gall0 if Gall0 else G0.copy()

    dictIndex = {contact_meta.index[i]: i for i in range(contact_meta.shape[0])}
    for ind in contact.index:
        t, i0, j0, Ci, Cj = contact.loc[ind]
        #t,i0,j0 = contact.loc[ind][['t','i','j']]
        if t > tmin and t < tfinal:
            i, j = dictIndex[i0], dictIndex[j0]

            ti = (int((t - tmin)//dt_G))*dt_G

            if Gt[ti].has_edge(i, j):
                Gt[ti][i][j]["weight"] = Gt[ti][i][j]["weight"] + 1.
            else:
                Gt[ti].add_edge(i, j, weight=1)

            if Gall.has_edge(i, j):
                Gall[i][j]["weight"] = Gall[i][j]["weight"] + 1.0
            else:
                Gall.add_edge(i, j, weight=1.0)

    return Gall, Gt


class GraphData:
    def __init__(self, dic):
        self.Gt = dic['Gt']
        self.dt_sec = dic['dt_sec']
        self.meta = dic['meta']
        self.TeachersDict = dic['TeachersDict']

        # self.sec0 = 86400
        # self.secf = 30880
        self.sec0 = dic['sec0']
        self.secf = dic['secf']

        self.ndays = len(self.Gt)
        self.keys = ['Gt', 'dt_sec', 'meta', 'TeachersDict']

        self.net_cycle = dic['net_cycle']
        self.type = dic['type']

    def __getitem__(self, k):
        if k in self.keys:
            if k == 'Gt':
                return self.Gt
            elif k == 'dt_sec':
                return self.dt_sec
            elif k == 'meta':
                return self.meta
            elif k == 'TeachersDict':
                return self.TeachersDict
        else:
            raise KeyError('{} is not a key of GraphData'.format(k))
#    def get_G_for_day(self,d):


class GraphDataSchool(GraphData):
    def __init__(self, dic):
        super().__init__(dic)
        #self.net_cycle = round(len(self.Gt) * 7/4 )

    def activitytimes(self, t):
        day = np.array(t % 7, int)
        second = np.array((t % 1)*(3600*24), int)-self.sec0
        cond_day = (day == 0) + (day == 1) + (day == 3) + (day == 4)
        cond_sec = (second > 0)*(second < self.secf)
        second20 = (second//self.dt_sec)*self.dt_sec
        week = np.array(t//7, int)
        dayactivity = day.copy()
        dayactivity[dayactivity > 2] -= 1
        day2 = (4 * week + dayactivity) % self.ndays
        return cond_day*cond_sec, second20, day2

    def workingday(self, t):
        day = np.array(t % 7, int)
        return (day == 0) + (day == 1) + (day == 3) + (day == 4)


class GraphDataOffice(GraphData):
    def __init__(self, dic):
        super().__init__(dic)
        self.list_week_days = [0, 1, 2, 3, 4]

        cycle_day = np.arange(self.net_cycle) % 7
        cond_day = np.logical_or.reduce(
            [cycle_day == d for d in self.list_week_days])
        self.index_to_day = (np.cumsum(cond_day)-1) % self.ndays
        self.index_to_day[~cond_day] = self.index_to_day[~cond_day] * 0 + 10000
        
    def day_to_Gt_ind(self, d):
        return self.index_to_day[d % self.net_cycle]

    def activitytimes(self, t):
        dt = t[1] - t[0]
        day = np.array(t % 7, int)
        second = np.array(((t % 1)*(3600*24)).round(), int)-self.sec0
        cond_day = np.logical_or.reduce(
            [day == d for d in self.list_week_days])
        cond_sec = (second >= 0)*(second <= self.secf)
        second20 = np.array((second/self.dt_sec).round(), int)*self.dt_sec
        day2 = np.array(dt * (np.cumsum(cond_day) - 1) +
                        t[0], int) % self.ndays
        day2[~cond_day] = day2[~cond_day] * 0 + 10000
        return cond_day*cond_sec, second20, day2



class GraphDataAllDays(GraphData):
    def __init__(self, dic):
        super().__init__(dic)
        self.list_week_days = [0, 1, 2, 3, 4, 5, 6]

        cycle_day = np.arange(self.net_cycle) % 7
        cond_day = np.logical_or.reduce(
            [cycle_day == d for d in self.list_week_days])
        self.index_to_day = (np.cumsum(cond_day)-1) % self.ndays
        self.index_to_day[~cond_day] = self.index_to_day[~cond_day] * 0 + 10000

    def day_to_Gt_ind(self, d):
        return self.index_to_day[d % self.net_cycle]

    def activitytimes(self, t):
        day = np.array(t % 7, int)
        second = np.array((t % 1)*(3600*24), int)-self.sec0
        cond_day = (day != 10000.2)
        cond_sec = (second > 0)*(second < self.secf)
        second20 = (second//self.dt_sec)*self.dt_sec
        week = np.array(t//7, int)
        dayactivity = day.copy()
        # dayactivity[dayactivity>2] -= 1
        day2 = (7 * week + dayactivity) % self.ndays
        return cond_day*cond_sec, second20, day2

    def workingday(self, t):
        return t == t


class GraphDataDays(GraphData):
    def __init__(self, dic, day_list):
        super().__init__(dic)
        self.changedays(day_list)

    def changedays(self, day_list):
        self.day_list = day_list
        self.n_workdays = len(day_list)

    def activitytimes(self, t):
        day = np.array(t % 7, int)
        second = np.array((t % 1)*(3600*24), int)-self.sec0
        cond_day = 0*day
        for D in self.day_list:
            cond_day += (day == D)
        cond_sec = (second > 0)*(second < self.secf)
        second20 = (second//self.dt_sec)*self.dt_sec
        week = np.array(t//7, int)
        dayactivity = day.copy()
        #dayactivity[dayactivity>2] -= 1
        day2 = (self.n_workdays * week + dayactivity) % self.ndays
        return cond_day*cond_sec, second20, day2


# %%
