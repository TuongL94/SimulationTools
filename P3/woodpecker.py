#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 20:42:44 2017
@author: Jonathan
"""
from woodpecker_model import init_woodpecker
from woodpecker_model import res
from state_event_woodpecker import state_event
from handle_event_woodpecker import handle_event
from assimulo.problem import Implicit_Problem
from assimulo.solvers import IDA
import numpy as N
import matplotlib.pyplot as P

t0 = 0.0
y0,yd0,switches0 = init_woodpecker()

model = Implicit_Problem(res,y0,yd0,t0,sw0=switches0)
model.state_events = state_event
model.handle_event = handle_event
#model.name = 'Woodpecker model with state events'

sim = IDA(model)
sim.algvar = [1, 1, 1, 0, 0, 0, 0, 0]
sim.suppress_alg = True
#sim.atol = [1e-5,1e-5,1e-5,1e15,1e15,1e15,1e15,1e15]

tfinal = 0.5
t,y,yd = sim.simulate(tfinal)

P.plot(t,y[:,1])
sim.print_event_data()

