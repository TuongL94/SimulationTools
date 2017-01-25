# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 20:42:44 2017
@author: Jonathan Astermark, Anders Hansson, Tuong Lam and Oskar Smedman
********************************
    Main program - Only need to run this!
********************************
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
tfinal = 2.0
y0,yd0,switches0 = init_woodpecker()

model = Implicit_Problem(res,y0,yd0,t0,sw0=switches0)
model.state_events = state_event
model.handle_event = handle_event

sim = IDA(model)
sim.algvar = [1, 1, 1, 0, 0, 0, 0, 0]
sim.suppress_alg = True
sim.atol = [1e-5,1e-5,1e-5,1e15,1e15,1e15,1e15,1e15]
t,y,yd = sim.simulate(tfinal)

"""
    Plotting

    y[:,0] - plots z (height)
    y[:,1] - plots sleeve angle
    y[:,2] - plots bird angle
"""
P.plot(t,y[:,0])
P.show()
#sim.print_event_data()
