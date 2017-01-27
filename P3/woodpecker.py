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
import matplotlib.pyplot as P
import numpy as N

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

#sim.print_event_data()

"""
    Plotting

    y[:,0] - plots z (height)
    y[:,1] - plots sleeve angle
    y[:,2] - plots bird angle
"""

# Plot z (height)
P.plot(t,y[:,0])
P.title('Woodpecker')
P.ylabel('Height z (m)')
P.xlabel('Time (s)')
P.savefig('assim_woodpecker_z.eps')
P.show()

# Plot phi_S
P.plot(t,y[:,1])
P.title('Woodpecker')
P.ylabel('Sleeve angle $\phi_S$ (rad)')
P.xlabel('Time (s)')
P.savefig('assim_woodpecker_phiS.eps')
P.show()

# Plot phi_B
P.plot(t,y[:,2])
P.title('Woodpecker')
P.ylabel('Sleeve angle $\phi_B$ (rad)')
P.xlabel('Time (s)')
P.savefig('assim_woodpecker_phiB.eps')
P.show()

# Plot all 3 variables in one figure using subplots
P.figure()
P.subplot(311)
P.title('Woodpecker')
P.plot(t,y[:,0])
ymin,ymax = P.ylim()
P.yticks(N.arange(ymin,ymax+(ymax-ymin)/5,(ymax-ymin)/5))
P.ylabel('z (m)')
P.subplot(312)
P.plot(t,y[:,1])
ymin,ymax = P.ylim()
P.yticks(N.arange(ymin,ymax+(ymax-ymin)/5,(ymax-ymin)/5))
P.ylabel('$\phi_S$ (rad)')
P.subplot(313)
P.plot(t,y[:,2])
P.ylabel('$\phi_B$ (rad)')
P.xlabel('Time (s)')
P.tight_layout()
P.savefig('assim_composite.eps')

# Plot single period in one figure
#start = 500
#stop = 1000
#subtime = t[start:stop]
#subint = N.zeros(N.shape(y[start:stop,0:3]))
#subint[:,0] = (y[start:stop,0]-N.min(y[start:stop,0])) / (N.max(y[start:stop,0])-N.min(y[start:stop,0]))
#subint[:,1] = (y[start:stop,1]-N.min(y[start:stop,1])) / (N.max(y[start:stop,1])-N.min(y[start:stop,1]))
#subint[:,2] = (y[start:stop,2]-N.min(y[start:stop,2])) / (N.max(y[start:stop,2])-N.min(y[start:stop,2]))
#P.figure()
#P.plot(subtime,subint)
#P.ylim( (-0.1,1.1) )
