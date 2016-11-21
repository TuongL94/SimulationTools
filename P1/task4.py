# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 12:39:07 2016
@author: Jonathan
"""
#from scipy import *
import pylab as P
import numpy as N
from assimulo.problem import Explicit_Problem  #Imports the problem formulation from Assimulo
from assimulo.solvers import CVode #Imports the solver CVode from Assimulo
#import matplotlib.pyplot as plt??

"""
Variable parameters
"""
#Spring constant
k=100

#Initial parameters
angle = 30 #0 is completely vertical, 90 is completely horizontal
stretch = .02

#Simulation parameters
tf=5.0
rtol = 1e-2 #default 1e-06
hmax = 0.0 #default 0 (infinity)
maxord = 1 #default 5
h0 = 0.0 #default 0.0

"""
Set up model
"""
#Determine initial conditions
length=1
lengthx = length*N.sin(N.pi*angle/180)
lengthy = -length*N.cos(N.pi*angle/180)
stretchx = stretch*N.sin(N.pi*angle/180)
stretchy = -stretch*N.cos(N.pi*angle/180)
y0=N.array([lengthx+stretchx, lengthy+stretchy, 0.0, 0.0])
t0=0.0

#Define rhs
def rhs(t,y):
    lamb = lambda y1,y2: k*(N.sqrt(y1**2+y2**2)-1) / N.sqrt(y1**2+y2**2)
    ydot=N.zeros(y.shape)
    ydot[0] = y[2]
    ydot[1] = y[3]
    ydot[2] = -y[0]*lamb(y[0],y[1])
    ydot[3] = -y[1]*lamb(y[0],y[1]) - 1

    return ydot

#Create Assimulo problem model
mod_pen = Explicit_Problem(rhs,y0,t0)
mod_pen.name = 'Elastic Pendulum with CVode'


"""
Run simulation
"""
#Create solver object
sim_pen = CVode(mod_pen)

#Simulation parameters
#atol = 1.e-6*ones(shape(y0))
atol = rtol*N.array([1, 1, 1, 1]) #default 1e-06
sim_pen.atol = atol
sim_pen.rtol = rtol
sim_pen.maxh = hmax
sim_pen.maxord = maxord
sim_pen.inith = h0

#Simulate
t,y=sim_pen.simulate(tf)


"""
Plot results
"""
#Plot
#P.plot(t,y)
#P.legend(['$x$','$y$','$\dot{x}$','$\dot{y}$'])
#P.title('Elastic pendulum with k = ' + str(k) + ', stretch = ' + str(stretch))
#P.xlabel('time')
#P.ylabel('state')
#show()
P.plot(y[:,0],y[:,1])
P.plot(y[0,0],y[0,1],'ro')
P.title('Elastic pendulum with k = ' + str(k) + ', stretch = ' + str(stretch)
        + '\nCVode with rtol = ' + str(rtol)
#        + ', atol = ' + str(atol)
            + ', maxord = ' + str(maxord))
P.xlabel('x')
P.ylabel('y')
P.show()
#filename = 'pend' + str(k) + '.eps'
P.savefig('pendtest.ps')
