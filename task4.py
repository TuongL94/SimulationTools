# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 12:39:07 2016
@author: Jonathan
"""
#from scipy import *
from pylab import * #apparently a bad idea?
#import numpy as N
from numpy import * #apparently a bad idea?
from assimulo.problem import Explicit_Problem  #Imports the problem formulation from Assimulo
from assimulo.solvers import CVode #Imports the solver CVode from Assimulo
#import matplotlib.pyplot as plt??

"""
Variable parameters
"""
#Spring constant
k=10

#Initial parameters
angle = 30 #0 is completely vertical, 90 is completely horizontal
stretch = 0.2

#Simulation length
tf=10.0


"""
Set up model
"""
#Determine initial conditions
length=1
lengthx = length*sin(pi*angle/180)
lengthy = length*cos(pi*angle/180)
stretchx = stretch*sin(pi*angle/180)
stretchy = stretch*cos(pi*angle/180)
y0=array([lengthx+stretchx, lengthy+stretchy, 0.0, 0.0])
t0=0.0

#Define rhs
def rhs(t,y):
    lamb = lambda y1,y2: k*(sqrt(y1**2+y2**2)-1) / sqrt(y1**2+y2**2)
    ydot=zeros(y.shape)
    ydot[0] = y[2]
    ydot[1] = y[3]
    ydot[2] = -y[0]*lamb(y[0],y[1])
    ydot[3] = -y[1]*lamb(y[0],y[1]) - 1

    return ydot

#Create Assimulo problem model
mod_pen = Explicit_Problem(rhs,y0,t0)
mod_pen.name = 'Elastic Pendelum with CVode'


"""
Run simulation
"""
#Create solver object
sim_pen = CVode(mod_pen)

#Simulation parameters
sim_pen.atol = 1.e-6 #default 1e-06
sim_pen.rtol = 1.e-6 #default 1e-06
sim_pen.maxord = 5 #default 5

#Simulate
t,y=sim_pen.simulate(tf)


"""
Plot results
"""
#Plot
plot(t,y)
legend(['x','y','ghi','jkl'])
title('plot')
xlabel('time')
ylabel('state')
#show()