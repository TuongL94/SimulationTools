# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 13:59:47 2016

@author: Anders
"""
import numpy as N
import pylab as P
from assimulo.problem import Explicit_Problem
from assimulo.solvers import CVode
import matplotlib.pyplot as plt


def rhs(t,y):
    '''
    The function rhs() defines the right hand side of the differential equation
    on the form y'=f(y,t).
    :input: t is the current time (scalar), y is the current function values (vector)
    :output: returns the left hand side of the differential equation
    '''
    k=10 #k is the spring constant
    yd=N.zeros((len(y),1)) #Allocate space
    F=force(k,y[0],y[1]) # The force from the spring. Defined so that the force is positive if it is towards the center of rotation
    
    #Equation system
    yd[0]=y[2]
    yd[1]=y[3]
    yd[2]=-y[0]*F
    yd[3]=-y[1]*F-1

    return yd
    
    
def force(k,y1,y2):
    '''
    The function force() calculates the force acting on the pendulum due to the spring
    with spring constant k.
    :input: k is the spring constant, y1 is the horizontal coordinate and y2 the vertical
    of the pendulum.
    :output: returns the magnitude of the force, defined as positive towards the center of rotation
    '''
    return k*((sqrt(y1**2+y2**2)-1)/(sqrt(y1**2+y2**2)))

#Initial values    
y0=N.array([0.0,1.2,0.0,0.0])
t0=0.0

#Creating the model (problem)
model = Explicit_Problem(rhs, y0, t0)
model.name = 'Spring Pendulum ODE'

#Creating a simulation object
sim = CVode(model)

#Setting end of simulation time
tfinal = 10.0

#Simulates and plots
t, y = sim.simulate(tfinal)
sim.plot()

#Plots the horizontal and vertical positions and velocities separately
plt.figure()
plt.subplot(411)
plt.plot(t,y[:,0])
plt.title('y1')
plt.subplot(412)
plt.plot(t,y[:,1])
plt.title('y2')
plt.subplot(413)
plt.plot(t,y[:,2])
plt.title('der(y1)')
plt.subplot(414)
plt.plot(t,y[:,3])
plt.title('der(y2)')




