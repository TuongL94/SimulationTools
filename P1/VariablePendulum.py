# -*- coding: utf-8 -*-
"""
Created on Sat Nov 12 11:37:28 2016

@author: Oskar
"""

# Till projekt 1, uppgift 3.

import numpy as np
from assimulo.solvers import CVode
from assimulo.problem import Explicit_Problem
from assimulo.exception import Explicit_ODE_Exception
import matplotlib.pyplot as plt
from bdf2 import BDF_2
from ee import EE




# Input:
solver_array=["BDF_2","EE"]# Testar olika lösare, input i textformat
k_array = 10.0**np.arange(2,6) # Testar olika fjäderkonstanter
stretch_array = [.02]  # Testar olika förspänningar; 0 = inte utsträckt


global k
global stretch
global solver



def rhs(t,y):
    def lamb(x1,x2): return k*(1 - 1/np.sqrt(x1**2 + x2**2))
    react = lamb(y[0],y[1])
    yd=np.array(range(4))*0.0
    
    yd[0] = y[2]
    yd[1] = y[3]
    yd[2] = -y[0]*react
    yd[3] = -y[1]*react-1

    #return np.array(range(4))
    return yd
    #return yd
    

def doSimulate():

        theta0 = 1.0 # radianer, startvinkel från lodrätt
        
        tfinal = 3
        
        title = 'k = {0}, stretch = {1}, {2}'.format(k,stretch,solver)
        
        r0 = 1.0 + stretch
        x0 = r0*np.sin(theta0)
        y0 = -r0*np.cos(theta0)
        t0 = 0.0
        
        y_init = np.array([x0,y0,0,0])
        
        ElasticSpring = Explicit_Problem(rhs,y_init,t0)
        if solver.lower() == "cvode":
            sim = CVode(ElasticSpring)
        elif solver.lower() =="bdf_2":
            sim = BDF_2(ElasticSpring)
        elif solver.lower() =="ee":
            sim = EE(ElasticSpring)
        else:
            sim == None
            raise ValueError('Expected "CVode", "EE" or "BDF_2"')
        sim.report_continuously=False
        
        npoints = 100*tfinal
        
        #t,y = sim.simulate(tfinal,npoints)
        try:
            t,y = sim.simulate(tfinal)
            xpos, ypos = y[:,0], y[:,1]
            #plt.plot(t,y[:,0:2])
            plt.plot(xpos,ypos)
            plt.plot(xpos[0],ypos[0],'or')
            plt.xlabel('y_1')
            plt.ylabel('y_2')
            plt.title(title)
            plt.axis('equal')
            plt.show()
        except Explicit_ODE_Exception:
            print("Solver did not converge within the specified no. of steps for {0}.".format(title))

        
#def multiSimulate(k_array,stretch_array,solver_array):
for next_solver in solver_array:
    solver = next_solver
    for next_k in k_array:
        k = next_k
        for next_stretch in stretch_array:
            stretch = next_stretch
            doSimulate()

#if __name__ == "__main__":
#    multiSimulate(k_array,stretch_array,solver_array)
