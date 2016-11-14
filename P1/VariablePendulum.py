# -*- coding: utf-8 -*-
"""
Created on Sat Nov 12 11:37:28 2016

@author: Oskar
"""

# Till projekt 1, uppgift 3.
# Använder nu bara cvode och inga andra lösare.

import numpy as np
from assimulo.solvers import CVode
from assimulo.problem import Explicit_Problem
import matplotlib.pyplot as plt

global k
global stretch


# Input längs ner.



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
        
        tfinal = 5
        
        r0 = 1.0 + stretch
        x0 = r0*np.sin(theta0)
        y0 = -r0*np.cos(theta0)
        t0 = 0.0
        
        y_init = np.array([x0,y0,0,0])
        
        ElasticSpring = Explicit_Problem(rhs,y_init,t0)
        sim = CVode(ElasticSpring)
        sim.report_continuously=False
        
        npoints = 100*tfinal
        
        t,y = sim.simulate(tfinal,npoints)
        
        xpos, ypos = y[:,0], y[:,1]
        #plt.plot(t,y[:,0:2])
        plt.plot(xpos,ypos)
        title = 'k = {0}, stretch = {1}'.format(k,stretch)
        plt.xlabel('y_1')
        plt.ylabel('y_2')
        plt.title(title)
        plt.axis('equal')
        plt.show()
        
        
def multiSimulate(k_array,stretch_array):
    for next_k in k_array:
        global k
        k = next_k
        for next_stretch in stretch_array:
            global stretch
            stretch = next_stretch
            doSimulate()

if __name__ == "__main__":
    # Input
    k_array = [1e3,1e5] # Testar olika fjäderkonstanter
    stretch_array = [0,.1]  # Testar olika förspänningar; 0 = inte utsträckt
    multiSimulate(k_array,stretch_array)
