# -*- coding: utf-8 -*-
"""
Created on Sun Nov 13 13:55:27 2016

@author: Anders
"""
from assimulo.explicit_ode import Explicit_ODE
from assimulo.ode import *
import numpy as np
import matplotlib.pyplot as mpl
import scipy.linalg as SL
import scipy
# from assimulo.solvers import CVode
import pdb


class BDF_3WithNewt(Explicit_ODE):
    """
    BDF-3 solver with Newton's method as corrector
    """
    tol=1.e-8     
    maxit=100     
    maxsteps=5000
    
    def __init__(self, problem):
        Explicit_ODE.__init__(self, problem) #Calls the base class
        
        #Solver options
        self.options["h"] = 0.001
        
        #Statistics
        self.statistics["nsteps"] = 0
        self.statistics["nfcns"] = 0
    
    def _set_h(self,h):
            self.options["h"] = float(h)

    def _get_h(self):
        return self.options["h"]
        
    h=property(_get_h,_set_h)
        
    def integrate(self, t, y, tf, opts):
        """
        _integrates (t,y) values until t > tf
        """
        h = self.options["h"]
        h = min(h, abs(tf-t))
        
        #Lists for storing the result
        tres = []
        yres = []

        
        for i in range(self.maxsteps):
            #pdb.set_trace()
            if t >= tf:
                break
            self.statistics["nsteps"] += 1
            
            if i==0:  # Two first steps using explicit euler
                t_np1,y_np1 = self.step_EE(t,y, h)
                t_nm1=0 #Temporarily
                y_nm1=0 #Temporarily
            elif i==1:
                t_np1,y_np1 = self.step_EE(t,y, h)
                
            else:   
                t_np1, y_np1 = self.step_BDF3WithNewt([t,t_nm1,t_nm2], [y,y_nm1,y_nm2], h)
            t,t_nm1,t_nm2=t_np1,t,t_nm1
            y,y_nm1,y_nm2=y_np1,y,y_nm1
            #pdb.set_trace()
            
            tres.append(t)
            yres.append(y.copy())
        
            h=min(self.h,np.abs(tf-t))
        else:
            raise Explicit_ODE_Exception('Final time not reached within maximum number of steps')
        
        return ID_PY_OK, tres, yres
    
    def step_EE(self, t, y, h):
        """
        This calculates the next step in the integration with explicit Euler.
        """
        self.statistics["nfcns"] += 1
        
        f = self.problem.rhs
        return t + h, y + h*f(t, y) 
        
    def step_BDF3WithNewt(self,T,Y, h):
        """
        BDF-3 with Newton Iteration and Zero order predictor
        
        alpha_0*y_np1+alpha_1*y_n+alpha_2*y_nm1+alpha_3*y_nm2=h f(t_np1,y_np1)
        alpha=[11/6,-3,3/2,-1/3]
        """
        alpha=[11./6,-3.,3./2,-1./3]
        f=self.problem.rhs
        
        t_n,t_nm1,t_nm2=T
        y_n,y_nm1,y_nm2=Y
        # predictor
        t_np1=t_n+h
        y_np1_i=y_n   # zero order predictor

        #Defines a function for Newton's method
        def F(y):
            return (-(alpha[1]*y_n+alpha[2]*y_nm1+alpha[3]*y_nm2)+h*f(t_np1,y))/alpha[0]-y
                    
        y_np1_ip1,dict,rand1,rand2=scipy.optimize.fsolve(F,y_np1_i,full_output=True)
        self.statistics["nfcns"] += dict['nfev']
            
        return t_np1, y_np1_ip1
        
             
        

            
    def print_statistics(self, verbose=NORMAL):
        self.log_message('Final Run Statistics            : {name} \n'.format(name=self.problem.name),        verbose)
        self.log_message(' Step-length                    : {stepsize} '.format(stepsize=self.options["h"]), verbose)
        self.log_message(' Number of Steps                : '+str(self.statistics["nsteps"]),          verbose)               
        self.log_message(' Number of Function Evaluations : '+str(self.statistics["nfcns"]),         verbose)
            
        self.log_message('\nSolver options:\n',                                    verbose)
        self.log_message(' Solver            : BDF3',                     verbose)
        self.log_message(' Solver type       : Newton\n',                      verbose)

