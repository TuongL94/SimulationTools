from assimulo.explicit_ode import Explicit_ODE
from assimulo.ode import *
from scipy.optimize import fsolve
import numpy as np
import matplotlib.pyplot as mpl
import scipy.linalg as SL


class BDF_4(Explicit_ODE):
    tol=1.e-8
    maxit=100
    maxsteps=100000

    def __init__(self, problem):
        Explicit_ODE.__init__(self, problem) #Calls the base class

        #Solver options
        self.options["h"] = 0.01

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

        #List/array for storing the 4 time and function values needed to compute new function values
        tprev = np.zeros(4)
        yprev = []

        for i in range(self.maxsteps):
            if i > 2:
                if tprev[-1] >= tf:
                    break
            else:
                if tprev[i] >= tf:
                    break
            self.statistics["nsteps"] += 1

            if i >= 0 and i < 3:
                if i == 0:
                    tprev[i] = t
                    yprev.append(y)
                t_new,y_new = self.step_EE(tprev[i],yprev[i], h)
                tprev[i+1] = t_new
                yprev.append(y_new)
                tres.append(tprev[i+1])
                yres.append(yprev[i+1])
                continue
            else:
                t_np1, y_np1 = self.step_BDF4([tprev[0],tprev[1],tprev[2],tprev[3]], [yprev[0],yprev[1],yprev[2],yprev[3]],h)
            tprev[0],tprev[1],tprev[2],tprev[3]= tprev[1],tprev[2],tprev[3],t_np1
            yprev[0],yprev[1],yprev[2],yprev[3]= yprev[1],yprev[2],yprev[3],y_np1

            tres.append(tprev[-1])
            yres.append(yprev[-1])

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

    def step_BDF4(self,T,Y, h):
        """
        BDF-4 with Newton iteration used for correction and zero order predictor

        alpha_0*y_np1+alpha_1*y_n+alpha_2*y_nm1+alpha_3*y_nm2+alpha_4*y_nm3=h f(t_np1,y_np1)
        alpha=[25/12,-4,3,-4/3,1/4]
        """
        alpha=[25./12.,-4.,3.,-4./3.,1./4.]
        f=self.problem.rhs

        # predictor
        t_np1=T[-1]+h
        y_np1_i=Y[-1]   # zero order predictor
        # corrector with Newton's iteration

            #self.statistics["nfcns"] += 1
        def g(x):
            return alpha[0]*x + alpha[1]*Y[3] + alpha[2]*Y[2] + alpha[3]*Y[1] + alpha[4]*Y[0]-h*f(t_np1,x)

        y_np1,dict,ier,mesg = fsolve(g,y_np1_i,full_output=True)
        self.statistics["nfcns"] += dict['nfev']
        return t_np1,y_np1

    def print_statistics(self, verbose=NORMAL):
        self.log_message('Final Run Statistics            : {name} \n'.format(name=self.problem.name),        verbose)
        self.log_message(' Step-length                    : {stepsize} '.format(stepsize=self.options["h"]), verbose)
        self.log_message(' Number of Steps                : '+str(self.statistics["nsteps"]),          verbose)
        self.log_message(' Number of Function Evaluations : '+str(self.statistics["nfcns"]),         verbose)

        self.log_message('\nSolver options:\n',                                    verbose)
        self.log_message(' Solver            : BDF4',                     verbose)
        self.log_message(' Solver type       : Fixed step\n',                      verbose)
