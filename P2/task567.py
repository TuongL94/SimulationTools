from scipy import *
from assimulo.problem import Implicit_Problem
from assimulo.solvers import IDA
from Project_2 import squeezer as sq
from Project_2 import squeezer2 as sq2
import matplotlib.pyplot as plt

t0 = 0
tfinal = 0.03

# Model 1, squeezer with index 3 constraints
y0_1,yp0_1 = sq.init_squeezer()
ncp_1 = 100
model_1 = Implicit_Problem(sq.squeezer, y0_1, yp0_1, t0)
model_1.name = 'Squeezer index 3 constraints'
sim_1 = IDA(model_1)
sim_1.suppress_alg = True
sim_1.algvar = array([1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0])
sim_1.atol = ones(20,)
for i in range(size(sim_1.atol)):
    if 0 <= i <= 6:
        sim_1.atol[i] = pow(10,-6)
    else:
        sim_1.atol[i] = pow(10,5)
t, y, yd = sim_1.simulate(tfinal, ncp_1)
plt.plot(t,y[:,14:shape(y)[1]]) # Plots the Lagrange multipliers
plt.show()
plt.plot(t,y[:,0:7]%(2*pi)) # Plots the position variables i.e beta,gamma etc, modulo 2pi.
plt.show()

# Model 2, squeezer with index 2 constraints
y0_2,yp0_2 = sq2.init_squeezer()
ncp_2 = 100
model_2 = Implicit_Problem(sq2.squeezer2, y0_2, yp0_2, t0)
model_2.name = 'Squeezer index 2 constraints'
sim_2 = IDA(model_2)
sim_2.suppress_alg = True
sim_2.algvar = array([1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0])
sim_2.atol = ones(20,)
for i in range(size(sim_2.atol)):
    if 0 <= i <= 6:
        sim_2.atol[i] = pow(10,-6)
    else:
        sim_2.atol[i] = pow(10,5)
t2, y2, yd2 = sim_2.simulate(tfinal, ncp_2)
plt.plot(t2,y2[:,14:shape(y2)[1]]) # Plots the Lagrange multipliers
plt.show()
plt.plot(t,y[:,0:7]%(2*pi)) # Plots the position variables i.e beta,gamma etc, modulo 2pi.
plt.show()
