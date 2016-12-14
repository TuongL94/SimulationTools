from numpy import *
from assimulo.problem import Explicit_Problem
from assimulo.solvers import Dopri5
from Project_2 import squeezer3 as sq
from Project_2 import modulo
import matplotlib.pyplot as plt

t0 = 0
tfinal = 0.03
y0,yp0 = sq.init_squeezer()
initial_conditions = y0
model = Explicit_Problem(sq.squeezer3,initial_conditions,t0)
sim = Dopri5(model)
t,y = sim.simulate(tfinal)
plt.plot(t,modulo.mod(y[:,0:7],2*pi))
plt.show()
plt.plot(t,y[:,14:shape(y)[1]])
plt.show()
