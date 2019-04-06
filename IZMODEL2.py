from pylab import * 
import math
import numpy as np
import matplotlib.pyplot as plt
import random as ra
from scipy.sparse import csr_matrix
from numpy.random import uniform
#http://www.mjrlab.org/wp-content/uploads/2014/05/network_python_tutorial2013.pdf
tmax = 1000
dt = 0.5

a = 0.02
b = 0.2
c = -65
d = 8
tau_s = 10
tr = [  200/dt, 700/dt  ]
rate_in = 2
n_in = 100
w_in = 0.07
W_in = w_in* np.ones( n_in )

T = math.ceil( tmax/ dt )
v = np.zeros ( T )
u = np.zeros ( T )

v[ 0 ] = -70
u[ 0 ] = -14

s_in = np.zeros ( n_in )
E_in = np.zeros ( n_in )
prate = dt * rate_in * 1e-3

print(" prate -> ", prate)
print(" n_in -> ", n_in)
#print("ran uniform -> " ,uniform(size=n_in))

h = np.ones(n_in)
lastsp = -np.infty * np.ones( n_in)

for t in np.arange(T-1):
    if t > tr[ 0 ] and t < tr [ 1 ]:
        p = uniform(size = n_in ) < prate;
       
    else:
        p = 0;

print(" p -> ", p)

s_in = ( 1 - dt / tau_s )* s_in + p
i = np.dot(W_in, s_in*E_in)
i -= np.dot(W_in, s_in)*v[ t ]

print(v[ t ])

if v[ t ] < 35:
    dv = (0.04*v[ t ] +5 ) *v[ t ]+140 -u[ t ]
    v[ t+1 ] = v[ t ] + ( dv + i) * dt
    du = a * (b * v[ t ] - u[ t ])
    u[ t + 1] = u [ t ] + dt*du
else :
    v[ t ] = 35
    v[ t + 1] = c
    u[ t + 1] = u[t]+d

figure()
tvec = arange(0, tmax, dt)
plot(tvec,v,'b', label = 'Voltage trace')
xlabel('Time [ms]')
ylabel('Membrane voltage [mV]')
title("A single qIF neuron with %d Poisson inputs"% n_in)
show()

