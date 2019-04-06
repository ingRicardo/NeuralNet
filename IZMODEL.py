import pylab
import math
import numpy as np
import matplotlib.pyplot as plt


#http://www.mjrlab.org/wp-content/uploads/2014/05/network_python_tutorial2013.pdf

# initialize parameters
tmax = 1000
dt = 0.5

# Neuron / network pars
a = 0.02 # RS , IB : 0. 0 2 , FS : 0. 1
b = 0.2 # RS , IB , FS : 0. 2
c = -65 # RS , FS : −65 IB : −55
d=8 # RS : 8 , IB : 4 , FS : 2

# Input pars
iapp=10
tr = [200/dt , 700/dt] # stm time

#reserve memory
T = math.ceil(tmax/dt)
v = np.zeros(T)
u = np.zeros(T)
v[0] = -70 #resting potential
u[0] = -14 # steady state

for t in range(T -1):
    if t > tr[0] and t < tr[1]:
        i = iapp
    else:
        i = 0

    if  ( v[ t ] < 35):
        # update ODE
        dv = ( 0.04 * v[ t ] + 5 ) * v[t] + 140 -u[ t ]
        v[ t +1 ] = v[ t ] + ( dv + i )* dt
        du = a * ( b*v[ t ]  - u [ t ] )
        u[ t + 1 ] = u [ t ] + dt * du
    else:
        # spike !
        v [ t ] = 35
        v[ t + 1 ] = c
        u[ t +1 ] = u[ t ] + d
        
plt.figure()
tvec  = np.arange(0, tmax/dt)
plt.plot( tvec, v, 'b', label = 'Voltage trace')
plt.xlabel( ' Time [ms] ' )
plt.ylabel( '  Membrane voltage [mV] ' )
plt.title ( " A single qIF neuron with current step input6 ")
plt.show()
    
