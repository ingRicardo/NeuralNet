from pylab import *
from math import ceil
from numpy.core.multiarray import zeros, arange


tmax = 1000
dt= 0.5

a = 0.02
b = 0.2
c = -65

d=8

iapp = 10
tr = [200/dt, 700/dt]
T = ceil(tmax/dt)
print T
v = zeros(int(T))
print v
u = zeros(int(T))
print u
v[0] = -70
u[0] = -14

for t in arange(T-1):
    if t> tr[0] and t<tr[1]:
        i= iapp
    else:
        i=0
        
    if v[int(t)] < 35:
        dv = (0.04 * v[int(t)] + 5) * v[int(t)]+140 - u[int(t)]
        v[int(t) + 1] = v[int(t)] + (dv + i) * dt
        du = a* (b * v[int(t)] - u[int(t)])
        u[int(t) + 1] = u[int(t)] + dt*du
    else:
        v[int(t)] = 35
        v[int(t) + 1]= c
        u[ int(t) + 1] = u[int(t)] + d

figure()
tvec = arange(0, tmax, dt)
plot(tvec, v, 'b', label = 'Voltage trace')
xlabel('Time [ms]')
ylabel('Membrane voltage [mV]')
title("A single qIF neuron with current step input 6")
show()