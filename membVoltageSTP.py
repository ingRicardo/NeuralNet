from pylab import *
from math import ceil
from numpy.random import uniform

# 1) initialize parameters
tmax = 1000
dt = 0.5

#1.1) Neuron / Network pars
a = 0.02
b = 0.2
c = -65
d = 8
tau_s = 10 # decay of synapses [ms]
tau_d = 500 #synaptic depression [ms]
std_u = 0.5 # STP parameter

# 1.2) Input pars
tr = ([int(200/dt), int(700/dt)])
rate_in = 10 #input rate
n_in = 1 # number of inputs
w_in = 0.03 # input weights
W_in = w_in * ones(n_in) #  vector

# 2) reserve memory 
T = ceil(tmax/dt)
v = zeros(T)
u = zeros(T)
v[0] = -70 # resting potential
u[0] = -14 # steady state
s_in = zeros(n_in) # synaptic variable
E_in = zeros(n_in) # rev potential
prate = dt*rate_in*1e-3 # abbrev
h = ones(n_in)
lastsp = -infty*ones(n_in)


# 3) for_loop over time
for t in arange(T-1):
    # 3.1) get input
    if  t>tr[0] and t<tr[1] :
        # NEW: get input Poisson spikes
     #   print(" prate -> ", prate)
      #  print("uniform in range -> ",uniform(size = n_in)*1e-3)
        p = (uniform(size = n_in)) < prate;
        
        #update synaptic depression
        tmp = exp(dt*( lastsp[p] - t)/ tau_d)
        h[p] = 1 - (1+(std_u -1) * h[p]) * tmp
        lastsp[p] = t
      
    else:

        p = 0; # no input
     
    # NEW: calculate input current
    s_in = (1 - dt/tau_s)*s_in + p*h
 
    i = dot(W_in, s_in*E_in)
    i -= dot(W_in, s_in)*v[t]
    
 
    if v[t]<35:
        # 3.2) update ODE
        dv = (0.04* v[t] + 5 ) * v[t]+140 - u[t]
        v[t +1 ] = v[t] + (dv + i) *dt
        du = a *(b*v[t] - u[t])
        u[t+1] = u[t] + dt*du
    else:
        # 3.3) spike !
        v[t] = 35
        v[t+1] = c
        u[t+1] = u[t] +d
   
  
# 4) plot voltage trace
figure()
tvec = arange(0, tmax, dt)

plot(tvec,v,'b',label='Voltage trace')
xlabel('Time [ms]')
ylabel('Membrane voltage [mV]')
title("""A single qIF neuron with %d Poisson inputs """% n_in)
show()

