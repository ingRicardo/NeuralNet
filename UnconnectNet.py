from pylab import *
from math import ceil
from numpy.random import uniform


# 1) initialize parameters
tmax = 1000
dt = 0.5

#1.1) Neuron / Network pars
n = 1000 # number of neurons
pinh = 0.2 #prob of inh neuron
inh = (uniform(size=n)<pinh) # whether inh.
exc = logical_not(inh)
a = inh.choose(0.02,0.1) # exc = 0.02, inh =0.1
b = 0.2
c = -65
d = inh.choose(8,2) # exc=8, inh = 2
tau_s = 10 # decay of synapses [ms]

#tau_d = 500 #synaptic depression [ms]
#std_u = 0.5 # STP parameter

# 1.2) Input pars
tr = ([int(200/dt), int(700/dt)])
rate_in = 2 #input rate
n_in = 100 # number of inputs
w_in = 0.07 # input weights
pconn_in = 0.1 #input conn prob.

C = uniform(size=(n, n_in))<pconn_in
W_in = C.choose(0, w_in) #  matrix

# 2) reserve memory 
T = ceil(tmax/dt)
print("T -> ",T)
v = zeros((T,n)) # now matrix
print("v -> ", v)
u = zeros((T,n)) # now matrix
v[0] = -70 # set 1st row
u[0] = -14 
s_in = zeros(n_in) 
E_in = zeros(n_in) 
prate = dt*rate_in*1e-3 # abbrev
#h = ones(n_in)
lastsp = -infty*ones(n_in)


# 3) for_loop over time
for t in arange(T-1):
    # 3.1) get input
    if  t>tr[0] and t<tr[1] :
        # NEW: get input Poisson spikes
        p = (uniform(size = n_in)) < prate;
        
    else:

        p = 0; # no input
     
    # NEW: calculate input current
    s_in = (1 - dt/tau_s)*s_in + p
    i = W_in.dot(s_in*E_in)
    i -= W_in.dot(s_in)*v[t]
    
    print("v[t] -> ",v[t])
    #NEW: handle all neurons
    fired = v[t]>=35
    
    print("fired -> ", fired)
  #  if fired:
        # 3.2) update ODE , simply update all
    dv = (0.04*v[t]+5)*v[t]+140-u[t]
    res = v[t] + (dv+i)*dt
    v[t+1] = v[t] + (dv+i)*dt
    du = a *(b*v[t] - u[t])
    u[t+1] = u[t] + dt*du
  #  else:
    # 3.3) spike !
    v[t][fired] = 35
    v[t+1][fired] = c
    u[t+1][fired] = u[t][fired] +d[fired]
   
  
# 4) plot voltage trace
# NEW: get spikes and plot
tspk, nspk =  nonzero(v == 35)
idx_i = in1d(nspk,nonzero(inh)[0]) # find inh
idx_e = logical_not(idx_i) # all others are exc

figure()
#tvec = arange(0, tmax, dt)
plot(tspk[idx_e]*dt,nspk[idx_e], 'k', label = 'Exc', markersize = 2)
plot(tspk[idx_i]*dt, nspk[idx_i],'r', label = 'inh', markersize = 2)
#plot(tvec,v,'b',label='Voltage trace')
xlabel('Time [ms]')
ylabel('Neuron number[\#]')
xlim((0,tmax))
title("""An unconnected network of %d qIF neurons """% n)
legend(loc='upper right')
show()

