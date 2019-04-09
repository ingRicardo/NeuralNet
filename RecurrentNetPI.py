from pylab import *
from math import ceil
from numpy.random import uniform
from scipy.sparse import csr_matrix
from scipy.linalg.special_matrices import circulant


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

#NEW recurrent parameter
width = pi/4 # half-width of the orientation tuning
w = 0.005 # average recurrent weight
pconn= 0.4 #set a bit higher
scaleEI = 2 # scale I->E
g_sc = 0.002 # scale of gamma
E = inh.choose(0,-85)
# NEW make weight matrix
W = zeros((n,n))
C = uniform(size = (n, n))
idx = nonzero(C<pconn) # sparse connectivity
W[idx] = gamma(w/g_sc, scale=g_sc, size=idx[0].size)
W[ix_(exc,inh)] *= scaleEI #submat indexing
theta = linspace(0, 2*pi, n) #NEW
R = circulant(cos(theta))>cos(width) #NEW
W[:,exc] = where(R[:,exc], W[:,exc],0)#NEW
W = csr_matrix(W) # make row sparse

# 1.2) Input pars
tr = ([int(200/dt), int(700/dt)])
rate_in = 2 #input rate
inwidth = pi/2
w_in = 0.07 # input weights
pconn_in = 0.2 #input conn prob.
n_in = 100 # number of inputs
C = uniform(size=(n, n_in))<pconn_in
W_in = C.choose(0, w_in) #  matrix
print("W_in -> ", W_in)
W_in[int(n/2),:]=0 #NEW

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
s= zeros(n) # rec synapses

#lastsp = -infty*ones(n_in)


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
    
    #NEW recurrent input
    s = (1 - dt/tau_s)*s +fired
    lsyn = W.dot(s*E) - W.dot(s)*v[t]
    i += lsyn # add to input vector
    
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
title("""A recurrent circular network of %d qIF neurons """% n)
legend(loc='upper right')
show()

