#%%
import numpy as np
import matplotlib.pyplot as plt
import time
import random
from IPython.display import display, Latex
# %%

#------------------------------------------------------------------
# The following equations are written in natural units 
# hcross = c = 1

# Action
def Action(x1,x2,x): 
    S = 0.5 * m * ( (x1 - x)**2 + (x2 - x)**2 + (omega**2) * (x**2) ) 
    return S  

def Action_total(path):
    S = 0
    for i in range(Ntau):
        S += 0.5*m*( (path[(i+1)%Ntau] - path[i] )**2 + (omega**2) * (path[i]**2) )
    return S

def MC_sweep(path0,Ntau,h):
    path_arr = []
    accept_rate = 0
    path = path0
    path_arr.append(list(path))

    index = np.arange(0,Ntau,1)
    random.shuffle(index)
    
    for i in range(Ntau):
        t = index[i]
        tmin = (t + Ntau - 1)% Ntau  # periodic boundary conditions
        tplu = (t + 1) % Ntau 
        x_new = path[t] + h * (random.random() - 0.5)
        S_old = Action(path[tplu],path[tmin],path[t])
        S_new = Action(path[tplu],path[tmin],x_new)
        if random.random() < np.exp( - (S_new - S_old) ):
            path[t] = x_new
            accept_rate += 1/Ntau
            path_arr.append(list(path))
    
    h = h * accept_rate/idrate

    return path, h


# %%

## Initial runs for the system to find suitable paths
Ntau = 1200# Number of time slices
path = []
h = 0.1
m = 1
omega = 1 
idrate = 0.8
dt = 1
time_arr = np.arange(0,Ntau,1)
#path0 = np.array([random.random() for i in range(Ntau)])
path0 = np.zeros(Ntau)
path_arr = []
path = path0
path_arr.append(list(path))
for i in range(100):
    path, h = MC_sweep(path,Ntau,h)
    #path_arr.append(list(path))
thermalized_path = path

#%%

N = 1000
Nsep = 1000
path_arr = []
start = time.time()
for i in range(N):
    path = thermalized_path
    for j in range(Nsep):
        path, h = MC_sweep(path,Ntau,h)
    path_arr.append(list(path))
path_arr_hist = np.array(path_arr).flatten()

end = time.time()
print("Time taken = ",end - start," seconds")
#%%
