#%%
import numpy as np
import matplotlib.pyplot as plt
import time
import random
# %%

#------------------------------------------------------------------
# The following equations are written in natural units 
# hcross = c = 1

# Action
def V_Harm(lamda,x):
    return 0.5* m* (omega**2) * (x**2)


def Action(V,x1,x2,x):
    S = 0.5 * m * ( (x1 - x)**2 + (x2 - x)**2) + V(lamda,x)  
    return S  

def Action_total(V,path):
    S = 0
    for i in range(Ntau):
        S += 0.5*m*( (path[(i + 1)%Ntau] - path[i] )**2 ) + V(lamda,path[i])
    return S


def MC_sweep(path0,Ntau,h,V):
    accept_rate = 0
    path = path0

    index = np.arange(0,Ntau,1)
    random.shuffle(index)
    
    for i in range(Ntau):
        t = index[i]
        tmin = (t + Ntau - 1)% Ntau  # periodic boundary conditions
        tplu = (t + 1) % Ntau 
        x_new = path[t] + h * (random.random() - 0.5)
        S_old = Action(V,path[tplu],path[tmin],path[t])
        S_new = Action(V,path[tplu],path[tmin],x_new)
        if random.random() < np.exp( - (S_new - S_old) ):
            path[t] = x_new
            accept_rate += 1/Ntau

    h = h * accept_rate/idrate

    return path, h


# %%

dt = 1
Ntau = int(120/dt) # Number of time slices
lamda = 2
path = []
h = 0.1
m = 1 * dt
omega = 1 * dt 
idrate = 0.8
time_arr = np.arange(0,Ntau,1)
#path0 = np.array([random.random() for i in range(Ntau)])

path0 = np.zeros(Ntau)
path_arr = []
path = path0
path_arr.append(list(path))

for i in range(100):
    path, h = MC_sweep(path,Ntau,h,V_Harm)
    path_arr.append(list(path))

thermalized_path = path
#%%

N = 1000
Nsep = 12
path_arr = []


for i in range(N):
    h = 0.1
    path = thermalized_path
    for j in range(Nsep):
        path, h = MC_sweep(path,Ntau,h,V_Harm)
    path_arr.append(list(path))


path_arr_hist = np.array(path_arr).flatten()
#%%

# Analytical solution
x = np.arange(-3,3,0.01)
psi0 = np.exp(- m*omega * x**2) * np.sqrt(m*omega/np.pi)

#%%
plt.hist(path_arr_hist,bins=500,density=True,histtype="step")
plt.plot(x,psi0)
plt.xlabel("x")
plt.ylabel(r"$|\psi_0|^2$")
plt.savefig("Harmonic_oscillator_ground_state_1000_runs.png",dpi=400)
plt.show()
plt.close()
#%%

plt.plot(path_arr[N-1], time_arr)
plt.xlabel("x")
plt.xlim((-4,4))
plt.ylabel("Time")
plt.savefig("Harmonic_oscillator_single_path.png",dpi=400)
plt.show()
plt.close()
# %%
