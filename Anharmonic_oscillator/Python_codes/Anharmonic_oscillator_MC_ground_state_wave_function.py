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

def V_Anharm(lamda,x):
    return lamda * (m**2) * (omega**3) * x**4 + 0.5* m* (omega**2) * (x**2)


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
    

def ground_state(path0,h,filename):
    path = path0

    for i in range(100):
        path, h = MC_sweep(path,Ntau,h,V_Anharm)

    thermalized_path = path

    path_arr = []

    for i in range(N):
        h = 0.1
        path = thermalized_path
        for j in range(Nsep):
            path, h = MC_sweep(path,Ntau,h,V_Anharm)
        path_arr.append(list(path))


    path_arr_hist = np.array(path_arr).flatten()

    # Analytical solution
    x = np.arange(min(path_arr_hist),max(path_arr_hist),0.01)
    psi0 = np.exp(- m*omega * x**2) * np.sqrt(m*omega/np.pi)

    
    plt.hist(path_arr_hist,bins=500,density=True,histtype="step")
    plt.plot(x,psi0)
    plt.xlabel("x")
    plt.ylabel(r"$|\psi_0|^2$")
    plt.savefig(filename+".png",dpi=400)
    plt.show()
    plt.close()



# %%
dt = 0.1
Ntau = int(120/dt) # Number of time slices
lamda = 0
h = 0.1
m = 1 * dt
omega = 1 * dt 
idrate = 0.8
time_arr = np.arange(0,Ntau,1)
N = 200
Nsep = 10
#path0 = [random.random() for i in 1:Ntau]
path0 = np.zeros(Ntau)

ground_state(path0,h,"Ground_state_Harmonic_oscillator_a")

#%%
dt = 1
Ntau = int(120/dt) # Number of time slices
lamda = 0
h = 0.1
m = 1 * dt
omega = 1 * dt 
time_arr = np.arange(0,Ntau,1)
path0 = np.zeros(Ntau)
N = 500
Nsep = 10
ground_state(path0,h,"Ground_state_Harmonic_oscillator_b")

#%%

dt = 0.1
Ntau = int(120/dt) # Number of time slices
lamda = 1e3
h = 0.1
m = 1 * dt
omega = 1 * dt 
time_arr = np.arange(0,Ntau,1)
path0 = np.zeros(Ntau)
N = 200
Nsep = 10
ground_state(path0,h,"Ground_state_Anharmonic_oscillator_c")

#%%

dt = 1
Ntau = int(120/dt) # Number of time slices
lamda = 1e3
h = 0.1
m = 1 * dt
omega = 1 * dt 
time_arr = np.arange(0,Ntau,1)
path0 = np.zeros(Ntau)
N = 500
Nsep = 10
ground_state(path0,h,"Ground_state_Anharmonic_oscillator_d")

# %%
