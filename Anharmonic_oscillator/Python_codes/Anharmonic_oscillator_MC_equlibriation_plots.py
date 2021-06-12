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
    
def equilibration(path0,h):
    
    path_arr = []
    path = path0
    path_arr.append(list(path))

    for i in range(100):
        path, h = MC_sweep(path,Ntau,h,V_Harm)
        path_arr.append(list(path))

    return path_arr, path

def expectation_plot(f,path_arr,filename,window_size,y_label):
    x_array = []
    for i in path_arr:
        s = 0
        for j in i:
            s += f(j)
        s = s/Ntau
        x_array.append(s)

    plt.plot(range(N), x_array)
    plt.xlabel("Monte-Carlo Sweeps")
    plt.ylabel(y_label)
    plt.savefig(filename+".svg")
    plt.show()
    plt.close()


def X2op(x):
    return (x)**2




# %%
dt = 0.2
Ntau = int(250/dt) # Number of time slices
lamda = 1
path = []
h = 0.1
m = 1 * dt
omega = 1 * dt 
idrate = 0.8
time_arr = np.arange(0,Ntau,1)
#path0 = [random.random() for i in 1:Ntau]

path0 = np.zeros(Ntau)

eq_arr_a, final_path = equilibration(path0,h)


N = len(eq_arr_a)

y_label = r"$\langle \hat{X}^2 \rangle$"
window_size = 2

expectation_plot(X2op,eq_arr_a,"Equilibration_plot_a",window_size,y_label)


#%%

dt = 1
Ntau = int(250/dt) # Number of time slices
lamda = 1
m = 1 * dt
omega = 1 * dt 
time_arr = np.arange(0,Ntau,1)
path0 = np.zeros(Ntau)

eq_arr_b, final_path = equilibration(path0,h)

expectation_plot(X2op,eq_arr_b,"Equilibration_plot_b",window_size,y_label)


#%%

dt = 0.2
Ntau = int(250/dt) # Number of time slices
lamda = 1e3
m = 1 * dt
omega = 1 * dt 
time_arr = np.arange(0,Ntau,1)
path0 = np.zeros(Ntau)

eq_arr_c, final_path = equilibration(path0,h)

expectation_plot(X2op,eq_arr_c,"Equilibration_plot_c",window_size,y_label)

#%%

dt = 1
Ntau = int(250/dt) # Number of time slices
lamda = 1e3
m = 1 * dt
omega = 1 * dt 
time_arr = np.arange(0,Ntau,1)
path0 = np.zeros(Ntau)

eq_arr_d, final_path = equilibration(path0,h)

expectation_plot(X2op,eq_arr_d,"Equilibration_plot_d",window_size,y_label)

#%%



