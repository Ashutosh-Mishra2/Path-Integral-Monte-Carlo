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

def V_Anharm(lamda,x,m,omega):
    return lamda * (m**2) * (omega**3) * x**4 + 0.5* m* (omega**2) * (x**2)


def V_Harm(lamda,x,m,omega):
    return 0.5* m* (omega**2) * (x**2)


def Action(V,x1,x2,x,m,omega,lamda):
    S = 0.5 * m * ( (x1 - x)**2 + (x2 - x)**2) + V(lamda,x,m,omega)  
    return S  

def Action_total(V,path,m,omega,lamda,Ntau):
    S = 0
    for i in range(Ntau):
        S += 0.5*m*( (path[(i + 1)%Ntau] - path[i] )**2 ) + V(lamda,path[i],m,omega)
    return S


def MC_sweep(path0,Ntau,h,V,m,omega,lamda):
    accept_rate = 0
    path = path0

    index = np.arange(0,Ntau,1)
    random.shuffle(index)
    
    for i in range(Ntau):
        t = index[i]
        tmin = (t + Ntau - 1)% Ntau  # periodic boundary conditions
        tplu = (t + 1) % Ntau 
        x_new = path[t] + h * (random.random() - 0.5)
        S_old = Action(V,path[tplu],path[tmin],path[t],m,omega,lamda)
        S_new = Action(V,path[tplu],path[tmin],x_new,m,omega,lamda)
        if random.random() < np.exp( - (S_new - S_old) ):
            path[t] = x_new
            accept_rate += 1/Ntau

    h = h * accept_rate/idrate

    return path, h    

def expectation_value(f,n,path_arr,Ntau,N,dt):
    s = 0
    for i in path_arr:
        for j in i:
            s += f(j,n,dt)
    s = s/(Ntau*N)
    return s


def Xop(x,n,dt):
    return (x*dt)**n


def ground_state_energy(path_arr,m,omega,lamda,Ntau,N,dt):
    m = m /dt
    omega = omega/dt
    x2 = expectation_value(Xop,2,path_arr,Ntau,N,dt)
    x4 = expectation_value(Xop,4,path_arr,Ntau,N,dt)
    return m*(omega**2)*x2 + 3*lamda*x4



def ground_state(path0,h,N,Ntau,m,omega,lamda,dt):
    path = path0

    for i in range(100):
        path, h = MC_sweep(path,Ntau,h,V_Harm,m,omega,lamda)

    thermalized_path = path

    path_arr = []

    for i in range(N):
        h = 0.1
        path = thermalized_path
        for j in range(Nsep):
            path, h = MC_sweep(path,Ntau,h,V_Harm,m,omega,lamda)
        path_arr.append(list(path))

    E0 = ground_state_energy(path_arr,m,omega,lamda,Ntau,N,dt)
    return E0


def ground_state_vs_dt_plot(dt_arr,N_arr,h_arr,lamda):
    E0_arr = []
    for j in range(len(N_arr)):
        dt = dt_arr[j]
        Ntau = int(250/dt)
        m = 1 * dt
        omega = 1 * dt
        time_arr = np.arange(0,Ntau,1)
        path0 = np.zeros(Ntau)
        E0_val = ground_state(path0,h_arr[j],N_arr[j],Ntau,m,omega,lamda,dt)
        E0_arr.append(E0_val)
    return E0_arr


# %%
dt_arr = [0.10,0.20,0.25,0.40,0.50,1.00]
N_arr = [500,500,500,500,500,500]
h_arr = [0.1,0.1,0.1,0.1,0.1,0.1]
Nsep = 100
lamda = 0 
idrate = 0.8

E0_arr = ground_state_vs_dt_plot(dt_arr,N_arr,h_arr,lamda)
print(E0_arr)
plt.scatter(dt_arr,E0_arr)
plt.plot(dt_arr,E0_arr)
plt.xscale("log")
plt.xlim((0.01,1))
plt.xlabel(r"$\delta t$")
plt.ylabel(r"$E_0$")
plt.savefig("Ground_state_energy_0.svg")
plt.show()
plt.close()
#%%
dt_arr = [0.05,0.10,0.20,0.25,0.40,0.50,1.00]
N_arr = [500,500,500,500,500,500,500]
h_arr = [0.1,0.1,0.1,0.1,0.1,0.1,0.1]
Nsep = 100
lamda = 1 
idrate = 0.8

E0_arr = ground_state_vs_dt_plot(dt_arr,N_arr,h_arr,lamda)
print(E0_arr)
plt.scatter(dt_arr,E0_arr)
plt.plot(dt_arr,E0_arr)
plt.xscale("log")
plt.xlim((0.01,1))
plt.xlabel(r"$\delta t$")
plt.ylabel(r"$E_0$")
plt.savefig("Ground_state_energy_1.svg")
plt.show()
plt.close()
#%%

dt_arr = [0.02,0.05,0.10,0.20,0.25,0.40,0.50,1.00]
N_arr = [300,500,500,500,500,500,500,500]
h_arr = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
Nsep = 100
lamda = 50
idrate = 0.8

E0_arr = ground_state_vs_dt_plot(dt_arr,N_arr,h_arr,lamda)
print(E0_arr)
plt.scatter(dt_arr,E0_arr)
plt.plot(dt_arr,E0_arr)
plt.xscale("log")
plt.xlim((0.01,1))
plt.xlabel(r"$\delta t$")
plt.ylabel(r"$E_0$")
plt.savefig("Ground_state_energy_50.svg")
plt.show()
plt.close()
#%%

dt_arr = [0.01,0.02,0.05,0.10,0.20,0.25,0.40,0.50,1.00]
N_arr = [300,300,500,500,500,500,500,500,500]
h_arr = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
Nsep = 100
lamda = 1000
idrate = 0.8

E0_arr = ground_state_vs_dt_plot(dt_arr,N_arr,h_arr,lamda)
print(E0_arr)
plt.scatter(dt_arr,E0_arr)
plt.plot(dt_arr,E0_arr)
plt.xscale("log")
plt.xlim((0.01,1))
plt.xlabel(r"$\delta t$")
plt.ylabel(r"$E_0$")
plt.savefig("Ground_state_energy_1000.svg")
plt.show()
plt.close()