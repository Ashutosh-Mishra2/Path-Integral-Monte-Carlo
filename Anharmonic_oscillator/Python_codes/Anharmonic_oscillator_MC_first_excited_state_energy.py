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


def Action(V,x1,x2,x,m,omega):
    S = 0.5 * m * ( (x1 - x)**2 + (x2 - x)**2) + V(lamda,x)  
    return S  

def Action_total(V,path,m,omega,lamda,Ntau):
    S = 0
    for i in range(Ntau):
        S += 0.5*m*( (path[(i + 1)%Ntau] - path[i] )**2 ) + V(lamda,path[i])
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
    

def expectation_value_corr(f,tau,n,path_arr,Ntau,N,dt):
    s0 = 0
    for i in path_arr:
        for j in i:
            s0 += f(j,n,dt)
    s0 = s0/(N*Ntau)
    return s0



def autocorrelation_function(f,tau,n,path_arr,Ntau,N,dt):
    s = 0
    for i in path_arr:
        for j in range(Ntau):
            s += f(i[j],i[(tau+j)%Ntau],n,dt)
    s = s/(Ntau*N)
    return s



def Xop(x,n,dt):
    return (x*dt)**n


def Xop_autocorr(x1,x2,n,dt):
    return (x1*dt)**n * (x2*dt)**n



def two_point_correlation(path_arr,m,omega,lamda,Ntau,N,dt,tau):
    m = m /dt
    omega = omega/dt
    x_autocorr = autocorrelation_function(Xop_autocorr,tau,1,path_arr,Ntau,N,dt)
    x_0  = expectation_value_corr(Xop,tau,1,path_arr,Ntau,N,dt)
    print(x_autocorr," ", x_0)
    return x_autocorr - (x_0**2)



def generate_paths(path0,h,N,Ntau,m,omega,lamda,dt):
    path = path0

    for i in range(100):
        path, h = MC_sweep(path,Ntau,h,V_Anharm,m,omega,lamda);
    

    thermalized_path = path

    path_arr = []

    for i in range(N):
        h = 0.1
        path = thermalized_path
        for j in range(Nsep):
            path, h = MC_sweep(path,Ntau,h,V_Harm)
        path_arr.append(list(path))
    
    return path_arr

def G2(n_arr,dt,lamda,N,h):
    G2_arr = []
    Ntau = int(2000/dt)
    m = 1 * dt
    omega = 1 * dt
    path0 = np.zeros(Ntau)
    path_arr = generate_paths(path0,h,N,Ntau,m,omega,lamda,dt)
    for n in n_arr:
        G2_val = two_point_correlation(path_arr,m,omega,lamda,Ntau,N,dt,n)
        G2_arr.append(G2_val)
    
    return G2_arr


def first_excited_state_energy(n_arr,dt_arr,lamda,N,h):
    E1 = []
    for dt in dt_arr:
        G2_temp = G2(n_arr,dt,lamda,N,h)
        G2_temp = np.log(G2_temp)
        slope = -(G2_temp[1] - G2_temp[0])/(n_arr[1] - n_arr[0])
        E1.append(slope)
    return E1

# %% 
N = 100
h = 0.1
n_arr = [1,3]
Nsep = 50
idrate = 0.8

lamda = 0 
dt_arr1 = [0.1, 0.2, 0.4, 0.5]
E1_arr_0 = first_excited_state_energy(n_arr,dt_arr1,lamda,N,h)

lamda = 1 
dt_arr2 = [0.1, 0.2, 0.4, 0.5]
E1_arr_1 = first_excited_state_energy(n_arr,dt_arr2,lamda,N,h)

lamda = 50 
dt_arr3 = [0.05,0.1, 0.2, 0.5]
E1_arr_50 = first_excited_state_energy(n_arr,dt_arr3,lamda,N,h)

lamda = 1000
n_arr = [1,2] 
dt_arr4 = [0.05,0.1, 0.2]
E1_arr_1000 = first_excited_state_energy(n_arr,dt_arr4,lamda,N,h)



plt.plot(dt_arr1, E1_arr_0 )
plt.xscale("log")
plt.xlabel(r"\delta t")
plt.ylabel(r"$E_1 - E_0$")
plt.savefig("first_excited_state_energy_0.svg")


plt.plot(dt_arr2, E1_arr_1 )
plt.xscale("log")
plt.xlabel(r"\delta t")
plt.ylabel(r"$E_1 - E_0$")
plt.savefig("first_excited_state_energy_1.svg")

plt.plot(dt_arr3, E1_arr_50 )
plt.xscale("log")
plt.xlabel(r"\delta t")
plt.ylabel(r"$E_1 - E_0$")
plt.savefig("first_excited_state_energy_50.svg")

plt.plot(dt_arr3, E1_arr_1000)
plt.xscale("log")
plt.xlabel(r"\delta t")
plt.ylabel(r"$E_1 - E_0$")
plt.savefig("first_excited_state_energy_1000.svg")
