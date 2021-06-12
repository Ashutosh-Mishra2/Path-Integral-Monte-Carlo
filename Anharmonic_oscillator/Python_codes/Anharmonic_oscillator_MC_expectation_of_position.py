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


# %%
dt = 1
Ntau = int(250/dt) # Number of time slices
lamda = 2
path = []
h = 0.1
m = 1 * dt
omega = 1 * dt 
idrate = 0.8
time_arr = np.arange(0,Ntau,1)
#path0 = [random.random() for i in 1:Ntau]

path0 = np.zeros(Ntau)
path_arr = []
path = path0
path_arr.append(list(path))

for i in range(100):
    path, h = MC_sweep(path,Ntau,h,V_Anharm);
    path_arr.append(list(path))

thermalized_path = path
#%%

N = 10000
Nsep = 50
path_arr = []


for i in range(N):
    h = 0.1
    path = thermalized_path
    for j in range(Nsep):
        path, h = MC_sweep(path,Ntau,h,V_Harm)
    path_arr.append(list(path))



#%%

# Calculate expectation of f(x)
def expectation_value(f,path_arr):
    s = 0
    for i in path_arr:
        for j in i:
            s += f(j)

    s = s/(Ntau*N)
    return s

def expectation_plot(f,path_arr,filename,window_size,y_label,expt):
    x_array = []
    for i in path_arr:
        s = 0
        for j in i:
            s += f(j)
        s = s/Ntau
        x_array.append(s)
    
    def moving_average(vs,n):
        return [sum(vs[i:(i+n-1)])/n for i in range(len(vs)-(n-1))]

    x_avg = moving_average(x_array,window_size)

    plt.scatter(range(window_size,N), x_array[window_size:],label = "MC")
    plt.plot(range(window_size-1,N), x_avg, label = "Moving average("+str(window_size)+")")
    plt.xlabel("Monte-Carlo Sweep")
    plt.ylabel(y_label)
    plt.title(y_label+" = "+str(expt))
    plt.savefig(filename+".png",dpi = 400)
    plt.show()
    plt.close()
#%%


def Xop(x):
    return x*dt

expt_x = expectation_value(Xop,path_arr)
print("Expectation value ⟨X⟩ =  ", expt_x)

y_label = r"$\langle \hat{X} \rangle$"
window_size = 1000

expectation_plot(Xop,path_arr,"Anarmonic_oscillator_expectation_of_X",window_size,y_label,expt_x)

#%%

def X2op(x):
    return (x*dt)**2


expt_x2 = expectation_value(X2op,path_arr)
print("Expectation value ⟨X²⟩ =  ", expt_x2)

y_label = r"$\langle \hat{X}**2 \rangle$"
window_size = 1000

expectation_plot(X2op,path_arr,"Anarmonic_oscillator_expectation_of_X2",window_size,y_label,expt_x2)

#%%


def X3op(x):
    return (x*dt)**3


expt_x3 = expectation_value(X3op,path_arr)
print("Expectation value ⟨X³⟩ =  ", expt_x3)

y_label = r"$\langle \hat{X}**3 \rangle$"
window_size = 1000

expectation_plot(X3op,path_arr,"Anarmonic_oscillator_expectation_of_X3",window_size,y_label,expt_x3)



#%%
def X4op(x):
    return (x*dt)**4

expt_x4 = expectation_value(X4op,path_arr)
print("Expectation value ⟨X⁴⟩ =  ", expt_x4)

y_label = r"$\langle \hat{X}**4 \rangle$"
window_size = 1000

expectation_plot(X4op,path_arr,"Anarmonic_oscillator_expectation_of_X4",window_size,y_label,expt_x4)
