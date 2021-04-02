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
    return lamda * x**4 + 0.5* m* (omega**2) * (x**2)

def V_Harm(x):
    return 0.5* m* (omega**2) * (x**2)

x = np.arange(-2,2,0.01)
omega = 1
m = 1
v_harm_arr = V_Harm(x)
v_anharm_arr = V_Anharm(1e3,x)
v_anharm_arr2 = V_Anharm(1,x)
v_anharm_arr3 = V_Anharm(10,x)



plt.plot(x,v_harm_arr,label="Harmonic Oscillator")
plt.plot(x,v_anharm_arr2,label="Anharmonic Oscillator (λ = 1)")
plt.plot(x,v_anharm_arr3,label="Anharmonic Oscillator (λ = 10)")
plt.plot(x,v_anharm_arr,label="Anharmonic Oscillator (λ = 1000)")
plt.legend()
plt.ylim((0,2))
plt.xlabel("x")
plt.ylabel("V(x)")
plt.savefig("potential_plot.svg")
plt.show()
plt.close()
# %%
