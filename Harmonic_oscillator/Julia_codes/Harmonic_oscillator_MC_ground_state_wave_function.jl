#%%
using Plots
pyplot()
using Random
using LaTeXStrings
# %%

#------------------------------------------------------------------
# The following equations are written in natural units 
# hcross = c = 1

# Action
function Action(x1,x2,x) 
    S = 0.5 * m * ( (x1 - x)^2 + (x2 - x)^2 + (omega^2) * (x^2) ) 
    return S
end  

function Action_total(path)
    S = 0
    for i in 1:Ntau
        S += 0.5*m*( (path[mod1(i + 1,Ntau)] - path[i] )^2 + (omega^2) * (path[i]^2) )
    return S
    end
end

function MC_sweep(path0,Ntau,h)
    #path_arr = []
    accept_rate = 0
    path = path0
    #append!(path_arr,path)

    index = collect(1:Ntau)
    shuffle!(index)
    
    for i in 1:Ntau
        t = index[i]
        tmin = mod1(t + Ntau - 1, Ntau)   # periodic boundary conditions
        tplu = t % Ntau + 1 
        x_new = path[t] + h * (rand() - 0.5)
        S_old = Action(path[tplu],path[tmin],path[t])
        S_new = Action(path[tplu],path[tmin],x_new)
        if rand() < exp( - (S_new - S_old) )
            path[t] = x_new
            accept_rate += 1/Ntau
            #append!(path_arr,path)
        end
    end
    
    h = h * accept_rate/idrate

    return path, h
end
    


# %%

Ntau = 120 # Number of time slices
path = []
h = 0.1
m = 1
omega = 1 
idrate = 0.8
dt = 1
time_arr = 1:Ntau
#path0 = [rand() for i in 1:Ntau]

path0 = zeros(Ntau)
path_arr = []
path = path0
append!(path_arr,path)

for i in 1:100
    path, h = MC_sweep(path,Ntau,h);
    #append!(path_arr,path)
end

thermalized_path = path
#%%

N = 10000
Nsep = 12
path_arr = []

lock = Threads.SpinLock()
Threads.@threads for i in 1:N
    h = 0.1
    path = thermalized_path
    for j in 1:Nsep
        path, h = MC_sweep(path,Ntau,h)
    end
    Threads.lock(lock)
    push!(path_arr,copy(path))
    Threads.unlock(lock)
end

path_arr_hist = collect(Iterators.flatten(path_arr))
histogram(path_arr_hist)
#%%

# Analytical solution
x = -3:0.01:3
psi0 = exp.(-x.^2)/sqrt(pi)

#%%
histogram(path_arr_hist,bins=1000,label = "MC",normed=true,grid=false)
plot!(x,psi0,label="Analytical",grid=false)
xlabel!("x")
ylabel!(L"$|\psi_0|^2$")
savefig("Harmonic_oscillator_wavefunction.svg")
#%%

#=  
S_arr = []
for i in path_arr
    append!(S_arr,Action_total(i))
end

plot(1:length(path_arr),S_arr)
# %%

#for i in path_arr
#    plt.plot(time_arr,i)
#end
# %%
#plot( path_arr[1], time_arr)
#plot!(path_arr[N-10], time_arr)
plot(path_arr[N], time_arr,grid = false, legend = false)
xlabel!("x")
xlims!((-4,4))
ylabel!("Time")
savefig("Harmonic_oscillator_single_path.svg") =#
# %%
