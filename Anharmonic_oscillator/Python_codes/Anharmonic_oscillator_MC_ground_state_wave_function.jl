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

function V_Anharm(lambda,x)
    return lambda * (m^2) * (omega^3) * x^4 + 0.5* m* (omega^2) * (x^2)
end

function V_Harm(lambda,x)
    return 0.5* m* (omega^2) * (x^2)
end

function Action(V,x1,x2,x)
    S = 0.5 * m * ( (x1 - x)^2 + (x2 - x)^2) + V(lambda,x)  
    return S
end  

function Action_total(V,path)
    S = 0
    for i in 1:Ntau
        S += 0.5*m*( (path[mod1(i + 1,Ntau)] - path[i] )^2 ) + V(lambda,path[i])
    return S
    end
end

function MC_sweep(path0,Ntau,h,V)
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
        S_old = Action(V,path[tplu],path[tmin],path[t])
        S_new = Action(V,path[tplu],path[tmin],x_new)
        if rand() < exp( - (S_new - S_old) )
            path[t] = x_new
            accept_rate += 1/Ntau
            #append!(path_arr,path)
        end
    end
    
    h = h * accept_rate/idrate

    return path, h
end
    

function ground_state(path0,h,filename)
    path = path0

    for i in 1:100
        path, h = MC_sweep(path,Ntau,h,V_Anharm);
    end

    thermalized_path = path

    path_arr = []

    lock = Threads.SpinLock()
    Threads.@threads for i in 1:N
        h = 0.1
        path = thermalized_path
        for j in 1:Nsep
            path, h = MC_sweep(path,Ntau,h,V_Anharm)
        end
        Threads.lock(lock)
        push!(path_arr,copy(path))
        Threads.unlock(lock)
    end

    path_arr_hist = collect(Iterators.flatten(path_arr))

    # Analytical solution
    x = minimum(path_arr_hist):0.01:maximum(path_arr_hist)
    psi0 = exp.(- m*omega * x.^2) * sqrt(m*omega/pi)

    
    p = histogram(path_arr_hist,bins=1000,label = "MC",normed=true,grid=false,histtype=:steps)
    #plot!(p, x,psi0,label="Analytical",grid=false)
    xlabel!(p,"x")
    ylabel!(p,L"$|\psi_0|^2$")
    savefig(p,filename*".svg")
    display(p) 
    
end



# %%
dt = 0.1
Ntau = Int(250/dt) # Number of time slices
lambda = 0
h = 0.1
m = 1 * dt
omega = 1 * dt 
idrate = 0.8
time_arr = 1:Ntau
N = 2000
Nsep = 100
#path0 = [rand() for i in 1:Ntau]
path0 = zeros(Ntau)

@time ground_state(path0,h,"Ground_state_Harmonic_oscillator_a")

#%%
dt = 1
Ntau = Int(250/dt) # Number of time slices
lambda = 0
h = 0.1
m = 1 * dt
omega = 1 * dt 
time_arr = 1:Ntau
path0 = zeros(Ntau)
N = 5000
Nsep = 100
@time ground_state(path0,h,"Ground_state_Harmonic_oscillator_b")

#%%

dt = 0.1
Ntau = Int(250/dt) # Number of time slices
lambda = 1e3
h = 0.1
m = 1 * dt
omega = 1 * dt 
time_arr = 1:Ntau
path0 = zeros(Ntau)
N = 2000
Nsep = 100
@time ground_state(path0,h,"Ground_state_Anharmonic_oscillator_c")

#%%

dt = 1
Ntau = Int(250/dt) # Number of time slices
lambda = 1e3
h = 0.1
m = 1 * dt
omega = 1 * dt 
time_arr = 1:Ntau
path0 = zeros(Ntau)
N = 5000
Nsep = 100
@time ground_state(path0,h,"Ground_state_Anharmonic_oscillator_d")
