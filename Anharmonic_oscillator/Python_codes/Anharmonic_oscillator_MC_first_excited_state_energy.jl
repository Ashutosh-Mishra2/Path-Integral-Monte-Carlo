using Plots: append!
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

function V_Anharm(lambda,x,m,omega)
    return lambda * (m^2) * (omega^3) * x^4 + 0.5* m* (omega^2) * (x^2)
end

function V_Harm(lambda,x,m,omega)
    return 0.5* m* (omega^2) * (x^2)
end

function Action(V,x1,x2,x,m,omega,lambda)
    S = 0.5 * m * ( (x1 - x)^2 + (x2 - x)^2) + V(lambda,x,m,omega)  
    return S
end  

function Action_total(V,path,m,omega,lambda)
    S = 0
    for i in 1:Ntau
        S += 0.5*m*( (path[mod1(i + 1,Ntau)] - path[i] )^2 ) + V(lambda,path[i],m,omega)
    return S
    end
end

function MC_sweep(path0,Ntau,h,V,m,omega,lambda)
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
        S_old = Action(V,path[tplu],path[tmin],path[t],m,omega,lambda)
        S_new = Action(V,path[tplu],path[tmin],x_new,m,omega,lambda)
        if rand() < exp( - (S_new - S_old) )
            path[t] = x_new
            accept_rate += 1/Ntau
            #append!(path_arr,path)
        end
    end
    
    h = h * accept_rate/idrate

    return path, h
end
    

function expectation_value_corr(f,tau,n,path_arr,Ntau,N,dt)
    s0 = 0
    for i in path_arr
        for j in i
            s0 += f(j,n,dt)
        end
    end
    s0 = s0/(N*Ntau)
    return s0
end


function autocorrelation_function(f,tau,n,path_arr,Ntau,N,dt)
    s = 0
    for i in path_arr
        for j in 1:Ntau
            s += f(i[j],i[mod1(tau+j,Ntau)],n,dt)
        end
    end
    s = s/(Ntau*N)
    return s
end


function Xop(x,n,dt)
    return (x*dt)^n
end

function Xop_autocorr(x1,x2,n,dt)
    return (x1*dt)^n * (x2*dt)^n
end


function two_point_correlation(path_arr,m,omega,lambda,Ntau,N,dt,tau)
    m = m /dt
    omega = omega/dt
    x_autocorr = autocorrelation_function(Xop_autocorr,tau,1,path_arr,Ntau,N,dt)
    x_0  = expectation_value_corr(Xop,tau,1,path_arr,Ntau,N,dt)
    println(x_autocorr," ", x_0)
    return x_autocorr - (x_0^2)
end


function generate_paths(path0,h,N,Ntau,m,omega,lambda,dt)
    path = path0

    for i in 1:100
        path, h = MC_sweep(path,Ntau,h,V_Anharm,m,omega,lambda);
    end

    thermalized_path = path

    path_arr = []

    lock = Threads.SpinLock()
    @time Threads.@threads for i in 1:N
        h = 0.1
        path = thermalized_path
        for j in 1:Nsep
            path, h = MC_sweep(path,Ntau,h,V_Anharm,m,omega,lambda)
        end
        Threads.lock(lock)
        push!(path_arr,copy(path))
        Threads.unlock(lock)
    end

    return path_arr
end

function G2(n_arr,dt,lambda,N,h)
    G2_arr = []
    Ntau = Int(floor(2000/dt))
    m = 1 * dt
    omega = 1 * dt
    path0 = zeros(Ntau)
    path_arr = generate_paths(path0,h,N,Ntau,m,omega,lambda,dt)
    for n in n_arr
        G2_val = two_point_correlation(path_arr,m,omega,lambda,Ntau,N,dt,n)
        append!(G2_arr,G2_val)
    end
    return G2_arr
end

function first_excited_state_energy(n_arr,dt_arr,lambda,N,h)
    E1 = []
    for dt in dt_arr
        G2_temp = G2(n_arr,dt,lambda,N,h)
        G2_temp = log.(G2_temp)
        slope = -(G2_temp[2] - G2_temp[1])/(n_arr[2] - n_arr[1])
        append!(E1,slope)
    end
    return E1
end

# %% 
N = 100
h = 0.1
n_arr = [1,3]
Nsep = 50
idrate = 0.8

lambda = 0 
dt_arr1 = [0.1, 0.2, 0.4, 0.5]
E1_arr_0 = first_excited_state_energy(n_arr,dt_arr1,lambda,N,h)

lambda = 1 
dt_arr2 = [0.1, 0.2, 0.4, 0.5]
E1_arr_1 = first_excited_state_energy(n_arr,dt_arr2,lambda,N,h)

lambda = 50 
dt_arr3 = [0.05,0.1, 0.2, 0.5]
E1_arr_50 = first_excited_state_energy(n_arr,dt_arr3,lambda,N,h)

lambda = 1000
n_arr = [1,2] 
dt_arr4 = [0.05,0.1, 0.2]
E1_arr_1000 = first_excited_state_energy(n_arr,dt_arr4,lambda,N,h)



plot(dt_arr1, E1_arr_0 , legend = false ,grids = false   ,xaxis = :log,  marker = ([:hex :d]))
xlabel!(L"\delta t")
ylabel!(L"$E_1 - E_0$")
savefig("first_excited_state_energy_0.svg")


plot(dt_arr2, E1_arr_1 , legend = false ,grids = false   ,xaxis = :log,  marker = ([:hex :d]))
xlabel!(L"\delta t")
ylabel!(L"$E_1 - E_0$")
savefig("first_excited_state_energy_1.svg")

plot(dt_arr3, E1_arr_50 , legend = false ,grids = false   ,xaxis = :log,  marker = ([:hex :d]))
xlabel!(L"\delta t")
ylabel!(L"$E_1 - E_0$")
savefig("first_excited_state_energy_50.svg")

plot(dt_arr4, E1_arr_1000 , legend = false ,grids = false   ,xaxis = :log,  marker = ([:hex :d]))
xlabel!(L"\delta t")
ylabel!(L"$E_1 - E_0$")
savefig("first_excited_state_energy_1000.svg")
