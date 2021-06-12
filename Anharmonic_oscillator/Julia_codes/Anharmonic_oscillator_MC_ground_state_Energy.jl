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
    

function expectation_value(f,n,path_arr,Ntau,N,dt)
    s = 0
    for i in path_arr
        for j in i
            s += f(j,n,dt)
        end
    end
    s = s/(Ntau*N)
    return s
end


function Xop(x,n,dt)
    return (x*dt)^n
end

function ground_state_energy(path_arr,m,omega,lambda,Ntau,N,dt)
    m = m /dt
    omega = omega/dt
    x2 = expectation_value(Xop,2,path_arr,Ntau,N,dt)
    x4 = expectation_value(Xop,4,path_arr,Ntau,N,dt)
    return m*(omega^2)*x2 + 3*lambda*x4
end


function ground_state(path0,h,N,Ntau,m,omega,lambda,dt)
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

    E0 = ground_state_energy(path_arr,m,omega,lambda,Ntau,N,dt)
    return E0
end


function ground_state_vs_dt_plot(dt_arr,N_arr,h_arr,lambda)
    E0_arr = []
    for j in 1:length(N_arr)
        dt = dt_arr[j]
        Ntau = Int(250/dt)
        m = 1 * dt
        omega = 1 * dt
        time_arr = 1:Ntau
        path0 = zeros(Ntau)
        E0_val = ground_state(path0,h_arr[j],N_arr[j],Ntau,m,omega,lambda,dt)
        append!(E0_arr,E0_val)
    end
    return E0_arr
end


# %%
dt_arr = [0.10,0.20,0.25,0.40,0.50,1.00]
N_arr = [500,500,500,500,500,500]
h_arr = [0.1,0.1,0.1,0.1,0.1,0.1]
Nsep = 100
lambda = 0 
idrate = 0.8

E0_arr = ground_state_vs_dt_plot(dt_arr,N_arr,h_arr,lambda)
println(E0_arr)
scatter(dt_arr,E0_arr,grid=false,legend=false,xaxis=:log)
plot!(dt_arr,E0_arr,grid=false,legend=false,xaxis=:log)
xlims!((0.01,1))
xlabel!(L"$\delta t$")
ylabel!(L"$E_0$")
savefig("Ground_state_energy_0.svg")

#%%
dt_arr = [0.05,0.10,0.20,0.25,0.40,0.50,1.00]
N_arr = [500,500,500,500,500,500,500]
h_arr = [0.1,0.1,0.1,0.1,0.1,0.1,0.1]
Nsep = 100
lambda = 1 
idrate = 0.8

E0_arr = ground_state_vs_dt_plot(dt_arr,N_arr,h_arr,lambda)
println(E0_arr)
scatter(dt_arr,E0_arr,grid=false,legend=false,xaxis=:log)
plot!(dt_arr,E0_arr,grid=false,legend=false,xaxis=:log)
xlims!((0.01,1))
xlabel!(L"$\delta t$")
ylabel!(L"$E_0$")
savefig("Ground_state_energy_1.svg")
#%%

dt_arr = [0.02,0.05,0.10,0.20,0.25,0.40,0.50,1.00]
N_arr = [300,500,500,500,500,500,500,500]
h_arr = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
Nsep = 100
lambda = 50
idrate = 0.8

E0_arr = ground_state_vs_dt_plot(dt_arr,N_arr,h_arr,lambda)
println(E0_arr)
scatter(dt_arr,E0_arr,grid=false,legend=false,xaxis=:log)
plot!(dt_arr,E0_arr,grid=false,legend=false,xaxis=:log)
xlims!((0.01,1))
xlabel!(L"$\delta t$")
ylabel!(L"$E_0$")
savefig("Ground_state_energy_50.svg")

#%%

dt_arr = [0.01,0.02,0.05,0.10,0.20,0.25,0.40,0.50,1.00]
N_arr = [300,300,500,500,500,500,500,500,500]
h_arr = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
Nsep = 100
lambda = 1000
idrate = 0.8

E0_arr = ground_state_vs_dt_plot(dt_arr,N_arr,h_arr,lambda)
println(E0_arr)
scatter(dt_arr,E0_arr,grid=false,legend=false,xaxis=:log)
plot!(dt_arr,E0_arr,grid=false,legend=false,xaxis=:log)
xlims!((0.01,1))
xlabel!(L"$\delta t$")
ylabel!(L"$E_0$")
savefig("Ground_state_energy_1000.svg")
