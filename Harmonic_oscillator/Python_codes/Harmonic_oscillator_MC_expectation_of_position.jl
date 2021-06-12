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
    


# %%
dt = 1
Ntau = Int(250/dt) # Number of time slices
lambda = 2
path = []
h = 0.1
m = 1 * dt
omega = 1 * dt 
idrate = 0.8
time_arr = 1:Ntau
#path0 = [rand() for i in 1:Ntau]

path0 = zeros(Ntau)
path_arr = []
path = path0
push!(path_arr,copy(path))

for i in 1:100
    path, h = MC_sweep(path,Ntau,h,V_Harm);
    push!(path_arr,copy(path))
end

thermalized_path = path
#%%

N = 10000
Nsep = 50
path_arr = []

lock = Threads.SpinLock()
Threads.@threads for i in 1:N
    h = 0.1
    path = thermalized_path
    for j in 1:Nsep
        path, h = MC_sweep(path,Ntau,h,V_Harm)
    end
    Threads.lock(lock)
    push!(path_arr,copy(path))
    Threads.unlock(lock)
end

#%%

# Calculate expectation of f(x)
function expectation_value(f,path_arr)
    s = 0
    for i in path_arr
        for j in i
            s += f(j)
        end
    end
    s = s/(Ntau*N)
    return s
end

function expectation_plot(f,path_arr,filename,window_size,y_label,expt)
    x_array = []
    for i in path_arr
        s = 0
        for j in i
            s += f(j)
        end
        s = s/Ntau
        append!(x_array,s)
    end

    moving_average(vs,n) = [sum(@view vs[i:(i+n-1)])/n for i in 1:(length(vs)-(n-1))]

    x_avg = moving_average(x_array,window_size)

    p = scatter(window_size:N, x_array[window_size:N],label = "MC", grid  = false,markerstrokestyle = :dot,markerstrokecolor = :black,markercolor = :black,markersize = 1, markeralpha = 0.6)
    plot!(p,window_size:N, x_avg, label = "Moving average("*string(window_size)*")", grid  = false, linecolor = :red)
    xlabel!(p,"Monte-Carlo Sweep")
    ylabel!(p,y_label)
    title!(y_label*" = "*string(expt))
    savefig(p,filename*".svg")
    display(p)
end

#%%

function Xop(x)
    return x*dt
end

expt_x = expectation_value(Xop,path_arr)
println("Expectation value ⟨X⟩ =  ", expt_x)

y_label = L"$\langle \hat{X} \rangle$"
window_size = 1000

expectation_plot(Xop,path_arr,"Harmonic_oscillator_expectation_of_X",window_size,y_label,expt_x)

#%%

function X2op(x)
    return (x*dt)^2
end

expt_x2 = expectation_value(X2op,path_arr)
println("Expectation value ⟨X²⟩ =  ", expt_x2)

y_label = L"$\langle \hat{X}^2 \rangle$"
window_size = 1000

expectation_plot(X2op,path_arr,"Harmonic_oscillator_expectation_of_X2",window_size,y_label,expt_x2)

#%%

function X3op(x)
    return (x*dt)^3
end

expt_x3 = expectation_value(X3op,path_arr)
println("Expectation value ⟨X³⟩ =  ", expt_x3)

y_label = L"$\langle \hat{X}^3 \rangle$"
window_size = 1000

expectation_plot(X3op,path_arr,"Harmonic_oscillator_expectation_of_X3",window_size,y_label,expt_x3)



#%%
function X4op(x)
    return (x*dt)^4
end

expt_x4 = expectation_value(X4op,path_arr)
println("Expectation value ⟨X⁴⟩ =  ", expt_x4)

y_label = L"$\langle \hat{X}^4 \rangle$"
window_size = 1000

expectation_plot(X4op,path_arr,"Harmonic_oscillator_expectation_of_X4",window_size,y_label,expt_x4)
