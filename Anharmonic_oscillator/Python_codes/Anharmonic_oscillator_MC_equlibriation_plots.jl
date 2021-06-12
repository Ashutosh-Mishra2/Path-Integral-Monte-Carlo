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
    
function equilibration(path0,h)
    
    path_arr = []
    path = path0
    push!(path_arr,copy(path))

    for i in 1:100
        path, h = MC_sweep(path,Ntau,h,V_Anharm);
        push!(path_arr,copy(path))
    end

    return path_arr, path
end

function expectation_plot(f,path_arr,filename,window_size,y_label)
    x_array = []
    for i in path_arr
        s = 0
        for j in i
            s += f(j)
        end
        s = s/Ntau
        append!(x_array,s)
    end

    p = plot(1:N, x_array[1:N],legend = false, grid  = false,markerstrokestyle = :dot,markerstrokecolor = :black,markercolor = :black,markersize = 1, markeralpha = 0.6)
    xlabel!(p,"Monte-Carlo Sweeps")
    ylabel!(p,y_label)
    savefig(p,filename*".svg")
    display(p)
end

function X2op(x)
    return (x)^2
end



# %%
dt = 0.2
Ntau = Int(250/dt) # Number of time slices
lambda = 1
path = []
h = 0.1
m = 1 * dt
omega = 1 * dt 
idrate = 0.8
time_arr = 1:Ntau
#path0 = [rand() for i in 1:Ntau]

path0 = zeros(Ntau)

eq_arr_a, final_path = equilibration(path0,h)


N = length(eq_arr_a)

y_label = L"$\langle \hat{X}^2 \rangle$"
window_size = 2

expectation_plot(X2op,eq_arr_a,"Equilibration_plot_a",window_size,y_label)


#%%

dt = 1
Ntau = Int(250/dt) # Number of time slices
lambda = 1
m = 1 * dt
omega = 1 * dt 
time_arr = 1:Ntau
path0 = zeros(Ntau)
eq_arr_b, final_path = equilibration(path0,h)

expectation_plot(X2op,eq_arr_b,"Equilibration_plot_b",window_size,y_label)


#%%

dt = 0.2
Ntau = Int(250/dt) # Number of time slices
lambda = 1e3
m = 1 * dt
omega = 1 * dt 
time_arr = 1:Ntau
path0 = zeros(Ntau)
eq_arr_c, final_path = equilibration(path0,h)

expectation_plot(X2op,eq_arr_c,"Equilibration_plot_c",window_size,y_label)

#%%

dt = 1
Ntau = Int(250/dt) # Number of time slices
lambda = 1e3
m = 1 * dt
omega = 1 * dt 
time_arr = 1:Ntau
path0 = zeros(Ntau)
eq_arr_d, final_path = equilibration(path0,h)

expectation_plot(X2op,eq_arr_d,"Equilibration_plot_d",window_size,y_label)

#%%



