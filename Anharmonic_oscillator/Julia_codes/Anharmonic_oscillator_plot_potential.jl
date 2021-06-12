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
    return lambda * x^4 + 0.5* m* (omega^2) * (x^2)
end

function V_Harm(x)
    return 0.5* m* (omega^2) * (x^2)
end

x = -2:0.01:2
omega = 1
m = 1
v_harm_arr = V_Harm.(x)
v_anharm_arr = V_Anharm.(1e3,x)
v_anharm_arr2 = V_Anharm.(1,x)
v_anharm_arr3 = V_Anharm.(10,x)



plot(x,v_harm_arr,label="Harmonic Oscillator",legend = :bottomleft,grid=false)
plot!(x,v_anharm_arr2,label="Anharmonic Oscillator (λ = 1)")
plot!(x,v_anharm_arr3,label="Anharmonic Oscillator (λ = 10)")
plot!(x,v_anharm_arr,label="Anharmonic Oscillator (λ = 1000)")
ylims!((0,2))
xlabel!("x")
ylabel!("V(x)")
savefig("potential_plot.svg")
