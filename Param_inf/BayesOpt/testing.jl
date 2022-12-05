using PlotlyJS, PyCall, DifferentialEquations, StaticArrays, BenchmarkTools, DataFrames, CSV, OrderedCollections, ScatteredInterpolation
include("/home/holliehindley/phd/rtc_models/Oct2022_model/rtc_model.jl")
include("/home/holliehindley/phd/rtc_models/sol_species_funcs.jl")
include("/home/holliehindley/phd/Param_inf/inf_setup.jl")
include("/home/holliehindley/phd/rtc_models/params_init_tspan.jl")


function rtc_bo_ω1(;ω_ab, ω_r)
    prob = ODEProblem(rtc_model_OD, (@SVector [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, rh_0, rd_0, rt_0, OD_0_wt3]), tspan2, (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, findmax(WT3)[1]]))
    solu = solve(prob, Rodas4(), saveat=t_2)
    obj = []
    df = DataFrame(solu)
    rename!(df, [:time, :rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt, :OD])
    OD = df[:, :OD]
    append!(obj, [abs2((i-j)/k) for (i,j,k) in zip(OD, WT3, WT3_std)])
    obj_final = sum(obj)
    return -(obj_final/(36))   
end

# in python writing the ranges to search for parameter value
py"""
param_range_ω = {'ω_ab': (0, 100), 'ω_r': (0, 100)}
"""

# import bayes_opt package from python
bayes_opt = pyimport("bayes_opt")

# setting the optimizer 
optimizer = bayes_opt.BayesianOptimization(f=rtc_bo_ω1, pbounds=py"param_range_ω", random_state=27, verbose=2) # verbose = 1 prints only when a maximum is observed (pink)

# timing the process and maximising the optimizer 
function timer()
    optimizer.maximize(init_points=2, n_iter=100, acq="ucb")#, xi=0.1) #kappa=2, 
end

@time timer()
print(optimizer.max)




# creating lists of values of parameters tried and errors for each one  
vals_ab, vals_r, errors, best_param_ab, best_param_r, best_error = results_two_param(optimizer)

# plot errors over iterations 
plot(range(10,length(errors)), errors[10:end], mode="markers+lines", Layout(xaxis_title="Iteration", yaxis_title="Error"))#, yaxis=(:log10, (1,Inf)))
plot(range(1,length(errors)), errors, mode="markers+lines", Layout(xaxis_title="Iteration", yaxis_title="Error"))#, yaxis=(:log10, (1,Inf)))
# png("error")

layout1 = Layout(xaxis_title="ω_ab", yaxis_title="ω_r", zaxis_title="Error")
plot(contour(z=errors, x=vals_ab, y=vals_r), layout1)


# 3D plot
x,y,z = plot_3D(vals_ab, vals_r, errors)

layout = Layout(scene=attr(xaxis_title="ω_ab", yaxis_title="ω_r", zaxis_title="Error"))
plot(surface(x=x,y=y, z=z), layout)



