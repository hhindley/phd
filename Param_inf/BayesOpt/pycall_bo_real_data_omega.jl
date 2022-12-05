using PlotlyJS, PyCall, DifferentialEquations, StaticArrays, BenchmarkTools, DataFrames, CSV, OrderedCollections, ScatteredInterpolation
include("/home/holliehindley/phd/rtc_models/Oct2022_model/rtc_model.jl")
include("/home/holliehindley/phd/rtc_models/sol_species_funcs.jl")
include("/home/holliehindley/phd/Param_inf/inf_setup.jl")
include("/home/holliehindley/phd/rtc_models/params_init_tspan.jl")


function rtc_bo_ω(;ω_ab, ω_r)
    obj_wt_ab, obj_wt_r = obj(rtc_model, initial, (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr]), tspan2, t_2, "mrna", WT1, WT1_std)

    # obj_gr_wt1 = obj_OD(rtc_model_OD, (@SVector [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, rh_0, rd_0, rt_0, OD_0_wt2]), (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, findmax(WT2)[1]]), tspan2, t_2, "OD", WT2, WT2_std)

    # obj_gr_wt2 = obj_OD(rtc_model_OD, (@SVector [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, rh_0, rd_0, rt_0, OD_0_wt3]), (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, findmax(WT3)[1]]), tspan2, t_2, "OD", WT3, WT3_std)

    # obj_gr_wt3 = obj_OD(rtc_model_OD, (@SVector [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, rh_0, rd_0, rt_0, OD_0_wt4]), (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, findmax(WT4)[1]]), tspan4, t_4, "OD", WT4, WT4_std)

    # return -(obj_wt_ab/36)
    # return -(obj_wt_r/36)
    # return -(obj_gr_wt1/36)
    # return -(obj_gr_wt2/36)
    # return -(obj_gr_wt3/(66^2))

    return -(obj_wt_ab/36 + obj_wt_r/36)
    # return -(obj_wt_ab/36 + obj_wt_r/36 + obj_gr_wt1/36)
    # return -(obj_wt_ab/36 + obj_wt_r/36 + obj_gr_wt1/36 + obj_gr_wt2/36)
    # return -(obj_wt_ab/36 + obj_wt_r/36 + obj_gr_wt1/36 + obj_gr_wt2/36 + obj_gr_wt3/(66^2))


end


# in python writing the ranges to search for parameter value
py"""
param_range_ω = {'ω_ab': (0, 1), 'ω_r': (40, 100)}
"""

# py"""
# param_range_ω = {'θtlr': (0, 100), 'g_max': (0, 100)}
# """

# import bayes_opt package from python
bayes_opt = pyimport("bayes_opt")

# setting the optimizer 
optimizer = bayes_opt.BayesianOptimization(f=rtc_bo_ω, pbounds=py"param_range_ω", random_state=27, verbose=2) # verbose = 1 prints only when a maximum is observed (pink)

# timing the process and maximising the optimizer 
function timer()
    optimizer.maximize(init_points=2, n_iter=50, acq="ei", xi=0.01) #kappa=2, xi = 0.0 (prefer exploitation), xi = 0.1 (prefer exploration)
end

@time timer()
print(optimizer.max)

# creating lists of values of parameters tried and errors for each one  
vals_ab, vals_r, errors, best_param_ab, best_param_r, best_error = results_two_param(optimizer)

# plot errors over iterations 
plot(range(1,length(errors)), errors, mode="markers+lines", Layout(xaxis_title="Iteration", yaxis_title="Error", yaxis_range=[0, 100]))#, yaxis=(:log10, (1,Inf)))

plot(range(10,length(errors)), errors[10:end], mode="markers+lines", Layout(xaxis_title="Iteration", yaxis_title="Error"))#, yaxis=(:log10, (1,Inf)))
plot(range(1,length(errors)), errors, mode="markers+lines", Layout(xaxis_title="Iteration", yaxis_title="Error"))#, yaxis=(:log10, (1,Inf)))
# png("error")

layout1 = Layout(xaxis_title="ω_ab", yaxis_title="ω_r", zaxis_title="Error")
plot(contour(z=errors, x=vals_ab, y=vals_r), layout1)


# 3D plot
x,y,z = plot_3D(vals_ab, vals_r, errors)

layout = Layout(scene=attr(xaxis_title="ω_ab", yaxis_title="ω_r", zaxis_title="Error"))
plot(surface(x=x,y=y, z=z), layout)





# will work for one param not two 
vals, errors, errors_ori, best_param, best_error, best_error_ori = results(optimizer)

# plot params over iterations
Plots.plot(range(1,length(optimizer.space.target)), vals, markershapes=[:circle], ylabel=("Param attempt"), xlabel=("Iteration"), legend=false, xticks = 0:1:length(optimizer.space.target))#, size=(1000,500))#, yaxis=(:log10, (1,Inf)))
hline!([4.14], labels= false)
# png("param_attempts")

# errors and params
Plots.scatter(vals, errors, ylabel=("Error"), xlabel=("Parameter"), legend = false)#,  yaxis=(:log10, (1,Inf)))
scatter!(best_param, best_error, color = "red", label = "", markershape=[:circle])
# png("error_vs_param")


# adujusted errors dont think this is a thing
Plots.plot(range(1,length(optimizer.space.target)), errors_ori, markershapes=[:circle], ylabel=("Adjusted error"), xlabel=("Iteration"), legend=false, xticks = 0:5:length(optimizer.space.target), size=(1000,500))#, yaxis=(:log10, (1,Inf)))
# png("adjusted errors")
# errors and params with original errors by sqrt and /15 (undoing sum?)
Plots.scatter(vals, errors_ori, ylabel=("Adjusted error"), xlabel=("Parameter"), legend = false)#,  yaxis=(:log10, (1,Inf)))
scatter!(best_param, best_error_ori, color = "red", label = "")
# png("adjusted_error_vs_param")




# trying different values of kappa
vals1, errors1, errors_ori1, best_param1, best_error1, best_error_ori1 = [], [], [], [], [], []
function plotting_κ(x, y)

    for i in collect(0:2:10)
        optimizer = bayes_opt.BayesianOptimization(f=rtc_bo, pbounds=py"param_range", random_state=27, verbose=2) 
        optimizer.maximize(init_points=2, n_iter=10, acq="ucb", kappa=i)
        vals, errors, errors_ori, best_param, best_error, best_error_ori = results(optimizer)
        push!(vals1, vals)
        push!(errors1, errors)
        push!(errors_ori1, errors_ori)
        push!(best_param1, best_param)
        push!(best_error1, best_error)
        push!(best_error_ori1, best_error_ori)
    end

    p = Plots.plot(layout=(3,2), size=(1000,800), xlabel="Iterations", ylabel="Error", xticks = 0:1:length(optimizer.space.target))
    for (i, j) in zip(1:6, collect(0:2:10))
            Plots.plot!(p, x, y[i], subplot=i, markershapes=[:circle], labels="κ = $j")
     end
    return display(p)
end
plotting_κ(1:12, errors1)
# png("kappa_error")

vals1, errors1, errors_ori1, best_param1, best_error1, best_error_ori1 = [], [], [], [], [], []
function plotting_κ2(x, y)

    for i in collect(0:2:10)
        optimizer = bayes_opt.BayesianOptimization(f=rtc_bo, pbounds=py"param_range", random_state=27, verbose=2) 
        optimizer.maximize(init_points=2, n_iter=10, acq="ucb", kappa=i)
        vals, errors, errors_ori, best_param, best_error, best_error_ori = results(optimizer)
        push!(vals1, vals)
        push!(errors1, errors)
        push!(errors_ori1, errors_ori)
        push!(best_param1, best_param)
        push!(best_error1, best_error)
        push!(best_error_ori1, best_error_ori)
    end

    p = Plots.plot(layout=(3,2), size=(1000,800), xlabel="Iterations", ylabel="Parameter estimate", xticks = 0:1:length(optimizer.space.target))
    for (i, j) in zip(1:6, collect(0:2:10))
            Plots.plot!(p, x, y[i], subplot=i, markershapes=[:circle], labels="κ = $j")
            hline!([4.14 4.14 4.14 4.14 4.14 4.14], labels= false)
     end
    return display(p)
end
plotting_κ2(1:12, vals1)
# png("kappa_param")

vals1, errors1, errors_ori1, best_param1, best_error1, best_error_ori1 = [], [], [], [], [], []
function plotting_κ1(x, y)

    for i in collect(0:2:10)
        optimizer = bayes_opt.BayesianOptimization(f=rtc_bo, pbounds=py"param_range", random_state=27, verbose=2) 
        optimizer.maximize(init_points=2, n_iter=10, acq="ucb", kappa=i)
        vals, errors, errors_ori, best_param, best_error, best_error_ori = results(optimizer)
        push!(vals1, vals)
        push!(errors1, errors)
        push!(errors_ori1, errors_ori)
        push!(best_param1, best_param)
        push!(best_error1, best_error)
        push!(best_error_ori1, best_error_ori)
    end

    p = Plots.plot(layout=(3,2), size=(1000,800), xlabel="Parameter estimate", ylabel="Error")
    for (i, j) in zip(1:6, collect(0:2:10))
            Plots.scatter!(p, x[i], y[i], subplot=i, markershapes=[:circle], labels="κ = $j", color=palette(:default)[1])
            scatter!(best_param1, best_error1, color = "red", label = "")
    end
    return display(p)
end
plotting_κ1(vals1, errors1)
# png("kappa_error_param")
