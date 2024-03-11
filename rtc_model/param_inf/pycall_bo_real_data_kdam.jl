using Plots, PyCall, DifferentialEquations, StaticArrays, BenchmarkTools, DataFrames, CSV, OrderedCollections
include("$PATHrtc_models/Oct2022_model/rtc_model.jl")
include("$PATHrtc_models/sol_species_funcs.jl")
include("$PATHParam_inf/inf_setup.jl")
include("$PATHrtc_models/params_init_tspan.jl")

# infer kdam for hpx conditions 
function rtc_bo_hpx(;kdam)
    obj_hpx_ab, obj_hpx_r = obj(rtc_model, initial, (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr]), tspan2, t_2, "mrna", hpx, hpx_std)
    
    obj_hpx_rtcon = obj_OD(rtc_model_OD, (@SVector [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, rh_0, rd_0, rt_0, OD_0_hpxon]), (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, findmax(hpx_rtcon)[1]]), tspan2, t_2, "OD", hpx_rtcon, hpx_rtcon_std)

    return -(obj_hpx_ab/36 + obj_hpx_r + obj_hpx_rtcon/36)
end

# infer kdam for colD conditions 
function rtc_bo_colD(;kdam)
    obj_wt_colD = obj_OD(rtc_model_OD, (@SVector [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, rh_0, rd_0, rt_0, OD_0_wtcolD]), (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, findmax(WT_colD)[1]]), tspan4, t_4, "OD", WT_colD, WT_colD_std)
    obj_nA_colD = obj_OD(rtc_model_knockout_OD, (@SVector [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, rh_0, rd_0, rt_0, OD_0_nAcolD]), (@SVector [L, c, kr, Vmax_init, Km_init, 0, ω_a, ω_r, θtscr, g_max, θtlr, km_a, km_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, findmax(nA_colD)[1]]), tspan4, t_4, "OD", nA_colD, nA_colD_std)
    obj_nB_colD = obj_OD(rtc_model_knockout_OD, (@SVector [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, rh_0, rd_0, rt_0, OD_0_nB_colD]), (@SVector [L, c, kr, Vmax_init, Km_init, ω_a, 0, ω_r, θtscr, g_max, θtlr, km_a, km_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, findmax(nB_colD)[1]]), tspan4, t_4, "OD", nB_colD, nB_colD_std)
    obj_nBB_colD = obj_OD(rtc_model_OD, (@SVector [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, rh_0, rd_0, rt_0, OD_0_nBBcolD]), (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, findmax(nB_B_colD)[1]]), tspan4, t_4, "OD", nB_B_colD, nB_B_colD_std)
    obj_nBBmut_colD = obj_OD(rtc_model_knockout_OD, (@SVector [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, rh_0, rd_0, rt_0, OD_0_bBBmutcolD]), (@SVector [L, c, kr, Vmax_init, Km_init, ω_a, 0, ω_r, θtscr, g_max, θtlr, km_a, km_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, findmax(nB_Bmut_colD)[1]]), tspan4, t_4, "OD", nB_Bmut_colD, nB_Bmut_colD_std)
    obj_nR_colD = obj_OD(rtc_model_OD, (@SVector [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, rh_0, rd_0, rt_0, OD_0_nRcolD]), (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, 0, θtscr, g_max, θtlr, km_a, km_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, findmax(nR_colD)[1]]), tspan4, t_4, "OD", nR_colD, nR_colD_std)
    return -(obj_wt_colD/(66^2) + obj_nA_colD/(66^2) + obj_nB_colD/(66^2) + obj_nBB_colD/(66^2) + obj_nBBmut_colD/(66^2) + obj_nR_colD/(66^2))
end

# in python writing the ranges to search for parameter value
py"""
param_range_dam = {'kdam': (0, 200)}
"""
# import bayes_opt package from python
bayes_opt = pyimport("bayes_opt")


function run_bo(func)
# setting the optimizer 
    optimizer = bayes_opt.BayesianOptimization(f=func, pbounds=py"param_range_dam", random_state=27, verbose=2) # verbose = 1 prints only when a maximum is observed (pink)
    function timer()
        optimizer.maximize(init_points=2, n_iter=10, acq="ucb", kappa=2)#, xi=0.0)
    end
    return @time timer()
end

hpx = run_bo(rtc_bo_hpx)

colD = run_bo(rtc_bo_colD)


print(optimizer.max)



# creating lists of values of parameters tried and errors for each one  

vals, errors, errors_ori, best_param, best_error, best_error_ori = results(optimizer)

# plot errors over iterations 
Plots.plot(range(1,length(optimizer.space.target)), errors, markershapes=[:circle], ylabel=("Error"), xlabel=("Iteration"), legend=false, xticks = 0:5:length(optimizer.space.target), size=(1000,500))#, yaxis=(:log10, (1,Inf)))
# png("error")
Plots.plot(range(1,length(optimizer.space.target)), errors_ori, markershapes=[:circle], ylabel=("Adjusted error"), xlabel=("Iteration"), legend=false, xticks = 0:5:length(optimizer.space.target), size=(1000,500))#, yaxis=(:log10, (1,Inf)))
# png("adjusted errors")

# plot params over iterations
Plots.plot(range(1,length(optimizer.space.target)), vals, markershapes=[:circle], ylabel=("Param attempt"), xlabel=("Iteration"), legend=false, xticks = 0:1:length(optimizer.space.target))#, size=(1000,500))#, yaxis=(:log10, (1,Inf)))
hline!([4.14], labels= false)
# png("param_attempts")

# errors and params
Plots.scatter(vals, errors, ylabel=("Error"), xlabel=("Parameter"), legend = false)#,  yaxis=(:log10, (1,Inf)))
scatter!(best_param, best_error, color = "red", label = "", markershape=[:circle])
# png("error_vs_param")

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
