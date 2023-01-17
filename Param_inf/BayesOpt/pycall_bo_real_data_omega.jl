using PlotlyJS, PyCall, DifferentialEquations, StaticArrays, DataInterpolations, LabelledArrays, BenchmarkTools, DataFrames, CSV, OrderedCollections, ScatteredInterpolation
include("/home/holliehindley/phd/rtc_models/Oct2022_model/rtc_model.jl")
include("/home/holliehindley/phd/rtc_models/sol_species_funcs.jl")
include("/home/holliehindley/phd/Param_inf/inf_setup.jl")
include("/home/holliehindley/phd/rtc_models/params_init_tspan.jl")
# include("/home/holliehindley/phd/rtc_models/Oct2022_model/units/millerunit_conversion.jl")
csv_mRNA_conc = DataFrame(CSV.File("/home/holliehindley/phd/data/mRNA_conc_uM.csv"))
mRNA_conc = csv_mRNA_conc[:,1]
mRNA_std = csv_mRNA_conc[:,2]

csv_wt = DataFrame(CSV.File("/home/holliehindley/phd/data/results_rtcOFF_grfit.csv"))
csv_wt = select!(csv_wt, Not(["log(OD)", "log(OD) error", "gr error", "od"]))
lam_wt, new_df_wt = extend_gr_curve(csv_wt)

# solu = sol_with_t(rtc_model1!, initial, [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, atp, na, nb, nr, lam_wt], tspan2, t_2);
# plotly_plot_sol(solu, "")

function rtc_bo_ω(;ω_ab, ω_r)
    # @show ω_ab, ω_r
    obj_wt_ab = obj(rtc_model1!, initial, ([L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, atp, na, nb, nr, lam_wt]), tspan2, t_2, "mrna", mRNA_conc, mRNA_std)
    return -(obj_wt_ab/36)
end


# in python writing the ranges to search for parameter value
py"""
param_range_ω = {'ω_ab': (0, 0.1), 'ω_r': (0, 0.1)}
"""


# import bayes_opt package from python
bayes_opt = pyimport("bayes_opt")

# setting the optimizer 
optimizer = bayes_opt.BayesianOptimization(f=rtc_bo_ω, pbounds=py"param_range_ω", random_state=27, verbose=2) # verbose = 1 prints only when a maximum is observed (pink)

# timing the process and maximising the optimizer 
function timer()
    optimizer.maximize(init_points=2, n_iter=100, acq="ei", xi=0.01) #kappa=2, xi = 0.0 (prefer exploitation), xi = 0.1 (prefer exploration)
end

@time timer()
print(optimizer.max)

function posterior(optimizer, x_obs, y_obs, grid)
    optimizer._gp.fit(x_obs, y_obs)
    mu, sigma = optimizer._gp.predict(grid, return_std=true)
    return mu, sigma
end

function plot_gp(optimizer, x, y)
    p = make_subplots(rows=2, cols=1, shared_xaxes=true, vertical_spacing=0.02)
    x_obs = Array([[res["params"]["ω_ab"]] for res in optimizer.res])
    y_obs = Array([[res["target"]] for res in optimizer.res])
    mu, sigma = posterior(optimizer, x_obs, y_obs, x)

    steps = length(optimizer.space)
    add_trace!(p, scatter(x=vec(x), y=y, name="Target"), row=1, col=1)
    add_trace!(p, scatter(x=reduce(vcat, x_obs), y=reduce(vcat, y_obs), line=attr(width=0), name="Observation", marker=attr(color="red")), row=1, col=1)
    add_trace!(p, scatter(x=vec(x), y=mu, line=attr(color="green"), name="Prediction"), row=1, col=1)
    
    add_trace!(p, scatter(x=vec(x), y=(@. mu+sigma), line=attr(width=0), showlegend=false), row=1, col=1)
    add_trace!(p, scatter(x=vec(x), y=(@. mu-sigma), line=attr(width=0), fill="tonexty", showlegend=true, name="95% confidence interval"), row=1, col=1)
    
    utility_function = bayes_opt.UtilityFunction(kind="ucb", kappa=5, xi=0)
    utility = utility_function.utility(x, optimizer._gp, 0)

    add_trace!(p, scatter(x=vec(x), y=utility, name="Utility function"), row=2, col=1)
    add_trace!(p, scatter(x=Float64[x[argmax(utility)]], y=Float64[maximum(utility)], marker=attr(symbol="star", color="yellow", size=8, line=attr(color="black", width=0.5)), name="Next best guess"), row=2, col=1)
    relayout!(p, title_text=("Gaussian Process and Utility Function After $steps Steps"))
    return p
    # return x_obs, y_obs
end


x = reshape(csv_wt."t", (6,1))
y = mRNA_conc
# plot(scatter(x=x, y=y))
plot_gp(optimizer, x, y)

Iterators.flatten(x_obs)
plot(scatter(x=x, y=y))

# creating lists of values of parameters tried and errors for each one  
vals_ab, vals_r, errors, best_param_ab, best_param_r, best_error = results_two_param(optimizer)

# plot errors over iterations 
plot(range(1,length(errors)), errors, mode="markers+lines", Layout(xaxis_title="Iteration", yaxis_title="Error", yaxis_range=[0, 100]))#, yaxis=(:log10, (1,Inf)))

plot(range(10,length(errors)), errors[10:end], mode="markers+lines", Layout(xaxis_title="Iteration", yaxis_title="Error"))#, yaxis=(:log10, (1,Inf)))
plot(range(1,length(errors)), errors, mode="markers+lines", Layout(xaxis_title="Iteration", yaxis_title="Error"))#, yaxis=(:log10, (1,Inf)))
# png("error")

layout1 = Layout(xaxis_title="ω_ab", yaxis_title="ω_r", zaxis_title="Error")
# plot(contour(z=errors, x=vals_ab, y=vals_r, layout1))





# 3D plot
x,y,z = plot_3D(vals_ab, vals_r, errors)

layout = Layout(scene=attr(xaxis_title="ω_ab", yaxis_title="ω_r", zaxis_title="Error"))
plot(surface(x=x,y=y, z=z), layout)
plot(contour(z=z, x=x, y=y, contours_start=0, contours_end=50))

# plotting how well the data fits 
data_plot = scatter(x=dfc[!,1]*60, y=WT1, name="data")
plot(data_plot)
solu_t = sol_with_t(rtc_model, initial, params, tspan2, t_2)
solu = sol(rtc_model, initial,params, tspan2)

plotly_plot_sol(solu_t)
plotly_plot_sol(solu)
plotly_plot_sol_withdata(solu_t)
plotly_plot_sol_withdata(solu)

rm_a = get_curve(solu_t, :rm_a); rm_r = get_curve(solu_t, :rm_r); 
rma_curve = scatter(x=solu_t.t, y=rm_a, name="mRNA RtcA")
rmr_curve = scatter(x=solu_t.t, y=rm_r, name="mRNA RtcR")
plot([data_plot, rma_curve, rmr_curve])

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
