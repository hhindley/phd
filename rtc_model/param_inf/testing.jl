using PlotlyJS, PyCall, DifferentialEquations, LabelledArrays, DataInterpolations, StaticArrays, BenchmarkTools, DataFrames, CSV, OrderedCollections, ScatteredInterpolation
include("$PATHrtc_models/Oct2022_model/rtc_model.jl")
include("$PATHrtc_models/sol_species_funcs.jl")
include("$PATHParam_inf/inf_setup.jl")
include("$PATHrtc_models/params_init_tspan.jl")


csv_wt = DataFrame(CSV.File("$PATHdata/results_rtcOFF_grfit.csv"))
csv_wt = select!(csv_wt, Not(["log(OD)", "log(OD) error", "gr error", "od"]))
lam_wt, new_df_wt = extend_gr_curve(csv_wt)

function rtc_bo_ω1(;ω_ab)
    obj_wt_ab, obj_wt_r = obj(rtc_model1!, initial, ([L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, atp, na, nb, nr, lam_wt]), tspan2, t_2, "mrna", WT1, WT1_std)
    return -(obj_wt_ab/36)
end

# in python writing the ranges to search for parameter value
py"""
param_range_ω = {'ω_ab': (0, 1)}
"""

# import bayes_opt package from python
bayes_opt = pyimport("bayes_opt")

# setting the optimizer 
optimizer = bayes_opt.BayesianOptimization(f=rtc_bo_ω1, pbounds=py"param_range_ω", random_state=27, verbose=2) # verbose = 1 prints only when a maximum is observed (pink)

# timing the process and maximising the optimizer 
function timer()
    optimizer.maximize(init_points=2, n_iter=0, acq="ucb")#, xi=0.1) #kappa=2, 
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
y = WT1
# plot(scatter(x=vec(x), y=y))
plot_gp(optimizer, x, y)

function more_iter(optimizer, x, y)
    optimizer.maximize(init_points=0, n_iter=1, acq="ucb")
    return plot_gp(optimizer,x,y)
end

more_iter(optimizer, x, y)





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



















function target(x)
    return @. exp(-(x-2)^2)+exp(-(x-6)^2/10)+1/(x^2+1)
end
function target1(;x)
    return @. exp(-(x-2)^2)+exp(-(x-6)^2/10)+1/(x^2+1)
end

x = collect(range(-2, 10, 10000))
y = target(x) 

# in python writing the ranges to search for parameter value
py"""
param_range = {'x': (-2, 10)}
"""

# import bayes_opt package from python
bayes_opt = pyimport("bayes_opt")

# setting the optimizer 
optimizer = bayes_opt.BayesianOptimization(f=target1, pbounds=py"param_range", random_state=27, verbose=2) # verbose = 1 prints only when a maximum is observed (pink)

optimizer.maximize(init_points=2, n_iter=0, kappa=5)


function posterior(optimizer, x_obs, y_obs, grid)
    optimizer._gp.fit(x_obs, y_obs)
    mu, sigma = optimizer._gp.predict(grid, return_std=true)
    return mu, sigma
end

function plot_gp(optimizer, x, y)
    p = make_subplots(rows=2, cols=1, shared_xaxes=true, vertical_spacing=0.02)
    x_obs = Array([[res["params"]["x"]] for res in optimizer.res])
    y_obs = Array([[res["target"]] for res in optimizer.res])
    mu, sigma = posterior(optimizer, x_obs, y_obs, x)

    steps = length(optimizer.space)
    add_trace!(p, scatter(x=vec(x), y=y, name="Target"), row=1, col=1)
    add_trace!(p, scatter(x=reduce(vcat, x_obs), y=reduce(vcat, y_obs), line=attr(width=0), name="Observation", marker=attr(color="red")), row=1, col=1)
    # add_trace!(p, scatter(x=((x_obs)[2]), y=y_obs[2], name="Observation", marker=attr(color="red"), showlegend=false), row=1, col=1)
    add_trace!(p, scatter(x=vec(x), y=mu, line=attr(color="green"), name="Prediction"), row=1, col=1)
    
    add_trace!(p, scatter(x=vec(x), y=(@. mu+1.96*sigma), line=attr(width=0), showlegend=false), row=1, col=1)
    add_trace!(p, scatter(x=vec(x), y=(@. mu-1.96*sigma), line=attr(width=0), fill="tonexty", showlegend=true, name="95% confidence interval"), row=1, col=1)
    
    utility_function = bayes_opt.UtilityFunction(kind="ucb", kappa=5, xi=0)
    utility = utility_function.utility(x, optimizer._gp, 0)

    add_trace!(p, scatter(x=vec(x), y=utility, name="Utility function"), row=2, col=1)
    add_trace!(p, scatter(x=Float64[x[argmax(utility)]], y=Float64[maximum(utility)], marker=attr(symbol="star", color="yellow", size=8, line=attr(color="black", width=0.5)), name="Next best guess"), row=2, col=1)
    relayout!(p, title_text=("Gaussian Process and Utility Function After $steps Steps"))
    return p
    # return x_obs, y_obs
end

x = reshape(x, (10000,1))
# y = WT1
# plot(scatter(x=vec(x), y=y))
plot_gp(optimizer, x, y)




optimizer.maximize(init_points=0, n_iter=1, kappa=5)#, xi=0.1) #kappa=2, 

plot_gp(optimizer, x, y)
# steps = length(optimizer.space)
# mu, utility, x_obs, y_obs = plot_gp(optimizer, x, y)
# p = make_subplots(rows=2, cols=1, shared_xaxes=true, vertical_spacing=0.02)
# add_trace!(p, scatter(x=vec(x), y=y), row=1, col=1)
# add_trace!(p, scatter(x=((x_obs)[1]), y=y_obs[1]), row=1, col=1)
# add_trace!(p, scatter(x=((x_obs)[2]), y=y_obs[2]), row=1, col=1)
# add_trace!(p, scatter(x=vec(x), y=mu), row=1, col=1)
# p

# add_trace!(p, scatter(x=vec(x), y=utility), row=2, col=1)
# p
# add_trace!(p, scatter(x=Float64[x[argmax(utility)]], y=Float64[maximum(utility)]), row=2, col=1)
# p

# relayout!(p, title_text=("Gaussian Process and Utility Function After $steps Steps"))
# p
# relayout!(p, xaxis_range[-2,10])
# Layout(xaxis_range=[-2, 10], yaxis_range=[0, maximum(utility)+0.5])


maximum(utility)

x[3446]
