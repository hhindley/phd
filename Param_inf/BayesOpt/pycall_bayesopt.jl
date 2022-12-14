using Plots, PyCall
include("/home/holliehindley/rtc_models_final_sept2022/julia_model.jl")

# set time span and how many time points to solve at 
tspan = (0, 100)
t = exp10.(range(-3,2,15))
pushfirst!(t, 0)

# get curve of species that we want to compare data/model to 
function get_curve(sol, species)
    df = DataFrame(sol)
    rename!(df, [:time, :rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt])
    species = df[:, species]
    return species
end

function get_rdrtca_rtrtcb(sol)
    df = DataFrame(sol)
    rename!(df, [:time, :rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt])
    rtca = df[:, :rtca]
    rtcb = df[:, :rtcb]
    rd = df[:, :rd]
    rt = df[:, :rt]
    rdrtca = rtca.*rd./(K1_tag.+rd.+K2_tag.*atp)
    rtrtcb = rtcb.*rt./(K1_rep.+rt.+K2_rep.*atp)
    return rdrtca, rtrtcb
end

# solving function
function sol(model, init, tspan, params, t)
    prob = ODEProblem(model, init, tspan, params)
    sol = solve(prob, Rodas4(), saveat=t)
    return sol
end

sol_syn = sol(rtc_model, init, tspan, params, t)

# plot solution 
Plots.plot(sol_syn[2:end], xaxis=(:log10, (1,Inf)), yaxis=(:log10, (1,Inf)))

rtca_syn = get_curve(sol_syn, :rtca)
rt_syn = get_curve(sol_syn, :rt)
rd_syn = get_curve(sol_syn, :rd)
rh_syn = get_curve(sol_syn, :rh)

rdrtca_syn, rtrtcb_syn = get_rdrtca_rtrtcb(sol_syn)

# objective function
function rtc_bo(;ω_ab,)#K1_tag, K2_tag, K1_rep, K2_rep, ω_ab)
    solu = sol(rtc_model, init, tspan, (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θ, max, thr, K1_tag, K2_tag, K1_rep, K2_rep, gr_c, d, krep, kdam, ktag, kdeg, kin, atp]), t)
    rtca = get_curve(solu, :rtca)
    # rm_a = get_curve(solu, :rm_a)
    obj = []
    for (i,j) in zip(rtca, rtca_syn)
        append!(obj, abs2(i-j))
    end
    return -sum(obj)
end

# in python writing the ranges to search for parameter value
py"""
param_range = {'ω_ab': (0, 50)}#,'K1_tag': (0.9, 1.1), 'K2_tag': (9.9, 10.1), 'K1_rep': (0.9, 1.1), 'K2_rep': (9.9, 10.1)}
"""

# import bayes_opt package from python
bayes_opt = pyimport("bayes_opt")

# setting the optimizer 
optimizer = bayes_opt.BayesianOptimization(f=rtc_bo, pbounds=py"param_range", random_state=27, verbose=2) # verbose = 1 prints only when a maximum is observed (pink)

# timing the process and maximising the optimizer 
function timer()
    optimizer.maximize(init_points=2, n_iter=30, acq="ucb", kappa=3)
end

@time timer()
print(optimizer_all.max)

# creating lists of values of parameters tried and errors for each one  
function results(optimizer)
    vals, errors = [], []
    for i in collect(1:length(optimizer.space.target))
        a = collect(values(optimizer.res)[i])
        append!(vals, collect(values(a[2][2])))
        append!(errors, -(collect(values(a[1][2]))))
    end

    # reversing sum of squares error calculation to get errors close to zero 
    errors_ori = sqrt.(errors/15)

    # best values 
    best_param = collect(values(collect(values(optimizer.max))[2]))
    best_error = -[collect(values(optimizer.max))[1]]

    # best error with sum of squares reversed 
    best_error_ori = sqrt.(best_error/15)
    return vals, errors, errors_ori, best_param, best_error, best_error_ori
end

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
