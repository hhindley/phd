using BayesianOptimization, GaussianProcesses, Distributions
include("/home/holliehindley/rtc_models_final_sept2022/julia_model.jl")

prob_syn = ODEProblem(rtc_model, init, tspan, params)
sol_syn = solve(prob_syn, Rodas4())

df_syn = DataFrame(sol_syn)
rename!(df_syn, [:time, :rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt])
rtca_syn = df_syn[:, :rtca]


trace1 = scatter(df, x=:time, y=:rm_a, name="rm_a")
trace2 = scatter(df, x=:time, y=:rtca, name="rtca")
trace3 = scatter(df, x=:time, y=:rm_b, name="rm_b")
trace4 = scatter(df, x=:time, y=:rtcb, name="rtcb")
trace5 = scatter(df, x=:time, y=:rm_r, name="rm_r")
trace6 = scatter(df, x=:time, y=:rtcr, name="rtcr")
trace7 = scatter(df, x=:time, y=:rh, name="rh")
trace8 = scatter(df, x=:time, y=:rd, name="rd")
trace9 = scatter(df, x=:time, y=:rt, name="rt")
l = Layout(xaxis_type="log", yaxis_type="log", yaxis_range=[0,5])
plot([trace1, trace2, trace3, trace4, trace5, trace6, trace7, trace8, trace9], l)



function rtc_bo(ω_ab)#K1_tag, K2_tag, K1_rep, K2_rep, ω_ab)
    prob = ODEProblem(rtc_model, init, tspan, params)
    sol = solve(prob, Rodas4())
    df = DataFrame(sol)
    rename!(df, [:time, :rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt])
    rtca = df[:, :rtca]
    obj = []
    for (i,j) in zip(rtca, rtca_syn)
        append!(obj, abs2(i-j))
    end
    return sum(obj)
end

# Choose as a model an elastic GP with input dimensions 2.
# The GP is called elastic, because data can be appended efficiently.
model = ElasticGPE(1,                            # 2 input dimensions
                #    mean = MeanConst(0.),         
                #    kernel = SEArd([0.], 5.),
                #    logNoise = 0.,
                   capacity = 3000)              # the initial capacity of the GP is 3000 samples.
# set_priors!(model.mean, [Normal(1, 2)])

# Optimize the hyperparameters of the GP using maximum a posteriori (MAP) estimates every 50 steps
modeloptimizer = MAPGPOptimizer(every = 50, #noisebounds = [-4, 3],       # bounds of the logNoise
                                kernbounds = [[-1], [10]],  # bounds of the 3 parameters GaussianProcesses.get_param_names(model.kernel)
                                maxeval = 40)
opt = BOpt(rtc_bo,
           model,
           UpperConfidenceBound(),                   # type of acquisition
           modeloptimizer,                        
           [1], [5],                     # lowerbounds, upperbounds         
           repetitions = 1,                          # evaluate the function for each input 5 times
           maxiterations = 100,                      # evaluate at 100 input positions
           sense = Min,                              # minimize the function
           acquisitionoptions = (method = :LD_LBFGS, # run optimization of acquisition function with NLopts :LD_LBFGS method
                                 restarts = 5,       # run the NLopt method from 5 random initial conditions each time.
                                 maxtime = 0.1,      # run the NLopt method for at most 0.1 second each time
                                 maxeval = 1000),    # run the NLopt methods for at most 1000 iterations (for other options see https://github.com/JuliaOpt/NLopt.jl)
            verbosity = Progress)

result = boptimize!(opt)











