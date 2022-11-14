using BayesOpt, DifferentialEquations, StaticArrays, Plots, BenchmarkTools

include("/home/holliehindley/rtc_models_final_sept2022/julia_model.jl")

config = ConfigParameters()         # initiates parameters 

config.n_iterations = 100
config.n_iter_relearn = 50
config.n_init_samples = 500     # how to change parameter values from original
config.noise = 1e-14

prob_syn = ODEProblem(rtc_model, init, tspan, params)
sol_syn = solve(prob_syn, Rodas4())

df = DataFrame(sol_syn)
rename!(df, [:time, :rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt])

rtca_syn = df[:, :rtca]

# obj = []
# for (i,j) in zip(rtca, rtca)
#     append!(obj, abs2(i-j))
# end
# obj
# sum(obj)

function rtc_bo(x)
    prob = ODEProblem(rtc_model, init, tspan, [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θ, max, thr, K1_tag, K2_tag, K1_rep, K2_rep, gr_c, d, krep, kdam, ktag, kdeg, kin, atp])
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

# x = range(-2,2, length=101)
# plot(x,@. sum(x^2))

lb =[0.9, 9.9, 0.9, 9.9, 4.1]; 
ub =[1.5, 10.5, 1.5, 10.5, 4.3]

optimizer, optimum = bayes_optimization(rtc_bo, lb, ub, config)