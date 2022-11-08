using Hyperopt, Plots
include("/home/holliehindley/rtc_models_final_sept2022/julia_model.jl")

prob_syn = ODEProblem(rtc_model, init, tspan, params)
sol_syn = solve(prob_syn, Rodas4())

df_syn = DataFrame(sol_syn)
rename!(df_syn, [:time, :rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt])
rtca_syn = df_syn[:, :rtca]

function rtc_bo(K1_tag, K2_tag, K1_rep, K2_rep, ω_ab)
    prob = ODEProblem(rtc_model, init, tspan, @SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θ, max, thr, K1_tag, K2_tag, K1_rep, K2_rep, gr_c, d, krep, kdam, ktag, kdeg, kin, atp])
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


ho = @hyperopt for i=1000,
    sampler = BOHB(), # This is default if none provided
    K1_tag = LinRange(0,4,1000),
    K2_tag = LinRange(8,15,1000),
    K1_rep = LinRange(0,4, 1000),
    K2_rep = LinRange(8,15, 1000),
    ω_ab = LinRange(1,7, 1000)
print(i, "\t", K1_tag, "\t", K2_tag, "\t", K1_rep, "   \t", K2_rep, "\t", ω_ab, "\t")
@show  rtc_bo(K1_tag, K2_tag, K1_rep, K2_rep, ω_ab)
end

Plots.plot(ho, size=(1000,800))

print(ho.history)
print(ho.results)
print(ho.candidates)
print(ho.params)
print(ho.minimum) # minumum of function rtc_bo

print(ho.minimizer) # best params

printmin(ho)