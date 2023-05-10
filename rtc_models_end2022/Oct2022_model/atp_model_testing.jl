using CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, PlotlyJS, Statistics
include("/home/holliehindley/phd/rtc_models/Oct2022_model/rtc_model.jl")
include("/home/holliehindley/phd/rtc_models/sol_species_funcs.jl")
include("/home/holliehindley/phd/rtc_models/params_init_tspan.jl")
include("/home/holliehindley/phd/Param_inf/inf_setup.jl")

# run model to steady state 
tspan = (0, 1e9)
params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr)
solu = sol(rtc_model, initial, tspan, params)
# p = plotly_plot_sol(solu, "log", "log")


rm_a_ss = get_ssval(solu, :rm_a)
rm_b_ss = get_ssval(solu, :rm_b)
rm_r_ss = get_ssval(solu, :rm_r)
rtca_ss = get_ssval(solu, :rtca)
rtcb_ss = get_ssval(solu, :rtcb)
rtcr_ss = get_ssval(solu, :rtcr)
rh_ss = get_ssval(solu, :rh)
rd_ss = get_ssval(solu, :rd)
rt_ss = get_ssval(solu, :rt)
atp_ss = 4000

d_atp = 0.005
tspan_atp = (0,1e9)
init_ss = [rm_a_ss, rtca_ss, rm_b_ss, rtcb_ss, rm_r_ss, rtcr_ss, rh_ss, rd_ss, rt_ss, atp_ss]
params_atp = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, na, nb, nr, d_atp] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :na, :nb, :nr, :d_atp)
# solu_atp = sol(rtc_model_atp!, init_ss, tspan, params_atp)


threshold = 255

condition(u,t,integrator) = any( (u[10] .< threshold) .& (u[10].>0))
# condition(u,t,integrator) = u[10] .< threshold


function affect!(integrator)
    integrator.u[10] = 0.000001 #(integrator.u .< threshold) .& (integrator.u .> 0)
    # @show integrator.u[10]
end

cb = DiscreteCallback(condition,affect!)

# d_atp = 0

solu_atp = solcb(rtc_model_atp!, init_ss, tspan, params_atp, cb)


atp = get_curve(solu_atp, :atp)
atp_end = solu_atp.t[findmin(atp)[2]]
p0 = plotly_plot_sol_atp(solu_atp, "log", "", "d_atp = $d_atp", true, atp_end)

# d_atp_range = 10 .^ range(log10(0.0001), log10(1))
# res = change_param_atp(d_atp_range, :d_atp, rtc_model_atp!, init_ss, all_species_atp)
# plot_change_param_sols_atp(d_atp_range, res, "d_atp", "log", "log")


# plots = []

# for i in d_atp_range
#     params_atp[:d_atp] = i
#     prob = ODEProblem(rtc_model_atp!, init_ss, tspan, params_atp, callback=cb)
#     solu = solve(prob, Rodas4())
#     if i == d_atp_range[1]
#         push!(plots, plotly_plot_sol_atp(solu, "log", "", "d_atp = $i", true))
#     else
#         push!(plots, plotly_plot_sol_atp(solu, "log", "", "d_atp = $i", false))
#     end
# end

# p = [p0 plots[1]; plots[2] plots[3]; plots[4] plots[5]]

