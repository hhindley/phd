using CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, PlotlyJS, Statistics
include("/home/holliehindley/phd/rtc_models/Oct2022_model/rtc_model.jl")
include("/home/holliehindley/phd/rtc_models/sol_species_funcs.jl")
include("/home/holliehindley/phd/rtc_models/params_init_tspan.jl")
include("/home/holliehindley/phd/Param_inf/inf_setup.jl")

# run model to steady state 
tspan = (0, 1e9)
params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr)
solu = sol(rtc_model, initial, tspan, params)
p = plotly_plot_sol(solu, "log", "log")


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

d_atp = 1
tspan_atp = (0,1e9)
init_ss = [rm_a_ss, rtca_ss, rm_b_ss, rtcb_ss, rm_r_ss, rtcr_ss, rh_ss, rd_ss, rt_ss, atp_ss]
params_atp = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, na, nb, nr, d_atp] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :na, :nb, :nr, :d_atp)
solu_atp = sol(rtc_model_atp!, init_ss, tspan, params_atp)

p = plotly_plot_sol_atp(solu_atp, "log", "log")


threshold = 0.1
condition(u,t,integrator) = any( (u[10] .< threshold) .& (u[10].>0))
function affect!(integrator)
    integrator.u[10] = 0 #(integrator.u .< threshold) .& (integrator.u .> 0)
    # @show integrator.u[10]
end

cb = DiscreteCallback(condition,affect!)
solu_atp = solcb(rtc_model_atp!, init_ss, tspan, params_atp, cb)

p = plotly_plot_sol_atp(solu_atp, "log", "log")
