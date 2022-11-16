using DifferentialEquations, StaticArrays, BenchmarkTools, DataFrames, Plots, OrderedCollections #, PlotlyJS

include("/home/holliehindley/phd/rtc_models/Oct2022_model/rtc_model.jl")
include("/home/holliehindley/phd/rtc_models/sol_species_funcs.jl")

L = 100; c = 0.01; kr = 10; Vmax_init = 5; Km_init = 55.829; ω_r = 4.14; ω_ab = 4.14;
θtscr = 20; g_max = 4; θtlr = 20; gr_c = 0.01; k_a = 10; k_b = 10; km = 5;
d = 0.01; krep = 10;  kdam = 0.05; ktag = 10; kdeg = 0.001; kin = 0.4; atp = 10;
rtca_0 = 1; rtcb_0 = 1;

params = @SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km, k_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp]

solu = sol(rtc_model, init, tspan, params)

Plots.plot(solu[2:end], ylabel="[species]", labels=["rm_a" "rtca" "rm_b" "rtcb" "rm_r" "rtcr" "rh" "rd" "rt"], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)))

dict_res = get_all_curves(solu, all_species)


Plots.plot(solu.t, [dict_res["rtca"],dict_res["rtcb"],dict_res["rtcr"]], xaxis=(:log10, (0.01,Inf)), labels=["RtcA" "RtcB" "RtcR"])
Plots.plot(solu.t, [dict_res["rm_a"], dict_res["rm_b"], dict_res["rm_r"]], xaxis=(:log10, (0.01,Inf)), labels=["rm_a" "rm_b" "rm_r"])
Plots.plot(solu.t, [dict_res["rh"], dict_res["rt"], dict_res["rd"]], xaxis=(:log10, (0.01,Inf)), labels=["Rt" "Rh" "Rd"])


tlr_el = g_max*atp/(θtlr+atp)
tlr = rh[end]*rm_a[end]*tlr_el

#need changing
r_tot = rh+rt+rd
Plots.plot(solu.t, [(@.rh/r_tot *100), (@.rt/r_tot *100), (@.rd/r_tot *100)], xaxis=(:log10, (0.01,Inf)), labels=["Rt" "Rh" "Rd"])


param_range = collect(0:1:100)
results = change_param(param_range, "ω_ab")
plot(param_range, results[1])





