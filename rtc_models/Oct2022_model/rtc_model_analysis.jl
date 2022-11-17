using DifferentialEquations, StaticArrays, BenchmarkTools, DataFrames, Plots, OrderedCollections #, PlotlyJS

include("/home/holliehindley/phd/rtc_models/Oct2022_model/rtc_model.jl")
include("/home/holliehindley/phd/rtc_models/sol_species_funcs.jl")



params = @SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km, k_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr]

solu = sol(rtc_model, init, tspan, params)

Plots.plot(solu[2:end], ylabel="[species]", labels=["rm_a" "rtca" "rm_b" "rtcb" "rm_r" "rtcr" "rh" "rd" "rt"], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)))

dict_res = get_all_curves(solu, all_species)

plot(solu.t, dict_res[:rtca])#, xaxis=(:log10, (0.01,Inf)))



Plots.plot(solu.t, [dict_res[:rtca],dict_res[:rtcb],dict_res[:rtcr]], xaxis=(:log10, (0.01,Inf)), labels=["RtcA" "RtcB" "RtcR"])
Plots.plot(solu.t, [dict_res["rm_a"], dict_res["rm_b"], dict_res["rm_r"]], xaxis=(:log10, (0.01,Inf)), labels=["rm_a" "rm_b" "rm_r"])
Plots.plot(solu.t, [dict_res["rh"], dict_res["rt"], dict_res["rd"]], xaxis=(:log10, (0.01,Inf)), labels=["Rt" "Rh" "Rd"])

plot(solu.t, dict_res[:rt])
plot(solu.t, dict_res[:rh])
plot(solu.t, dict_res[:rd])

r_tot = reduce(vcat, (dict_res[:rh]+dict_res[:rd]+dict_res[:rt]))
Plots.plot(solu.t, [((@.dict_res[:rh][1]/r_tot) *100), (@.dict_res[:rt][1]/r_tot *100), (@.dict_res[:rh][1]/r_tot *100)], xaxis=(:log10, (0.01,Inf)), labels=["Rt" "Rh" "Rd"])



param_range = collect(0:1:1000)
results = change_param(param_range, "kdam")
plot(param_range, results[:rh])


