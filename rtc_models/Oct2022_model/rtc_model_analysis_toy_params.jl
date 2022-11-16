using DifferentialEquations, StaticArrays, BenchmarkTools, DataFrames, Plots, OrderedCollections #, PlotlyJS

include("/home/holliehindley/phd/rtc_models/Oct2022_model/rtc_model.jl")
include("/home/holliehindley/phd/rtc_models/sol_species_funcs.jl")

L = 100; c = 0.01; kr = 10; Vmax_init = 5; Km_init = 55.829; ω_r = 4.14; ω_ab = 4.14;
θtscr = 20; g_max = 4; θtlr = 20; gr_c = 0.01; k_a = 10; k_b = 10; km = 5;
d = 0.01; krep = 10;  kdam = 0.05; ktag = 10; kdeg = 0.001; kin = 0.4; atp = 10;
rtca_0 = 1; rtcb_0 = 1;

solu = sol(rtc_model, init, tspan, params)
Plots.plot(solu[2:end], ylabel="[species]", labels=["rm_a" "rtca" "rm_b" "rtcb" "rm_r" "rtcr" "rh" "rd" "rt"], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)))

for (i,j) in zip(results, species)
    push!(i, get_curve(solu, j))
end

dict_res = OrderedDict("rm_a"=>[], "rtca"=>[], "rm_b"=>[], "rtcb"=>[], "rm_r"=>[], "rtcr"=>[], "rh"=>[], "rd"=>[], "rt"=>[])

results[1]

Plots.plot(solu.t, [results[2],results[4],results[6]], xaxis=(:log10, (0.01,Inf)), labels=["RtcA" "RtcB" "RtcR"])
Plots.plot(solu.t, [rm_a, rm_b, rm_r], xaxis=(:log10, (0.01,Inf)), labels=["rm_a" "rm_b" "rm_r"])
Plots.plot(solu.t, [rt,rh,rd], xaxis=(:log10, (0.01,Inf)), labels=["Rt" "Rh" "Rd"])




rm_a = get_curve(solu, :rm_a); rm_b = get_curve(solu, :rm_b); rm_r = get_curve(solu, :rm_r)
rtca = get_curve(solu, :rtca); rtcb = get_curve(solu, :rtcb); rtcr = get_curve(solu, :rtcr)
rh = get_curve(solu, :rh); rd = get_curve(solu, :rd); rt = get_curve(solu, :rt)

tlr_el = g_max*atp/(θtlr+atp)
tlr = rh[end]*rm_a[end]*tlr_el

Plots.plot(solu.t, [rtca,rtcb,rtcr], xaxis=(:log10, (0.01,Inf)), labels=["RtcA" "RtcB" "RtcR"])
Plots.plot(solu.t, [rm_a, rm_b, rm_r], xaxis=(:log10, (0.01,Inf)), labels=["rm_a" "rm_b" "rm_r"])
Plots.plot(solu.t, [rt,rh,rd], xaxis=(:log10, (0.01,Inf)), labels=["Rt" "Rh" "Rd"])

r_tot = rh+rt+rd
Plots.plot(solu.t, [(@.rh/r_tot *100), (@.rt/r_tot *100), (@.rd/r_tot *100)], xaxis=(:log10, (0.01,Inf)), labels=["Rt" "Rh" "Rd"])


param_range = collect(0:1:100)
results = change_param(param_range, "ω_ab")
plot(param_range, results[1])





