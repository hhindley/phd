using DifferentialEquations, StaticArrays, BenchmarkTools, DataFrames, OrderedCollections, Plots #, PlotlyJS

include("/home/holliehindley/phd/rtc_models/Oct2022_model/rtc_model.jl")
include("/home/holliehindley/phd/rtc_models/sol_species_funcs.jl")

prob = ODEProblem(rtc_model_density, init_den, tspan, param_vector[1])
solu1 = solve(prob, Rodas4())#, abstol=1e-15, reltol=1e-12) doesn't solve when run at timespans more than 2100 in length

param_vector = @SVector [values(param_dict)]
plot(solu1[2:end], ylabel="[species]", labels=["rm_a" "rtca" "rm_b" "rtcb" "rm_r" "rtcr" "rh" "rd" "rt" "den"],  xaxis=(:log10, (1,Inf)))
den = get_curve(solu1, :den)
plot(solu1.t, den)
rh = get_curve(solu1, :rh)
plot(solu1.t, rh)
lam = gr_c*rh
plot(solu1.t, lam)


param_vector = @SVector [values(param_dict)]

solu1 = sol(rtc_model_density, init_den, tspan, param_vector[1])

solu = @time(sol(rtc_model, init, tspan, param_vector[1]))
get_curve(solu, :rh)

plot(solu[2:end], ylabel="[species]", labels=["rm_a" "rtca" "rm_b" "rtcb" "rm_r" "rtcr" "rh" "rd" "rt"],  xaxis=(:log10, (1,Inf)))
# savefig(plot(solu[2:end], ylabel="[species]", labels=["rm_a" "rtca" "rm_b" "rtcb" "rm_r" "rtcr" "rh" "rd" "rt"],  xaxis=(:log10, (1,Inf))), "rtc_plot.svg") #, palette=:seaborn_bright)

dict_res = get_all_curves(solu, all_species)

# transcription+translation of rtcb
# alpha = dict_res[:rt][1][end]/kr 
# fa = (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6)
# ra = fa*dict_res[:rtcr][1][end]
# tscr = (ra*Vmax_init*atp/(Km_init+atp))*(ω_ab*atp/(θtscr+atp))
# tlr = (1/408)*dict_res[:rh][1][end]*dict_res[:rm_b][1][end]*(g_max*atp/(θtlr+atp))
# tot = tscr+tlr

Plots.plot(solu.t, [dict_res[:rtca],dict_res[:rtcb],dict_res[:rtcr]], xaxis=(:log10, (0.01,Inf)), labels=["RtcA" "RtcB" "RtcR"])
Plots.plot(solu.t, [dict_res[:rm_a], dict_res[:rm_b], dict_res[:rm_r]], xaxis=(:log10, (0.01,Inf)), labels=["rm_a" "rm_b" "rm_r"])
Plots.plot(solu.t, [dict_res[:rh], dict_res[:rt], dict_res[:rd]], xaxis=(:log10, (0.01,Inf)), labels=["Rh" "Rt" "Rd"])

plot(solu.t, dict_res[:rt])
plot(solu.t, dict_res[:rh])
plot(solu.t, dict_res[:rd])

r_tot = reduce(vcat, (dict_res[:rh]+dict_res[:rd]+dict_res[:rt]))
Plots.plot(solu.t, [((@.dict_res[:rh][1]/r_tot) *100), (@.dict_res[:rt][1]/r_tot *100), (@.dict_res[:rd][1]/r_tot *100)], xaxis=(:log10, (0.01,Inf)), labels=["Rh" "Rt" "Rd"])


# vary param
kdeg_range = collect(0:0.01:0.2)
results_kdeg = change_param(kdeg_range, "kdeg")
plot(kdeg_range, results_kdeg[:rd], legend=false, xlabel="kdeg", ylabel="[Rd]")

kdam_range = collect(0:0.1:10)
results_kdam = change_param(kdam_range, "kdam")
plot(kdam_range, [results_kdam[:rh],results_kdam[:rd]], labels=["Rh" "Rd"], xlabel="kdam", ylabel="[Species]")

