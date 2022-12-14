using DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, DataFrames, OrderedCollections, PlotlyJS

include("/home/holliehindley/phd/rtc_models/Oct2022_model/rtc_model.jl")
include("/home/holliehindley/phd/rtc_models/sol_species_funcs.jl")
include("/home/holliehindley/phd/rtc_models/params_init_tspan.jl")

solu = @time(sol(rtc_model, initial, params, tspan))
plotly_plot_sol(solu)

get_curve(solu, :rh)

solu_OD = sol(rtc_model_OD, init_OD, params, tspan)
plotly_plot_sol_OD(solu_OD)


dict_res = get_all_curves(solu, all_species)

Plots.plot(solu.t, [dict_res[:rm_a], dict_res[:rtca]], labels=["rm_a" "RtcA"])
Plots.plot(solu.t, [dict_res[:rm_b], dict_res[:rtcb]], labels=["rm_b" "RtcB"])
Plots.plot(solu.t, [dict_res[:rm_b], dict_res[:rtcb], dict_res[:rm_a], dict_res[:rtca]], labels=["rm_b" "RtcB" "rm_a" "RtcA"])
Plots.plot(solu.t, [dict_res[:rm_r], dict_res[:rtcr]], labels=["rm_r" "RtcR"])
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
results_kdeg = change_param(kdeg_range, "kdeg", rtc_model_OD, init_OD, all_species_OD, param_dict_OD)
plot(kdeg_range, results_kdeg[:rd], Layout(xaxis_title="kdeg", yaxis_title="rd"))


kdam_range = collect(0:0.1:10)
results_kdam = change_param(kdam_range, "kdam", rtc_model_OD, init_OD, all_species_OD, param_dict_OD)
plot(kdam_range, results_kdam[:OD], Layout(xaxis_title="kdam", yaxis_title="OD"))
plot(kdam_range, results_kdam[:rm_a], Layout(xaxis_title="kdam", yaxis_title="rm_a"))

ω_ab_range = collect(1:0.1:2)
results_ωab = change_param(ω_ab_range, "ω_ab", rtc_model_OD, init_OD, all_species_OD, param_dict_OD)
plot(ω_ab_range, results_ωab[:OD], Layout(xaxis_title="ω_ab", yaxis_title="OD"))
plot(ω_ab_range, results_ωab[:rm_a], Layout(xaxis_title="ω_ab", yaxis_title="rm_a"))


ω_r_range = collect(0:1:100)
results_ωr = change_param(ω_r_range, "ω_r", rtc_model_OD, init_OD, all_species_OD, param_dict_OD)
plot(ω_r_range, results_ωr[:OD], Layout(xaxis_title="ω_r", yaxis_title="OD"))
plot(ω_r_range, results_ωr[:rm_a], Layout(xaxis_title="ω_r", yaxis_title="rm_a"))
