using DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, DataFrames, OrderedCollections, PlotlyJS, Statistics

include("/home/holliehindley/phd/rtc_models/Oct2022_model/rtc_model.jl")
include("/home/holliehindley/phd/rtc_models/sol_species_funcs.jl")
include("/home/holliehindley/phd/rtc_models/params_init_tspan.jl")

solu = @time(sol(rtc_model, initial, tspan, params))
plotly_plot_sol(solu)

rtca = get_curve(solu, :rtca)
c = rtca[length(rtca)-3:length(rtca)]
std(c)

params = [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr]
all_res = []
ω_ab_range = collect(0:0.1:5)
for i in ω_ab_range
    params[7] = i
    dict_res = OrderedDict(name => [] for name in all_species)
    for val in ω_ab_range  
        params[6] = val
        solu = sol(rtc_model, initial, tspan, params)
        for (i,j) in zip(values(dict_res), all_species)
            push!(i, get_ssval(solu, j))
        end
    end
    push!(all_res, dict_res)

end

rtcas = []
for i in (1:length(ω_ab_range))
    push!(rtcas, all_res[i][:rtca])
end

values(rtcas)
vec = []
for i in (1:length(ω_ab_range))
    append!(vec, values(rtcas[i]))
end
vec = reshape(vec, (length(ω_ab_range),length(ω_ab_range)))


plot(contour(x=ω_ab_range, y=ω_ab_range, z=vec, colorbar=attr(title="RtcA")), Layout(xaxis_title="ω_ab", yaxis_title="ω_r", title="RtcA conc at last time point"))








params = [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr]
all_res = []
ω_ab_range = collect(0:0.1:5)
for i in ω_ab_range
    params[7] = i
    dict_res = OrderedDict(name => [] for name in all_species)
    for val in ω_ab_range  
        params[6] = val
        solu = sol(rtc_model, initial, tspan, params)
        for (i,j) in zip(values(dict_res), all_species)
            push!(i, check_get_ssval(solu, j))
        end
    end
    push!(all_res, dict_res)

end

rtcas = []
for i in (1:length(ω_ab_range))
    push!(rtcas, all_res[i][:rtca])
end

values(rtcas)
vec = []
for i in (1:length(ω_ab_range))
    append!(vec, values(rtcas[i]))
end
vec = reshape(vec, (length(ω_ab_range),length(ω_ab_range)))


plot(contour(x=ω_ab_range, y=ω_ab_range, z=vec, colorbar=attr(title="RtcA")), Layout(xaxis_title="ω_ab", yaxis_title="ω_r", title="RtcA conc at last time point"))#, xaxis_type="log", yaxis_type="log"))












ntspan[2]




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


param_dict = OrderedDict("L"=>L, "c"=>c, "kr"=>kr, "Vmax_init"=>Vmax_init, "Km_init"=>Km_init, "ω_ab"=>ω_ab, "ω_r"=>ω_r, "θtscr"=>θtscr, "g_max"=>g_max, "θtlr"=>θtlr, "km_a"=>km_a, "km_b"=>km_b, "gr_c"=>gr_c, "d"=>d, "krep"=>krep, "kdam"=>kdam, "ktag"=>ktag, "kdeg"=>kdeg, "kin"=>kin, "atp"=>atp, "na"=>na, "nb"=>nb, "nr"=>nr, "k"=>k)

gr_c_range = collect(0:0.001:0.001)
dict_res = OrderedDict(name => [] for name in all_species)
for val in gr_c_range  
    param_dict["gr_c"] = val
    param = values(param_dict)
    solu = sol(rtc_model, initial, tspan, param)
    for (i,j) in zip(values(dict_res), all_species)
        push!(i, get_ssval(solu, j))
    end
end

plot(scatter(x=gr_c_range, y=dict_res[:rtca]))

dict_res[:rtca][1]