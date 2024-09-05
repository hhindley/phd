using StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed, TableOperations, JSON, Query, FindFirstFunctions, CategoricalArrays
using Parameters, LabelledArrays, BenchmarkTools
using Revise, LinearAlgebra, Printf, ModelingToolkit, OrderedCollections


using InteractiveViz, GLMakie, Parameters, CSV, DataFrames, DifferentialEquations, LabelledArrays, BenchmarkTools, Revise, LinearAlgebra, Printf, ModelingToolkit, OrderedCollections, Colors

include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/switching_funcs.jl"))

include(joinpath(homedir(), "phd/general_funcs/all_model_funcs.jl"))
include(joinpath(homedir(), "phd/general_funcs/solving.jl"))
include(joinpath(homedir(), "phd/rtc_model/parameters/rtc_params.jl"))
include(joinpath(homedir(), "phd/rtc_model/parameters/rtc_params_molecs.jl"))

include(joinpath(homedir(), "phd/rtc_model/models/rtc_orig.jl"))
include(joinpath(homedir(), "phd/rtc_model/functions/bf_funcs/bf_funcs.jl"))

mount_path, folders, folders_dict = load_file_structure("kdam_testing/keyvals2")
folders_dict = Dict(filter(pair -> pair.first in [6,7,8,9,10,11,12,13,14,15], folders_dict))
# folders_dict = Dict(filter(pair -> pair.first in [6], folders_dict))
dict_times, dict_kdamvals, dict_titles, dict_results, dict_reacts, dict_props, dict_counts, dict_hists = load_data(mount_path, folders, folders_dict, reacts=false, props=false, hists=false)

# plot one result
folder = 6; index = 1; species = "rtca"; num_plots = 1;
f = plot_results("plot_results", dict_results[folder][index], num_plots, folders_dict[folder], titles=[dict_titles[folder][index]], species="$species", xlabel="time", ylabel="$species", size=(800,650), tosave=false, conc=true)
f_rib = plot_results("plot_results", dict_results[folder][index], 1, folders_dict[folder], titles=[dict_titles[folder][index]], species=[:rh, :rd, :rt], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false, conc=false)
f_mrna = plot_results("plot_results", dict_results[folder][index], 1, folders_dict[folder], species=[:rm_a, :rm_b, :rm_r], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false, conc=true)
f_prot = plot_results("plot_results", dict_results[folder][index], 1, folders_dict[folder], titles=[dict_titles[folder][index]], species=[:rtca, :rtcb, :rtcr], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false, conc=false)

display(GLMakie.Screen(), f)

# display(GLMakie.Screen(), f_all)
display(GLMakie.Screen(), f_rib)
display(GLMakie.Screen(), f_mrna)
display(GLMakie.Screen(), f_prot)

DataInspector(f_rib)
DataInspector(f_mrna)
DataInspector(f_prot)


zoom = dict_results[folder][2][1:100000,:]

# two versions of determining the switching rate
# if rtca and rtcb are both above a certain threshold
# or if rtca and rtcb are both increasing then we are in the on state, if they are not changing or are decreasing then we are in the off state 
threshold=5
df, start_times, stop_times = switch_times(dict_results[7][2], threshold=threshold)
switch_rate_on, switch_rate_off = calc_switch_rate(df, start_times, stop_times)


# plot switch times
f = Figure()
ax = Axis(f[1,1], title="Threshold for on state = $threshold, switch rate = $(round(switch_rate_on, digits=6))", xlabel="time", ylabel="molecule number")
lines!(ax, df.time, df.rtca, label="rtca")
lines!(ax, df.time, df.rtcb, label="rtcb")
scatter!(ax, start_times, ones(length(start_times)), markersize=10, label="start on")
scatter!(ax, stop_times, ones(length(stop_times)), markersize=10, label="stop on")
axislegend()
DataInspector(f)
display(GLMakie.Screen(), f)    

# plot on states
fdf = filter(row -> row.on == 1, df)
f_rtca = Figure()
ax2 = Axis(f_rtca[1,1], title="Threshold for on state = $threshold, time on = $(round(time_on, digits=2)), switch rate = $(round(switch_rate, digits=6))", xlabel="time", ylabel="molecule number")
lines!(ax2, df.time, df.rtca, label="rtca")
scatter!(ax2, fdf.time, fdf.rtca, markersize=5, color=:pink, label="rtca on")
axislegend()
DataInspector(f_rtca)
display(GLMakie.Screen(), f_rtca) 

f_rtcb = Figure()
ax3 = Axis(f_rtcb[1,1], title="Threshold for on state = $threshold, time on = $(round(time_on, digits=2)), switch rate = $(round(switch_rate, digits=6))", xlabel="time", ylabel="molecule number")
lines!(ax3, df.time, df.rtcb, label="rtcb")
scatter!(ax3, fdf.time, fdf.rtcb, markersize=5, color=:pink, label="rtcb on")
axislegend()
DataInspector(f_rtcb)
display(GLMakie.Screen(), f_rtcb) 


switch_vals_on = []
switch_vals_off = []
for i in eachindex(dict_results)
    println(i)
    switch_val_on, switch_val_off = get_all_switch_rates(i, threshold=threshold)
    push!(switch_vals_on, switch_val_on)
    push!(switch_vals_off, switch_val_off)
end


f_switch = Figure()
ax = Axis(f_switch[1,1], xlabel="kdam", ylabel="switch rate", title="threshold method", yscale=log10)
[lines!(ax, dict_kdamvals[i][:kdam], switch_vals_on[i], label="on→off $i") for i in eachindex(dict_kdamvals)]
[lines!(ax, dict_kdamvals[i][:kdam], switch_vals_off[i], label="off→on $i") for i in eachindex(dict_kdamvals)]

axislegend(position=:rc)
display(GLMakie.Screen(), f_switch)

f = Figure()
ax = Axis(f[1,1], xlabel="kdam", ylabel="switch rate", title="difference method", yscale=log10)
lines!(ax, dict_kdamvals[folder][:kdam], switch_vals_on1, label="on→off")
lines!(ax, dict_kdamvals[folder][:kdam], switch_vals_off1, label="off→on")
axislegend()


fracs_on = []
fracs_off = []
for i in eachindex(switch_vals_on)
    println(i)
    frac_on, frac_off = calc_frac_times(switch_vals_on[i], switch_vals_off[i])
    push!(fracs_on, frac_on)
    push!(fracs_off, frac_off)
end

f_frac = Figure()
ax = Axis(f_frac[1,1], xlabel="kdam", ylabel="fraction of time in state", title="threshold method")
[lines!(ax, dict_kdamvals[i][:kdam], fracs_on[i], label="on state $i") for i in eachindex(dict_kdamvals)]
[lines!(ax, dict_kdamvals[i][:kdam], fracs_off[i], label="off state $i") for i in eachindex(dict_kdamvals)]
axislegend(position=:rc)
display(GLMakie.Screen(), f_frac)


# threshold as unstable steady state line
br = get_br_molec(rtc_model, ssvals_rtc_molec, params_rtc_molec, 1.5)
bf = bf_point_df(br)
df = create_br_df(br)
kdam1 = findall(x->x==bf.kdam[1],df.kdam)[1]
kdam2 = findall(x->x==bf.kdam[2],df.kdam)[1]
unstab_rtca = df[!,:rtca][kdam1:kdam2]
unstab_rtcb = df[!,:rtca][kdam1:kdam2]
unstab_kdam = df[!,:kdam][kdam1:kdam2]

min_bs = ceil(minimum(unstab_kdam), digits=1)
max_bs = floor(maximum(unstab_kdam), digits=1)

min_kdam_ind = find_closest_index(dict_kdamvals[6][:kdam], min_bs)
max_kdam_ind = find_closest_index(dict_kdamvals[6][:kdam], max_bs)

function find_closest_index(array, value)
    differences = abs.(array .- value)
    closest_index = argmin(differences)
    return closest_index
end

thresholds=[]
for i in dict_kdamvals[6][:kdam][min_kdam_ind:max_kdam_ind]
    println(i)
    index = find_closest_index(unstab_kdam, i)
    push!(thresholds, unstab_rtca[index])
end

start_vals = fill(thresholds[1], length(dict_kdamvals[6][:kdam][1:min_kdam_ind-1]))
end_vals = fill(thresholds[end], length(dict_kdamvals[6][:kdam][max_kdam_ind+1:end]))

pushfirst!(thresholds, start_vals...)
push!(thresholds, end_vals...)

thresholds




thresholds_rtca, thresholds_rtcb = get_unstab_threshold_array(6)


lines(df.kdam, df.rtca)
lines!(df.kdam, df.rtcb)
lines!(dict_kdamvals[6][:kdam], thresholds_rtca)
lines!(dict_kdamvals[6][:kdam], thresholds_rtcb)

# average conc for on and off states
threshold = 5 
species = :rh
means_on = get_av_conc_state(dict_results[6][1], threshold, species, on=true)
means_off = get_av_conc_state(dict_results[6][1], threshold, species, on=false)

means_on = get_av_conc_state(dict_results[6][1], thresholds_rtca[1], species, on=true)
means_off = get_av_conc_state(dict_results[6][1], thresholds_rtca[1], species, on=false)

# using singular threshold
species_ons = []
species_offs = []
for i in eachindex(dict_results)
    println(i)
    species_on, species_off = get_all_av_conc(i, 5, species)
    push!(species_ons, species_on)
    push!(species_offs, species_off)
end

cmap = :rainbow_bgyr_35_85_c72_n256
f = Figure()
ax = Axis(f[1,1], xlabel="Damage rate (min-1)", ylabel="[$species] (μM)")
[scatter!(ax, kdams, species_ons[i], color=fracs_on[i], colorrange=(0,1), colormap=cmap) for i in eachindex(species_ons)]
[scatter!(ax, kdams, species_offs[i], color=fracs_off[i], colorrange=(0,1), colormap=cmap) for i in eachindex(species_offs)]
Colorbar(f[1, 2], limits = (0,1), colormap = cmap, label="fraction of time in state")

[thresholds_rtca, thresholds_rtcb]

# using dynamic threshold
species_ons_dt = []
species_offs_dt = []
for i in eachindex(dict_results)
    println(i)
    species_on, species_off = get_all_av_conc(i, thresholds_rtca, species)
    push!(species_ons_dt, species_on)
    push!(species_offs_dt, species_off)
end













df = CSV.read("/Users/s2257179/Desktop/res.csv", DataFrame)
kdams = CSV.read("/Users/s2257179/Desktop/kdam.csv", DataFrame)

convert_to_vector(str) = eval(Meta.parse(str))

kdams = kdams.kdam
rtca_ons = [convert_to_vector(df.rtca_ons[i]) for i in 1:nrow(df)]
rtcb_ons = [convert_to_vector(df.rtcb_ons[i]) for i in 1:nrow(df)]
rtca_offs = [convert_to_vector(df.rtca_offs[i]) for i in 1:nrow(df)]
rtcb_offs = [convert_to_vector(df.rtcb_offs[i]) for i in 1:nrow(df)]
fracs_on = [convert_to_vector(df.fracs_on[i]) for i in 1:nrow(df)]
fracs_off = [convert_to_vector(df.fracs_off[i]) for i in 1:nrow(df)]


rtca_conc = [rtca_ons; rtca_offs]
rtcb_conc = [rtcb_ons; rtcb_offs]
fracs = [fracs_on; fracs_off]

cmap = :rainbow_bgyr_35_85_c72_n256
f = Figure()
ax = Axis(f[1,1], xlabel="Damage rate (min-1)", ylabel="[RtcA] (μM)")
[scatter!(ax, kdams, rtca_conc[i], color=fracs[i], colorrange=(0,1), colormap=cmap) for i in eachindex(rtca_conc)]
Colorbar(f[1, 2], limits = (0,1), colormap = cmap, label="fraction of time in state")

f = Figure()
ax = Axis(f[1,1], xlabel="Damage rate (min-1)", ylabel="[RtcB] (μM)")
[scatter!(ax, kdams, rtcb_conc[i], color=fracs[i], colorrange=(0,1), colormap=cmap) for i in eachindex(rtcb_conc)]
Colorbar(f[1, 2], limits = (0,1), colormap = cmap, label="fraction of time in state")

df = DataFrame(:rtca_ons=>rtca_ons, :rtcb_ons=>rtcb_ons, :rtca_offs=>rtca_offs, :rtcb_offs=>rtcb_offs, :fracs_on=>fracs_on, :fracs_off=>fracs_off)
df_kdam_vals = DataFrame(:kdam=>dict_kdamvals[6][:kdam])

CSV.write("/Users/s2257179/Desktop/res.csv", df)
CSV.write("/Users/s2257179/Desktop/kdam.csv", df_kdam_vals)



























# this file so far does three things
    # calculates the switching rate
    # calculates the fraction of time in each state
    # calculates the average concentration in each state
# for calculating switching rate
threshold=5
df, start_times, stop_times = switch_times(dict_results[7][2], threshold=threshold)
switch_rate_on, switch_rate_off = calc_switch_rate(df, start_times, stop_times)

# switch times does most of the same things as the get_av_conc_state function

function switch_times(res; threshold=nothing)
    df = res[!,[:time, :rtca, :rtcb]]
    condition = (df.rtca .> threshold) .| (df.rtcb .> threshold)
    
    df = insertcols!(df, 1, :on => condition)
    condition_changes = diff(Int.(condition))
    start_indices = findall(x -> x == 1, condition_changes)
    stop_indices = findall(x -> x == -1, condition_changes)
    start_times = []
    for i in start_indices
        push!(start_times, df.time[i])
    end
    stop_times = []
    for i in stop_indices
        push!(stop_times, df.time[i])
    end
    if start_times[1] > stop_times[1] && start_times[end] > stop_times[end]
        # println("start_times[1] > stop_times[1] && start_times[end] > stop_times[end]")
        pushfirst!(start_times, df.time[1])
        push!(stop_times, df.time[end])
    elseif start_times[end] > stop_times[end]
        # println("start_times[end] > stop_times[end]")
        push!(stop_times, df.time[end])
    elseif start_times[1] > stop_times[1] 
        # println("start_times[1] > stop_times[1] ")
        pushfirst!(start_times, df.time[1])
    end
    return df, start_times, stop_times
end