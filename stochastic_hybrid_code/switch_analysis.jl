using InteractiveViz, GLMakie, StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed, TableOperations, JSON, Query, FindFirstFunctions, CategoricalArrays, Colors

include(joinpath(homedir(), "phd/stochastic_hybrid_code/analysis_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_funcs.jl"))

mount_path = "/Users/s2257179/stoch_files/kdam_testing/"
# mount_path = "/Users/s2257179/stoch_files/hysteresis/"
# mount_path = "/Users/s2257179/stoch_files/threshold_testing/"

all_items = readdir(mount_path)
folders = [item for item in all_items if isdir(joinpath(mount_path, item)) && !occursin("DS", item)]
folders_dict = Dict(i => folder for (i, folder) in enumerate(folders))

folders_dict = Dict(filter(pair -> pair.first in [6], folders_dict))

dict_times, dict_kdamvals, dict_titles, dict_results, dict_reacts, dict_props, dict_counts, dict_hists = setup_dicts(folders_dict)

for i in eachindex(folders_dict)
    println(i)
    dict_times[i], dict_kdamvals[i], dict_titles[i], dict_results[i], dict_reacts[i], dict_props[i] = LoadDataVars(folders[i]);
    dict_hists[i] = load_hist_files(joinpath(mount_path, folders_dict[i], "hists"))
    dict_counts[i] = prod_tot_count(dict_reacts[i])
end

dict_kdamvals[6][:kdam]
# plot one result
folder = 6; index = 2; species = "rtca"; num_plots = 1;
f = plot_results("plot_results", dict_results[folder][index], num_plots, folders_dict[folder], titles=[dict_titles[folder][index]], species="$species", xlabel="time", ylabel="$species", size=(800,650), tosave=false)
f_rib = plot_results("plot_results", dict_results[folder][index], 1, folders_dict[folder], titles=[dict_titles[folder][index]], species=[:rh, :rd, :rt], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false, conc=true)
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



# if you want to plot with log scale then need to remove the zeros and replace with a very small value
epsilon = 1e-5
dict_results[folder][index] = replace_zeros_in_dataframe(dict_results[folder][index])
nodam = plot_results("plot_results", dict_results[folder][index], 1, folders_dict[folder], species=all_species, titles=[dict_titles[folder][index]], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false, yscale=identity, conc=false)
nodam_conc = plot_results("plot_results", dict_results[folder][index], 1, folders_dict[folder], species=all_species, titles=[dict_titles[folder][index]], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false, yscale=identity, conc=true)

display(GLMakie.Screen(), nodam)
display(GLMakie.Screen(), nodam_conc)

DataInspector(nodam)
DataInspector(nodam_conc)



zoom = dict_results[folder][2][1:100000,:]

function get_switch_times(res)
    df = res[!,[:time, :rtca, :rtcb]]
    df = insertcols!(df, 1, :index => 1:nrow(df))
    # fdf = filter(row -> row.rtca > 2 && row.rtcb > 2, df)
    fdf = df[(df.rtca .> 5) .& (df.rtcb .> 5), :]

    condition = diff(fdf.time) .> 1
    indices = findall(x -> x == 1, condition)

    conditions = [(any(df.rtca[fdf.index[i]:fdf.index[i+1]] .== 0) && any(df.rtcb[fdf.index[i]:fdf.index[i+1]] .== 0)) for i in indices]

    filtered_indices = [i for (i, cond) in zip(indices, conditions) if cond]

    start_times = [fdf.time[1]]
    push!(start_times, (df.time[fdf.index[filtered_indices .+ 1]]...))
    stop_times = [fdf.time[end]]
    pushfirst!(stop_times, (df.time[fdf.index[filtered_indices]]...))

    return start_times, stop_times
end

start_times, stop_times = get_switch_times(zoom)

time_on = sum(stop_times.-start_times)  

switch_rate = 1/time_on


@elapsed start_times, stop_times = get_switch_times(dict_results[folder][3])
@elapsed start_times, stop_times = get_switch_times(dict_results[folder][5])

f = Figure()
ax = Axis(f[1,1])
lines!(ax, df.time, df.rtca)
scatter!(ax, start_times, ones(length(start_times)), markersize=10)
scatter!(ax, stop_times, ones(length(stop_times)), markersize=10)


df = zoom[!,[:time, :rtca, :rtcb]]
# df = insertcols!(df, 1, :index => 1:nrow(df))
condition = (df.rtca .> 1.99) .& (df.rtcb .> 1.99)
df = insertcols!(df, 1, :on => condition)

condition_changes = diff(Int.(condition))
true_indices = findall(x -> x == 1, condition_changes)
if condition[1]
    true_indices = [1; true_indices .+ 1]
else
    true_indices = true_indices .+ 1
end

on_index = []
for index in true_indices
    if index + 10 <= nrow(df)
        rtca_changes = diff(df.rtca[index:index+10])
        rtcb_changes = diff(df.rtcb[index:index+10])
        if any(rtca_changes .!= 0) && any(rtcb_changes .!= 0)
            push!(on_index, index)
        end
        # println(any(rtca_changes .!= 0) && any(rtcb_changes .!= 0))
    end
end
on_index
s1 = []
for i in on_index
    push!(s1, df.time[i])
end
s1

f = Figure()
ax = Axis(f[1,1])
lines!(ax, zoom.time, zoom.rtca)
scatter!(ax, s1, ones(length(s1)), markersize=10)
scatter!(ax, stop_times, ones(length(stop_times)), markersize=10)

df[630:649,:]
df[648,:]





# two versions of determining the switching rate
# if rtca and rtcb are both above a certain threshold
threshold = 10
df = dict_results[folder][2][!,[:time, :rtca, :rtcb]]
condition = (df.rtca .> threshold) .& (df.rtcb .> threshold)
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
f = Figure()
ax = Axis(f[1,1])
lines!(ax, df.time, df.rtca)
scatter!(ax, start_times, ones(length(start_times)), markersize=10)
scatter!(ax, stop_times, ones(length(stop_times)), markersize=10)


start_times
stop_times
stop_times = stop_times[2:end]

time_on = sum(stop_times.-start_times)
switch_rate = 1/time_on

# or if rtca and rtcb are both increasing then we are in the on state, if they are not changing or are decreasing then we are in the off state 








start_dict = Dict(i => [] for i in 1:length(dict_results[folder]))
stop_dict = Dict(i => [] for i in 1:length(dict_results[folder]))
for i in eachindex(dict_results[folder])
    println(i)
    start_times, stop_times = get_switch_times(dict_results[folder][i])
    start_dict[i] = start_times
    stop_dict[i] = stop_times
end

f = Figure()
ax = Axis(f[1,1])
lines!(ax, df.time, df.rtca)
scatter!(ax, start_times, ones(length(start_times)), markersize=10)
scatter!(ax, stop_times, ones(length(stop_times)), markersize=10)


switch_rates = Dict{Int, Any}()
for i in eachindex(start_dict)
    time_on = sum(stop_dict[i].-start_dict[i])
    switch_rates[i] = 1/time_on
end

sorted_keys = sort(collect(keys(switch_rates)))
switch_vals = [switch_rates[key] for key in sorted_keys]

f = Figure()
ax = Axis(f[1,1])
lines!(ax, dict_kdamvals[folder][:kdam], switch_vals)


