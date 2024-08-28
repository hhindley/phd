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

# two versions of determining the switching rate
# if rtca and rtcb are both above a certain threshold
# or if rtca and rtcb are both increasing then we are in the on state, if they are not changing or are decreasing then we are in the off state 
function switch_times(res; threshold=nothing)
    df = res[!,[:time, :rtca, :rtcb]]
    if threshold !== nothing
        condition = (df.rtca .> threshold) .& (df.rtcb .> threshold)
    else
        diff_rtca = diff(df.rtca)
        diff_rtcb = diff(df.rtcb)
        condition = (diff_rtca .> 1e-10) .| (diff_rtcb .> 1e-10)
        pushfirst!(condition, condition[1])
    end
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
    return df, start_times, stop_times
end
threshold=10
df, start_times, stop_times = switch_times(dict_results[folder][2], threshold=threshold)
df1, start_times1, stop_times1 = switch_times(dict_results[folder][2])

time_on = sum(stop_times[2:end].-start_times)
time_off = df.time[end] - time_on
switch_rate_on = 1/time_on
switch_rate_off = 1/time_off

time_on1 = sum(stop_times1.-start_times1)
time_off1 = df1.time[end] - time_on1
switch_rate1 = 1/time_on1
switch_rate_off1 = 1/time_off1

# plot switch times
f = Figure()
ax = Axis(f[1,1], title="Threshold for on state = $threshold, time on = $(round(time_on, digits=2)), switch rate = $(round(switch_rate, digits=6))", xlabel="time", ylabel="molecule number")
lines!(ax, df.time, df.rtca, label="rtca")
lines!(ax, df.time, df.rtcb, label="rtcb")
scatter!(ax, start_times, ones(length(start_times)), markersize=10, label="start on")
scatter!(ax, stop_times, ones(length(stop_times)), markersize=10, label="stop on")
axislegend()
DataInspector(f)
display(GLMakie.Screen(), f)    

f1 = Figure()
ax1 = Axis(f1[1,1], title="Threshold for on state = $threshold, time on = $(round(time_on, digits=2)), switch rate = $(round(switch_rate, digits=6))", xlabel="time", ylabel="molecule number")
lines!(ax1, df1.time, df1.rtca, label="rtca")
lines!(ax1, df1.time, df1.rtcb, label="rtcb")
scatter!(ax1, start_times1, ones(length(start_times1)), markersize=10, label="start on")
scatter!(ax1, stop_times1, ones(length(stop_times1)), markersize=10, label="stop on")
axislegend()
DataInspector(f1)
display(GLMakie.Screen(), f1)

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

fdf1 = filter(row -> row.on == 1, df1)
f_rtca1 = Figure()
ax4 = Axis(f_rtca1[1,1], title="On = protein increasing, time on = $(round(time_on1, digits=2)), switch rate = $(round(switch_rate1, digits=6))")
lines!(ax4, df1.time, df1.rtca, label="rtca")
scatter!(ax4, fdf1.time, fdf1.rtca, markersize=5, color=:pink, label="rtca on")
axislegend()
DataInspector(f_rtca1)
display(GLMakie.Screen(), f_rtca1) 

f_rtcb1 = Figure()
ax5 = Axis(f_rtcb1[1,1], title="On = protein increasing, time on = $(round(time_on1, digits=2)), switch rate = $(round(switch_rate1, digits=6))")
lines!(ax5, df1.time, df1.rtcb, label="rtcb")
scatter!(ax5, fdf1.time, fdf1.rtcb, markersize=5, color=:pink, label="rtcb on")
axislegend()
DataInspector(f_rtcb1)
display(GLMakie.Screen(), f_rtcb1) 


function get_all_switch_rates(folder; threshold=nothing)
    dfs = Dict{Any, Any}()  
    start_dict = Dict(i => [] for i in 1:length(dict_results[folder]))
    stop_dict = Dict(i => [] for i in 1:length(dict_results[folder]))
    for i in eachindex(dict_results[folder])
        println(i)
        df, start_times, stop_times = switch_times(dict_results[folder][i], threshold=threshold)
        dfs[i] = df
        start_dict[i] = start_times
        stop_dict[i] = stop_times
    end

    switch_rates_on = Dict{Int, Any}()
    switch_rates_off = Dict{Int, Any}()
    for i in eachindex(start_dict)
        if length(stop_dict[i]) !== length(start_dict[i])
            time_on = sum(stop_dict[i][2:end].-start_dict[i])
        else
            time_on = sum(stop_dict[i].-start_dict[i])
        end
        switch_rates_on[i] = 1/time_on
        time_off = dfs[i].time[end] - time_on
        switch_rates_off[i] = 1/time_off

    end

    sorted_keys = sort(collect(keys(switch_rates)))
    switch_vals_on = [switch_rates_on[key] for key in sorted_keys]
    switch_vals_off = [switch_rates_off[key] for key in sorted_keys]

    return switch_vals_on, switch_vals_off
end

switch_vals_on, switch_vals_off = get_all_switch_rates(6, threshold=threshold)
switch_vals_on1, switch_vals_off1 = get_all_switch_rates(6)

f = Figure()
ax = Axis(f[1,1], xlabel="kdam", ylabel="switch rate")
lines!(ax, dict_kdamvals[folder][:kdam], switch_vals_on, label="threshold method")
lines!(ax, dict_kdamvals[folder][:kdam], switch_vals_on1, label="difference method")
axislegend()

f = Figure()
ax = Axis(f[1,1], xlabel="kdam", ylabel="switch rate")
lines!(ax, dict_kdamvals[folder][:kdam], switch_vals_off, label="threshold method")
lines!(ax, dict_kdamvals[folder][:kdam], switch_vals_off1, label="difference method")
axislegend()

f = Figure()
ax = Axis(f[1,1], xlabel="kdam", ylabel="switch rate")
lines!(ax, dict_kdamvals[folder][:kdam], switch_vals_on, label="on→off")
lines!(ax, dict_kdamvals[folder][:kdam], switch_vals_off, label="off→on")
axislegend()

f = Figure()
ax = Axis(f[1,1], xlabel="kdam", ylabel="switch rate")
lines!(ax, dict_kdamvals[folder][:kdam], switch_vals_on1, label="on→off")
lines!(ax, dict_kdamvals[folder][:kdam], switch_vals_off1, label="off→on")
axislegend()
