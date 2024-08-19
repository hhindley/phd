using StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed, TableOperations, JSON, Query, FindFirstFunctions, CategoricalArrays, Colors

# using PlotlyJS
using InteractiveViz, WGLMakie

include(joinpath(homedir(), "phd/stochastic_hybrid_code/analysis_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_funcs.jl"))

mount_path = "/Users/s2257179/stoch_files/kdam_testing/"

all_items = readdir(mount_path)
folders = [item for item in all_items if isdir(joinpath(mount_path, item)) && !occursin("DS", item)]
folders_dict = Dict(i => folder for (i, folder) in enumerate(folders))

folders_dict = Dict(filter(pair -> pair.first in [7, 12, 8, 13, 10, 11, 9, 14], folders_dict))

dict_times, dict_kdamvals, dict_titles, dict_results, dict_reacts, dict_props, dict_counts, dict_hists = setup_dicts(folders_dict)

for i in eachindex(folders_dict)
    println(i)
    dict_times[i], dict_kdamvals[i], dict_titles[i], dict_results[i], dict_reacts[i], dict_props[i] = LoadDataVars(folders[i]);
    dict_hists[i] = load_hist_files(joinpath(mount_path, folders_dict[i], "hists"))
    dict_counts[i] = prod_tot_count(dict_reacts[i])
end

dict_results

dict_plot_times, dict_plot_counts, dict_plot_results, dict_plot_hists, dict_stoch_reacts, dict_plot_props = setup_plot_dicts()

for i in eachindex(folders_dict)
    println(i)
    dict_plot_times[i] = plot_times(dict_times[i], "$(folders_dict[i])", folder=folders_dict[i])    
    dict_plot_counts[i] = plot_totstochcount(dict_kdamvals[i][:kdam], dict_counts[i], "$(folders_dict[i])", folder=folders_dict[i])
    dict_stoch_reacts[i] = plot_results("plot_stoch_reacts", dict_reacts[i], length(dict_kdamvals[i][:kdam]), folders_dict[i], xlabel="reaction", ylabel="count", titles=dict_titles[i], size=(1000,650), tosave=true)
end

# to get all plots (probably not needed often)
for specie in all_species
    println(specie)
    for i in eachindex(folders_dict)
        println(i)
        dict_plot_results[i, specie] = plot_results("plot_results", dict_results[i], length(dict_kdamvals[i][:kdam]), folders_dict[i], species=specie, xlabel="time", ylabel="$specie", titles=dict_titles[i], size=(1000,650), tosave=true);
        dict_plot_hists[i,specie] = plot_results("plot_hists", dict_hists[i], length(dict_kdamvals[i][:kdam]), folders_dict[i], species=specie, xlabel="$specie", ylabel="frequency", titles=dict_titles[i], hidelabels=[true, true], linkaxes=true, size=(1000,650), tosave=true);
    end
end

(dict_plot_hists[1, :rh])
(dict_plot_hists[1, :rd])
(dict_plot_hists[1, :rm_a])

# all propensities 
for folder in eachindex(folders_dict)
    println(folder)
    for i in 1:length(dict_props[folder])
        println(dict_kdamvals[folder][:kdam][i])
        dict_plot_props[folder, dict_kdamvals[folder][:kdam][i]] = plot_prop(dict_results[folder], dict_props[folder], i, "kdam_$(dict_kdamvals[folder][:kdam][i]), thresh_$(dict_kdamvals[folder][:threshold][i])", dict_kdamvals[folder][:threshold], folders_dict[folder], maxval=500, tosave=true, size=(800,650))
    end
end

# making the hists have less bins to see if we see anything different 
agg_hists = Dict{Int64, DataFrame}()
for i in eachindex(dict_hists[1]["rm_a"])
    agg_hists[i] = hists_with_less_bins(dict_hists[1]["rm_a"][i], 10)
end
agg_hists
dict_hists
agg_hists[5]

f, ax = plot_hist(agg_hists[8], maxval=1300)

display(f)

keys(dict_plot_props)
display(dict_plot_props[1, 0.005])

# plot all hists on top of eachother 
for specie in all_species
    plot_hists_overlay(1, "$specie", 2, last=[7], folder=folders_dict[1], maxval=2e4)
end


# average molecule numbers over whole simulation
mean(dict_results[1][1].rm_a)
mean(dict_results[1][2].rm_a)
mean(dict_results[1][3].rm_a)
mean(dict_results[1][4].rm_a)
mean(dict_results[1][5].rm_a)
mean(dict_results[1][6].rm_a)
mean(dict_results[1][7].rm_a)
mean(dict_results[1][8].rm_a)
mean(dict_results[1][9].rm_a)
mean(dict_results[1][10].rm_a)
mean(dict_results[1][11].rm_a)
mean(dict_results[1][12].rm_a)
mean(dict_results[1][13].rm_a)
mean(dict_results[1][14].rm_a)
mean(dict_results[1][15].rm_a)
