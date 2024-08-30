using StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed, TableOperations, JSON, Query, FindFirstFunctions, CategoricalArrays, Colors

# using PlotlyJS
using InteractiveViz, WGLMakie

# include(joinpath(homedir(), "phd/stochastic_hybrid_code/analysis_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_funcs.jl"))


mount_path = "/Users/s2257179/stoch_files/kdam_testing/"
all_items = readdir(mount_path)
folders = [item for item in all_items if isdir(joinpath(mount_path, item)) && !occursin("DS", item)]
folders_dict = Dict(i => folder for (i, folder) in enumerate(folders))

folders_dict = Dict(filter(pair -> pair.first in [7], folders_dict))

dict_times, dict_kdamvals, dict_titles, dict_results, dict_reacts, dict_props, dict_counts, dict_hists = setup_dicts(folders_dict)

for i in eachindex(folders_dict)
    println(i)
    dict_times[i], dict_kdamvals[i], dict_titles[i], dict_results[i], dict_reacts[i], dict_props[i] = LoadDataVars(folders[i]);
    dict_hists[i] = load_hist_files(joinpath(mount_path, folders_dict[i], "hists"))
    dict_counts[i] = prod_tot_count(dict_reacts[i])
end

# all propensities 
dict_plot_times, dict_plot_counts, dict_plot_results, dict_plot_hists, dict_stoch_reacts, dict_plot_props = setup_plot_dicts()
for folder in eachindex(folders_dict)
    println(folder)
    for i in 1:length(dict_props[folder])
        println(dict_kdamvals[folder][:kdam][i])
        if dict_kdamvals[folder][:kdam][i] == 0.0
            continue
        end
        dict_plot_props[folder, dict_kdamvals[folder][:kdam][i]] = plot_prop(dict_results[folder], dict_props[folder], i, "kdam_$(dict_kdamvals[folder][:kdam][i]), thresh_150", set_thresh=150, folders_dict[folder], maxval=500, tosave=true, size=(800,650))
    end
end

# individual
plot_prop(dict_results[folder], dict_props[folder], 2, "kdam_$(dict_kdamvals[folder][:kdam][2]), thresh_150", set_thresh=150, folders_dict[folder], maxval=500, tosave=false, size=(800,650))