using StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed, TableOperations, JSON, Query, FindFirstFunctions, CategoricalArrays, Colors
using InteractiveViz, GLMakie

# include(joinpath(homedir(), "phd/stochastic_hybrid_code/analysis_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_funcs.jl"))


mount_path, folders, folders_dict = load_file_structure("kdam_testing")
folders_dict = Dict(filter(pair -> pair.first in [7], folders_dict))
dict_times, dict_kdamvals, dict_titles, dict_results, dict_reacts, dict_props, dict_counts, dict_hists = load_data(mount_path, folders, folders_dict)


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

dict_plot_results[6, :rtca]

(dict_plot_hists[1, :rh])
(dict_plot_hists[1, :rd])
(dict_plot_hists[1, :rm_a])

# all propensities 
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
