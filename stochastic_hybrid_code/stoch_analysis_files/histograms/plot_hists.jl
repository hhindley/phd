using StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed, TableOperations, JSON, Query, FindFirstFunctions, CategoricalArrays
using InteractiveViz, GLMakie

include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_funcs.jl"))

mount_path, folders, folders_dict = load_file_structure("kdam_testing")
folders_dict = Dict(filter(pair -> pair.first in [7], folders_dict))
dict_times, dict_kdamvals, dict_titles, dict_results, dict_reacts, dict_props, dict_counts, dict_hists = load_data(mount_path, folders, folders_dict)

folder = 7; index = 2; specie = :rh;
# to plot one at a time
p = plot_results("plot_hists", dict_hists[folder]["$specie"][index], 1, folders_dict[folder], species=specie, xlabel="$specie", ylabel="frequency", titles=[dict_titles[folder][index]], hidelabels=[true, true], linkaxes=true, size=(1000,650), tosave=false);
display(GLMakie.Screen(), p)


# plot all from one folder 
p = plot_results("plot_hists", dict_hists[folder], length(dict_kdamvals[folder][:kdam]), folders_dict[folder], species=specie, xlabel="$specie", ylabel="frequency", titles=dict_titles[folder], hidelabels=[true, true], linkaxes=true, size=(1000,650), tosave=false);
display(GLMakie.Screen(), p)



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

