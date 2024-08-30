using StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed, TableOperations, JSON, Query, FindFirstFunctions, CategoricalArrays
using InteractiveViz, GLMakie

include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_funcs.jl"))


mount_path, folders, folders_dict = load_file_structure("kdam_testing")
folders_dict = Dict(filter(pair -> pair.first in [7], folders_dict))
dict_times, dict_kdamvals, dict_titles, dict_results, dict_reacts, dict_props, dict_counts, dict_hists = load_data(mount_path, folders, folders_dict)

folder = 6; index = 1; species = "rtca"; num_plots = 1;
f = plot_results("plot_results", dict_results[folder][index], num_plots, folders_dict[folder], titles=[dict_titles[folder][index]], species="$species", xlabel="time", ylabel="$species", size=(800,650), tosave=false, conc=true)
f_rib = plot_results("plot_results", dict_results[folder][index], 1, folders_dict[folder], titles=[dict_titles[folder][index]], species=[:rh, :rd, :rt], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false, conc=false)
f_mrna = plot_results("plot_results", dict_results[folder][index], 1, folders_dict[folder], species=[:rm_a, :rm_b, :rm_r], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false, conc=true)
f_prot = plot_results("plot_results", dict_results[folder][index], 1, folders_dict[folder], titles=[dict_titles[folder][index]], species=[:rtca, :rtcb, :rtcr], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false, conc=false)

display(GLMakie.Screen(), f.plot)
DataInspector(f)
