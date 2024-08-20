using InteractiveViz, GLMakie, StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed, TableOperations, JSON, Query, FindFirstFunctions, CategoricalArrays, Colors

include(joinpath(homedir(), "phd/stochastic_hybrid_code/analysis_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_funcs.jl"))

mount_path = "/Users/s2257179/stoch_files/kdam_testing/"

all_items = readdir(mount_path)
folders = [item for item in all_items if isdir(joinpath(mount_path, item)) && !occursin("DS", item)]
folders_dict = Dict(i => folder for (i, folder) in enumerate(folders))

folders_dict = Dict(filter(pair -> pair.first in [4, 8], folders_dict))

dict_times, dict_kdamvals, dict_titles, dict_results, dict_reacts, dict_props, dict_counts, dict_hists = setup_dicts(folders_dict)

for i in eachindex(folders_dict)
    println(i)
    dict_times[i], dict_kdamvals[i], dict_titles[i], dict_results[i], dict_reacts[i], dict_props[i] = LoadDataVars(folders[i]);
    dict_hists[i] = load_hist_files(joinpath(mount_path, folders_dict[i], "hists"))
    dict_counts[i] = prod_tot_count(dict_reacts[i])
end

# plot one result
folder = 8; index = 4; species = "rtca"; num_plots = 1;
f = plot_results("plot_results", dict_results[folder][index], num_plots, folders_dict[folder], titles=[dict_titles[folder][index]], species="$species", xlabel="time", ylabel="$species", size=(1000,650), tosave=false)

zoom = dict_results[folder][index][149100:162000,:]

f_all = plot_results("plot_results", zoom, 1, folders_dict[folder], species=all_species, xlabel="time", ylabel="specie", size=(1000,650), tosave=false, linkaxes=false)
f_rib = plot_results("plot_results", zoom, 1, folders_dict[folder], species=[:rh, :rd, :rt], xlabel="time", ylabel="specie", size=(1000,650), tosave=false, linkaxes=false)
f_mrna = plot_results("plot_results", zoom, 1, folders_dict[folder], species=[:rm_a, :rm_b, :rm_r], xlabel="time", ylabel="specie", size=(1000,650), tosave=false, linkaxes=false)
f_prot = plot_results("plot_results", zoom, 1, folders_dict[folder], species=[:rtca, :rtcb, :rtcr], xlabel="time", ylabel="specie", size=(1000,650), tosave=false, linkaxes=false)

display(GLMakie.Screen(), f_all)
display(GLMakie.Screen(), f_rib)
display(GLMakie.Screen(), f_mrna)
display(GLMakie.Screen(), f_prot)

DataInspector(f_rib)
DataInspector(f_mrna)
DataInspector(f_prot)

prop = plot_prop(dict_results[folder], dict_props[folder], index, "kdam_$(dict_kdamvals[folder][:kdam][index]), thresh_150", set_thresh=150, folders_dict[folder], maxval=500, tosave=false, size=(800,650))
DataInspector(prop)
