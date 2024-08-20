using InteractiveViz, GLMakie, StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed, TableOperations, JSON, Query, FindFirstFunctions, CategoricalArrays, Colors

include(joinpath(homedir(), "phd/stochastic_hybrid_code/analysis_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_funcs.jl"))

mount_path = "/Users/s2257179/stoch_files/kdam_testing/"
mount_path = "/Users/s2257179/stoch_files/hysteresis/"

all_items = readdir(mount_path)
folders = [item for item in all_items if isdir(joinpath(mount_path, item)) && !occursin("DS", item)]
folders_dict = Dict(i => folder for (i, folder) in enumerate(folders))

folders_dict = Dict(filter(pair -> pair.first in [2], folders_dict))

dict_times, dict_kdamvals, dict_titles, dict_results, dict_reacts, dict_props, dict_counts, dict_hists = setup_dicts(folders_dict)

for i in eachindex(folders_dict)
    println(i)
    dict_times[i], dict_kdamvals[i], dict_titles[i], dict_results[i], dict_reacts[i], dict_props[i] = LoadDataVars(folders[i]);
    dict_hists[i] = load_hist_files(joinpath(mount_path, folders_dict[i], "hists"))
    dict_counts[i] = prod_tot_count(dict_reacts[i])
end

# plot one result
folder = 2; index = 10; species = "rh"; num_plots = 1;
f = plot_results("plot_results", dict_results[folder][index], num_plots, folders_dict[folder], titles=[dict_titles[folder][index]], species="$species", xlabel="time", ylabel="$species", size=(800,650), tosave=false)
f_rib = plot_results("plot_results", dict_results[folder][index], 1, folders_dict[folder], species=[:rh, :rd, :rt], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false)
f_mrna = plot_results("plot_results", dict_results[folder][index], 1, folders_dict[folder], species=[:rm_a, :rm_b, :rm_r], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false)
f_prot = plot_results("plot_results", dict_results[folder][index], 1, folders_dict[folder], species=[:rtca, :rtcb, :rtcr], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false)


get_column(zoom, :rh)./get_column(zoom, :volume)
zoom = dict_results[folder][index][149100:162000,:]
f = plot_results("plot_results", zoom, num_plots, folders_dict[folder], species="$species", xlabel="time", ylabel="$species", size=(800,650), tosave=false)

# f_all = plot_results("plot_results", zoom, 1, folders_dict[folder], species=all_species, xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false)
f_rib = plot_results("plot_results", zoom, 1, folders_dict[folder], species=[:rh, :rd, :rt], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false)
f_mrna = plot_results("plot_results", zoom, 1, folders_dict[folder], species=[:rm_a, :rm_b, :rm_r], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false)
f_prot = plot_results("plot_results", zoom, 1, folders_dict[folder], species=[:rtca, :rtcb, :rtcr], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false)

# display(GLMakie.Screen(), f_all)
display(GLMakie.Screen(), f_rib)
display(GLMakie.Screen(), f_mrna)
display(GLMakie.Screen(), f_prot)

DataInspector(f_rib)
DataInspector(f_mrna)
DataInspector(f_prot)

zoom_props = [dict_props[folder][index][s][149100:162000] for s in eachindex(react_names[1:end-1])]

prop = plot_prop(zoom, zoom_props, index, "kdam_$(dict_kdamvals[folder][:kdam][index]), thresh_150", zoom=true, set_thresh=150, folders_dict[folder], maxval=500, tosave=false, size=(800,650))
DataInspector(prop)



using StatsBase

specie=:rd
data = dict_results[folder][index][!,specie]

lag_range = range(1,200000,length=100)
lags = Int.(round.(collect(lag_range)))
autocors = autocor(data, lags)

f = Figure()
ax = Axis(f[1, 1], xlabel = "lags", ylabel = "autocorrelation",
    title = "Autocorrelation of $specie")
pa = lines!(lag_range, autocors)
DataInspector(pa)


dict_kdamvals[folder][:kdam]

f_rib_on = plot_results("plot_results", dict_results[folder][3], 1, folders_dict[folder], species=[:rh, :rd, :rt], titles=[dict_titles[folder][3]], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false)
f_mrna_on = plot_results("plot_results", dict_results[folder][3], 1, folders_dict[folder], species=[:rm_a, :rm_b, :rm_r], titles=[dict_titles[folder][3]], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false)
f_prot_on = plot_results("plot_results", dict_results[folder][3], 1, folders_dict[folder], species=[:rtca, :rtcb, :rtcr], titles=[dict_titles[folder][3]], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false)

f_rib_off = plot_results("plot_results", dict_results[folder][15], 1, folders_dict[folder], species=[:rh, :rd, :rt], titles=[dict_titles[folder][15]], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false)
f_mrna_off = plot_results("plot_results", dict_results[folder][15], 1, folders_dict[folder], species=[:rm_a, :rm_b, :rm_r], titles=[dict_titles[folder][15]], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false)
f_prot_off = plot_results("plot_results", dict_results[folder][15], 1, folders_dict[folder], species=[:rtca, :rtcb, :rtcr], titles=[dict_titles[folder][15]], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false)

f_rib_off = plot_results("plot_results", dict_results[2][16], 1, folders_dict[folder], species=[:rh, :rd, :rt], titles=[dict_titles[folder][16]], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false)
f_mrna_off = plot_results("plot_results", dict_results[2][16], 1, folders_dict[folder], species=[:rm_a, :rm_b, :rm_r], titles=[dict_titles[folder][16]], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false)
f_prot_off = plot_results("plot_results", dict_results[2][16], 1, folders_dict[folder], species=[:rtca, :rtcb, :rtcr], titles=[dict_titles[folder][16]], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false)

display(GLMakie.Screen(), f_rib_on)
display(GLMakie.Screen(), f_mrna_on)
display(GLMakie.Screen(), f_prot_on)

display(GLMakie.Screen(), f_rib_off)
display(GLMakie.Screen(), f_mrna_off)
display(GLMakie.Screen(), f_prot_off)

DataInspector(f_rib_off)
DataInspector(f_mrna_off)
DataInspector(f_prot_off)
