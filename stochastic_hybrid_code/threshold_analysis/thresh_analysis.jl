using StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed, TableOperations, JSON, Query, FindFirstFunctions, CategoricalArrays, Colors

# using PlotlyJS
using InteractiveViz, WGLMakie

include(joinpath(homedir(), "phd/stochastic_hybrid_code/analysis_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_funcs.jl"))

mount_path = "/Users/s2257179/stoch_files/threshold_testing/"

all_items = readdir(mount_path)
folders = [item for item in all_items if isdir(joinpath(mount_path, item)) && !occursin("DS", item)]
folders_dict = Dict(i => folder for (i, folder) in enumerate(folders))

delete!(folders_dict, 12)

folders_dict = Dict(filter(pair -> pair.first in [1, 8], folders_dict))

dict_times, dict_threshvals, dict_titles, dict_results, dict_reacts, dict_props, dict_counts, dict_hists = setup_dicts(folders_dict)

for i in eachindex(folders_dict)
    println(i)
    dict_times[i], dict_threshvals[i], dict_titles[i], dict_results[i], dict_reacts[i], dict_props[i] = LoadDataVars(folders[i]);
    dict_hists[i] = load_hist_files(joinpath(mount_path, folders_dict[i], "hists"))
    # calculate total counts
    dict_counts[i] = prod_tot_count(dict_reacts[i])
end


dict_plot_times, dict_plot_counts, dict_plot_results, dict_plot_hists, dict_stoch_reacts, dict_plot_props = setup_plot_dicts()

for i in eachindex(folders_dict)
    println(i)
    # plot times
    dict_plot_times[i] = plot_times(dict_times[i], "$(folders_dict[i])", folder=folders_dict[i])
    
    dict_plot_counts[i] = plot_totstochcount(dict_threshvals[i], dict_counts[i], "$(folders_dict[i])", folder=folders_dict[i])

    dict_stoch_reacts[i] = plot_results("plot_stoch_reacts", dict_reacts[i], length(dict_threshvals[i]), folders_dict[i], xlabel="reaction", ylabel="count", titles=dict_titles[i], size=(1000,650), tosave=true)

end



i = 1; specie=:rtca
dict_plot_results[i, specie] = plot_results("plot_results", dict_results[i], length(dict_threshvals[i]), folders_dict[i], species=specie, xlabel="time", ylabel="$specie", titles=dict_titles[i], size=(1000,650))#, folder=folders_dict[i]);
dict_plot_hists[i,specie] = plot_results("plot_hists", dict_hists[i], length(dict_threshvals[i]), folders_dict[i], species=specie, xlabel="$specie", ylabel="frequency", titles=dict_titles[i], hidelabels=[false, false], linkaxes=false, size=(1000,650))#, folder=folders_dict[i]);

display(dict_plot_times[6])
display(dict_plot_results[i, :rtca])
display(dict_plot_counts[6])
display(dict_stoch_reacts[6])



# to get all plots (probably not needed often)
for specie in all_species
    println(specie)
    for i in eachindex(folders_dict)
        println(i)
        dict_plot_results[i, specie] = plot_results("plot_results", dict_results[i], length(dict_threshvals[i]), folders_dict[i], species=specie, xlabel="time", ylabel="$specie", titles=dict_titles[i], size=(1000,650), tosave=true);
        dict_plot_hists[i,specie] = plot_results("plot_hists", dict_hists[i], length(dict_threshvals[i]), folders_dict[i], species=specie, xlabel="$specie", ylabel="frequency", titles=dict_titles[i], hidelabels=[false, false], linkaxes=false, size=(1000,650), tosave=true);
    end
end

# all propensities 
for folder in eachindex(folders_dict)
    println(folder)
    for i in eachindex(dict_props[folder])
        println(dict_threshvals[folder][i])
        dict_plot_props[folder, dict_threshvals[folder][i]] = plot_prop(dict_results[folder], dict_props[folder], i, "threshold_$(dict_threshvals[folder][i])", dict_threshvals[folder], folders_dict[folder], maxval=500, tosave=true, size=(800,650))
    end
end

folder = 1; i = 1
plot_prop(dict_results[folder], dict_props[folder], i, "threshold_$(dict_threshvals[folder][i])", dict_threshvals[folder], folders_dict[folder], maxval=500, tosave=true, size=(800,650))



f1 = plot_prop(dict_results[i], dict_props[i], 4, "threshold_$(dict_threshvals[i][4])", dict_threshvals[i], folders_dict[i])




display(f1)
f3_3 = plot_prop(df_results3, df_props3, 2, "threshold_$(threshold_vals3[2])", threshold_vals3)


df_results11[1].rm_a

h2 = hist(df_results11[2].rtca[100000:200000])
h3 = hist!(df_results11[2].rtca[200000:300000])
h4 = hist!(df_results11[2].rtca[300000:400000])
h5 = hist!(df_results11[2].rtca[400000:500000])
h6 = hist!(df_results11[2].rtca[500000:600000])
h7 = hist!(df_results11[2].rtca[700000:800000])
h8 = hist!(df_results11[2].rtca[800000:end])



10000*log(2)/lam_val

df_results11[2][1:end,:rm_a]