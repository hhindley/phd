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

folders_dict = Dict(filter(pair -> pair.first in [7,8,10,11], folders_dict))

dict_times, dict_kdamvals, dict_titles, dict_results, dict_reacts, dict_props, dict_counts, dict_hists = setup_dicts(folders_dict)

for i in eachindex(folders_dict)
    println(i)
    dict_times[i], dict_kdamvals[i], dict_titles[i], dict_results[i], dict_reacts[i], dict_props[i] = LoadDataVars(folders[i]);
    dict_hists[i] = load_hist_files(joinpath(mount_path, folders_dict[i], "hists"))
    dict_counts[i] = prod_tot_count(dict_reacts[i])
end

dict_plot_times, dict_plot_counts, dict_plot_results, dict_plot_hists, dict_stoch_reacts, dict_plot_props = setup_plot_dicts()


ind_005 = [[11, 2], [10, 5], [7, 1]]
ind_01 = [[11, 3], [10, 7], [7, 2]]
ind_05 = [[11, 4], [10, 8], [7, 4]]
ind_1 = [[11, 5], [10, 9], [7, 6], [8, 2], [8, 3]]
ind_5 = [[11, 13], [10, 10], [7, 8], [8, 5]]

dict_reacts_new = Dict("0.005" => [dict_reacts[ind_005[1][1]][ind_005[1][2]], dict_reacts[ind_005[2][1]][ind_005[2][2]], dict_reacts[ind_005[3][1]][ind_005[3][2]]], 
                        "0.01" => [dict_reacts[ind_01[1][1]][ind_01[1][2]], dict_reacts[ind_01[2][1]][ind_01[2][2]], dict_reacts[ind_01[3][1]][ind_01[3][2]]], 
                        "0.05" => [dict_reacts[ind_05[1][1]][ind_05[1][2]], dict_reacts[ind_05[2][1]][ind_05[2][2]], dict_reacts[ind_05[3][1]][ind_05[3][2]]], 
                        "0.1" => [dict_reacts[ind_1[1][1]][ind_1[1][2]], dict_reacts[ind_1[2][1]][ind_1[2][2]], dict_reacts[ind_1[3][1]][ind_1[3][2]], dict_reacts[ind_1[4][1]][ind_1[4][2]], dict_reacts[ind_1[5][1]][ind_1[5][2]]], 
                        "0.5" => [dict_reacts[ind_5[1][1]][ind_5[1][2]], dict_reacts[ind_5[2][1]][ind_5[2][2]], dict_reacts[ind_5[3][1]][ind_5[3][2]], dict_reacts[ind_5[4][1]][ind_5[4][2]]]
)

dict_kdamvals_new = Dict("0.005" => Dict(:threshold=>[dict_kdamvals[ind_005[1][1]][:threshold][ind_005[1][2]], dict_kdamvals[ind_005[2][1]][:threshold][ind_005[2][2]], dict_kdamvals[ind_005[3][1]][:threshold][ind_005[3][2]]], :kdam=>[dict_kdamvals[ind_005[1][1]][:kdam][ind_005[1][2]],dict_kdamvals[ind_005[2][1]][:kdam][ind_005[2][2]],dict_kdamvals[ind_005[3][1]][:kdam][ind_005[3][2]]]),
                        "0.01" => Dict(:threshold=>[dict_kdamvals[ind_01[1][1]][:threshold][ind_01[1][2]], dict_kdamvals[ind_01[2][1]][:threshold][ind_01[2][2]], dict_kdamvals[ind_01[3][1]][:threshold][ind_01[3][2]]], :kdam=>[dict_kdamvals[ind_01[1][1]][:kdam][ind_01[1][2]], dict_kdamvals[ind_01[2][1]][:kdam][ind_01[2][2]], dict_kdamvals[ind_01[3][1]][:kdam][ind_01[3][2]]]),
                        "0.05" => Dict(:threshold=>[dict_kdamvals[ind_05[1][1]][:threshold][ind_05[1][2]],dict_kdamvals[ind_05[2][1]][:threshold][ind_05[2][2]],dict_kdamvals[ind_05[3][1]][:threshold][ind_05[3][2]]], :kdam=>[dict_kdamvals[ind_05[1][1]][:kdam][ind_05[1][2]],dict_kdamvals[ind_05[2][1]][:kdam][ind_05[2][2]],dict_kdamvals[ind_05[3][1]][:kdam][ind_05[3][2]]]),
                        "0.1" => Dict(:threshold=>[dict_kdamvals[ind_1[1][1]][:threshold][ind_1[1][2]],dict_kdamvals[ind_1[2][1]][:threshold][ind_1[2][2]],dict_kdamvals[ind_1[3][1]][:threshold][ind_1[3][2]],dict_kdamvals[ind_1[4][1]][:threshold][ind_1[4][2]],dict_kdamvals[ind_1[5][1]][:threshold][ind_1[5][2]]], :kdam=>[dict_kdamvals[ind_1[1][1]][:kdam][ind_1[1][2]],dict_kdamvals[ind_1[2][1]][:kdam][ind_1[2][2]],dict_kdamvals[ind_1[3][1]][:kdam][ind_1[3][2]],dict_kdamvals[ind_1[4][1]][:kdam][ind_1[4][2]],dict_kdamvals[ind_1[5][1]][:kdam][ind_1[5][2]]]),
                        "0.5" => Dict(:threshold=>[dict_kdamvals[ind_5[1][1]][:threshold][ind_5[1][2]],dict_kdamvals[ind_5[2][1]][:threshold][ind_5[2][2]],dict_kdamvals[ind_5[3][1]][:threshold][ind_5[3][2]],dict_kdamvals[ind_5[4][1]][:threshold][ind_5[4][2]]], :kdam=>[dict_kdamvals[ind_5[1][1]][:kdam][ind_5[1][2]],dict_kdamvals[ind_5[2][1]][:kdam][ind_5[2][2]],dict_kdamvals[ind_5[3][1]][:kdam][ind_5[3][2]],dict_kdamvals[ind_5[4][1]][:kdam][ind_5[4][2]]])
)

dict_titles_new = Dict("0.005" => [dict_titles[ind_005[1][1]][ind_005[1][2]], dict_titles[ind_005[2][1]][ind_005[2][2]], dict_titles[ind_005[3][1]][ind_005[3][2]]], 
                        "0.01" => [dict_titles[ind_01[1][1]][ind_01[1][2]], dict_titles[ind_01[2][1]][ind_01[2][2]], dict_titles[ind_01[3][1]][ind_01[3][2]]], 
                        "0.05" => [dict_titles[ind_05[1][1]][ind_05[1][2]], dict_titles[ind_05[2][1]][ind_05[2][2]], dict_titles[ind_05[3][1]][ind_05[3][2]]], 
                        "0.1" => [dict_titles[ind_1[1][1]][ind_1[1][2]], dict_titles[ind_1[2][1]][ind_1[2][2]], dict_titles[ind_1[3][1]][ind_1[3][2]], dict_titles[ind_1[4][1]][ind_1[4][2]], dict_titles[ind_1[5][1]][ind_1[5][2]]], 
                        "0.5" => [dict_titles[ind_5[1][1]][ind_5[1][2]], dict_titles[ind_5[2][1]][ind_5[2][2]], dict_titles[ind_5[3][1]][ind_5[3][2]], dict_titles[ind_5[4][1]][ind_5[4][2]]]
)

folders_dict_new = Dict("0.005" => "kdam_0.005", "0.01" => "kdam_0.01", "0.05" => "kdam_0.05", "0.1" => "kdam_0.1", "0.5" => "kdam_0.5")

dict_counts_new = Dict("0.005" => [dict_counts[ind_005[1][1]][ind_005[1][2]], dict_counts[ind_005[2][1]][ind_005[2][2]], dict_counts[ind_005[3][1]][ind_005[3][2]]], 
                        "0.01" => [dict_counts[ind_01[1][1]][ind_01[1][2]], dict_counts[ind_01[2][1]][ind_01[2][2]], dict_counts[ind_01[3][1]][ind_01[3][2]]], 
                        "0.05" => [dict_counts[ind_05[1][1]][ind_05[1][2]], dict_counts[ind_05[2][1]][ind_05[2][2]], dict_counts[ind_05[3][1]][ind_05[3][2]]], 
                        "0.1" => [dict_counts[ind_1[1][1]][ind_1[1][2]], dict_counts[ind_1[2][1]][ind_1[2][2]], dict_counts[ind_1[3][1]][ind_1[3][2]], dict_counts[ind_1[4][1]][ind_1[4][2]], dict_counts[ind_1[5][1]][ind_1[5][2]]], 
                        "0.5" => [dict_counts[ind_5[1][1]][ind_5[1][2]], dict_counts[ind_5[2][1]][ind_5[2][2]], dict_counts[ind_5[3][1]][ind_5[3][2]], dict_counts[ind_5[4][1]][ind_5[4][2]]]
)

dict_results_new = Dict("0.005" => [dict_results[ind_005[1][1]][ind_005[1][2]], dict_results[ind_005[2][1]][ind_005[2][2]], dict_results[ind_005[3][1]][ind_005[3][2]]], 
                        "0.01" => [dict_results[ind_01[1][1]][ind_01[1][2]], dict_results[ind_01[2][1]][ind_01[2][2]], dict_results[ind_01[3][1]][ind_01[3][2]]], 
                        "0.05" => [dict_results[ind_05[1][1]][ind_05[1][2]], dict_results[ind_05[2][1]][ind_05[2][2]], dict_results[ind_05[3][1]][ind_05[3][2]]], 
                        "0.1" => [dict_results[ind_1[1][1]][ind_1[1][2]], dict_results[ind_1[2][1]][ind_1[2][2]], dict_results[ind_1[3][1]][ind_1[3][2]], dict_results[ind_1[4][1]][ind_1[4][2]], dict_results[ind_1[5][1]][ind_1[5][2]]], 
                        "0.5" => [dict_results[ind_5[1][1]][ind_5[1][2]], dict_results[ind_5[2][1]][ind_5[2][2]], dict_results[ind_5[3][1]][ind_5[3][2]], dict_results[ind_5[4][1]][ind_5[4][2]]]
)

for i in eachindex(dict_reacts_new)
    plot_results("plot_stoch_reacts", dict_reacts_new[i], length(dict_kdamvals_new[i][:kdam]), folders_dict_new[i], xlabel="reaction", ylabel="count", titles=dict_titles_new[i], size=(1000,650), tosave=true)
    plot_totstochcount(dict_kdamvals_new[i][:kdam], dict_counts_new[i], "$(folders_dict_new[i])", folder=folders_dict_new[i])
end

specie = :rtca
for i in eachindex(dict_reacts_new)
    println(i)
    plot_results("plot_results", dict_results_new[i], length(dict_kdamvals_new[i][:kdam]), folders_dict_new[i], species=specie, xlabel="time", ylabel="$specie", titles=dict_titles_new[i], size=(1000,650), tosave=true);
    # dict_plot_hists[i,specie] = plot_results("plot_hists", dict_hists[i], length(dict_kdamvals[i][:kdam]), folders_dict[i], species=specie, xlabel="$specie", ylabel="frequency", titles=dict_titles[i], hidelabels=[true, true], linkaxes=true, size=(1000,650), tosave=true);
end
