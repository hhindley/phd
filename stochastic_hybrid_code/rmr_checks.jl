using InteractiveViz, GLMakie, StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed, TableOperations, JSON, Query, FindFirstFunctions, CategoricalArrays, Colors

include(joinpath(homedir(), "phd/stochastic_hybrid_code/analysis_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_funcs.jl"))

mount_path = "/Users/s2257179/stoch_files/kdam_testing/"
mount_path2 = "/Users/s2257179/stoch_files/hysteresis/"

all_items = readdir(mount_path)
folders = [item for item in all_items if isdir(joinpath(mount_path, item)) && !occursin("DS", item)]
folders_dict = Dict(i => folder for (i, folder) in enumerate(folders))
folders_dict = Dict(filter(pair -> pair.first in [1,2,3,4,5,6,8,9], folders_dict))

all_items2 = readdir(mount_path2)
folders2 = [item for item in all_items2 if isdir(joinpath(mount_path2, item)) && !occursin("DS", item)]
folders_dict2 = Dict(i => folder for (i, folder) in enumerate(folders2))

dict_times, dict_kdamvals, dict_titles, dict_results, dict_reacts, dict_props, dict_counts, dict_hists = setup_dicts(folders_dict)

for i in eachindex(folders_dict)
    println(i)
    dict_times[i], dict_kdamvals[i], dict_titles[i], dict_results[i], dict_reacts[i], dict_props[i] = LoadDataVars(folders[i]);
    dict_hists[i] = load_hist_files(joinpath(mount_path, folders_dict[i], "hists"))
    dict_counts[i] = prod_tot_count(dict_reacts[i])
end

dict_times2, dict_kdamvals2, dict_titles2, dict_results2, dict_reacts2, dict_props2, dict_counts2, dict_hists2 = setup_dicts(folders_dict2)

for i in eachindex(folders_dict2)
    println(i)
    dict_times2[i], dict_kdamvals2[i], dict_titles2[i], dict_results2[i], dict_reacts2[i], dict_props2[i] = LoadDataVars(folders2[i]);
    dict_hists2[i] = load_hist_files(joinpath(mount_path2, folders_dict2[i], "hists"))
    dict_counts2[i] = prod_tot_count(dict_reacts2[i])
end

function get_counts(dict_reacts, dict_kdamvals)
    txr_counts = []
    for i in eachindex(dict_reacts)
        for j in eachindex(dict_reacts[i])
            count = filter(row -> row.reaction == :tscr_r, dict_reacts[i][j]).count
            if :threshold in keys(dict_kdamvals[i])
                thresh = dict_kdamvals[i][:threshold][j]
            else 
                thresh = 150
            end
            push!(txr_counts, (count, dict_kdamvals[i][:kdam][j], thresh))   
        end
    end

    counts = [i[1] for i in txr_counts]
    counts = vcat(counts...)
    kdamvals = [i[2] for i in txr_counts]
    thresh = [i[3] for i in txr_counts]

    return counts, kdamvals, thresh
end

counts, kdamvals, thresh = get_counts(dict_reacts, dict_kdamvals)
counts2, kdamvals2, thresh2 = get_counts(dict_reacts2, dict_kdamvals2)

f = Figure()
ax = Axis(f[1, 1], xlabel="kdam", ylabel="tscr_r count")
plot!(ax, kdamvals, counts, color=thresh, colormap=:rainbow)
plot!(ax, kdamvals2, counts2, color=thresh2, colormap=:rainbow)
Colorbar(f[1,2], limits=(minimum(thresh), maximum(thresh)), colormap=:rainbow, label="threshold")