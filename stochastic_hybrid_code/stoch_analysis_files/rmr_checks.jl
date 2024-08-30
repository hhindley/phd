using InteractiveViz, GLMakie, StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed, TableOperations, JSON, Query, FindFirstFunctions, CategoricalArrays, Colors

include(joinpath(homedir(), "phd/stochastic_hybrid_code/analysis_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_funcs.jl"))


mount_path, folders, folders_dict = load_file_structure("kdam_testing")
folders_dict = Dict(filter(pair -> pair.first in [7], folders_dict))
dict_times, dict_kdamvals, dict_titles, dict_results, dict_reacts, dict_props, dict_counts, dict_hists = load_data(mount_path, folders, folders_dict)

mount_path, folders, folders_dict = load_file_structure("hysteresis")
folders_dict = Dict(filter(pair -> pair.first in [7], folders_dict))
dict_times, dict_kdamvals, dict_titles, dict_results, dict_reacts, dict_props, dict_counts, dict_hists = load_data(mount_path, folders, folders_dict)


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