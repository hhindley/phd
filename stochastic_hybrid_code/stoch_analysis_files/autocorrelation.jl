using InteractiveViz, GLMakie, StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed, TableOperations, JSON, Query, FindFirstFunctions, CategoricalArrays, Colors

include(joinpath(homedir(), "phd/stochastic_hybrid_code/analysis_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_funcs.jl"))

mount_path, folders, folders_dict = load_file_structure("kdam_testing")
folders_dict = Dict(filter(pair -> pair.first in [7], folders_dict))
dict_times, dict_kdamvals, dict_titles, dict_results, dict_reacts, dict_props, dict_counts, dict_hists = load_data(mount_path, folders, folders_dict)

# autocorrelation
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