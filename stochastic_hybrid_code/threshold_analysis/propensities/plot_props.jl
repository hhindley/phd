using StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed, TableOperations, JSON, Query, FindFirstFunctions, CategoricalArrays

# using PlotlyJS
using InteractiveViz, GLMakie

include(joinpath(homedir(), "phd/stochastic_hybrid_code/analysis_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_funcs.jl"))

df_props = load_files("/Users/s2257179/stoch_files/thresh_test_arrow_files_14_06/props")


plot_props(df_results);

# f = Figure()
# ax = Axis(f[1,1])
# for r in names(df_props[1])
#     iscatter!(ax, df_results[1].time, df_props[1][:,r])
# end
# ylims!(ax, 0,maximum([maximum(i) for i in eachcol(df_props[1])]))



