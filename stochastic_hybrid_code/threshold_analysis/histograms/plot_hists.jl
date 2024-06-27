using StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed, TableOperations, JSON, Query, FindFirstFunctions, CategoricalArrays

include(joinpath(homedir(), "phd/stochastic_hybrid_code/analysis_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_funcs.jl"))


# plotting hist_freq
# bins = [i[1] for i in hist_freq.bin]
# push!(bins, hist_freq.bin[end][2])
# bin_c = (bins[1:end-1] .+ bins[2:end]) ./ 2
# f = Figure()
# ax = Axis(f[1,1], yscale=log10)
# barplot!(ax, bin_c, hist_freq.freq, width=diff(bins), gap=0)

