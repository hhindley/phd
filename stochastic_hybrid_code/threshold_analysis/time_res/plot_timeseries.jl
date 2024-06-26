using StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed, TableOperations, JSON, Query, FindFirstFunctions, CategoricalArrays

# using PlotlyJS
using InteractiveViz, GLMakie

include(joinpath(homedir(), "phd/stochastic_hybrid_code/analysis_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_funcs.jl"))

df_results = load_files("/Users/s2257179/stoch_files/thresh_test_arrow_files_14_06/results")

threshold_vals0 = range(10,310,length=20)
threshold_vals1 = range(246,310,length=5)
threshold_vals = [threshold_vals0[1:15]; threshold_vals1]

# times2 = DataFrame(CSV.File("/Users/s2257179/stoch_files/thresh_times_last5.csv"))

# plotBIG(times2.threshold, times2.time/60/60, "threshold", "time (hours)")

plot_subplots(df_results, :rm_a, threshold_vals);
plot_subplots(df_results, :rm_b, threshold_vals);
plot_subplots(df_results, :rm_r, threshold_vals);

plot_subplots(df_results, :rtca, threshold_vals)
plot_subplots(df_results, :rtcb, threshold_vals);
plot_subplots(df_results, :rtcr, threshold_vals);

plot_subplots(df_results, :rh, threshold_vals);
plot_subplots(df_results, :rd, threshold_vals);
plot_subplots(df_results, :rt, threshold_vals);