using StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase
using PlotlyJS
# using InteractiveViz, GLMakie

include("/home/hollie_hindley/Documents/phd/stochastic_hybrid_code/analysis_funcs.jl")

threshold_vals = range(10,310,length=20)
df_times = DataFrame(CSV.File("/home/hollie_hindley/Documents/stochastic_hybrid/thresh_times.csv")) # times for kdam_test1

plot(scatter(x=df_times.threshold, y=df_times.time/60/60))


dfs = [DataFrame(CSV.File(joinpath("/home/hollie_hindley/Documents/stochastic_hybrid/thresh_test", file), header=["event", "time", "rm_a", "rtca", "rm_b", "rtcb", "rm_r", "rtcr", "rh", "rd", "rt", "volume", "totprop"])) for file in readdir("/home/hollie_hindley/Documents/stochastic_hybrid/thresh_test")]

props, df_rs, df_ps = df_sort_all(dfs)

plot(scatter(x=df_ps[1].time, y=df_ps[1].rtca))



dfs1 = load_arrow_files("/home/hollie_hindley/Documents/stochastic_hybrid/thresh_test_arrow")


isequal(dfs, dfs1)