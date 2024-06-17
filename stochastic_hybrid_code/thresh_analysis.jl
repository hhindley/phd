using StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed, TableOperations, JSON, Query
# using PlotlyJS
using InteractiveViz, GLMakie

include(joinpath(homedir(), "phd/stochastic_hybrid_code/analysis_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))

@elapsed df_props, df_reacts, df_res = loadsort_all_arrow_files("/Users/s2257179/stoch_files/thresh_test_arrow_files_14_06/")

function plotBIG(x, y)
    xvals = range(minimum(x), maximum(x), length=length(x))
    return ilines(xvals, y)
end

threshold_vals = range(10,310,length=20)
df_times = DataFrame(CSV.File(joinpath(homedir(), "/Users/s2257179/stoch_files/thresh_times.csv"))) # times for kdam_test1

ilines(threshold_vals, df_times.time/60/60)


plotBIG(df_results.time, df_results.rm_a)



loadsort_all_arrow_files("/Users/s2257179/stoch_files/thresh_test_arrow_files_14_06")



df_props = Arrow.Table("/Users/s2257179/phd/stochastic_hybrid_code/test/props/test.arrow") |> DataFrame
df_react = Arrow.Table("/Users/s2257179/phd/stochastic_hybrid_code/test/reacts/test.arrow") |> DataFrame
df_res = Arrow.Table("/Users/s2257179/phd/stochastic_hybrid_code/test/results/test.arrow") |> DataFrame

df = Arrow.Table("/Users/s2257179/phd/stochastic_hybrid_code/test/test.arrow") |> DataFrame