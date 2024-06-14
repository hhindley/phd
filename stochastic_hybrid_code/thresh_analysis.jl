using StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed, TableOperations, JSON
# using PlotlyJS
using InteractiveViz, GLMakie

include(joinpath(homedir(), "phd/stochastic_hybrid_code/analysis_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))

function plotBIG(x, y)
    xvals = range(minimum(x), maximum(x), length=length(x))
    return ilines(xvals, y)
end
homedir()
@elapsed threshold_vals = range(10,310,length=20)
df_times = DataFrame(CSV.File(joinpath(homedir(), "/Users/s2257179/stoch_files/thresh_times.csv"))) # times for kdam_test1

# scatter(rand(4))
ilines(threshold_vals, df_times.time/60/60)
plot(scatter(x=df_times.threshold, y=df_times.time/60/60), Layout(xaxis_title="threshold", yaxis_title="time (hours)"))

loadandsort_arrow_file("/Users/s2257179/stoch_files/thresh_test_arrow/")

df_props, df_reacts, df_res = loadsort_all_arrow_files("/Users/s2257179/stoch_files/thresh_test_arrow/")

ilines(threshold_vals, df_res[1].rm_a)


@elapsed df_prop, df_react, df_results = loadandsort_arrow_file("/Users/s2257179/stoch_files/thresh_test_arrow/thresh_262.63157894736844.arrow")

@elapsed df_prop, df_react, df_results = loadandsort_arrow_file("/Users/s2257179/stoch_files/thresh_test_arrow/thresh_104.73684210526316.arrow")

plotBIG(df_results.time, df_results.rm_a)




@elapsed df = Arrow.Table("/Users/s2257179/stoch_files/thresh_test_arrow/thresh_10.0.arrow") |> DataFrame

@elapsed df_r = @view df[findall(row -> length(row[:event]) < 3 && row[:event][1] != 0, eachrow(df)), :]
@elapsed df_p = @view df[findall(row -> length(row[:event]) > 3, eachrow(df)), :]
@time react_names = [:tscr_ab, :tscr_r, :tlr_a, :tlr_b, :tlr_r, :Vinflux, :Vdam, :Vtag, :Vrep, :deg_rd, :deg_rma, :deg_rmb, :deg_rmr, :V]
@time props = str_to_arr.(df_p.event)
@time df_prop = DataFrame(transpose(hcat(props...)), react_names)
@time df_r.event = [parse.(Float64, subarray) for subarray in df_r.event]
@time df_react = combine(groupby(df_r, :event), nrow => :count)[2:12,:]
@time insertcols!(df_react, :reaction => react_names[Int64.(df_react.event)])


@time df_p = get_prop_cols(df)


@elapsed df = Arrow.Table("/Users/s2257179/stoch_files/thresh_test_arrow/thresh_10.0.arrow") |> DataFrame
@elapsed df_p = @view df[findall(row -> length(row[:event]) > 3, eachrow(df)), :]


 
react_names = [:tscr_ab, :tscr_r, :tlr_a, :tlr_b, :tlr_r, :Vinflux, :Vdam, :Vtag, :Vrep, :deg_rd, :deg_rma, :deg_rmb, :deg_rmr, :V]

@elapsed atab = Arrow.Table("/Users/s2257179/phd/stochastic_hybrid_code/test.arrow")
@elapsed df_p = atab |> TableOperations.filter(x -> length(x.event) > 3) |> DataFrame
@elapsed df_r = atab |> TableOperations.filter(x -> length(x.event) < 3 && x.event[1] !=0) |> DataFrame
@elapsed df_prop = DataFrame(transpose(hcat(df_p.event...)), react_names)
@elapsed df_react = combine(groupby(df_r, :event), nrow => :count)[2:11,:]
@elapsed df_react.event = reduce(vcat, map(v -> round.(Int64, v), collect.(df_react.event)))
insertcols!(df_react, :reaction => react_names[df_react.event])


@elapsed loadandsort_arrow_file("/Users/s2257179/phd/stochastic_hybrid_code/test.arrow")