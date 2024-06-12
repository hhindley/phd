using StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed
# using PlotlyJS
using InteractiveViz, GLMakie

include(joinpath(homedir(), "phd/stochastic_hybrid_code/analysis_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))

@elapsed threshold_vals = range(10,310,length=20)
df_times = DataFrame(CSV.File(joinpath(homedir(), "Documents/stochastic_hybrid/thresh_times.csv"))) # times for kdam_test1

# scatter(rand(4))
# ilines(threshold_vals, df_times.time/60/60)
plot(scatter(x=df_times.threshold, y=df_times.time/60/60), Layout(xaxis_title="threshold", yaxis_title="time (hours)"))

loadandsort_arrow_file("/Users/s2257179/stoch_files/thresh_test_arrow/")

df_props, df_reacts, df_res = loadsort_all_arrow_files("/Users/s2257179/stoch_files/thresh_test_arrow/")

ilines(threshold_vals, df_res[1].rm_a)


@elapsed df_prop, df_react, df_results = loadandsort_arrow_file("/Users/s2257179/stoch_files/thresh_test_arrow/thresh_104.73684210526316.arrow")

@elapsed df_prop1, df_react1, df_results1 = new_func_arrow("/Users/s2257179/stoch_files/thresh_test_arrow/thresh_104.73684210526316.arrow")



@elapsed df = Arrow.Table("/Users/s2257179/stoch_files/thresh_test_arrow/thresh_10.0.arrow") |> DataFrame

df_r = get_reacts(df)
df_r.event = [parse.(Float64, subarray) for subarray in df_r.event]
df_react = combine(groupby(df_r, :event), nrow => :count)[2:12,:]


@elapsed df_r = get_reacts(df)
@elapsed df_p = get_prop_cols(df)

react_names = [:tscr_ab, :tscr_r, :tlr_a, :tlr_b, :tlr_r, :Vinflux, :Vdam, :Vtag, :Vrep, :deg_rd, :deg_rma, :deg_rmb, :deg_rmr, :V]

@elapsed props = str_to_arr.(df_p.event)

@elapsed df_prop = DataFrame(transpose(hcat(props...)), react_names)

@elapsed df_results = select!(df_p, Not(:event))


@elapsed df_r.event = [parse.(Float64, subarray) for subarray in df_r.event]

@elapsed df_react = combine(groupby(df_r, :event), nrow => :count)[2:12,:]
insertcols!(df_react, :reaction => react_names[Int64.(df_react.event)])






df_prop
df_react

df_results

plotBIG(df_results.time, df_results.rm_a)









@elapsed df_p = df[findall(row -> length(row[:event]) > 3, eachrow(df)), :]

@elapsed df_r = df[findall(row -> length(row[:event]) < 3 && row[:event][1] != 0, eachrow(df)), :]

react_names = [:tscr_ab, :tscr_r, :tlr_a, :tlr_b, :tlr_r, :Vinflux, :Vdam, :Vtag, :Vrep, :deg_rd, :deg_rma, :deg_rmb, :deg_rmr, :V]

@elapsed rep = replace.(df_p.event, r"[\[\]\(Any)]" => "")
@elapsed df_prop = split.(rep, ",")


# @elapsed df_prop = split.(replace.(df_p.event, r"[\[\]\(Any)]" => ""), ",")
@elapsed df_prop = permutedims(mapcols(x -> parse.(Float64, x), DataFrame(df_prop, :auto)))
rename!(df_prop, react_names)

df_results = select(df_p, Not(:event))
df_r.event = [parse.(Float64, subarray) for subarray in df_r.event]

df_react = combine(groupby(df_r, :event), nrow => :count)[2:12,:]
insertcols!(df_react, :reaction => react_names[Int64.(df_react.event)])

a = df_results.time
a = range(minimum(df_results.time), maximum(df_results.time), length=length(df_results.time))
ilines(a, df_results.rm_a)


function plotBIG(x, y)
    xvals = range(minimum(x), maximum(x), length=length(x))
    return ilines(xvals, y)
end



