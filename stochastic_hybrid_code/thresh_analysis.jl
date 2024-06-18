using StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed, TableOperations, JSON, Query
# using PlotlyJS
# using InteractiveViz, GLMakie

include(joinpath(homedir(), "phd/stochastic_hybrid_code/analysis_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))

# @elapsed df_props, df_reacts, df_res = loadsort_all_arrow_files("/Users/s2257179/stoch_files/thresh_test_arrow_files_14_06/")

# function plotBIG(x, y)
#     xvals = range(minimum(x), maximum(x), length=length(x))
#     return ilines(xvals, y)
# end

# threshold_vals = range(10,310,length=20)
# df_times = DataFrame(CSV.File(joinpath(homedir(), "/Users/s2257179/stoch_files/thresh_times.csv"))) # times for kdam_test1

# ilines(threshold_vals, df_times.time/60/60)


# plotBIG(df_results.time, df_results.rm_a)



# loadsort_all_arrow_files("/home/hollie_hindley/Documents/stochastic_hybrid/thresh_test_arrow_files_14_06")



# arrow_conv("/home/hollie_hindley/Documents/stochastic_hybrid/test", "/home/hollie_hindley/Documents/stochastic_hybrid/test_arrow")

# df = CSV.File("/home/hollie_hindley/Documents/stochastic_hybrid/thresh_test_last5/thresh_246.0.dat") |> DataFrame

# for chunk in CSV.Chunks("/home/hollie_hindley/Documents/stochastic_hybrid/test/test.dat"; header=["event", "time", "rm_a", "rtca", "rm_b", "rtcb", "rm_r", "rtcr", "rh", "rd", "rt", "volume", "totprop", "missing"])
#     chunk.event = Array{Float64}.(JSON.parse.(chunk.event))
#     # Arrow.Writer("/home/hollie_hindley/Documents/stochastic_hybrid/test/test_arrow_1.arrow", chunk)
# end

# # df_props = Arrow.Table("/Users/s2257179/phd/stochastic_hybrid_code/test/props/test.arrow") |> DataFrame
# # df_react = Arrow.Table("/Users/s2257179/phd/stochastic_hybrid_code/test/reacts/test.arrow") |> DataFrame
# # df_res = Arrow.Table("/Users/s2257179/phd/stochastic_hybrid_code/test/results/test.arrow") |> DataFrame

# # df = Arrow.Table("/Users/s2257179/phd/stochastic_hybrid_code/test/test.arrow") |> DataFrame
# open(Arrow.Writer, "/home/hollie_hindley/Documents/stochastic_hybrid/test/test_arrow.arrow") do writer
# for chunk in CSV.Chunks("/home/hollie_hindley/Documents/stochastic_hybrid/test/test.dat"; header=["event", "time", "rm_a", "rtca", "rm_b", "rtcb", "rm_r", "rtcr", "rh", "rd", "rt", "volume", "totprop", "missing"], drop=["missing"])
#     df_chunk = DataFrame(chunk) 
#     df_chunk.event = Array{Float64}.(JSON.parse.(df_chunk.event))
#     Arrow.write(writer, df_chunk)
# end
# end


# open(Arrow.Writer, joinpath(arrow_folder_path, splitext(basename(file))[1] * ".arrow")) do writer
#     for chunk in CSV.Chunks(joinpath(folder_path, file); header=["event", "time", "rm_a", "rtca", "rm_b", "rtcb", "rm_r", "rtcr", "rh", "rd", "rt", "volume", "totprop", "missing"], drop=["missing"])
#         df_chunk = DataFrame(chunk) 
#         df_chunk.event = Array{Float64}.(JSON.parse.(df_chunk.event))
#         Arrow.write(writer, df_chunk)
#     end
# end

# arrow_conv("/home/hollie_hindley/Documents/stochastic_hybrid/test", "/home/hollie_hindley/Documents/stochastic_hybrid/test_arrow")






arrow_conv("/home/hollie_hindley/Documents/stochastic_hybrid/thresh_test_last5", "/home/hollie_hindley/Documents/stochastic_hybrid/thresh_test_last5_arrow")

# df_res = Arrow.Table("/home/hollie_hindley/Documents/stochastic_hybrid/test_arrow/results/test.arrow") |> DataFrame
# df_props = Arrow.Table("/home/hollie_hindley/Documents/stochastic_hybrid/test_arrow/props/test.arrow") |> DataFrame
# df_reacts = Arrow.Table("/home/hollie_hindley/Documents/stochastic_hybrid/test_arrow/reacts/test.arrow") |> DataFrame

# arrow_conv("/home/hollie_hindley/Documents/stochastic_hybrid/test", "/home/hollie_hindley/Documents/stochastic_hybrid/test_arrow")
# df1, df2, df3 = loadandsort_arrow_file("/home/hollie_hindley/Documents/stochastic_hybrid/test_arrow/test.arrow")

# isequal(df_reacts, df2)

# df_reacts
# df2

# df_reacts = DataFrame[]
# react_names=[:tscr_ab, :tscr_r, :tlr_a, :tlr_b, :tlr_r, :Vinflux, :Vdam, :Vtag, :Vrep, :deg_rd, :deg_rma, :deg_rmb, :deg_rmr, :V]
# for chunk in CSV.Chunks("/home/hollie_hindley/Documents/stochastic_hybrid/test/test.dat"; header=["event", "time", "rm_a", "rtca", "rm_b", "rtcb", "rm_r", "rtcr", "rh", "rd", "rt", "volume", "totprop", "missing"], drop=["missing"])
#     df_chunk = DataFrame(chunk) 
#     df_chunk.event = Array{Float64}.(JSON.parse.(df_chunk.event))

#     dfr_chunk = filter(row -> length(row.event) < 3 && row.event[1] != 0, df_chunk)
#     dfp_chunk = filter(row -> length(row.event) > 3, df_chunk)

#     df_prop_chunk = DataFrame(transpose(hcat(dfp_chunk.event...)), react_names)
#     df_react_chunk = combine(groupby(dfr_chunk, :event), nrow => :count)
#     if !isempty(df_react_chunk)
#         df_react_chunk.event = reduce(vcat, map(v -> round.(Int64, v), collect.(df_react_chunk.event)))
#         insertcols!(df_react_chunk, :reaction => react_names[df_react_chunk.event])
#         push!(df_reacts, df_react_chunk)
#     end
    
# end
# df_react = combine(groupby(vcat(df_reacts...), [:reaction, :event]), :count => sum => :count)

