using StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed
using PlotlyJS
# using InteractiveViz, WGLMakie

include(joinpath(homedir(), "phd/stochastic_hybrid_code/analysis_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))

threshold_vals = range(10,310,length=20)
df_times = DataFrame(CSV.File(joinpath(homedir(), "Documents/stochastic_hybrid/thresh_times.csv"))) # times for kdam_test1

# scatter(rand(4))
# ilines(threshold_vals, df_times.time/60/60)
plot(scatter(x=df_times.threshold, y=df_times.time/60/60), Layout(xaxis_title="threshold", yaxis_title="time (hours)"))


df_props, df_reacts, df_res = loadsort_all_arrow_files("/home/hollie_hindley/Documents/stochastic_hybrid/thresh_test_arrow")

# df_props

# @elapsed props, df_rs, df_ps = df_sort_all(dfs)

# plot(scatter(x=df_ps[1].time, y=df_ps[1].rtca))


# df_props, df_r, df_p = df_sort(dfs[1])

# df_p = filter!(row -> length(row[:event]) > 1, df)

# df = Arrow.Table("/home/hollie_hindley/Documents/stochastic_hybrid/thresh_test_arrow/thresh_10.0.arrow") |> DataFrame
# ls = [length(i) for i in df.event]
# unique(ls)
# if 100 in ls
#     println("yes")
# end
# sort(unique(ls))

# length(df.event)
# df_p = Arrow.Table("/home/hollie_hindley/Documents/stochastic_hybrid/thresh_test_arrow/thresh_10.0.arrow") |> DataFrame |> filter(row -> length(row[:event]) > 3)
# df_r = Arrow.Table("/home/hollie_hindley/Documents/stochastic_hybrid/thresh_test_arrow/thresh_10.0.arrow") |> DataFrame |> filter(row -> length(row[:event]) < 3 && row[:event][1] != 0)

# props = split.(replace.(df_p.event, r"[\[\]\(Any)]" => ""), ",")
# props = permutedims(mapcols(x -> parse.(Float64, x), DataFrame(props, :auto)))
# rename!(props, react_names)


# df_p = select(df_p, Not(:event))
# df_r.event = [parse.(Float64, subarray) for subarray in df_r.event]



# files = readdir("/home/hollie_hindley/Documents/stochastic_hybrid/thresh_test_arrow")[1:5]
# files = files[1:5]
# @distributed for (i, file) in collect(enumerate(files))
#     println(i, file)
# end


# props = df_p.event
# props = [split(replace(i, r"[\[\]\(Any)]" => ""), ",") for i in props]
# props = [parse.(Float64, subarray) for subarray in props]


# df_props = DataFrame([name => Vector{Float64}(undef, length(props)) for name in react_names])

# for (i, prop) in enumerate(props)
#     df_props[i, :] = prop
# end

# df_props.V


# df_r


# df_react = combine(groupby(df_r, :event), nrow => :count)[2:12,:]
# insertcols!(df_react, :reaction => react_names[Int64.(df_react.event)])






# @elapsed df_prop, df_react, df_res = loadandsort_arrow_file("/home/hollie_hindley/Documents/stochastic_hybrid/thresh_test_arrow/thresh_10.0.arrow")



# df1, df2 = Arrow.Table("/home/hollie_hindley/Documents/stochastic_hybrid/thresh_test_arrow/thresh_10.0.arrow") |> DataFrame |> filter(row -> length(row[:event]) > 3)
# df_r = Arrow.Table("/home/hollie_hindley/Documents/stochastic_hybrid/thresh_test_arrow/thresh_10.0.arrow") |> DataFrame |> filter(row -> length(row[:event]) < 3 && row[:event][1] != 0)

# df1


# df2

# df = Arrow.Table("/home/hollie_hindley/Documents/stochastic_hybrid/thresh_test_arrow/thresh_10.0.arrow") |> DataFrame

# df1 = filter(row -> length(row[:event]) > 3, df)
# df2 = filter(row -> length(row[:event]) < 3 && row[:event][1] != 0, df)