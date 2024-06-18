using StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed, TableOperations, JSON, Query
# using PlotlyJS
using InteractiveViz, GLMakie

include(joinpath(homedir(), "phd/stochastic_hybrid_code/analysis_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_funcs.jl"))

df_results = load_files("/Users/s2257179/stoch_files/thresh_test_arrow_files_14_06/results")
df_props = load_files("/Users/s2257179/stoch_files/thresh_test_arrow_files_14_06/props")
df_reacts = load_files("/Users/s2257179/stoch_files/thresh_test_arrow_files_14_06/reacts")

threshold_vals0 = range(10,310,length=20)
threshold_vals1 = range(246,310,length=5)
threshold_vals = [threshold_vals0[1:15]; threshold_vals1]

times2 = DataFrame(CSV.File("/Users/s2257179/stoch_files/thresh_times_last5.csv"))

plotBIG(times2.threshold, times2.time/60/60, "threshold", "time (hours)")

for i in 1:length(df_results)
    display(GLMakie.Screen(), plotBIG(df_results[i].time, df_results[i].rm_a, xtitle="time", ytitle="rm_a", title="threshold = $(threshold_vals[i])"))
end

plotBIG(df_results[1].time, df_results[1].rtca, xtitle="time", ytitle="rtca")
# display(GLMakie.Screen(), plotBIG(df_results[1].time, df_results[1].rm_a, xtitle="time", ytitle="rtca"))

df_results[1]

plot_subplots(df_results, :rm_a, threshold_vals)

f = Figure()
for i in 1:length(df_results)
    xvals = range(minimum(data[i].time), maximum(data[i].time), length=length(data[i].time))
    ax = Axis(f[1, i], xlabel = "time", ylabel = "$species", title="$(threshold_vals[i])")
    ilines!(f[1,i], xvals, data[i][:,species])
end

f = Figure()
xvals = range(minimum(df_results[1].time), maximum(df_results[1].time), length=length(df_results[1].time))
ax1 = Axis(f[1, 1], xlabel = "", ylabel = "", title="")
ax2 = Axis(f[2,1], xlabel = "", ylabel = "", title="")
hidexdecorations!(ax1)
linkaxes!(ax1, ax2)

ilines!(f[1,1], xvals, df_results[1][:,:rm_a])
ilines!(f[2,1], range(minimum(df_results[2].time), maximum(df_results[2].time), length=length(df_results[2].time)), df_results[2][:,:rm_a])

plotBIG(df_results[1].time, df_results[1][:,:rm_a])



f = Figure()
total_rows = 5
total_columns = 4

for j in 1:total_columns
    for i in 1:total_rows
        data_ind = i + total_rows * (j - 1)
        println(data_ind)
        ax = Axis(f[i, j], xlabel = "time", ylabel = "rm_a", title="threshold: $(threshold_vals[data_ind])")
        ilines!(f[i,j], makiex(df_results[data_ind].time), df_results[data_ind][:,:rm_a])

        # Hide x-axis decorations for axes not in the bottom row
        if i != total_rows
            hidexdecorations!(ax, grid=false)
        end

        # Hide y-axis decorations for axes not in the first column
        if j > 1
            hideydecorations!(ax, grid=false)
        end
    end
end

# delete the last threshold results from the source data 


plot_subplots(df_results, :rm_a, threshold_vals)




