using StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed, TableOperations, JSON, Query, FindFirstFunctions, CategoricalArrays

# using PlotlyJS
using InteractiveViz, GLMakie

include(joinpath(homedir(), "phd/stochastic_hybrid_code/analysis_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_funcs.jl"))

df_results0 = load_files("/Users/s2257179/stoch_files/thresh_test_arrow_files_14_06/results")
df_results = load_files("/Users/s2257179/stoch_files/thresh_1906_final_files/results")

threshold_vals1 = range(10,310,length=20)
threshold_vals2 = range(246,310,length=5)
threshold_vals0 = [threshold_vals1[1:15]; threshold_vals2]

times2 = DataFrame(CSV.File("/Users/s2257179/stoch_files/thresh_times_1906.csv"))

times_plot = plotBIG(times2.threshold, times2.time/60/60, xtitle="threshold", ytitle="time (hours)", title="")
save("/Users/s2257179/phd/stochastic_hybrid_code/thresh_plots2/run_times.png", times_plot)

# to plot one at a time
f = Figure()
ax = Axis(f[1,1], yscale=identity)
plot_timeres(get_column(df_results[15], :time), get_column(df_results[15], :rm_a), f[1,1])




titles = ["threshold: $(round(threshold_vals1[i], digits=2))" for i in eachindex(threshold_vals1)]

f_rma = plot_results("plot_results", df_results, :rm_a, 5, 4, xlabel="time", ylabel="rm_a", titles=titles, folder="thresh_plots2")
f_rmb = plot_results("plot_results", df_results, :rm_b, 5, 4, xlabel="time", ylabel="rm_b", titles=titles, folder="thresh_plots2")
f_rmr = plot_results("plot_results", df_results, :rm_r, 5, 4, xlabel="time", ylabel="rm_r", titles=titles, folder="thresh_plots2")

f_rtca = plot_results("plot_results", df_results, :rtca, 5, 4, xlabel="time", ylabel="rtca", titles=titles, folder="thresh_plots2")
f_rtcb = plot_results("plot_results", df_results, :rtcb, 5, 4, xlabel="time", ylabel="rtcb", titles=titles, folder="thresh_plots2")
f_rtcr = plot_results("plot_results", df_results, :rtcr, 5, 4, xlabel="time", ylabel="rtcr", titles=titles, folder="thresh_plots2")

f_rh = plot_results("plot_results", df_results, :rh, 5, 4, xlabel="time", ylabel="rh", titles=titles, folder="thresh_plots2")
f_rd = plot_results("plot_results", df_results, :rd, 5, 4, xlabel="time", ylabel="rd", titles=titles, folder="thresh_plots2")
f_rt = plot_results("plot_results", df_results, :rt, 5, 4, xlabel="time", ylabel="rt", titles=titles, folder="thresh_plots2")


display(GLMakie.Screen(), f_rtcb.plot)
display(GLMakie.Screen(), f_rtca.plot)

titles0 = ["threshold: $(round(threshold_vals1[i], digits=2))" for i in eachindex(threshold_vals1)]

f_rma0 = plot_results("plot_results", df_results0, :rm_a, 5, 4, xlabel="time", ylabel="rm_a", titles=titles0, folder="thresh_plots")
f_rmb0 = plot_results("plot_results", df_results0, :rm_b, 5, 4, xlabel="time", ylabel="rm_b", titles=titles0, folder="thresh_plots")
f_rmr0 = plot_results("plot_results", df_results0, :rm_r, 5, 4, xlabel="time", ylabel="rm_r", titles=titles0, folder="thresh_plots")

f_rtca0 = plot_results("plot_results", df_results0, :rtca, 5, 4, xlabel="time", ylabel="rtca", titles=titles0, folder="thresh_plots")
f_rtcb0 = plot_results("plot_results", df_results0, :rtcb, 5, 4, xlabel="time", ylabel="rtcb", titles=titles0, folder="thresh_plots")
f_rtcr0 = plot_results("plot_results", df_results0, :rtcr, 5, 4, xlabel="time", ylabel="rtcr", titles=titles0, folder="thresh_plots")

f_rh0 = plot_results("plot_results", df_results0, :rh, 5, 4, xlabel="time", ylabel="rh", titles=titles0, folder="thresh_plots")
f_rd0 = plot_results("plot_results", df_results0, :rd, 5, 4, xlabel="time", ylabel="rd", titles=titles0, folder="thresh_plots")
f_rt0 = plot_results("plot_results", df_results0, :rt, 5, 4, xlabel="time", ylabel="rt", titles=titles0, folder="thresh_plots")


