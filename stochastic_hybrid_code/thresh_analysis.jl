using StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed, TableOperations, JSON, Query
# using PlotlyJS
using InteractiveViz, GLMakie

include(joinpath(homedir(), "phd/stochastic_hybrid_code/analysis_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_funcs.jl"))

df_results = load_files("/Users/s2257179/stoch_files/thresh_test_arrow_files_14_06/results")
df_props = load_files("/Users/s2257179/stoch_files/thresh_test_arrow_files_14_06/props")
df_reacts = load_files("/Users/s2257179/stoch_files/thresh_test_arrow_files_14_06/reacts")
df_reacts = [DataFrame(df_reacts[i]) for i in eachindex(df_reacts)]
react_names=[:tscr_ab, :tscr_r, :tlr_a, :tlr_b, :tlr_r, :Vinflux, :Vdam, :Vtag, :Vrep, :deg_rd, :deg_rma, :deg_rmb, :deg_rmr, :V]
for i in eachindex(df_reacts)
    for num in 1:13
        if num ∉ df_reacts[i].event
            push!(df_reacts[i], (event=num, count=0, reaction=react_names[num]))
        end
    end
end

# for i in eachindex(df_reacts)
#     df_reacts[i].count .+= 1e-3  # Adjust counts for each dataframe
# end

threshold_vals0 = range(10,310,length=20)
threshold_vals1 = range(246,310,length=5)
threshold_vals = [threshold_vals0[1:15]; threshold_vals1]

# times2 = DataFrame(CSV.File("/Users/s2257179/stoch_files/thresh_times_last5.csv"))

# plotBIG(times2.threshold, times2.time/60/60, "threshold", "time (hours)")

plot_subplots(df_results, :rm_a, threshold_vals)
plot_subplots(df_results, :rm_b, threshold_vals)
plot_subplots(df_results, :rm_r, threshold_vals)

plot_subplots(df_results, :rtca, threshold_vals)
plot_subplots(df_results, :rtcb, threshold_vals)
plot_subplots(df_results, :rtcr, threshold_vals)

plot_subplots(df_results, :rh, threshold_vals)
plot_subplots(df_results, :rd, threshold_vals)
plot_subplots(df_results, :rt, threshold_vals)



using Colors
plotBIG(df_results[1].time, df_results[1].rtca)

df_reacts[1]

f = Figure()
ax = Axis(f[1,1],yscale=log10)
barplot!(ax, 1:length(collect(df_reacts[11].event)), df_reacts[11].count)
# ylims!(ax, 1, 10000000)

f = Figure()
ax = Axis(f[1,1],yscale=log10)
barplot!(ax, 1:length(collect(df_reacts[1].event)), df_reacts[1].count)

df_reacts[11]


plot_subplots_hists(df_reacts, threshold_vals)

tot_counts = Int64[]
for i in eachindex(df_reacts)
    push!(tot_counts,sum(df_reacts[i].count))
end


f = Figure()
ax = Axis(f[1,1],xlabel="threshold", ylabel="total stochastic reaction count")
barplot!(threshold_vals, tot_counts)