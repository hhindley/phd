using StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed, TableOperations, JSON, Query, FindFirstFunctions, CategoricalArrays, Colors

# using PlotlyJS
using InteractiveViz, WGLMakie

include(joinpath(homedir(), "phd/stochastic_hybrid_code/analysis_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_funcs.jl"))

mount_path = "/Users/s2257179/stoch_files"
folder = "run_individually_0407_final_files"
hist_path = joinpath(mount_path, "hists_run_individually_0407")
filepath = joinpath(mount_path, folder)

df_reacts = load_files(joinpath(filepath, "reacts"))
df_reacts = create_df_reacts(df_reacts)


df_times = CSV.File("/Users/s2257179/stoch_files/test_0307_5_values_from_original_range_times.csv") |> DataFrame

f=Figure()
ax=Axis(f[1,1],xlabel="threshold", ylabel="time (hours)")
lines!(ax,df_times.threshold, df_times.time/60/60, label="no rounding")

threshold_vals = df_times.threshold

function prod_tot_count(df_reacts)
    tot_counts = Int64[]
    for i in eachindex(df_reacts)
        push!(tot_counts,sum(df_reacts[i].count))
    end
    return tot_counts
end

tot_counts = prod_tot_count(df_reacts)


f = Figure()
ax = Axis(f[1,1],xlabel="threshold", ylabel="total stochastic reaction count")
barplot!(threshold_vals, tot_counts)



df_reacts
f = Figure()
ax = Axis(f[1,1],xlabel="reactions", ylabel="reaction count", title="thresh_val", xticks=(1:13, react_names_str), xticklabelrotation=45)
barplot!(df_reacts[1].event, df_reacts[1].count)

titles = ["threshold: $(round(threshold_vals[i], digits=2))" for i in eachindex(threshold_vals)]

f = plot_results("plot_stoch_reacts", df_reacts, 5, xlabel="reaction", ylabel="count", titles=titles, hidelabels=[true, true], linkaxes=true, species="stoch_reacts")#, folder="thresh_plots")
display(f)


dfs = Dict{Symbol, DataFrame}()
for rn in react_names[1:end-1]
    df = build_reaction_count_df(df_reacts, rn, threshold_vals)
    if !haskey(dfs, rn)
        dfs[rn] = df
    end
end

titles_reacts = ["$i" for i in react_names[1:end-1]]

p1 = plot_results("plot_individual_reacts", dfs, 12, titles=titles_reacts, xlabel="threshold", ylabel="count", hidelabels=[true, false], linkaxes=false, folder="thresh_plots", species="individual_reacts")




df_results = load_files(joinpath(filepath, "results"))

f_rtca = plot_results("plot_results", df_results, 5, species=:rtca, xlabel="time", ylabel="rtca", titles=titles, folder="testing_figs")


df_props = load_files(joinpath(filepath, "props"), dataframe=false)

f1 = plot_prop(df_results, df_props, 1, "thresh_plots2/props", "threshold_$(threshold_vals[1])", threshold_vals, 31856.296174439733)
f2 = plot_prop(df_results, df_props, 2, "thresh_plots2/props", "threshold_$(threshold_vals[2])", threshold_vals, 31856.296174439733)
f3 = plot_prop(df_results, df_props, 3, "thresh_plots2/props", "threshold_$(threshold_vals[3])", threshold_vals, 31856.296174439733)
f4 = plot_prop(df_results, df_props, 4, "thresh_plots2/props", "threshold_$(threshold_vals[4])", threshold_vals, 31856.296174439733)
f5 = plot_prop(df_results, df_props, 5, "thresh_plots2/props", "threshold_$(threshold_vals[5])", threshold_vals, 31856.296174439733)

display(f1)
display(f2)


f = plot_props(df_results, df_props, 5, threshold_vals)


dfs = load_hist_files(hist_path)

titles = ["threshold: $(round(threshold_vals[i], digits=2))" for i in eachindex(threshold_vals)]


rma_hist = plot_results("plot_hists", dfs, 5, species=:rh, xlabel="rh conc", ylabel="frequency", titles=titles)#, folder="thresh_plots2/hists")

