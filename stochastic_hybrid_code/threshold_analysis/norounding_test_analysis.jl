using StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed, TableOperations, JSON, Query, FindFirstFunctions, CategoricalArrays, Colors

# using PlotlyJS
using InteractiveViz, WGLMakie

include(joinpath(homedir(), "phd/stochastic_hybrid_code/analysis_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_funcs.jl"))

function create_df_reacts(df_reacts)
    df_reacts = [DataFrame(df_reacts[i]) for i in eachindex(df_reacts)]
    for i in eachindex(df_reacts)
        for num in 1:13
            if num ∉ df_reacts[i].event
                push!(df_reacts[i], (event=num, count=0, reaction=react_names[num]))
            end
        end
    end
    return df_reacts
end

df_reacts = load_files("/Users/s2257179/stoch_files/thresh_0107_noround_final_files/reacts")
df_reacts = create_df_reacts(df_reacts)

df_reacts_round = load_files("/Users/s2257179/stoch_files/thresh_0107_final_files/reacts")
df_reacts_round = create_df_reacts(df_reacts_round)



df_times = CSV.File("/Users/s2257179/stoch_files/thresh_times_0107_noround.csv") |> DataFrame
df_times_round = CSV.File("/Users/s2257179/stoch_files/thresh_times_0107.csv") |> DataFrame

f=Figure()
ax=Axis(f[1,1],xlabel="threshold", ylabel="time (hours)")
lines!(ax,df_times.threshold, df_times.time/60/60, label="no rounding")
lines!(ax,df_times_round.threshold, df_times_round.time/60/60, label="rounding")
axislegend(ax)

threshold_vals = df_times.threshold

function prod_tot_count(df_reacts)
    tot_counts = Int64[]
    for i in eachindex(df_reacts)
        push!(tot_counts,sum(df_reacts[i].count))
    end
    return tot_counts
end

tot_counts = prod_tot_count(df_reacts)
tot_counts_round = prod_tot_count(df_reacts_round)
colors=cgrad(:tab10)
combined_df = DataFrame(threshold=repeat(threshold_vals, 2),
                        counts=vcat(tot_counts, tot_counts_round),
                        group=repeat([1, 2], inner=length(threshold_vals)))
f = Figure()
ax = Axis(f[1,1],xlabel="threshold", ylabel="total stochastic reaction count")
barplot!(combined_df.threshold, combined_df.counts, dodge=combined_df.group, color=colors[combined_df.group])
Legend(f[1,2],[PolyElement(polycolor = colors[i]) for i in unique(combined_df.group)], ["no rounding", "rounding"])



df_reacts
f = Figure()
ax = Axis(f[1,1],xlabel="reactions", ylabel="reaction count", title="thresh_val", xticks=(1:13, react_names_str), xticklabelrotation=45)
barplot!(df_reacts[1].event, df_reacts[1].count)

titles = ["threshold: $(round(threshold_vals[i], digits=2))" for i in eachindex(threshold_vals)]

f = plot_results("plot_stoch_reacts", df_reacts, 5, xlabel="reaction", ylabel="count", titles=titles, hidelabels=[true, true], linkaxes=true, species="stoch_reacts")#, folder="thresh_plots")
display(f)
f1 = plot_results("plot_stoch_reacts", df_reacts_round, 5, xlabel="reaction", ylabel="count", titles=titles, hidelabels=[true, true], linkaxes=true, species="stoch_reacts")#, folder="thresh_plots")
# display(GLMakie.Screen(), f.plot)


dfs = Dict{Symbol, DataFrame}()
for rn in react_names[1:end-1]
    df = build_reaction_count_df(df_reacts, rn, threshold_vals)
    if !haskey(dfs, rn)
        dfs[rn] = df
    end
end

dfs_round = Dict{Symbol, DataFrame}()
for rn in react_names[1:end-1]
    df = build_reaction_count_df(df_reacts_round, rn, threshold_vals)
    if !haskey(dfs_round, rn)
        dfs_round[rn] = df
    end
end

f = Figure()
ax = Axis(f[1,1],xlabel="threshold", ylabel="count")
barplot!(dfs[:Vtag].threshold, dfs[:Vtag].count)
dfs

titles_reacts = ["$i" for i in react_names[1:end-1]]

p1 = plot_results("plot_individual_reacts", dfs, 12, titles=titles_reacts, xlabel="threshold", ylabel="count", hidelabels=[true, false], linkaxes=false, folder="thresh_plots", species="individual_reacts")
p2 = plot_results("plot_individual_reacts", dfs_round, 12, titles=titles_reacts, xlabel="threshold", ylabel="count", hidelabels=[true, false], linkaxes=false, folder="thresh_plots2", species="individual_reacts")




df_results = load_files("/Users/s2257179/stoch_files/thresh_0107_noround_final_files/results")
df_results_round = load_files("/Users/s2257179/stoch_files/thresh_0107_final_files/results")

f_rtca = plot_results("plot_results", df_results, 5, species=:rtca, xlabel="time", ylabel="rtca", titles=titles, folder="testing_figs")
f_rtca1 = plot_results("plot_results", df_results_round, 5, species=:rtca, xlabel="time", ylabel="rtca", titles=titles, folder="thresh_plots2")


df_props = load_files("/Users/s2257179/stoch_files/thresh_0107_noround_final_files/props", dataframe=false)
df_props_round = load_files("/Users/s2257179/stoch_files/thresh_0107_final_files/props", dataframe=false)

f1 = plot_prop(df_results, df_props, 1, "thresh_plots2/props", "threshold_$(threshold_vals[1])", threshold_vals, 31856.296174439733)
