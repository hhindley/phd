using StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed, TableOperations, JSON, Query, FindFirstFunctions, CategoricalArrays

# using PlotlyJS
using InteractiveViz, WGLMakie

include(joinpath(homedir(), "phd/stochastic_hybrid_code/analysis_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_funcs.jl"))

df_reacts = load_files("/Users/s2257179/stoch_files/thresh_0107_noround_final_files/reacts")
df_reacts = [DataFrame(df_reacts[i]) for i in eachindex(df_reacts)]
for i in eachindex(df_reacts)
    for num in 1:13
        if num ∉ df_reacts[i].event
            push!(df_reacts[i], (event=num, count=0, reaction=react_names[num]))
        end
    end
end


df_times = CSV.File("/Users/s2257179/stoch_files/thresh_times_0107_noround.csv") |> DataFrame
lines(df_times.threshold, df_times.time/60/60)

threshold_vals1 = range(10,310,length=20)
threshold_vals2 = range(246,310,length=5)
threshold_vals0 = [threshold_vals1[1:15]; threshold_vals2]



# f = Figure();
# ax = Axis(f[1,1],yscale=log10);
# barplot!(ax, 1:length(collect(df_reacts[16].event)), df_reacts[16].count, bottom=1);
# display(GLMakie.Screen(), f)


tot_counts0 = Int64[]
for i in eachindex(df_reacts0)
    push!(tot_counts0,sum(df_reacts0[i].count))
end

tot_counts = Int64[]
for i in eachindex(df_reacts)
    push!(tot_counts,sum(df_reacts[i].count))
end

f = Figure()
ax = Axis(f[1,1],xlabel="threshold", ylabel="total stochastic reaction count")
barplot!(threshold_vals0, tot_counts0)
save("/Users/s2257179/phd/stochastic_hybrid_code/thresh_plots/total_counts.png", f)

f = Figure()
ax = Axis(f[1,1],xlabel="threshold", ylabel="total stochastic reaction count")
barplot!(threshold_vals1, tot_counts)
save("/Users/s2257179/phd/stochastic_hybrid_code/thresh_plots2/total_counts.png", f)

df_reacts
react_names_str = [string(i) for i in react_names]
f = Figure(size=(1450,900))
ax = Axis(f[1,1],xlabel="reactions", ylabel="reaction count", title="thresh_val", xticks=(1:14, react_names_str), xticklabelrotation=45)
barplot!(df_reacts[1].event, df_reacts[1].count)

titles = ["threshold: $(round(threshold_vals1[i], digits=2))" for i in eachindex(threshold_vals1)]
titles0 = ["threshold: $(round(threshold_vals1[i], digits=2))" for i in eachindex(threshold_vals1)]

f = plot_results("plot_stoch_reacts", df_reacts0, 5, 4, xlabel="reaction", ylabel="count", titles=titles, hidelabels=[true, false], linkaxes=false, species="stoch_reacts", folder="thresh_plots")
f = plot_results("plot_stoch_reacts", df_reacts, 5, 4, xlabel="reaction", ylabel="count", titles=titles, hidelabels=[true, false], linkaxes=false, species="stoch_reacts", folder="thresh_plots2")

display(GLMakie.Screen(), f.plot)


dfs = Dict{Symbol, DataFrame}()
for rn in react_names[1:end-1]
    df = build_reaction_count_df(df_reacts, rn, threshold_vals1)
    if !haskey(dfs, rn)
        dfs[rn] = df
    end
end

dfs0 = Dict{Symbol, DataFrame}()
for rn in react_names[1:end-1]
    df = build_reaction_count_df(df_reacts0, rn, threshold_vals0)
    if !haskey(dfs0, rn)
        dfs0[rn] = df
    end
end

f = Figure()
ax = Axis(f[1,1],xlabel="threshold", ylabel="count")
barplot!(dfs[:Vtag].threshold, dfs[:Vtag].count)
dfs

titles_reacts = ["$i" for i in react_names[1:end-1]]
push!(titles_reacts, "empty")
push!(titles_reacts, "empty")

f = create_subplots("plot_individual_reacts", 5, 3, titles=titles_reacts)
f = add_subplots(f, "plot_individual_reacts", dfs, 5, 3; species=:count, linkaxes=false)

p1 = plot_results("plot_individual_reacts", dfs0, 5, 3, titles=titles_reacts, xlabel="threshold", ylabel="count", hidelabels=[true, false], linkaxes=false, folder="thresh_plots", species="individual_reacts")
p2 = plot_results("plot_individual_reacts", dfs, 5, 3, titles=titles_reacts, xlabel="threshold", ylabel="count", hidelabels=[true, false], linkaxes=false, folder="thresh_plots2", species="individual_reacts")

display(GLMakie.Screen(), p1.plot)