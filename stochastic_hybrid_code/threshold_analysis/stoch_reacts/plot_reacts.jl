using StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed, TableOperations, JSON, Query, FindFirstFunctions, CategoricalArrays

# using PlotlyJS
using InteractiveViz, GLMakie

include(joinpath(homedir(), "phd/stochastic_hybrid_code/analysis_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_funcs.jl"))

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


# f = Figure();
# ax = Axis(f[1,1],yscale=log10);
# barplot!(ax, 1:length(collect(df_reacts[16].event)), df_reacts[16].count, bottom=1);
# display(GLMakie.Screen(), f)



# tot_counts = Int64[]
# for i in eachindex(df_reacts)
#     push!(tot_counts,sum(df_reacts[i].count))
# end


# f = Figure()
# ax = Axis(f[1,1],xlabel="threshold", ylabel="total stochastic reaction count")
# barplot!(threshold_vals, tot_counts)

# save("/Users/s2257179/phd/stochastic_hybrid_code/thresh_plots/total_counts.png", f)

