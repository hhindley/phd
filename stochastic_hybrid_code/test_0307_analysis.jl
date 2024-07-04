using StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed, TableOperations, JSON, Query, FindFirstFunctions, CategoricalArrays, Colors

# using PlotlyJS
using InteractiveViz, WGLMakie

include(joinpath(homedir(), "phd/stochastic_hybrid_code/analysis_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_funcs.jl"))

mount_path = "/Users/s2257179/stoch_files"
folder = "test_0307_5_values_from_original_range_final_files"
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












df = Arrow.Table("/Users/s2257179/stoch_files/test_0307_5_values_from_original_range_final_files/results/thresh_183.68421052631578.arrow") |> DataFrame

hist_freq = DataFrame(CSV.File("/Users/s2257179/Desktop/test_rh.csv"))
hist_freq1 = DataFrame(CSV.File("/Users/s2257179/Desktop/rh_NEW.csv"))

hist_freq.bin = Array{Float64}.(JSON.parse.(hist_freq.bin))
hist_freq1.bin = Array{Float64}.(JSON.parse.(hist_freq1.bin))



grouped_df, hist_df, bin_edges = get_grouped_df(df, :rh)
counts=[]
for df in grouped_df
    push!(counts, length(df.t))
end

counts


f = Figure(size=(1000,650))
ax = Axis(f[1,1], title="my data")
mydat = plot_hist(hist_freq, f[1,1])
ax = Axis(f[1,2], title="my data2")
mydat = plot_hist(hist_freq1, f[1,2])


ax = Axis(f[1,2], title="makie")
f_lib = hist!(f[1,2],df.rh, bins=50)

bins = [i[1] for i in hist_freq.bin]
push!(bins, hist_freq.bin[end][2])
bin_c = (bins[1:end-1] .+ bins[2:end]) ./ 2
ax = Axis(f[1,3], title="my data counts")
barplot!(f[1,3], bin_c, counts, width=diff(bins), gap=0)

linkaxes!(f.content...)



