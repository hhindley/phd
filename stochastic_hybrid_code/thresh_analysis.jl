using StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed, TableOperations, JSON, Query, FindFirstFunctions

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

plot_subplots(df_results, :rm_a, threshold_vals);
plot_subplots(df_results, :rm_b, threshold_vals);
plot_subplots(df_results, :rm_r, threshold_vals);

plot_subplots(df_results, :rtca, threshold_vals)
plot_subplots(df_results, :rtcb, threshold_vals);
plot_subplots(df_results, :rtcr, threshold_vals);

plot_subplots(df_results, :rh, threshold_vals);
plot_subplots(df_results, :rd, threshold_vals);
plot_subplots(df_results, :rt, threshold_vals);

plot_props(df_results);


# f = Figure();
# ax = Axis(f[1,1],yscale=log10);
# barplot!(ax, 1:length(collect(df_reacts[16].event)), df_reacts[16].count, bottom=1);
# display(GLMakie.Screen(), f)


plot_subplots_hists(df_reacts, threshold_vals)

# tot_counts = Int64[]
# for i in eachindex(df_reacts)
#     push!(tot_counts,sum(df_reacts[i].count))
# end


# f = Figure()
# ax = Axis(f[1,1],xlabel="threshold", ylabel="total stochastic reaction count")
# barplot!(threshold_vals, tot_counts)

# save("/Users/s2257179/phd/stochastic_hybrid_code/thresh_plots/total_counts.png", f)



# f = Figure()
# ax = Axis(f[1,1])
# for r in names(df_props[1])
#     iscatter!(ax, df_results[1].time, df_props[1][:,r])
# end
# ylims!(ax, 0,maximum([maximum(i) for i in eachcol(df_props[1])]))






hist_df = DataFrame(rtca = sort(df_results[1].rtca))
length(hist_df.rtca) ÷ 20

new_col = repeat(1:20, inner=length(hist_df.rtca) ÷ 20)
addons = repeat([20], inner=length(hist_df.rtca) % 20)

hist_df.group = vcat(new_col, addons)

hist_df

groupby(hist_df, :group)

hist_df.rtca

using CategoricalArrays
bin_edges = 0:1000:maximum(hist_df.rtca) + 10

maximum(hist_df.rtca)
bin_indices = cut(hist_df.rtca, bin_edges)

function custom_cut(values, edges)
    bins = zeros(Int, length(values))
    for (i, value) in enumerate(values)
        for (j, edge) in enumerate(edges)
            if j == length(edges) || (value >= edges[j] && value < edges[j+1])
                bins[i] = j
                break
            end
        end
    end
    return bins
end

# Apply the custom cut function to assign each rtca value to a bin
bin_indices = custom_cut(hist_df.rtca, bin_edges)

bin_counts = countmap(bin_indices)

bin_counts_full = [get(bin_counts, i, 0) for i in 1:length(bin_edges)-1]

fig = Figure()
ax = Axis(fig[1, 1], xlabel = "RTCA Range", ylabel = "Count", title = "RTCA Histogram")
barplot!(ax, bin_edges[1:end-1], bin_counts_full)

# Add the bin_indices as a column to your DataFrame
hist_df.bin = bin_indices

length(bin_edges)-1


total_time = df_results[1].time[end]
hist_df = DataFrame(t = df_results[1].time, s = df_results[1][:,:rtca])
bin_edges = range(minimum(hist_df.s),maximum(hist_df.s), length=51)
hist_times = hist_df.t
freqs = Float64[]
for bin in 1:(length(bin_edges)-1)
    println("$bin, $(bin+1)")
    num_in_range = get_bins(hist_df, bin_edges, bin)
    println("got num in range")
    time_in_state = get_timeInState(num_in_range, hist_df)
    println("got time in state")
    push!(freqs, (length(num_in_range.s)*tot_time_in_state)/total_time)
end

function get_bins(hist_df, bin_edges, bin)
    lb = bin_edges[bin]
    ub = bin_edges[bin+1]
    # num_in_range = filter(x -> lb <= x.s < ub, hist_df)
    allowmissing!(hist_df)
    num_in_range = subset(hist_df, AsTable(:s) => (@. x-> lb <= x.s < ub))
    return num_in_range
end

function get_timeInState(num_in_range, hist_df)
    tot_time_in_state = 0
    hist_times = hist_df.t  
    for row in eachrow(num_in_range)
        print(row)
        index = FindFirstFunctions.findfirstequal(row.t, hist_times)
        if index != 1
            prev_time = hist_times[index-1]
            time_in_state = row.t - prev_time
            tot_time_in_state += time_in_state
        end
    end
    return tot_time_in_state
end

num_in_range = get_bins(hist_df, bin_edges, 1)
@elapsed tot = get_timeInState(num_in_range, hist_df)

