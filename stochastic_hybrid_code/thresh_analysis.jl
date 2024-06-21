using StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed, TableOperations, JSON, Query, FindFirstFunctions, CategoricalArrays

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


df = Arrow.Table("/Users/s2257179/phd/stochastic_hybrid_code/thresh_10.0.arrow") |> DataFrame

function get_grouped_df(data, specie)
    total_time = df.time[end]
    hist_df = DataFrame(t = data.time, s = data[:,specie])
    bin_edges = range(minimum(hist_df.s),maximum(hist_df.s)+0.01, length=51)

    hist_df.bin = cut(hist_df.s, bin_edges)
    hist_df.bin_index = levelcode.(hist_df.bin)
    hist_df.actual_index = 1:length(hist_df.s)

    grouped_df = groupby(hist_df, :bin_index)
    return grouped_df, hist_df
end

grouped_df, hist_df = get_grouped_df(df, :rh)

function get_time_in_state(grouped_df)
    in_sequence = false
    stops = []
    starts = []
    # println("looping through indexes")
    for i in 1:length(grouped_df.actual_index) - 1
        # println("checking if increasing by 1")
        if grouped_df.actual_index[i] + 1 == grouped_df.actual_index[i + 1]
            if !in_sequence
                # Sequence starts increasing by 1
                in_sequence = true
                # println("storing the start index")
                push!(starts, i)  # Store the start index
            end
        else
            if in_sequence
                # println("sequence stops increasing by 1")
                # Sequence stops increasing by 1
                in_sequence = false
                # println("storing the stop index")
                push!(stops, i)  # Store the stop index
            end
        end
    end

    # Handle the case where the sequence ends with the last element
    # println("sequence ends with last element case")
    if in_sequence
        push!(stops, length(grouped_df.actual_index))
    end

    # println("calculating total time in each state!")
    tot_time_in_state = 0
    for i in eachindex(stops)
        tot_time_in_state += grouped_df.t[stops[i]]-grouped_df.t[starts[i]]
    end

    println("finished!")
    return tot_time_in_state
end



tot_time = get_time_in_state(grouped_df[20])


all_histos=[]
for i in [:rm_a, :rm_b, :rm_r, :rtca, :rtcb, :rtcr, :rh, :rt, :rd]
    grouped_df, hist_df = get_grouped_df(df, i)
    hist_freq = DataFrame(bin = unique(hist_df.bin), freq = zeros(length(unique(hist_df.bin))))
    total_time = hist_df.t[end]
    for (df,i) in ProgressBar(zip(grouped_df, 1:length(grouped_df)))
        tot_time_in_state = get_time_in_state(df)
        hist_freq[i, :freq] = (length(df.t)*tot_time_in_state)/total_time
    end
    push!(all_histos, hist_freq)
end


hist_freq.bin = string.(hist_freq.bin)

hist_freq
f = Figure()
ax = Axis(f[1,1], yscale=log10)
barplot!(ax, bins, hist_freq.freq)

bin_edges = range(minimum(hist_df.s),maximum(hist_df.s)+0.01, length=51)
bins = range(bin_edges[1], bin_edges[end], length=50)