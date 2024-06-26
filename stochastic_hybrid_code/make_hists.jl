using StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed, TableOperations, JSON, Query, FindFirstFunctions, CategoricalArrays

include(joinpath(homedir(), "phd/stochastic_hybrid_code/analysis_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_funcs.jl"))


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
function get_grouped_df(data, specie)
    hist_df = DataFrame(t = data.time, s = data[:,specie])
    bin_edges = range(minimum(hist_df.s),maximum(hist_df.s)+0.01, length=51)

    hist_df.bin = cut(hist_df.s, bin_edges)
    hist_df.bin_index = levelcode.(hist_df.bin)
    hist_df.actual_index = 1:length(hist_df.s)

    all_bins_df = DataFrame(bin_index = 1:50)
    full_hist_df = leftjoin(all_bins_df, hist_df, on = :bin_index)
        
    grouped_df = groupby(full_hist_df, :bin_index)
    # grouped_df = groupby(hist_df, :bin_index)
    return grouped_df, hist_df, bin_edges
end
function make_hist_freq(bin_edges)
    full_bins=Array[]
    for i in eachindex(collect(bin_edges))
        if i != 51
            push!(full_bins, [collect(bin_edges)[i], collect(bin_edges)[i+1]])
        end
    end

    hist_freq = DataFrame([full_bins], [:bin])
    z1 = DataFrame(freq = zeros(1:50))
    
    return hcat(hist_freq, z1)
end

# grouped_df, hist_df, bin_edges = get_grouped_df(df, :rh)
# hist_freq = make_hist_freq(bin_edges)
# total_time = hist_df.t[end]
# for (df,i) in ProgressBar(zip(grouped_df, 1:length(grouped_df)))
#     tot_time_in_state = get_time_in_state(df)
#     hist_freq[i, :freq] = (length(df.t)*tot_time_in_state)/total_time
# end
# CSV.write("/Users/s2257179/phd/stochastic_hybrid_code/hists/thresh10/rh.csv", hist_freq)

folderpath = "/home/hollie_hindley/Documents/stochastic_hybrid/thresh_test_arrow_files_14_06/results"
files = readdir(folderpath)[1:2]
for file in files
    filepath = joinpath(folderpath, file)
    res = Arrow.Table(filepath) |> DataFrame
    if !isdir("/home/hollie_hindley/Documents/stochastic_hybrid/hists/$(file[1:end-6])")
        mkdir("/home/hollie_hindley/Documents/stochastic_hybrid/hists/$(file[1:end-6])")
    end
    for i in [:rm_a, :rm_b, :rm_r, :rtca, :rtcb, :rtcr, :rh, :rt, :rd]
        grouped_df, hist_df, bin_edges = get_grouped_df(res, i)
        hist_freq = make_hist_freq(bin_edges)
        total_time = hist_df.t[end]
        for (df,i) in ProgressBar(zip(grouped_df, 1:length(grouped_df)))
            tot_time_in_state = get_time_in_state(df)
            hist_freq[i, :freq] = (length(df.t)*tot_time_in_state)/total_time
        end
        # push!(all_histos, hist_freq)
        CSV.write("/home/hollie_hindley/Documents/stochastic_hybrid/hists/$(file[1:end-6])/$i.csv", hist_freq)
    end
end
