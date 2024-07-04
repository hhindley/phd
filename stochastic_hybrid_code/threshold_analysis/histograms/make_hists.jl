using StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed, TableOperations, JSON, Query, FindFirstFunctions, CategoricalArrays, Colors

include(joinpath(homedir(), "phd/stochastic_hybrid_code/analysis_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_funcs.jl"))


function get_time_in_state(grouped_df)
    in_sequence = false
    stops = Int[]
    starts = Int[]
    index_length = length(grouped_df.actual_index)
    # println("looping through indexes")
    for i in 1:index_length - 1
        current_index = grouped_df.actual_index[i]
        next_index = grouped_df.actual_index[i + 1]
        # println("checking if increasing by 1")
        if current_index + 1 == next_index
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
        push!(stops, index_length)
    end

    # println("calculating total time in each state!")
    tot_time_in_state = 0
    for i in eachindex(stops)
        tot_time_in_state += grouped_df.t[stops[i]]-grouped_df.t[starts[i]]
    end

    # println("finished!")
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

grouped_df, hist_df, bin_edges = get_grouped_df(df, :rh)
hist_freq = make_hist_freq(bin_edges)
total_time = hist_df.t[end]
for (df,i) in ProgressBar(zip(grouped_df, 1:length(grouped_df)))
    tot_time_in_state = get_time_in_state(df)
    hist_freq[i, :freq] = (length(df.t)*tot_time_in_state)/total_time
end
# CSV.write("/Users/s2257179/phd/stochastic_hybrid_code/hists/thresh10/rh.csv", hist_freq)

using Plots
histogram(df.rh[1:end], xlims=(800,3e3), bins=50)


bins = [i[1] for i in hist_freq.bin]
push!(bins, hist_freq.bin[end][2])
bin_c = (bins[1:end-1] .+ bins[2:end]) ./ 2
bar(bin_c, hist_freq.freq)


df = Arrow.Table("/home/hollie_hindley/Documents/stochastic_hybrid/test_0307_5_values_from_original_range_final_files/results/thresh_183.68421052631578.arrow") |> DataFrame
@elapsed grouped_df, hist_df, bin_edges = get_grouped_df(df, :rh)
@elapsed hist_freq = make_hist_freq(bin_edges)
total_time = hist_df.t[end]
Threads.nthreads()
@elapsed Threads.@threads for i in 1:length(grouped_df)
    df = grouped_df[i]
    tot_time_in_state = get_time_in_state(df)
    hist_freq[i, :freq] = (length(df.t)*tot_time_in_state)/total_time
end
length(grouped_df[3].t)
lengths = [(length(grouped_df[i].t)) for i in 1:length(grouped_df)]
findmax(lengths)

df = grouped_df[1]
function testing(df)
    in_sequence = false
    stops = Int[]
    starts = Int[]
    index_length = length(df.actual_index)
    for i in 1:index_length - 1
        current_index = df.actual_index[i]
        next_index = df.actual_index[i + 1]
        if current_index + 1 == next_index
            if !in_sequence
                in_sequence = true
                push!(starts, i)  
            end
        else
            if in_sequence

                in_sequence = false
                push!(stops, i) 
            end
        end
    end
    return stops, starts
end
@elapsed stops, starts = testing(df)

# Handle the case where the sequence ends with the last element
# println("sequence ends with last element case")

function testing_vectorized(df)
    # Step 1: Calculate differences between consecutive indices
    diffs = diff(df.actual_index)
    
    # Step 2: Identify where the difference is exactly 1 (true for consecutive indices)
    in_sequence = diffs .== 1
    
    # Step 3: Find starts - where in_sequence switches from false to true
    starts = findall(diff([false; in_sequence]) .== 1)
    
    # Step 4: Find stops - where in_sequence switches from true to false
    stops = findall(diff([in_sequence; false]) .== -1)
    
    return stops, starts
end
@elapsed stops_v, starts_v = testing_vectorized(df)

diffs = diff(df.actual_index)
in_sequence = diffs .== 1
starts = findall(diff([false; in_sequence]) .== 1)
stops = findall(diff([in_sequence; false]) .== -1)

hist_freq


ex = DataFrame(x=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4],y=[1,2,3,4,5,6,7,8,9,15,16,17,18,19])
d = diff(ex.y)
in_seq = d .== 1
sta = findall(diff([false; in_seq]) .== 1)
sto = findall(diff([in_seq; false]) .== -1)



stops
stops_v

@elapsed if in_sequence
    push!(stops, index_length)
end

# println("calculating total time in each state!")
tot_time_in_state = 0
@elapsed for i in eachindex(stops)
    tot_time_in_state += grouped_df[3].t[stops[i]]-grouped_df[3].t[starts[i]]
end


# orig_hist_freq
# orig_hist_freq = hist_freq

CSV.write("/home/hollie_hindley/phd/stochastic_hybrid_code/threshold_analysis/histograms/test_rh.csv", hist_freq)






# folderpath = "/home/hollie_hindley/Documents/stochastic_hybrid/thresh_test_arrow_files_14_06/results"
# files = readdir(folderpath)
# for file in files
#     filepath = joinpath(folderpath, file)
#     res = Arrow.Table(filepath) |> DataFrame
#     if !isdir("/home/hollie_hindley/Documents/stochastic_hybrid/hists/$(file[1:end-6])")
#         mkdir("/home/hollie_hindley/Documents/stochastic_hybrid/hists/$(file[1:end-6])")
#     end
#     for i in [:rm_a, :rm_b, :rm_r, :rtca, :rtcb, :rtcr, :rh, :rt, :rd]
#         grouped_df, hist_df, bin_edges = get_grouped_df(res, i)
#         hist_freq = make_hist_freq(bin_edges)
#         total_time = hist_df.t[end]
#         for (df,i) in ProgressBar(zip(grouped_df, 1:length(grouped_df)))
#             tot_time_in_state = get_time_in_state(df)
#             hist_freq[i, :freq] = (length(df.t)*tot_time_in_state)/total_time
#         end
#         # push!(all_histos, hist_freq)
#         CSV.write("/home/hollie_hindley/Documents/stochastic_hybrid/hists/$(file[1:end-6])/$i.csv", hist_freq)
#     end
# end
