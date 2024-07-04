using StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed, TableOperations, JSON, Query, FindFirstFunctions, CategoricalArrays, Colors

include(joinpath(homedir(), "phd/stochastic_hybrid_code/analysis_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_funcs.jl"))


function vec_time_in_state(df)
    d = diff(df.actual_index)
    in_seq = d .== 1
    starts_v = findall(diff([false; in_seq]) .== 1)
    stops_v = findall(diff([in_seq; false]) .== -1) .+ 1
    tot_time_in_state = sum(df.t[stops_v] .- df.t[starts_v])
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

# df = Arrow.Table("/Users/s2257179/stoch_files/test_0307_5_values_from_original_range_final_files/results/thresh_183.68421052631578.arrow") |> DataFrame

# grouped_df, hist_df, bin_edges = get_grouped_df(df, :rh)
# hist_freq = make_hist_freq(bin_edges)
# total_time = hist_df.t[end]
# @elapsed for i in 1:length(grouped_df)
#     df = grouped_df[i]
#     tot_time_in_state = vec_time_in_state(df)
#     hist_freq[i, :freq] = (length(df.t)*tot_time_in_state)/total_time
# end

function create_histogram_files(folder_to_convert, folder_to_store_hists)
    mainpath = "/home/hollie_hindley/Documents/stochastic_hybrid"
    folder = folder_to_convert # change this for different folders 
    folderpath = joinpath(joinpath(mainpath,folder), "results")
    folder_to_store_hists = folder_to_store_hists # change this for different folders
    files = readdir(folderpath)
    for file in files
        filepath = joinpath(folderpath, file)
        res = Arrow.Table(filepath) |> DataFrame
        if !isdir(joinpath(mainpath, folder_to_store_hists))
            mkdir(joinpath(mainpath, folder_to_store_hists))
        end
        if !isdir(joinpath(joinpath(mainpath, folder_to_store_hists), "$(file[1:end-6])"))
            mkdir(joinpath(joinpath(mainpath, folder_to_store_hists), "$(file[1:end-6])"))
        end
        for i in [:rm_a, :rm_b, :rm_r, :rtca, :rtcb, :rtcr, :rh, :rt, :rd]
            grouped_df, hist_df, bin_edges = get_grouped_df(res, i)
            hist_freq = make_hist_freq(bin_edges)
            total_time = hist_df.t[end]
            for i in 1:length(grouped_df)
                df = grouped_df[i]
                tot_time_in_state = vec_time_in_state(df)
                hist_freq[i, :freq] = (length(df.t)*tot_time_in_state)/total_time
            end
            CSV.write(joinpath(joinpath(joinpath(mainpath, folder_to_store_hists), "$(file[1:end-6])"), "$i.csv"), hist_freq)
        end
    end
end

create_histogram_files("run_individually_0407_final_files", "hists_run_individually_0407")