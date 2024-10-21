using InteractiveViz, WGLMakie, StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed, TableOperations, JSON, Query, FindFirstFunctions, CategoricalArrays, Colors
using JLD2, InteractiveViz, GLMakie, Statistics, DataFrames, ColorSchemes, KernelDensity, Arrow, StatsBase

include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/stoch_analysis_files/data_processing_funcs.jl"))
kdams = [0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5]


# loading data (works best on server when in uni)
mount_path, folders, folders_dict = load_file_structure("kdam_testing")
folders_dict = Dict(filter(pair -> pair.first in [7], folders_dict))
dict_times, dict_kdamvals, dict_titles, dict_results, dict_reacts, dict_props, dict_counts, dict_hists = load_data(mount_path, folders, folders_dict)


# plot one result
folder = 8; index = 4; species = "rm_a"; num_plots = 1;
f = plot_results("plot_results", dict_results[folder][index], num_plots, folders_dict[folder], titles=[dict_titles[folder][index]], species="$species", xlabel="time", ylabel="$species", size=(1000,650), tosave=false)

# testing if at steady state
# set range for sample size for histogram comparisons
tests = collect(range(1, length(dict_results[1][1].rh), step=75000))
fig = Figure()
ax=Axis(fig[1,1])
for i in 1:length(tests)-1
    println(i)
    h = hist!(ax, dict_results[1][1].rh[tests[end-i]:tests[end-(i-1)]], bins=20)
    text!(ax, "$i", position = Point2f(1e4, 1800), align = (:center, :center), fontsize = 24, color = :red)
    sleep(1)
end # stop loop at point where histogram looks different 
n = 7 # the number where the histogram looks different
tests[end]-tests[end-(n-1)] # will be the length of the ss sample size 
100*(450000/length(dict_results[11][1].rh)) # percentage of the total simulation that the ss sample size is
ss_sample = dict_results[11][1].rh[tests[end-6]:end] # the steady state sample
mean(ss_sample) # mean of the steady state sample
# in this example the steady state region was 15% of the total simulation, so for now will take 10% as the steady state sample in hysteresis analysis


# loading data out of server

type_kdam = "low_kdam"
@load "/Users/s2257179/Desktop/saved_variables/$type_kdam/$(type_kdam)_stops.jld2" df_lengths df_stops 
df_lengths_low = df_lengths; df_stops_low = df_stops
df_rtca_low = DataFrame(Arrow.Table("/Users/s2257179/Desktop/saved_variables/$type_kdam/$(type_kdam)_rtca.arrow"))
df_times_low = DataFrame(Arrow.Table("/Users/s2257179/Desktop/saved_variables/$type_kdam/$(type_kdam)_times.arrow"))

res_low, times_res_low = remove_missing(df_rtca_low, df_times_low)

sims, times = split_into_simulations(res_low, times_res_low)

# testing if at steady state
# set range for sample size for histogram comparisons
tests = collect(range(1, length(sims[1][0.0]), step=10000))
fig = Figure()
ax=Axis(fig[1,1])
for i in 1:length(tests)-1
    println(i)
    h = hist!(ax, sims[1][0.0][tests[end-i]:tests[end-(i-1)]], bins=20)
    text!(ax, "$i", position = Point2f(1e4, 1800), align = (:center, :center), fontsize = 24, color = :red)
    sleep(1)
end # stop loop at point where histogram looks different 
n = 7 # the number where the histogram looks different
tests[end]-tests[end-(n-1)] # will be the length of the ss sample size 
100*(450000/length(dict_results[11][1].rh)) # percentage of the total simulation that the ss sample size is
ss_sample = dict_results[11][1].rh[tests[end-6]:end] # the steady state sample
mean(ss_sample) # mean of the steady state sample
# in this example the steady state region was 15% of the total simulation, so for now will take 10% as the steady state sample in hysteresis analysis

