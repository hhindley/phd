using InteractiveViz, WGLMakie, StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed, TableOperations, JSON, Query, FindFirstFunctions, CategoricalArrays, Colors

include(joinpath(homedir(), "phd/stochastic_hybrid_code/analysis_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_funcs.jl"))

mount_path = "/Users/s2257179/stoch_files/kdam_testing/"

all_items = readdir(mount_path)
folders = [item for item in all_items if isdir(joinpath(mount_path, item)) && !occursin("DS", item)]
folders_dict = Dict(i => folder for (i, folder) in enumerate(folders))

folders_dict = Dict(filter(pair -> pair.first in [11], folders_dict))

dict_times, dict_kdamvals, dict_titles, dict_results, dict_reacts, dict_props, dict_counts, dict_hists = setup_dicts(folders_dict)

dict_titles[11][5]
# plot one result
folder = 12; index = 1; species = "rtca"; num_plots = 1;
f = plot_results("plot_results", dict_results[folder][index], num_plots, folders_dict[folder], titles=[dict_titles[folder][index]], species="$species", xlabel="time", ylabel="$species", size=(1000,650), tosave=false)

# testing if at steady state
# set range for sample size for histogram comparisons
tests = collect(range(1, length(dict_results[11][1].rh), step=75000))
fig = Figure()
ax=Axis(fig[1,1])
for i in 1:length(tests)-1
    println(i)
    h = hist!(ax, dict_results[11][1].rh[tests[end-i]:tests[end-(i-1)]], bins=20)
    text!(ax, "$i", position = Point2f(1e4, 1800), align = (:center, :center), fontsize = 24, color = :red)
    sleep(1)
end # stop loop at point where histogram looks different 
n = 7 # the number where the histogram looks different
tests[end]-tests[end-(n-1)] # will be the length of the ss sample size 
100*(450000/length(dict_results[11][1].rh)) # percentage of the total simulation that the ss sample size is
ss_sample = dict_results[11][1].rh[tests[end-6]:end] # the steady state sample
mean(ss_sample) # mean of the steady state sample
# in this example the steady state region was 15% of the total simulation, so for now will take 10% as the steady state sample in hysteresis analysis
