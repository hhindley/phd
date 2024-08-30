using StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed, TableOperations, JSON, Query, FindFirstFunctions, CategoricalArrays, Colors

# using PlotlyJS
using InteractiveViz, GLMakie

# include(joinpath(homedir(), "phd/stochastic_hybrid_code/analysis_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_funcs.jl"))


mount_path = "/Users/s2257179/stoch_files/kdam_testing/"

all_items = readdir(mount_path)
folders = [item for item in all_items if isdir(joinpath(mount_path, item)) && !occursin("DS", item)]
folders_dict = Dict(i => folder for (i, folder) in enumerate(folders))

folders_dict = Dict(filter(pair -> pair.first in [7], folders_dict))

dict_times, dict_kdamvals, dict_titles, dict_results, dict_reacts, dict_props, dict_counts, dict_hists = setup_dicts(folders_dict)

for i in eachindex(folders_dict)
    println(i)
    dict_times[i], dict_kdamvals[i], dict_titles[i], dict_results[i], dict_reacts[i], dict_props[i] = LoadDataVars(folders[i]);
    dict_hists[i] = load_hist_files(joinpath(mount_path, folders_dict[i], "hists"))
    dict_counts[i] = prod_tot_count(dict_reacts[i])
end


dict_plot_times, dict_plot_counts, dict_plot_results, dict_plot_hists, dict_stoch_reacts, dict_plot_props = setup_plot_dicts()
# plotting all
for i in eachindex(folders_dict)
    println(i)
    dict_plot_times[i] = plot_times(dict_times[i], "$(folders_dict[i])", folder=folders_dict[i])    
    dict_plot_counts[i] = plot_totstochcount(dict_kdamvals[i][:kdam], dict_counts[i], "$(folders_dict[i])", folder=folders_dict[i])
    dict_stoch_reacts[i] = plot_results("plot_stoch_reacts", dict_reacts[i], length(dict_kdamvals[i][:kdam]), folders_dict[i], xlabel="reaction", ylabel="count", titles=dict_titles[i], size=(1000,650), tosave=true)
end

# plotting individual 
folder = 7; index = 2;
plot_results("plot_stoch_reacts", dict_reacts[folder][index], 1, folders_dict[folder], xlabel="reaction", ylabel="count", titles=[dict_titles[folder][index]], size=(1000,650), tosave=false)


# fraction of stoch reacts for each condition
function calc_frac_stoch_react(folder)
    for df in eachindex(dict_reacts[folder])
        fraction = [i/sum(dict_reacts[folder][df].count) for i in dict_reacts[folder][1].count]
        dict_reacts[folder][df] = insertcols!(dict_reacts[folder][df], 4, :fraction => fraction)
    end
end

calc_frac_stoch_react(6)
calc_frac_stoch_react(7)
calc_frac_stoch_react(8)
calc_frac_stoch_react(9)

colors = Makie.wong_colors()
labels = ["kdam = $i" for i in dict_kdamvals[6][:kdam]]
events = [event for df in dict_reacts[6] for event in df.event]
fractions = [fraction for df in dict_reacts[6] for fraction in df.fraction]
groups = repeat(Int.(collect(range(1, length(dict_reacts[6]), length=length(dict_reacts[6])))), inner=length(dict_reacts[6][1].event))
elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]

f = Figure()
ax = Axis(f[1,1], xticks=(1:13, react_names_str), xticklabelrotation=45)
barplot!(ax, events, fractions, dodge=groups, color=colors[groups])
Legend(f[1,2], elements, labels)
