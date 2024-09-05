# using StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed, TableOperations, JSON, Query, FindFirstFunctions, CategoricalArrays
# using Parameters, LabelledArrays, BenchmarkTools
# using Revise, LinearAlgebra, Printf, ModelingToolkit, OrderedCollections
using InteractiveViz, GLMakie, Parameters, CSV, DataFrames, DifferentialEquations, LabelledArrays, BenchmarkTools, Revise, LinearAlgebra, Printf, ModelingToolkit, OrderedCollections, Colors, JLD2

include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/switching_funcs.jl"))

include(joinpath(homedir(), "phd/general_funcs/all_model_funcs.jl"))
include(joinpath(homedir(), "phd/general_funcs/solving.jl"))
include(joinpath(homedir(), "phd/rtc_model/parameters/rtc_params.jl"))
include(joinpath(homedir(), "phd/rtc_model/parameters/rtc_params_molecs.jl"))

include(joinpath(homedir(), "phd/rtc_model/models/rtc_orig.jl"))
include(joinpath(homedir(), "phd/rtc_model/functions/bf_funcs/bf_funcs.jl"))

mount_path, folders, folders_dict = load_file_structure("kdam_testing/keyvals2")
folders_dict
folders_dict = Dict(filter(pair -> pair.first in [11,12,13,14,15,16,17,18,19,20], folders_dict))
# folders_dict = Dict(filter(pair -> pair.first in [6], folders_dict))
dict_times, dict_kdamvals, dict_titles, dict_results, dict_reacts, dict_props, dict_counts, dict_hists = load_data(mount_path, folders, folders_dict, reacts=false, props=false, hists=false)

@save "/Users/s2257179/phd/stochastic_hybrid_code/saved_data/keyvals2_raw_data.jld2" folders_dict dict_times dict_kdamvals dict_titles dict_results dict_reacts dict_props dict_counts dict_hists

# plot one result
folder = 6; index = 1; species = "rtca"; num_plots = 1;
f = plot_results("plot_results", dict_results[folder][index], num_plots, folders_dict[folder], titles=[dict_titles[folder][index]], species="$species", xlabel="time", ylabel="$species", size=(800,650), tosave=false, conc=true)
f_rib = plot_results("plot_results", dict_results[folder][index], 1, folders_dict[folder], titles=[dict_titles[folder][index]], species=[:rh, :rd, :rt], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false, conc=false)
f_mrna = plot_results("plot_results", dict_results[folder][index], 1, folders_dict[folder], species=[:rm_a, :rm_b, :rm_r], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false, conc=true)
f_prot = plot_results("plot_results", dict_results[folder][index], 1, folders_dict[folder], titles=[dict_titles[folder][index]], species=[:rtca, :rtcb, :rtcr], xlabel="time", ylabel="specie", size=(800,650), tosave=false, linkaxes=false, conc=false)

# single threshold 
threshold = 5
# thresholds as unstable steady state line 
thresholds_rtca, thresholds_rtcb = get_unstab_threshold_array(collect(keys(folders_dict))[1]) # argument just has to be any folder number to get kdam vals
thresholds_rtca[1]
thresholds_rtcb[1]

# for calculating one kdam val in one folder 
folder = 18; kdam = 1; species = :rh; # [:rtca, :rtcb]; 
condition = set_condition(dict_results[folder][kdam], threshold_rtca=threshold) # single threshold 
condition = set_condition(dict_results[folder][kdam], threshold_rtca=thresholds_rtca[kdam], threshold_rtcb=thresholds_rtcb[kdam]) # different thresholds for different kdam vals 
start_indices, stop_indices = start_stop_indices(dict_results[folder][kdam], condition)
start_times, stop_times = switch_times(dict_results[folder][kdam], start_indices, stop_indices)
switch_rate_on, switch_rate_off = calc_switch_rate(dict_results[folder][kdam], start_times, stop_times)
frac_on, frac_off = calc_frac_times(switch_rate_on, switch_rate_off)
species_mean_on = calc_av_state_conc(dict_results[folder][kdam], species, start_indices, stop_indices)
species_mean_off = calc_av_state_conc(dict_results[folder][kdam], species, start_indices, stop_indices, on=false)

# plot switch times
f = Figure()
ax = Axis(f[1,1], title="Threshold for on state = $threshold, switch rate = $(round(switch_rate_on, digits=6))", xlabel="time", ylabel="molecule number")
lines!(ax, dict_results[folder][kdam].time, dict_results[folder][kdam].rh, label="rtca")
# lines!(ax, dict_results[folder][kdam].time, dict_results[folder][kdam].rtcb, label="rtcb")
scatter!(ax, start_times, ones(length(start_times)), markersize=10, label="start on")
scatter!(ax, stop_times, ones(length(stop_times)), markersize=10, label="stop on")
lines!(ax, [0, dict_results[folder][kdam].time[end]], [thresholds_rtca[1], thresholds_rtca[1]], color=:red, label="threshold")
axislegend()
DataInspector(f)
display(GLMakie.Screen(), f)  


# calculating for all kdam vals in one folder 
# start_indices_f, stop_indices_f = folder_indices(folder, 5) # single threshold
@elapsed start_indices_f, stop_indices_f = folder_indices(folder, thresholds_rtca, threshold_rtcb=thresholds_rtcb) # different thresholds for different kdam vals
@elapsed switch_rates_on_f, switch_rates_off_f, fracs_on_f, fracs_off_f = folder_switchrates_fracs(folder, start_indices_f, stop_indices_f)
@elapsed species_mean_on_f, species_mean_off_f = folder_concs(folder, species, start_indices_f, stop_indices_f)

species_mean_on_f[0.0]
species_mean_off_f[0.0]

# all folders in folders_dict
all_start_indices, all_stop_indices = all_indices(folders_dict, threshold) # single threshold 
all_switch_rates_on, all_switch_rates_off, all_fracs_on, all_fracs_off = all_switchrates_fracs(folders_dict, all_start_indices, all_stop_indices)
all_species_mean_on, all_species_mean_off = all_concs(folders_dict, [:rtca, :rh], all_start_indices, all_stop_indices)

all_start_indices_bs, all_stop_indices_bs = all_indices(folders_dict, thresholds_rtca, threshold_rtcb=thresholds_rtcb) 
all_switch_rates_on_bs, all_switch_rates_off_bs, all_fracs_on_bs, all_fracs_off_bs = all_switchrates_fracs(folders_dict, all_start_indices_bs, all_stop_indices_bs)
all_species_mean_on_bs, all_species_mean_off_bs = all_concs(folders_dict, [:rtca, :rh], all_start_indices_bs, all_stop_indices_bs)

@save "/Users/s2257179/phd/stochastic_hybrid_code/saved_data/data.jld2" folders_dict all_start_indices all_stop_indices all_switch_rates_on all_switch_rates_off all_fracs_on all_fracs_off all_species_mean_on all_species_mean_off all_start_indices_bs all_stop_indices_bs all_switch_rates_on_bs all_switch_rates_off_bs all_fracs_on_bs all_fracs_off_bs all_species_mean_on_bs all_species_mean_off_bs
