using StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Statistics, Arrow, FilePathsBase, Distributed, TableOperations, JSON, Query, FindFirstFunctions, CategoricalArrays
using Parameters, LabelledArrays, BenchmarkTools
using Revise, LinearAlgebra, Printf, ModelingToolkit, OrderedCollections


using InteractiveViz, GLMakie, Parameters, CSV, DataFrames, DifferentialEquations, LabelledArrays, BenchmarkTools, Revise, LinearAlgebra, Printf, ModelingToolkit, OrderedCollections, Colors

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
folders_dict = Dict(filter(pair -> pair.first in [6,7,8,9,10,11,12,13,14,15], folders_dict))
# folders_dict = Dict(filter(pair -> pair.first in [6], folders_dict))
dict_times, dict_kdamvals, dict_titles, dict_results, dict_reacts, dict_props, dict_counts, dict_hists = load_data(mount_path, folders, folders_dict, reacts=false, props=false, hists=false)

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

# for calculating one kdam val in one folder 
folder = 7; kdam = 20; species = [:rtca, :rtcb]; 
condition = set_condition(dict_results[folder][kdam], threshold) # single threshold 
condition1 = set_condition(dict_results[folder][kdam], threshold_rtca=thresholds_rtca[kdam], threshold_rtcb=thresholds_rtcb[kdam]) # different thresholds for different kdam vals 
start_indices, stop_indices = start_stop_indices(dict_results[folder][kdam], condition)
start_times, stop_times = switch_times(dict_results[folder][kdam], start_indices, stop_indices)
switch_rate_on, switch_rate_off = calc_switch_rate(dict_results[folder][kdam], start_times, stop_times)
frac_on, frac_off = calc_frac_times(switch_rate_on, switch_rate_off)
species_mean_on = calc_av_state_conc(dict_results[folder][kdam], species, start_indices, stop_indices)
species_mean_off = calc_av_state_conc(dict_results[folder][kdam], species, start_indices, stop_indices, on=false)

# plot switch times
f = Figure()
ax = Axis(f[1,1], title="Threshold for on state = $threshold, switch rate = $(round(switch_rate_on, digits=6))", xlabel="time", ylabel="molecule number")
lines!(ax, dict_results[folder][kdam].time, dict_results[folder][kdam].rtca, label="rtca")
lines!(ax, dict_results[folder][kdam].time, dict_results[folder][kdam].rtcb, label="rtcb")
scatter!(ax, start_times, ones(length(start_times)), markersize=10, label="start on")
scatter!(ax, stop_times, ones(length(stop_times)), markersize=10, label="stop on")
axislegend()
DataInspector(f)
display(GLMakie.Screen(), f)  


# calculating for all kdam vals in one folder 
start_indices_f, stop_indices_f = folder_indices(folder, 5) # single threshold
start_indices_f2, stop_indices_f2 = folder_indices(folder, thresholds_rtca, threshold_rtcb=thresholds_rtcb) # different thresholds for different kdam vals
switch_rates_on_f, switch_rates_off_f, fracs_on_f, fracs_off_f = folder_switchrates_fracs(folder, start_indices_f, stop_indices_f)
species_mean_on_f, species_mean_off_f = folder_concs(folder, species, start_indices_f, stop_indices_f)



# all folders in folders_dict
all_start_indices, all_stop_indices = all_indices(folders_dict, threshold) # single threshold 
all_start_indices, all_stop_indices = all_indices(folders_dict, thresholds_rtca, threshold_rtcb=thresholds_rtcb) 
all_switch_rates_on, all_switch_rates_off, all_fracs_on, all_fracs_off = all_switchrates_fracs(folders_dict, all_start_indices, all_stop_indices)
all_species_mean_on, all_species_mean_off = all_concs(folders_dict, species, all_start_indices, all_stop_indices)

# sorting results for plotting
kdams = sort(collect(keys(all_switch_rates_on[collect(keys(folders_dict))[1]])))

# plotting all results
f_switch = Figure()
ax = Axis(f_switch[1,1], xlabel="kdam", ylabel="switch rate", yscale=log10)
[lines!(ax, kdams, [all_switch_rates_on[i][key] for key in kdams], label="on→off $i") for i in eachindex(all_switch_rates_on)]
[lines!(ax, kdams, [all_switch_rates_off[i][key] for key in kdams], label="on→off $i") for i in eachindex(all_switch_rates_on)]
axislegend(position=:rc)
display(GLMakie.Screen(), f_switch)

f_frac = Figure()
ax = Axis(f_frac[1,1], xlabel="kdam", ylabel="fraction of time in state", title="threshold method")
[lines!(ax, kdams, [all_fracs_on[i][key] for key in kdams], label="on state $i") for i in eachindex(all_fracs_on)]
[lines!(ax, kdams, [all_fracs_off[i][key] for key in kdams], label="off state $i") for i in eachindex(all_fracs_off)]
axislegend(position=:rc)
display(GLMakie.Screen(), f_frac)

specie_plot = 2
cmap = :rainbow_bgyr_35_85_c72_n256
f = Figure()
ax = Axis(f[1,1], xlabel="Damage rate (min-1)", ylabel="[$(species[specie_plot])] (μM)")
[scatter!(ax, kdams, [all_species_mean_on[i][key][specie_plot] for key in kdams], color=[all_fracs_on[i][key] for key in kdams], colorrange=(0,1), colormap=cmap) for i in eachindex(all_species_mean_on)]
[scatter!(ax, kdams, [all_species_mean_off[i][key][specie_plot] for key in kdams], color=[all_fracs_off[i][key] for key in kdams], colorrange=(0,1), colormap=cmap) for i in eachindex(all_species_mean_off)]
Colorbar(f[1, 2], limits = (0,1), colormap = cmap, label="fraction of time in state")