# using Parameters, DataFrames, DifferentialEquations, Revise, LinearAlgebra, Printf, ModelingToolkit, JLD2

# using Parameters, DataFrames, DifferentialEquations, Revise, ModelingToolkit, JLD2

# using ModelingToolkit, DifferentialEquations, Colors, OrderedCollections, Printf, LabelledArrays, BifurcationKit, Revise


using Parameters, CSV, DataFrames, DifferentialEquations, LabelledArrays, BenchmarkTools
using Revise, LinearAlgebra, Printf, ModelingToolkit, OrderedCollections, Colors, JLD2, Statistics

include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/file_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/plotting_funcs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/switching_funcs.jl"))

include(joinpath(homedir(), "phd/general_funcs/all_model_funcs.jl"))
include(joinpath(homedir(), "phd/general_funcs/solving.jl"))
include(joinpath(homedir(), "phd/rtc_model/parameters/rtc_params.jl"))
include(joinpath(homedir(), "phd/rtc_model/parameters/rtc_params_molecs.jl"))

include(joinpath(homedir(), "phd/rtc_model/models/rtc_orig.jl"))
include(joinpath(homedir(), "phd/rtc_model/functions/bf_funcs/bf_funcs.jl"))


mount_path, folders, folders_dict = load_file_structure("kdam_testing/keyvals2", server=true)
folders_dict
folders_dict = Dict(filter(pair -> pair.first in [2,4,6,8,10,11,12,13,14,15,16,17,18,19,20], folders_dict))
# folders_dict = Dict(filter(pair -> pair.first in [6], folders_dict))
dict_times, dict_kdamvals, dict_titles, dict_results, dict_reacts, dict_props, dict_counts, dict_hists = load_data(mount_path, folders, folders_dict, reacts=false, props=false, hists=false)

@save "/home/hollie_hindley/Documents/stochastic_hybrid/saved_data/keyvals2_raw_data.jld2" folders_dict dict_times dict_kdamvals dict_titles dict_results dict_reacts dict_props dict_counts dict_hists

thresholds_rtca, thresholds_rtcb = get_unstab_threshold_array(collect(keys(folders_dict))[1]) # argument just has to be any folder number to get kdam vals

# all folders in folders_dict
all_start_indices2, all_stop_indices2 = all_indices(folders_dict, 2) # single threshold 
all_switch_rates_on2, all_switch_rates_off2, all_fracs_on2, all_fracs_off2 = all_switchrates_fracs(folders_dict, all_start_indices2, all_stop_indices2)
all_species_mean_on2, all_species_mean_off2 = all_concs(folders_dict, [:rtca, :rm_a, :rh], all_start_indices2, all_stop_indices2)

all_start_indices5, all_stop_indices5 = all_indices(folders_dict, 5) # single threshold 
all_switch_rates_on5, all_switch_rates_off5, all_fracs_on5, all_fracs_off5 = all_switchrates_fracs(folders_dict, all_start_indices5, all_stop_indices5)
all_species_mean_on5, all_species_mean_off5 = all_concs(folders_dict, [:rtca, :rm_a, :rh], all_start_indices5, all_stop_indices5)

all_start_indices10, all_stop_indices10 = all_indices(folders_dict, 10) # single threshold 
all_switch_rates_on10, all_switch_rates_off10, all_fracs_on10, all_fracs_off10 = all_switchrates_fracs(folders_dict, all_start_indices10, all_stop_indices10)
all_species_mean_on10, all_species_mean_off10 = all_concs(folders_dict, [:rtca, :rm_a, :rh], all_start_indices10, all_stop_indices10)

all_start_indices_bs, all_stop_indices_bs = all_indices(folders_dict, thresholds_rtca, threshold_rtcb=thresholds_rtcb) 
all_switch_rates_on_bs, all_switch_rates_off_bs, all_fracs_on_bs, all_fracs_off_bs = all_switchrates_fracs(folders_dict, all_start_indices_bs, all_stop_indices_bs)
all_species_mean_on_bs, all_species_mean_off_bs = all_concs(folders_dict, [:rtca, :rm_a, :rh], all_start_indices_bs, all_stop_indices_bs)

@save "/home/hollie_hindley/Documents/stochastic_hybrid/saved_data/data_thresh_2.jld2" all_start_indices2 all_stop_indices2 all_switch_rates_on2 all_switch_rates_off2 all_fracs_on2 all_fracs_off2 all_species_mean_on2 all_species_mean_off2
@save "/home/hollie_hindley/Documents/stochastic_hybrid/saved_data/data_thresh_5.jld2" all_start_indices5 all_stop_indices5 all_switch_rates_on5 all_switch_rates_off5 all_fracs_on5 all_fracs_off5 all_species_mean_on5 all_species_mean_off5
@save "/home/hollie_hindley/Documents/stochastic_hybrid/saved_data/data_thresh_10.jld2" all_start_indices10 all_stop_indices10 all_switch_rates_on10 all_switch_rates_off10 all_fracs_on10 all_fracs_off10 all_species_mean_on10 all_species_mean_off10
@save "/home/hollie_hindley/Documents/stochastic_hybrid/saved_data/data_thresh_bs.jld2" all_start_indices_bs all_stop_indices_bs all_switch_rates_on_bs all_switch_rates_off_bs all_fracs_on_bs all_fracs_off_bs all_species_mean_on_bs all_species_mean_off_bs thresholds_rtca thresholds_rtcb

