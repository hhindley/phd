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


# mount_path, folders, folders_dict = load_file_structure("kdam_testing/keyvals2_low_kdam", server=true)
# mount_path, folders, folders_dict = load_file_structure("hysteresis/low", server=true)
mount_path, folders, folders_dict = load_file_structure("kdam_testing/high_kdam", server=true)
folders_dict
# folders

# mount_path = "/home/hollie_hindley/Documents/stochastic_hybrid/kdam_testing/high_kdam"
# all_items = readdir("/home/hollie_hindley/Documents/stochastic_hybrid/kdam_testing/high_kdam")
# folders = [item for item in all_items if !occursin("DS", item)]
# data_folders = [readdir(joinpath("/home/hollie_hindley/Documents/stochastic_hybrid/kdam_testing/high_kdam", folders[i])) for i in eachindex(folders)]
# all_folders = [joinpath(folders[folder], readdir(joinpath("/home/hollie_hindley/Documents/stochastic_hybrid/kdam_testing/high_kdam", folders[folder]))[i]) for folder in eachindex(folders) for i in eachindex(data_folders[folder])]
# folders_dict = Dict(i => folder for (i, folder) in enumerate(all_folders))

# folder = folders_dict[1]
# filepath = joinpath(mount_path, folder)
# times_file = (replace(folder, "final_files" => "") * "times.csv")[8:end]
# if isfile(joinpath(filepath, times_file))
#     df_times = CSV.File(joinpath(filepath, times_file)) |> DataFrame
# else
#     df_times = []#CSV.File(timefilepath) |> DataFrame
# end
# joinpath(filepath, times_file)



# LoadDataVars(folders_dict[1])

# joinpath("/home/hollie_hindley/Documents/stochastic_hybrid/kdam_testing/high_kdam", folders_dict[1][1])



dict_times, dict_kdamvals, dict_titles, dict_results, dict_reacts, dict_props, dict_counts, dict_hists = load_data(mount_path, folders, folders_dict, reacts=false, props=false, hists=false)

# @save "/home/hollie_hindley/Documents/stochastic_hybrid/saved_variables/keyvals2_raw_data.jld2" folders_dict dict_times dict_kdamvals dict_titles dict_results dict_reacts dict_props dict_counts dict_hists

thresholds_rtca, thresholds_rtcb = get_unstab_threshold_array(collect(keys(folders_dict))[1]) # argument just has to be any folder number to get kdam vals

# check threshold for using unstable region of bistable curve 
# using PlotlyJS
# plot(scatter(x=dict_kdamvals[1][:kdam], y=thresholds_rtca))

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

indices = Dict("start"=>Dict("2"=>all_start_indices2, "5"=>all_start_indices5, "10"=>all_start_indices10, "bs"=>all_start_indices_bs),
                    "stop"=>Dict("2"=>all_stop_indices2, "5"=>all_stop_indices5, "10"=>all_stop_indices10, "bs"=>all_stop_indices_bs))

switch_rates = Dict("on"=>Dict("2"=>all_switch_rates_on2, "5"=>all_switch_rates_on5, "10"=>all_switch_rates_on10, "bs"=>all_switch_rates_on_bs),
                    "off"=>Dict("2"=>all_switch_rates_off2, "5"=>all_switch_rates_off5, "10"=>all_switch_rates_off10, "bs"=>all_switch_rates_off_bs))

fracs = Dict("on"=>Dict("2"=>all_fracs_on2, "5"=>all_fracs_on5, "10"=>all_fracs_on10, "bs"=>all_fracs_on_bs),
            "off"=>Dict("2"=>all_fracs_off2, "5"=>all_fracs_off5, "10"=>all_fracs_off10, "bs"=>all_fracs_off_bs))

species_mean = Dict("on"=>Dict("2"=>all_species_mean_on2, "5"=>all_species_mean_on5, "10"=>all_species_mean_on10, "bs"=>all_species_mean_on_bs),
                    "off"=>Dict("2"=>all_species_mean_off2, "5"=>all_species_mean_off5, "10"=>all_species_mean_off10, "bs"=>all_species_mean_off_bs))

thresholds_bs = Dict("rtca"=>thresholds_rtca, "rtcb"=>thresholds_rtcb)

@save "/home/hollie_hindley/Documents/stochastic_hybrid/saved_variables/high_kdam.jld2" indices switch_rates fracs species_mean thresholds_bs
# @save "/home/hollie_hindley/Documents/stochastic_hybrid/saved_variables/low_kdam.jld2" indices switch_rates fracs species_mean thresholds_bs