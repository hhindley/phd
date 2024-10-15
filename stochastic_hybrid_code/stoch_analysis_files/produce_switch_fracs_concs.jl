using Parameters, CSV, DataFrames, DifferentialEquations, LabelledArrays, BenchmarkTools
using Revise, LinearAlgebra, Printf, ModelingToolkit, OrderedCollections, Colors, JLD2, Statistics, DuckDB

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
mount_path, folders, folders_dict = load_file_structure("kdam_testing/low_kdam", server=true)
folders_dict

dict_times, dict_kdamvals, dict_titles, dict_results, dict_reacts, dict_props, dict_counts, dict_hists = load_data(mount_path, folders, folders_dict, reacts=false, props=false, hists=false)

kdam_res = Dict{String, Any}()
kdam_range = [0.0,0.02,0.04,0.06,0.08,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5]
# saving only RtcA results for use on laptop out of server and vpn 
for ind in 1:20
    kdam_res["$(kdam_range[ind])"] = vcat([dict_results[i][ind].rtca for i in eachindex(dict_results)]...)
end
max_length = maximum(length.(values(kdam_res)))
for key in keys(kdam_res)
    col_length = length(kdam_res[key])
    if col_length < max_length
        kdam_res[key] = vcat(kdam_res[key], fill(missing, max_length - col_length))
    end
end
kdam_res
df = DataFrame(kdam_res)

@save "/home/hollie_hindley/Documents/stochastic_hybrid/saved_variables/low_kdam_rtca.jld2" df



size_in_bytes = Base.summarysize(df)

# Convert the size to gigabytes
size_in_gb = size_in_bytes / (1024^3)


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