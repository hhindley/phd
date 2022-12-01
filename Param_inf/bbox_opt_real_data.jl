using BlackBoxOptim, Plots, PyCall, DifferentialEquations, StaticArrays, BenchmarkTools, DataFrames, OrderedCollections
include("/home/holliehindley/phd/rtc_models/Oct2022_model/rtc_model.jl")
include("/home/holliehindley/phd/rtc_models/sol_species_funcs.jl")
include("/home/holliehindley/phd/Param_inf/inf_setup.jl")
include("/home/holliehindley/phd/rtc_models/params_init_tspan.jl")



x=[1,1]
param_dict_bbopt = copy(param_dict) # work out a way to include this in choosing params function so can be used on the below line and dont have to alter the exact dictionary
param_dict["ω_ab"] = x[1]; param_dict["ω_r"] = x[2];
param_dict


function rtc_bo_ω(x)
    obj_wt = compare_data_and_sol(rtc_model, tspan2, t_2, "mrnas", WT1, WT1_std)
    return (sum(obj_wt))
end

res = bboptimize(rtc_bo_ω; SearchRange = [(0, 10.0), (0, 10.0)], NumDimensions = 2, MaxSteps = 3500)

best_candidate(res)
best_fitness(res)