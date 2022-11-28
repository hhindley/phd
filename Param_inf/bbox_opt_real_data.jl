using BlackBoxOptim, Plots, PyCall, DifferentialEquations, StaticArrays, BenchmarkTools, DataFrames, OrderedCollections
include("/home/holliehindley/phd/rtc_models/Oct2022_model/rtc_model.jl")
include("/home/holliehindley/phd/rtc_models/sol_species_funcs.jl")
include("/home/holliehindley/phd/Param_inf/inf_setup.jl")
include("/home/holliehindley/phd/rtc_models/params_init_tspan.jl")



x=[1,1]
param_dict_new = copy(param_dict)
param_dict_new["ω_ab"] = x[1]; param_dict_new["ω_r"] = x[2]

function rtc_bo_ω(x)
    obj_wt = compare_data_and_sol(rtc_model, init, tspan_promoter, param_dict_new, t_promoter, "mrnas", WT1, WT1_std)
    return (sum(obj_wt))
end

res = bboptimize(rtc_bo_ω; SearchRange = [(0, 10.0), (0, 10.0)], NumDimensions = 2, MaxSteps = 3500)

best_candidate(res)
best_fitness(res)