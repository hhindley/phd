using BlackBoxOptim, Plots, PyCall, DifferentialEquations, StaticArrays, BenchmarkTools, DataFrames, OrderedCollections
include("/home/holliehindley/phd/rtc_models/Oct2022_model/rtc_model.jl")
include("/home/holliehindley/phd/rtc_models/sol_species_funcs.jl")
include("/home/holliehindley/phd/Param_inf/inf_setup.jl")


# tspan = (0, 100)
# t = exp10.(range(-3,2,15))
# pushfirst!(t, 0)

# sol_syn = sol_with_t(rtc_model, init, tspan, params, t)
# rtca_syn = get_curve(sol_syn, :rtca)
# rtcr_syn = get_curve(sol_syn, :rtcr)
# rd_syn = get_curve(sol_syn, :rd)
# rh_syn = get_curve(sol_syn, :rh)

# function rtc_bo1(x)
#     solu = sol_with_t(rtc_model, init, tspan_promoter, (@SVector [L, c, kr, Vmax_init, Km_init, x[1], x[2], θtscr, g_max, θtlr, km_a, km_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr]), t_promoter)
#     rtca = get_curve(solu, :rtca)
#     rtcr = get_curve(solu, :rtcr)
#     obj = []
#     for (i,j) in zip(rtca, rtca_syn)
#         append!(obj, abs2(i-j))
#     end
#     for (i,j) in zip(rtcr, rtcr_syn)
#         append!(obj, abs2(i-j))
#     end
#     return sum(obj)
# end

res = bboptimize(rtc_bo1; SearchRange = [(0, 10.0), (0, 10.0)], NumDimensions = 2, MaxSteps = 3500)

x=[1,1]
param_dict_new = copy(param_dict)
param_dict_new["ω_ab"] = x[1]; param_dict_new["ω_r"] = x[2]

function rtc_bo_ω(x)
    obj_wt = compare_data_and_sol(rtc_model, tspan_promoter, param_dict_new, t_promoter, "mrnas", WT1, WT1_std)
    return (sum(obj_wt))
end

res = bboptimize(rtc_bo_ω; SearchRange = [(0, 10.0), (0, 10.0)], NumDimensions = 2, MaxSteps = 3500)

best_candidate(res)
best_fitness(res)