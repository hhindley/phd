using BlackBoxOptim, Plots, PyCall, DifferentialEquations, StaticArrays, BenchmarkTools, DataFrames, OrderedCollections
include("/home/holliehindley/phd/rtc_models/Oct2022_model/rtc_model.jl")
include("/home/holliehindley/phd/rtc_models/sol_species_funcs.jl")
include("/home/holliehindley/phd/Param_inf/inf_setup.jl")
include("/home/holliehindley/phd/rtc_models/params_init_tspan.jl")


function rtc_bo_ω(x)
    obj_wt_ab, obj_wt_r = obj(rtc_model, initial, (@SVector [L, c, kr, Vmax_init, Km_init, x[1], x[2], θtscr, g_max, θtlr, km_a, km_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr]), tspan2, t_2, "mrna", WT1, WT1_std)

    return (obj_wt_ab/36 + obj_wt_r/36)
end

res = bboptimize(rtc_bo_ω; SearchRange = [(0, 1.0), (40.0, 100.0)], Method = :adaptive_de_rand_1_bin_radiuslimited, MaxSteps = 50000)


print("ω_ab = ", best_candidate(res)[1])
print("ω_r = ", best_candidate(res)[2])
print("error = ", best_fitness(res))

