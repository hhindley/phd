using BlackBoxOptim, PyCall, DifferentialEquations, StaticArrays, BenchmarkTools, DataFrames, OrderedCollections
include("$PATHrtc_models/Oct2022_model/rtc_model.jl")
include("$PATHrtc_models/sol_species_funcs.jl")
include("$PATHParam_inf/inf_setup.jl")
include("$PATHrtc_models/params_init_tspan.jl")

csv_mRNA_conc = DataFrame(CSV.File("$PATHdata/paper_conv.csv"))
mRNA_conc = csv_mRNA_conc[:,1]
mRNA_std = csv_mRNA_conc[:,2]

csv_wt = DataFrame(CSV.File("$PATHdata/results_rtcOFF_grfit.csv"))
csv_wt = select!(csv_wt, Not(["log(OD)", "log(OD) error", "gr error", "od"]))
lam_wt, new_df_wt = extend_gr_curve(csv_wt)

function rtc_bo_ω(x)
    obj_wt_ab = obj(rtc_model1!, initial, ([L, c, kr, Vmax_init, Km_init, x[1], x[2], θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, atp, na, nb, nr, lam_wt]), tspan2, t_2, "mrna", mRNA_conc, mRNA_std)

    return (obj_wt_ab/36)
end

res = bboptimize(rtc_bo_ω; SearchRange = [(0, 0.1), (0, 0.1)], Method = :adaptive_de_rand_1_bin_radiuslimited, MaxSteps = 50000)


print("ω_ab = ", best_candidate(res)[1])
print("ω_r = ", best_candidate(res)[2])
print("error = ", best_fitness(res))

