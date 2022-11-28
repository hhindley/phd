using Plots, PyCall, DifferentialEquations, StaticArrays, BenchmarkTools, DataFrames, CSV, OrderedCollections
include("/home/holliehindley/phd/rtc_models/Oct2022_model/rtc_model.jl")
include("/home/holliehindley/phd/rtc_models/sol_species_funcs.jl")
include("/home/holliehindley/phd/Param_inf/inf_setup.jl")
include("/home/holliehindley/phd/rtc_models/params_init_tspan.jl")


function rtc_bo_ω(;ω_ab, ω_r)
     obj_wt = compare_data_and_sol(rtc_model, init, tspan2, param_dict, t_2, "mrnas", WT1, WT1_std)

     obj_gr_wt1 = compare_data_and_sol(rtc_model_density, init_den, tspan2, param_dict, t_2, "den", WT2, WT2_std)
     obj_gr_wt2 = compare_data_and_sol(rtc_model_density, init_den, tspan2, param_dict, t_2, "den", WT3, WT3_std)
     obj_gr_wt3 = compare_data_and_sol(rtc_model_density, init_den, tspan4, param_dict, t_4, "den", WT4, WT4_std)

     return -sum(sum(obj) for obj in (obj_wt, obj_gr_wt1, obj_gr_wt2, obj_gr_wt3))
 end


 function rtc_bo_dam(;kdam)
     obj_wt = compare_data_and_sol(rtc_model, tspan, t_2, "mrnas", WT1, param_dict)

     return -(sum(obj_wt))
 end






# objective function for inferring ω_ab/ω_r
function rtc_bo(;ω_ab, ω_r)
   # ω_ab - promoter expression
    # obj_wt = []
    # solu_wt = sol_with_t(rtc_model, init, tspan, (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km, k_b, gr_c, d, krep, 0, ktag, kdeg, kin, atp, na, nb, nr]), t_2) # wt so damage at 0
    # rm_a_wt = get_curve(solu_wt, :rm_a)
    # rm_b_wt = get_curve(solu_wt, :rm_b)
    # append!(obj_wt, [abs2(i-j) for (i,j) in zip(rm_a_wt, WT1)])
    # append!(obj_wt, [abs2(i-j) for (i,j) in zip(rm_b_wt, WT1)])

    obj_wt = compare_data_and_sol(rtc_model, tspan, t_2, ω_ab, 0, "mrnas", WT1)

    # obj_hpx = []
    # solu_hpx = sol_with_t(rtc_model, init, tspan, (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km, k_b, gr_c, d, krep, 2, ktag, kdeg, kin, atp, na, nb, nr]), t_2) # hpx so damage higher 
    # rm_a_hpx = get_curve(solu_hpx, :rm_a)
    # rm_b_hpx = get_curve(solu_hpx, :rm_b)
    # append!(obj_hpx, [abs2(i-j) for (i,j) in zip(rm_a_hpx, hpx)])
    # append!(obj_hpx, [abs2(i-j) for (i,j) in zip(rm_b_hpx, hpx)])
    
    obj_hpx = compare_data_and_sol(rtc_model, tspan, t_2, ω_ab, 2, "mrnas", hpx)


    # growth data - compare to rh
    #wt 
    # obj_hpxWT = [] 
    # solu_hpx1WT = sol_with_t(rtc_model, init, tspan, (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km, k_b, gr_c, d, krep, 0, ktag, kdeg, kin, atp, na, nb, nr]), t_2)
    # rh_hpx1WT = get_curve(solu_hpx1WT, :rh)
    # append!(obj_hpx1WT, [abs2(i-j) for (i,j) in zip(rh_hpx1WT, WT2)])
    obj_hpxWT = compare_data_and_sol(rtc_model, tspan, t_2, ω_ab, 0, "rh", WT2)

    # solu_hpx2WT = sol_with_t(rtc_model, init, tspan, (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km, k_b, gr_c, d, krep, 0, ktag, kdeg, kin, atp, na, nb, nr]), t_2)
    # rh_hpx2WT = get_curve(solu_hpx2WT, :rh)
    # append!(obj_hpxWT, [abs2(i-j) for (i,j) in zip(rh_hpx2WT, WT3)])
    obj_hpx2WT = compare_data_and_sol(rtc_model, tspan, t_2, ω_ab, 0, "rh", WT3)

    # #hpx
    # obj_hpx1 = []
    # solu_hpx1 = sol_with_t(rtc_model, init, tspan, (@SVector [L, c, kr, Vmax_init, Km_init, 0, ω_r, θtscr, g_max, θtlr, km, k_b, gr_c, d, krep, 2, ktag, kdeg, kin, atp, na, nb, nr]), t_2)
    # rh_hpx1 = get_curve(solu_hpx1, :rh)
    # append!(obj_hpx1, [abs2(i-j) for (i,j) in zip(rh_hpx1, hpx_rtcoff)]) # rtc off - gene expression off - set w_ab to zero - damage still present
    obj_hpx1 = compare_data_and_sol(rtc_model, tspan, t_2, 0, 2, "rh", hpx_rtcoff)

    # obj_hpx2 = []
    # solu_hpx2 = sol_with_t(rtc_model, init, tspan, (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km, k_b, gr_c, d, krep, 2, ktag, kdeg, kin, atp, na, nb, nr]), t_2)
    # rh_hpx2 = get_curve(solu_hpx2, :rh)
    # append!(obj_hpx2, [abs2(i-j) for (i,j) in zip(rh_hpx2, hpx_rtcon)]) # rtc on and damage present 
    obj_hpx2 = compare_data_and_sol(rtc_model, tspan, t_2, ω_ab, 2, "rh", hpx_rtcon)

    # # colD
    # # wt colD
    # obj_colDWT = []
    # solu_colDWT = sol_with_t(rtc_model, init, tspan, (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km, k_b, gr_c, d, krep, 0, ktag, kdeg, kin, atp, na, nb, nr]), t_4) # wt no colD so damage zero
    # rh_colDWT = get_curve(solu_colDWT, :rh)
    # append!(obj_colDWT, [abs2(i-j) for (i,j) in zip(rh_colDWT, WT4)])
    obj_colDWT = compare_data_and_sol(rtc_model, tspan2, t_4, ω_ab, 0, "rh", WT4)

    # obj_colD = []
    # solu_colD = sol_with_t(rtc_model, init, tspan, (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km, k_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr]), t_4) # wt with colD so need to induce the system 
    # rh_colD = get_curve(solu_colD, :rh)
    # append!(obj_colD, [abs2(i-j) for (i,j) in zip(rh_colD, WT_colD)])
    obj_colD = compare_data_and_sol(rtc_model, tspan2, t_4, ω_ab, 0, "rh", WT_colD) # how to induce the system with colD

    # # rtca 
    # obj_a = []
    # solu_a = sol_with_t(rtc_model, init, tspan, (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km, k_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr]), t_4) # wt no colD but no rtcA so 
    # rh_a = get_curve(solu_a, :rh)
    # append!(obj_a, [abs2(i-j) for (i,j) in zip(rh_a, nA)])
    obj_a = compare_data_and_sol(rtc_model, tspan2, t_4, ω_ab, 0, "rh", nA) # use different model 

    # obj_colD_a = []
    # solu_colD_a = sol_with_t(rtc_model, init, tspan, (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km, k_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr]), t_4) # wt no colD but no rtcA so 
    # rh_colD_a = get_curve(solu_colD_a, :rh)
    # append!(obj_colD_a, [abs2(i-j) for (i,j) in zip(rh_colD_a, nA_colD)])
    obj_colD_a = compare_data_and_sol(rtc_model, tspan2, t_4, ω_ab, 0, "rh", nA_colD) # how to induce the system with colD, use different model

    # # rtcb
    # obj_b = []
    # solu_b = sol_with_t(rtc_model, init, tspan, (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km, k_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr]), t_4) # wt no colD but no rtcA so 
    # rh_b = get_curve(solu_b, :rh)
    # append!(obj_b, [abs2(i-j) for (i,j) in zip(rh_b, nB)])
    obj_b = compare_data_and_sol(rtc_model, tspan2, t_4, ω_ab, 0, "rh", nB) # use different model

    # obj_colD_b = []
    # solu_colD_b = sol_with_t(rtc_model, init, tspan, (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km, k_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr]), t_4) # wt no colD but no rtcA so 
    # rh_colD_b = get_curve(solu_colD_b, :rh)
    # append!(obj_colD_b, [abs2(i-j) for (i,j) in zip(rh_colD_b, nB_colD)])
    obj_colD_b = compare_data_and_sol(rtc_model, tspan2, t_4, ω_ab, 0, "rh", nB_colD) # how to induce the system with colD, use different model

    # obj_bcomp = []
    # solu_bcomp = sol_with_t(rtc_model, init, tspan, (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km, k_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr]), t_4) # wt with colD so need to induce the system 
    # rh_bcomp = get_curve(solu_bcomp, :rh)
    # append!(obj_bcomp, [abs2(i-j) for (i,j) in zip(rh_bcomp, nB_B)])
    obj_bcomp = compare_data_and_sol(rtc_model, tspan2, t_4, ω_ab, 0, "rh", nB_B) # use different model

    # obj_colD_bcomp = []
    # solu_colD_bcomp = sol_with_t(rtc_model, init, tspan, (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km, k_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr]), t_4) # wt with colD so need to induce the system 
    # rh_colD_bcomp = get_curve(solu_colD_bcomp, :rh)
    # append!(obj_colD_bcomp, [abs2(i-j) for (i,j) in zip(rh_colD_bcomp, nB_B_colD)])
    obj_colD_bcomp = compare_data_and_sol(rtc_model, tspan2, t_4, ω_ab, 0, "rh", nB_B_colD) # how to induce the system with colD, use different model

    # obj_colD_bmut = []
    # solu_colD_bmut = sol_with_t(rtc_model, init, tspan, (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km, k_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr]), t_4) # wt with colD so need to induce the system 
    # rh_colD_bmut = get_curve(solu_colD_bmut, :rh)
    # append!(obj_colD_bmut, [abs2(i-j) for (i,j) in zip(rh_colD_bmut, nB_Bmut)])
    obj_colD_bmut = compare_data_and_sol(rtc_model, tspan2, t_4, ω_ab, 0, "rh", nB_Bmut) # how to induce the system with colD, use different model

    # obj_bmut = []
    # solu_bmut = sol_with_t(rtc_model, init, tspan, (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km, k_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr]), t_4) # wt with colD so need to induce the system 
    # rh_bmut = get_curve(solu_bmut, :rh)
    # append!(obj_bmut, [abs2(i-j) for (i,j) in zip(rh_bmut, nB_Bmut_colD)])
    obj_bmut = compare_data_and_sol(rtc_model, tspan2, t_4, ω_ab, 0, "rh", nB_Bmut_colD) # how to induce the system with colD, use different model

    # # rtcr 
    # obj_r = []
    # solu_r = sol_with_t(rtc_model, init, tspan, (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km, k_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr]), t_4) # wt no colD but no rtcA so 
    # rh_r = get_curve(solu_r, :rh)
    # append!(obj_r, [abs2(i-j) for (i,j) in zip(rh_r, nR)])
    obj_r = compare_data_and_sol(rtc_model, tspan2, t_4, ω_ab, 0, "rh", nR) # how to induce the system with colD, use different model

    # obj_colD_r = []
    # solu_colD_r = sol_with_t(rtc_model, init, tspan, (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km, k_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr]), t_4) # wt no colD but no rtcA so 
    # rh_colD_r = get_curve(solu_colD_r, :rh)
    # append!(obj_colD_r, [abs2(i-j) for (i,j) in zip(rh_colD_r, nR_colD)])
    obj_colD_r = compare_data_and_sol(rtc_model, tspan2, t_4, ω_ab, 0, "rh", nR_colD) # how to induce the system with colD, use different model


    return -(sum(obj_wt)+sum(obj_hpx)+sum(obj_hpxWT)+sum(obj_hpx2WT)+sum(obj_hpx1)+sum(obj_hpx2)+sum(obj_colDWT)+sum(obj_colD)+sum(obj_a)+sum(obj_colD_a)+sum(obj_b)+sum(obj_colD_b)+sum(obj_bcomp)+sum(obj_colD_bcomp)+sum(obj_colD_bmut)+sum(obj_bmut)+sum(obj_r)+sum(obj_colD_r))
end

# objective function for inferring kdam
function rtc_bo(;kdam)
    # ω_ab - promoter expression
    # obj_wt = []
    # solu_wt = sol_with_t(rtc_model, init, tspan, (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km, k_b, gr_c, d, krep, 0, ktag, kdeg, kin, atp, na, nb, nr]), t_2) # wt so damage at 0
    # rm_a_wt = get_curve(solu_wt, :rm_a)
    # rm_b_wt = get_curve(solu_wt, :rm_b)
    # append!(obj_wt, [abs2(i-j) for (i,j) in zip(rm_a_wt, WT1)])
    # append!(obj_wt, [abs2(i-j) for (i,j) in zip(rm_b_wt, WT1)])

    obj_wt = compare_data_and_sol(rtc_model, tspan, t_2, ω_ab, 0, "mrnas", WT1)

    # obj_hpx = []
    # solu_hpx = sol_with_t(rtc_model, init, tspan, (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km, k_b, gr_c, d, krep, 2, ktag, kdeg, kin, atp, na, nb, nr]), t_2) # hpx so damage higher 
    # rm_a_hpx = get_curve(solu_hpx, :rm_a)
    # rm_b_hpx = get_curve(solu_hpx, :rm_b)
    # append!(obj_hpx, [abs2(i-j) for (i,j) in zip(rm_a_hpx, hpx)])
    # append!(obj_hpx, [abs2(i-j) for (i,j) in zip(rm_b_hpx, hpx)])
    
    obj_hpx = compare_data_and_sol(rtc_model, tspan, t_2, ω_ab, 2, "mrnas", hpx)


    # growth data - compare to rh
    #wt 
    # obj_hpxWT = [] 
    # solu_hpx1WT = sol_with_t(rtc_model, init, tspan, (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km, k_b, gr_c, d, krep, 0, ktag, kdeg, kin, atp, na, nb, nr]), t_2)
    # rh_hpx1WT = get_curve(solu_hpx1WT, :rh)
    # append!(obj_hpx1WT, [abs2(i-j) for (i,j) in zip(rh_hpx1WT, WT2)])
    obj_hpxWT = compare_data_and_sol(rtc_model, tspan, t_2, ω_ab, 0, "rh", WT2)

    # solu_hpx2WT = sol_with_t(rtc_model, init, tspan, (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km, k_b, gr_c, d, krep, 0, ktag, kdeg, kin, atp, na, nb, nr]), t_2)
    # rh_hpx2WT = get_curve(solu_hpx2WT, :rh)
    # append!(obj_hpxWT, [abs2(i-j) for (i,j) in zip(rh_hpx2WT, WT3)])
    obj_hpx2WT = compare_data_and_sol(rtc_model, tspan, t_2, ω_ab, 0, "rh", WT3)

    # #hpx
    # obj_hpx1 = []
    # solu_hpx1 = sol_with_t(rtc_model, init, tspan, (@SVector [L, c, kr, Vmax_init, Km_init, 0, ω_r, θtscr, g_max, θtlr, km, k_b, gr_c, d, krep, 2, ktag, kdeg, kin, atp, na, nb, nr]), t_2)
    # rh_hpx1 = get_curve(solu_hpx1, :rh)
    # append!(obj_hpx1, [abs2(i-j) for (i,j) in zip(rh_hpx1, hpx_rtcoff)]) # rtc off - gene expression off - set w_ab to zero - damage still present
    obj_hpx1 = compare_data_and_sol(rtc_model, tspan, t_2, 0, 2, "rh", hpx_rtcoff)

    # obj_hpx2 = []
    # solu_hpx2 = sol_with_t(rtc_model, init, tspan, (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km, k_b, gr_c, d, krep, 2, ktag, kdeg, kin, atp, na, nb, nr]), t_2)
    # rh_hpx2 = get_curve(solu_hpx2, :rh)
    # append!(obj_hpx2, [abs2(i-j) for (i,j) in zip(rh_hpx2, hpx_rtcon)]) # rtc on and damage present 
    obj_hpx2 = compare_data_and_sol(rtc_model, tspan, t_2, ω_ab, 2, "rh", hpx_rtcon)

    # # colD
    # # wt colD
    # obj_colDWT = []
    # solu_colDWT = sol_with_t(rtc_model, init, tspan, (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km, k_b, gr_c, d, krep, 0, ktag, kdeg, kin, atp, na, nb, nr]), t_4) # wt no colD so damage zero
    # rh_colDWT = get_curve(solu_colDWT, :rh)
    # append!(obj_colDWT, [abs2(i-j) for (i,j) in zip(rh_colDWT, WT4)])
    obj_colDWT = compare_data_and_sol(rtc_model, tspan2, t_4, ω_ab, 0, "rh", WT4)

    # obj_colD = []
    # solu_colD = sol_with_t(rtc_model, init, tspan, (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km, k_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr]), t_4) # wt with colD so need to induce the system 
    # rh_colD = get_curve(solu_colD, :rh)
    # append!(obj_colD, [abs2(i-j) for (i,j) in zip(rh_colD, WT_colD)])
    obj_colD = compare_data_and_sol(rtc_model, tspan2, t_4, ω_ab, 0, "rh", WT_colD) # how to induce the system with colD

    # # rtca 
    # obj_a = []
    # solu_a = sol_with_t(rtc_model, init, tspan, (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km, k_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr]), t_4) # wt no colD but no rtcA so 
    # rh_a = get_curve(solu_a, :rh)
    # append!(obj_a, [abs2(i-j) for (i,j) in zip(rh_a, nA)])
    obj_a = compare_data_and_sol(rtc_model, tspan2, t_4, ω_ab, 0, "rh", nA) # use different model 

    # obj_colD_a = []
    # solu_colD_a = sol_with_t(rtc_model, init, tspan, (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km, k_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr]), t_4) # wt no colD but no rtcA so 
    # rh_colD_a = get_curve(solu_colD_a, :rh)
    # append!(obj_colD_a, [abs2(i-j) for (i,j) in zip(rh_colD_a, nA_colD)])
    obj_colD_a = compare_data_and_sol(rtc_model, tspan2, t_4, ω_ab, 0, "rh", nA_colD) # how to induce the system with colD, use different model

    # # rtcb
    # obj_b = []
    # solu_b = sol_with_t(rtc_model, init, tspan, (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km, k_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr]), t_4) # wt no colD but no rtcA so 
    # rh_b = get_curve(solu_b, :rh)
    # append!(obj_b, [abs2(i-j) for (i,j) in zip(rh_b, nB)])
    obj_b = compare_data_and_sol(rtc_model, tspan2, t_4, ω_ab, 0, "rh", nB) # use different model

    # obj_colD_b = []
    # solu_colD_b = sol_with_t(rtc_model, init, tspan, (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km, k_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr]), t_4) # wt no colD but no rtcA so 
    # rh_colD_b = get_curve(solu_colD_b, :rh)
    # append!(obj_colD_b, [abs2(i-j) for (i,j) in zip(rh_colD_b, nB_colD)])
    obj_colD_b = compare_data_and_sol(rtc_model, tspan2, t_4, ω_ab, 0, "rh", nB_colD) # how to induce the system with colD, use different model

    # obj_bcomp = []
    # solu_bcomp = sol_with_t(rtc_model, init, tspan, (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km, k_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr]), t_4) # wt with colD so need to induce the system 
    # rh_bcomp = get_curve(solu_bcomp, :rh)
    # append!(obj_bcomp, [abs2(i-j) for (i,j) in zip(rh_bcomp, nB_B)])
    obj_bcomp = compare_data_and_sol(rtc_model, tspan2, t_4, ω_ab, 0, "rh", nB_B) # use different model

    # obj_colD_bcomp = []
    # solu_colD_bcomp = sol_with_t(rtc_model, init, tspan, (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km, k_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr]), t_4) # wt with colD so need to induce the system 
    # rh_colD_bcomp = get_curve(solu_colD_bcomp, :rh)
    # append!(obj_colD_bcomp, [abs2(i-j) for (i,j) in zip(rh_colD_bcomp, nB_B_colD)])
    obj_colD_bcomp = compare_data_and_sol(rtc_model, tspan2, t_4, ω_ab, 0, "rh", nB_B_colD) # how to induce the system with colD, use different model

    # obj_colD_bmut = []
    # solu_colD_bmut = sol_with_t(rtc_model, init, tspan, (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km, k_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr]), t_4) # wt with colD so need to induce the system 
    # rh_colD_bmut = get_curve(solu_colD_bmut, :rh)
    # append!(obj_colD_bmut, [abs2(i-j) for (i,j) in zip(rh_colD_bmut, nB_Bmut)])
    obj_colD_bmut = compare_data_and_sol(rtc_model, tspan2, t_4, ω_ab, 0, "rh", nB_Bmut) # how to induce the system with colD, use different model

    # obj_bmut = []
    # solu_bmut = sol_with_t(rtc_model, init, tspan, (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km, k_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr]), t_4) # wt with colD so need to induce the system 
    # rh_bmut = get_curve(solu_bmut, :rh)
    # append!(obj_bmut, [abs2(i-j) for (i,j) in zip(rh_bmut, nB_Bmut_colD)])
    obj_bmut = compare_data_and_sol(rtc_model, tspan2, t_4, ω_ab, 0, "rh", nB_Bmut_colD) # how to induce the system with colD, use different model

    # # rtcr 
    # obj_r = []
    # solu_r = sol_with_t(rtc_model, init, tspan, (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km, k_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr]), t_4) # wt no colD but no rtcA so 
    # rh_r = get_curve(solu_r, :rh)
    # append!(obj_r, [abs2(i-j) for (i,j) in zip(rh_r, nR)])
    obj_r = compare_data_and_sol(rtc_model, tspan2, t_4, ω_ab, 0, "rh", nR) # how to induce the system with colD, use different model

    # obj_colD_r = []
    # solu_colD_r = sol_with_t(rtc_model, init, tspan, (@SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km, k_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr]), t_4) # wt no colD but no rtcA so 
    # rh_colD_r = get_curve(solu_colD_r, :rh)
    # append!(obj_colD_r, [abs2(i-j) for (i,j) in zip(rh_colD_r, nR_colD)])
    obj_colD_r = compare_data_and_sol(rtc_model, tspan2, t_4, ω_ab, 0, "rh", nR_colD) # how to induce the system with colD, use different model


    return -(sum(obj_wt)+sum(obj_hpx)+sum(obj_hpxWT)+sum(obj_hpx2WT)+sum(obj_hpx1)+sum(obj_hpx2)+sum(obj_colDWT)+sum(obj_colD)+sum(obj_a)+sum(obj_colD_a)+sum(obj_b)+sum(obj_colD_b)+sum(obj_bcomp)+sum(obj_colD_bcomp)+sum(obj_colD_bmut)+sum(obj_bmut)+sum(obj_r)+sum(obj_colD_r))
end

# in python writing the ranges to search for parameter value
py"""
param_range_ω = {'ω_ab': (0, 10), 'ω_r': (0, 10)}#, 'kdam': (0, 1)}
"""

py"""
param_range_dam = {'kdam': (0, 1)}
"""
# import bayes_opt package from python
bayes_opt = pyimport("bayes_opt")

# setting the optimizer 
optimizer = bayes_opt.BayesianOptimization(f=rtc_bo_ω, pbounds=py"param_range_ω", random_state=27, verbose=2) # verbose = 1 prints only when a maximum is observed (pink)

# timing the process and maximising the optimizer 
function timer()
    optimizer.maximize(init_points=2, n_iter=100, acq="ucb", kappa=2)
end

@time timer()
print(optimizer.max)




# creating lists of values of parameters tried and errors for each one  
function results(optimizer)
    vals, errors = [], []
    for i in collect(1:length(optimizer.space.target))
        a = collect(values(optimizer.res)[i])
        append!(vals, collect(values(a[2][2])))
        append!(errors, -(collect(values(a[1][2]))))
    end

    # reversing sum of squares error calculation to get errors close to zero 
    errors_ori = sqrt.(errors/15)

    # best values 
    best_param = collect(values(collect(values(optimizer.max))[2]))
    best_error = -[collect(values(optimizer.max))[1]]

    # best error with sum of squares reversed 
    best_error_ori = sqrt.(best_error/15)
    return vals, errors, errors_ori, best_param, best_error, best_error_ori
end

vals, errors, errors_ori, best_param, best_error, best_error_ori = results(optimizer)

# plot errors over iterations 
Plots.plot(range(1,length(optimizer.space.target)), errors, markershapes=[:circle], ylabel=("Error"), xlabel=("Iteration"), legend=false, xticks = 0:5:length(optimizer.space.target), size=(1000,500))#, yaxis=(:log10, (1,Inf)))
# png("error")
Plots.plot(range(1,length(optimizer.space.target)), errors_ori, markershapes=[:circle], ylabel=("Adjusted error"), xlabel=("Iteration"), legend=false, xticks = 0:5:length(optimizer.space.target), size=(1000,500))#, yaxis=(:log10, (1,Inf)))
# png("adjusted errors")

# plot params over iterations
Plots.plot(range(1,length(optimizer.space.target)), vals, markershapes=[:circle], ylabel=("Param attempt"), xlabel=("Iteration"), legend=false, xticks = 0:1:length(optimizer.space.target))#, size=(1000,500))#, yaxis=(:log10, (1,Inf)))
hline!([4.14], labels= false)
# png("param_attempts")

# errors and params
Plots.scatter(vals, errors, ylabel=("Error"), xlabel=("Parameter"), legend = false)#,  yaxis=(:log10, (1,Inf)))
scatter!(best_param, best_error, color = "red", label = "", markershape=[:circle])
# png("error_vs_param")

# errors and params with original errors by sqrt and /15 (undoing sum?)
Plots.scatter(vals, errors_ori, ylabel=("Adjusted error"), xlabel=("Parameter"), legend = false)#,  yaxis=(:log10, (1,Inf)))
scatter!(best_param, best_error_ori, color = "red", label = "")
# png("adjusted_error_vs_param")




# trying different values of kappa
vals1, errors1, errors_ori1, best_param1, best_error1, best_error_ori1 = [], [], [], [], [], []
function plotting_κ(x, y)

    for i in collect(0:2:10)
        optimizer = bayes_opt.BayesianOptimization(f=rtc_bo, pbounds=py"param_range", random_state=27, verbose=2) 
        optimizer.maximize(init_points=2, n_iter=10, acq="ucb", kappa=i)
        vals, errors, errors_ori, best_param, best_error, best_error_ori = results(optimizer)
        push!(vals1, vals)
        push!(errors1, errors)
        push!(errors_ori1, errors_ori)
        push!(best_param1, best_param)
        push!(best_error1, best_error)
        push!(best_error_ori1, best_error_ori)
    end

    p = Plots.plot(layout=(3,2), size=(1000,800), xlabel="Iterations", ylabel="Error", xticks = 0:1:length(optimizer.space.target))
    for (i, j) in zip(1:6, collect(0:2:10))
            Plots.plot!(p, x, y[i], subplot=i, markershapes=[:circle], labels="κ = $j")
     end
    return display(p)
end
plotting_κ(1:12, errors1)
# png("kappa_error")

vals1, errors1, errors_ori1, best_param1, best_error1, best_error_ori1 = [], [], [], [], [], []
function plotting_κ2(x, y)

    for i in collect(0:2:10)
        optimizer = bayes_opt.BayesianOptimization(f=rtc_bo, pbounds=py"param_range", random_state=27, verbose=2) 
        optimizer.maximize(init_points=2, n_iter=10, acq="ucb", kappa=i)
        vals, errors, errors_ori, best_param, best_error, best_error_ori = results(optimizer)
        push!(vals1, vals)
        push!(errors1, errors)
        push!(errors_ori1, errors_ori)
        push!(best_param1, best_param)
        push!(best_error1, best_error)
        push!(best_error_ori1, best_error_ori)
    end

    p = Plots.plot(layout=(3,2), size=(1000,800), xlabel="Iterations", ylabel="Parameter estimate", xticks = 0:1:length(optimizer.space.target))
    for (i, j) in zip(1:6, collect(0:2:10))
            Plots.plot!(p, x, y[i], subplot=i, markershapes=[:circle], labels="κ = $j")
            hline!([4.14 4.14 4.14 4.14 4.14 4.14], labels= false)
     end
    return display(p)
end
plotting_κ2(1:12, vals1)
# png("kappa_param")

vals1, errors1, errors_ori1, best_param1, best_error1, best_error_ori1 = [], [], [], [], [], []
function plotting_κ1(x, y)

    for i in collect(0:2:10)
        optimizer = bayes_opt.BayesianOptimization(f=rtc_bo, pbounds=py"param_range", random_state=27, verbose=2) 
        optimizer.maximize(init_points=2, n_iter=10, acq="ucb", kappa=i)
        vals, errors, errors_ori, best_param, best_error, best_error_ori = results(optimizer)
        push!(vals1, vals)
        push!(errors1, errors)
        push!(errors_ori1, errors_ori)
        push!(best_param1, best_param)
        push!(best_error1, best_error)
        push!(best_error_ori1, best_error_ori)
    end

    p = Plots.plot(layout=(3,2), size=(1000,800), xlabel="Parameter estimate", ylabel="Error")
    for (i, j) in zip(1:6, collect(0:2:10))
            Plots.scatter!(p, x[i], y[i], subplot=i, markershapes=[:circle], labels="κ = $j", color=palette(:default)[1])
            scatter!(best_param1, best_error1, color = "red", label = "")
    end
    return display(p)
end
plotting_κ1(vals1, errors1)
# png("kappa_error_param")
