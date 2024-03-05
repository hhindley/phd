include("/home/holliehindley/phd/growth_model/parameters/gm_uM_parameters.jl")
include("/home/holliehindley/phd/rtc_model/parameters/rtc_params.jl")


koff_val = 0.06 # 1/min # inferred by elena # chloramphenicol unbinding rate 
kdam_p_val = 0.3 # represents percent damaged here 

ω_ab_val_comb = ω_ab_val*2
ω_r_val_comb = ω_r_val*10


# params_comb = @LArray [d, kr, L, c, wr_uM, we_uM, we_uM, wq_uM, ω_ab, ω_r, thetar_uM, thetax_uM, Kq_uM, nq, Vmax_init, Km_init, nr, nx, nR, nA, nB, gmax, Kgamma_uM, M_uM, kb_uM, ku, abx, kon, koff, vt, vm, s0_uM, Kt_uM, Km_uM, km_a, km_b, krep, kdam_p, ktag, kdeg, k_diss, ns] (:d, :kr, :L, :c, :w_rh, :w_t, :w_m, :w_q, :w_BA, :w_R, :θ_rh, :θ_nr, :Kq, :hq, :Vmax_init, :Km_init, :nrh, :nx, :nR, :nA, :nB, :gmax, :Kgamma, :M, :kb, :ku, :abx, :kon, :koff, :vt, :vm, :s0, :Kt, :Km, :km_a, :km_b, :krep, :kdam, :ktag, :kdeg, :k_diss, :ns)


# m_rh_0 = 0;
# m_t_0 = 0;
# m_m_0 = 0;
# m_q_0 = 0;
# m_R_0 = 0;
# m_A_0 = 0;
# m_B_0 = 0;
# c_rh_0 = 0;
# c_t_0 = 0;
# c_m_0 = 0;
# c_q_0 = 0;
# c_R_0 = 0;
# c_A_0 = 0;
# c_B_0 = 0;
# et_0 = 0;
# em_0 = 0;
# q_0 = 0;
# R_0 = 0;
# A_0 = 0;
# B_0 = 0;
# rh_0_uM = 11.29*sf #0.0166;
# rt_0 = 0 ;
# rd_0 = 0 ;
# z_rh_0 = 0;
# z_t_0 = 0;
# z_m_0 = 0;
# z_q_0 = 0;
# z_R_0 = 0;
# z_A_0 = 0;
# z_B_0 = 0;
# si_0 = 0;
# a_0_uM = 1000*sf# 1.661

# init_comb = [combined_model.c_rh=>0.0,combined_model.em=>0.0,combined_model.c_q=>0.0,combined_model.c_t=>0.0,combined_model.et=>0.0,combined_model.c_m=>0.0,combined_model.m_t=>0.0,combined_model.m_m=>0.0,combined_model.q=>0.0,combined_model.si=>0.0,combined_model.m_q=>0.0,combined_model.m_rh=>0.0,combined_model.rh=>0.0166,combined_model.a=>1.66,combined_model.z_m=>0.0,combined_model.z_q=>0.0,combined_model.z_t=>0.0,combined_model.z_m=>0.0,combined_model.m_R=>0.0,combined_model.m_A=>0.0,combined_model.m_B=>0.0,combined_model.]
# comb_init_v2 = [m_rh_0, m_t_0, m_m_0, m_q_0, m_R_0, m_A_0, m_B_0, c_rh_0, c_t_0, c_m_0, c_q_0, c_R_0, c_A_0, c_B_0, et_0, em_0, q_0, R_0, A_0, B_0, rh_0_uM, rt_0, rd_0, z_rh_0, z_t_0, z_m_0, z_q_0, z_R_0, z_A_0, z_B_0, 0,0,0,0,0,0,0,si_0, a_0_uM]

# species_comb = [:m_rh, :m_t, :m_m, :m_q, :m_R, :m_A, :m_B, :c_rh, :c_t, :c_m, :c_q, :c_R, :c_A, :c_B, :et, :em, :q, :R, :A, :B, :rh, :rt, :rd, :z_rh, :z_t, :z_m, :z_q, :z_R, :z_A, :z_B, :si, :a]