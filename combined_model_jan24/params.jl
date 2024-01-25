include("/home/holliehindley/phd/growth_model_2021/growth_model/uM_parameters.jl")
include("/home/holliehindley/phd/may23_rtc/rtc_parameters/params.jl")


abx = 0 # uM # abx concentration 
kon = k_cm # 1/min*uM # chloramphenicol binding rate 
koff = 0.06 # 1/min # inferred by elena # chloramphenicol unbinding rate 


comb_params = @LArray [dm, kr, L, c, wr_uM, we_uM, we_uM, wq_uM, ω_ab, ω_r, thetar_uM, thetax_uM, Kq_uM, nq, Vmax_init, Km_init, nr, nx, nR, nA, nB, gmax, Kgamma_uM, M_uM, kb_uM, ku, abx, kon, koff, vt, vm, s0_uM, Kt_uM, Km_uM, km_a, km_b, krep, kdam, ktag, kdeg, ns] (:d, :kr, :L, :c, :w_rh, :w_t, :w_m, :w_q, :w_BA, :w_R, :θ_rh, :θ_nr, :Kq, :hq, :Vmax_init, :Km_init, :nrh, :nt, :nm, :nq, :nR, :nA, :nB, :gmax, :Kgamma, :M, :kb, :ku, :abx, :kon, :koff, :vt, :vm, :s0, :Kt, :Km, :km_a, :km_b, :krep, :kdam, :ktag, :kdeg, :ns)

m_rh_0 = 0;
m_t_0 = 0;
m_m_0 = 0;
m_q_0 = 0;
m_R_0 = 0;
m_A_0 = 0;
m_B_0 = 0;
c_rh_0 = 0;
c_t_0 = 0;
c_m_0 = 0;
c_q_0 = 0;
c_R_0 = 0;
c_A_0 = 0;
c_B_0 = 0;
et_0 = 0;
em_0 = 0;
q_0 = 0;
R_0 = 0;
A_0 = 0;
B_0 = 0;
rh_0_uM = 11.29*sf #0.0166;
rt_0 = 0 ;
rd_0 = 0 ;
z_rh_0 = 0;
z_t_0 = 0;
z_m_0 = 0;
z_q_0 = 0;
z_R_0 = 0;
z_A_0 = 0;
z_B_0 = 0;
si_0 = 0;
a_0_uM = 1000*sf; # 1.661

comb_init = [m_rh_0, m_t_0, m_m_0, m_q_0, m_R_0, m_A_0, m_B_0, c_rh_0, c_t_0, c_m_0, c_q_0, c_R_0, c_A_0, c_B_0, et_0, em_0, q_0, R_0, A_0, B_0, rh_0_uM, rt_0, rd_0, z_rh_0, z_t_0, z_m_0, z_q_0, z_R_0, z_A_0, z_B_0, si_0, a_0_uM]