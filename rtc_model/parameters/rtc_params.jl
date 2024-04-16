
L_val = 10.; #10 # unitless 
c_val = 0.001; # unitless
kr_val = 0.125; # uM
Vmax_init_val = 39.51; # min-1 
Km_init_val = 250.; # uM
θtscr_val = 160.01; # uM
θtlr_val = 255.73; # uM
# k_b = 17.7; 
nA_val = 338.; # aa molec-1
nB_val = 408.; # aa molec-1
nR_val = (532.)*6; # aa molec-1
d_val = 0.2; # min-1 
krep_val = 137.; # uM-1 min-1 ???
ktag_val = 9780.;#0.1; # uM-1 min-1 ???
atp_val = 3000.;#4000;#2500; # uM
km_a_val = 20.; # uM
km_b_val = 16.; # uM
g_max_val = 1260.; # aa min-1 molec-1
kdeg_val = 0.001; # min-1
kin_val = 0.022/100 #2.381 # uM molec aa-1 ??? 
ω_ab_val = 1.e-5; #0.056/1e2#4#0.093; #0.0828304057748932;#4; # unitless 
ω_r_val = 1.e-6; #0.01/1e4 #0.0019*6 #70.53; #0.0019*6#79.43865871861044; #0.0019*6;  # uM min-1
# ω_a = 4; 
# ω_b = 4;
kdam_val =  0.0#0.000147;#0.05; # min-1
# k = 2; # carrying capacity - changes depending on the data?
lam_val = 0.014; # min-1
sf = 1e6/(6.022e23*1e-15) # uM molec-1
kb_gm_val = 1/sf # min-1 uM-1
ku_gm_val = 1. # min-1
kc_val = 0.6 # kc (uM-1) is the binding constant (kb/ku) (uM-1) for ribosomes to mRNAs - kc/sf = uM-1 
k_diss_val = 0.006 # min-1

k_inhib1_val = 1 #0.05
k_inhib2_val = 0.0025
inhib_val = 0.1
k_inhib1a_val = 0.3
k_inhib1b_val = 0.5

# k_inhib_vals = [0.1, 0.15, 0.2]

k_inhib_vals = [k_inhib1a_val, k_inhib1b_val, k_inhib1_val]


tspan = (0, 1e9);

# params_rtc = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, nA, nB, nR, lam, kc, k_diss] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam, :kc, :k_diss)
# params_inhib = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, nA, nB, nR, lam, kc, k_inhib1, k_inhib2, inhib] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam, :kc, :k_inhib1, :k_inhib2, :inhib)


# species_rtc = [:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt]
# species_inhib = [:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt, :rtc_i]

# solu = sol(rtc_model, init_rtc, tspan, params_rtc)
# ssvals_rtc = get_all_ssvals(solu, species_rtc)



# params_bf = (L = L, c = c, kr = kr, Vmax_init = Vmax_init, Km_init = Km_init,
# ω_ab = ω_ab, ω_r = ω_r, θtscr = θtscr, g_max = g_max, θtlr = θtlr, km_a = km_a, km_b = km_b, 
# d = d, krep = krep, kdam = kdam, ktag = ktag, kdeg = kdeg, kin = kin, atp = atp,
# na = nA, nb = nB, nr = nR, lam = lam, kc = kc, k_diss = k_diss) 

# params_bf_inhib = (L = L, c = c, kr = kr, Vmax_init = Vmax_init, Km_init = Km_init,
# ω_ab = ω_ab, ω_r = ω_r, θtscr = θtscr, g_max = g_max, θtlr = θtlr, km_a = km_a, km_b = km_b, 
# d = d, krep = krep, kdam = kdam, ktag = ktag, kdeg = kdeg, kin = kin, atp = atp,
# na = nA, nb = nB, nr = nR, lam = lam, kc = kc, k_diss = k_diss,
# k_inhib1=k_inhib1, k_inhib2=k_inhib2, inhib=inhib) 

#atp = 3578.9473684210525

# old bf params for first paper draft before finding error in model 
# @consts begin
#     L = 10; #10 
#     c = 0.001; 
#     kr = 0.125; 
#     Vmax_init = 39.51; 
#     Km_init = 250; 
#     θtscr = 160.01;  
#     θtlr = 255.73; 
#     # k_b = 17.7; 
#     na = 338; 
#     nb = 408; 
#     nr = 532*6;
#     d = 0.2; 
#     krep = 137; 
#     ktag = 9780;#0.1; 
#     # atp = 4000;#2500; 
#     km_a = 20; 
#     km_b = 16;
#     g_max = 2.0923; 
#     gr_c = 0.0008856; # 0.000599; 
#     kdeg = 0.001; 
#     # kin = 0.054; #2.381 
#     ω_ab = 4#4#0.093; #0.0828304057748932;#4; 
#     ω_r = 0.0019*6#2e-7 #0.0019*6 #70.53; #0.0019*6#79.43865871861044; #0.0019*6;  
#     ω_a = 4; 
#     ω_b = 4;
#     # kdam =  0.#0.000147;#0.05; 
#     k = 2; # carrying capacity - changes depending on the data?
#     # lam = 0.033;

#     # rtca_0 = 0#0.00894; 
#     # rtcb_0 = 0#0.0216; 
#     # rh_0 = 11.29; #69.56; #69.4
#     # rtcr_0 = 0# 0.0131 #0.04; # 8.67e-3; # change this based on keeping steady state level the whole time course (levels shouldn't really change)
#     # rm_a_0 = 0; 
#     # rm_b_0 = 0; 
#     # rm_r_0 = 0#0.0131#0.04 # 0; 
#     # rd_0 = 0; 
#     # rt_0 = 0;
# end

# tspan = (0,1e9)

# params2 = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
# θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
# krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
# kdeg = 0.001, kin = 0.022222222, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
# kdam =  0.01, lam = 0.014)

# params_ = @LArray [L, c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 0.5, ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
