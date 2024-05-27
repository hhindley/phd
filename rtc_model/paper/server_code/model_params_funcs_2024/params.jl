
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

k_inhib1_val = 1
k_inhib2_val = 0.0025
inhib_val = 0.1
k_inhib1a_val = 0.3
k_inhib1b_val = 0.5

k_inhib_vals = [k_inhib1a_val, k_inhib1b_val, k_inhib1_val]


tspan = (0, 1e9);
