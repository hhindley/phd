using OrderedCollections, StaticArrays
L = 10; 
c = 0.001; 
kr = 0.125; 
Vmax_init = 39.51; 
Km_init = 250; 
θtscr = 160.01;  
θtlr = 255.73; 
k_b = 17.7; 
na = 338; 
nb = 408; 
nr = 532*6;
d = 0.2; 
krep = 137; 
ktag = 0.1; 
atp = 2500; 
km_a = 20; 
km_b = 16;
g_max = 2.0923; 
gr_c = 0.0008856; # 0.000599; 
kdeg = 0.001; 
kin = 0.054; #2.381 
ω_ab = 0.01#4#0.093; #0.0828304057748932;#4; 
ω_r = 0.0019*6 #70.53; #0.0019*6#79.43865871861044; #0.0019*6;  
ω_a = 4; 
ω_b = 4;
kdam = 0.000147;#0.05; 
k = 2; # carrying capacity - changes depending on the data?

rtca_0 = 0.00894; 
rtcb_0 = 0.0216; 
rh_0 = 11.29; #69.56; #69.4
rtcr_0 = 0.04; # 8.67e-3; # change this based on keeping steady state level the whole time course (levels shouldn't really change)
rm_a_0 = 0; 
rm_b_0 = 0; 
rm_r_0 = 0.04 # 0; 
rd_0 = 0; 
rt_0 = 0;

OD_0 = 0.01;
OD_0_wt2 = 0.109;
OD_0_wt3 = 0;
OD_0_wt4 = 0.06;
OD_0_hpxon = 0;
OD_0_wtcolD = 0.0623;
OD_0_nAcolD = 0.0577;
OD_0_nB_colD = 0.061;
OD_0_nBBcolD = 0.0597;
OD_0_bBBmutcolD = 0.055;
OD_0_nRcolD = 0.061;

# param_dict = OrderedDict("L"=>L, "c"=>c, "kr"=>kr, "Vmax_init"=>Vmax_init, "Km_init"=>Km_init, "ω_ab"=>ω_ab, "ω_r"=>ω_r, "θtscr"=>θtscr, "g_max"=>g_max, "θtlr"=>θtlr, "km_a"=>km_a, "km_b"=>km_b, "gr_c"=>gr_c, "d"=>d, "krep"=>krep, "kdam"=>kdam, "ktag"=>ktag, "kdeg"=>kdeg, "kin"=>kin, "atp"=>atp, "na"=>na, "nb"=>nb, "nr"=>nr);
# param_dict_ko = OrderedDict("L"=>L, "c"=>c, "kr"=>kr, "Vmax_init"=>Vmax_init, "Km_init"=>Km_init, "ω_a"=>ω_a, "ω_b"=>ω_b, "ω_r"=>ω_r, "θtscr"=>θtscr, "g_max"=>g_max, "θtlr"=>θtlr, "km_a"=>km_a, "km_b"=>km_b, "gr_c"=>gr_c, "d"=>d, "krep"=>krep, "kdam"=>kdam, "ktag"=>ktag, "kdeg"=>kdeg, "kin"=>kin, "atp"=>atp, "na"=>na, "nb"=>nb, "nr"=>nr, "k"=>k);
# param_dict_OD =  OrderedDict("L"=>L, "c"=>c, "kr"=>kr, "Vmax_init"=>Vmax_init, "Km_init"=>Km_init, "ω_ab"=>ω_ab, "ω_r"=>ω_r, "θtscr"=>θtscr, "g_max"=>g_max, "θtlr"=>θtlr, "km_a"=>km_a, "km_b"=>km_b, "gr_c"=>gr_c, "d"=>d, "krep"=>krep, "kdam"=>kdam, "ktag"=>ktag, "kdeg"=>kdeg, "kin"=>kin, "atp"=>atp, "na"=>na, "nb"=>nb, "nr"=>nr, "k"=>k);
# params = (@SVector [values(param_dict)])[1]
initial = @SVector [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, rh_0, rd_0, rt_0];
init_OD = @SVector [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, rh_0, rd_0, rt_0, OD_0];

# # params_la = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)

tspan = (0, 1e9);