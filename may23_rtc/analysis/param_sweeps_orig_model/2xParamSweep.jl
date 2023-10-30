using Parameters, LabelledArrays, StaticArrays, CSV, DataFrames, DifferentialEquations, BenchmarkTools, OrderedCollections, DataInterpolations, PlotlyJS, Statistics
include("/home/holliehindley/phd/may23_rtc/functions/solving.jl"); include("/home/holliehindley/phd/may23_rtc/functions/plotting.jl"); include("/home/holliehindley/phd/may23_rtc/functions/sweep_params.jl"); include("/home/holliehindley/phd/may23_rtc/models/rtc_orig.jl"); include("/home/holliehindley/phd/may23_rtc/models/atp_lam_kin_t.jl"); include("/home/holliehindley/phd/may23_rtc/analysis/t_param_setup.jl");



# include("/home/holliehindley/phd/may23_rtc/functions/solving.jl"); include("/home/holliehindley/phd/may23_rtc/functions/plotting.jl")
# using DifferentialEquations, PlotlyJS, DataFrames, Parameters
# params = (L = 10, c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250, θtscr = 160.01, θtlr = 255.73, na = 338, nb = 408, nr = 532*6, d = 0.2, krep = 137, ktag = 9780, atp = 4000, km_a = 20, km_b = 16, g_max = 2.0923, gr_c = 0.0008856, kdeg = 0.001, kin = 0.054, ω_ab = 4, ω_r = 2e-7, ω_a = 4, ω_b = 4, kdam =  0.01, k = 2, lam = 0.033)
# initial = [0, 0, 0, 0, 0, 0, 11.29, 0, 0]
# tspan=(0,1e9)
# solu = sol(rtc_model, initial, tspan, params)
# plotly_plot_sol(solu, "log", "", "")






@consts begin   
    L = 10; #10 
    c = 0.001; 
    kr = 0.125; 
    Vmax_init = 39.51; 
    Km_init = 250; 
    θtscr = 160.01;  
    θtlr = 255.73; 
    # k_b = 17.7; 
    na = 338; 
    nb = 408; 
    nr = 532*6;
    d = 0.2; 
    krep = 137; 
    ktag = 9780;#0.1; 
    atp = 4000;#2500; 
    km_a = 20; 
    km_b = 16;
    g_max = 2.0923; 
    gr_c = 0.0008856; # 0.000599; 
    kdeg = 0.001; 
    kin = 0.054; #2.381 
    ω_ab = 4#4#0.093; #0.0828304057748932;#4; 
    ω_r = 2e-7 #0.0019*6 #70.53; #0.0019*6#79.43865871861044; #0.0019*6;  
    ω_a = 4; 
    ω_b = 4;
    kdam =  0.01#0.000147;#0.05; 
    k = 2; # carrying capacity - changes depending on the data?
    lam = 0.033;

    rtca_0 = 0#0.00894; 
    rtcb_0 = 0#0.0216; 
    rh_0 = 11.29; #69.56; #69.4
    rtcr_0 = 0# 0.0131 #0.04; # 8.67e-3; # change this based on keeping steady state level the whole time course (levels shouldn't really change)
    rm_a_0 = 0; 
    rm_b_0 = 0; 
    rm_r_0 = 0#0.0131#0.04 # 0; 
    rd_0 = 0; 
    rt_0 = 0;

    # tspan = (0, 1e9);
end
params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
initial = @SVector [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, rh_0, rd_0, rt_0]


L_range = 10 .^(range(1,stop=3,length=101))
c_range = 10 .^(range(-4,stop=0,length=101))
wab_range = 10 .^range(-8,stop=-1,length=101)
wr_range = 10 .^(range(-7,stop=-2,length=101))
atp_range = range(500,stop=5000,length=101)
kin_range = range(0,stop=0.2,length=101)
lam_range = range(0.001,stop=0.04,length=101)
kdam_range = 10 .^ range(-4,stop=0,length=101)

solu = sol(rtc_model, initial, tspan, params)
ss_init = ss_init_vals(solu)

all_species = [:rm_a, :rh, :rd, :rt]

#L and c
save_2x_plots(:L, :c, L_range, c_range, "L_c", ss_init, "log", "log") 


#L and wab 
save_2x_plots(:L, :ω_ab, L_range, wab_range, "L_wab", "log", "log")


#L and wr
save_2x_plots(:L, :ω_r, L_range, wr_range, "L_wr", "log", "log")


#L and atp
save_2x_plots(:L, :atp, L_range, atp_range, "L_atp", "log", "")


#L and kin 
save_2x_plots(:L, :kin, L_range, kin_range, "L_kin", "log", "")


#L and lam
save_2x_plots(:L, :lam, L_range, lam_range, "L_lam", "log", "")
 


#L and kdam 
save_2x_plots(:L, :kdam, L_range, kdam_range, "L_kdam", "log", "log")


#c and wab 
save_2x_plots(:c, :ω_ab, c_range, wab_range, "c_wab", "log", "log")


#c and wr 
save_2x_plots(:c, :ω_r, c_range, wr_range, "c_wr", "log", "log")


#c and atp
save_2x_plots(:c, :atp, c_range, atp_range, "c_atp", "log", "")

#c and kin 
save_2x_plots(:c, :kin, c_range, kin_range, "c_kin", "log", "")


#c and lam 
save_2x_plots(:c, :lam, c_range, lam_range, "c_lam", "log", "")

#c and kdam 
save_2x_plots(:c, :kdam, c_range, kdam_range, "c_kdam", "log", "log")



#wab and wr 
save_2x_plots(:ω_ab, :ω_r, wab_range, wr_range, "wab_wr", ss_init,"log", "log")

#wab and atp 
save_2x_plots(:ω_ab, :atp, wab_range, atp_range, "wab_atp", ss_init,"log", "")

#wab and kin
save_2x_plots(:ω_ab, :kin, wab_range, kin_range, "wab_kin", ss_init,"log", "")

#wab and lam 
save_2x_plots(:ω_ab, :lam, wab_range, lam_range, "wab_lam", ss_init,"log", "")

#wab and kdam 
save_2x_plots(:ω_ab, :kdam, wab_range, kdam_range, "wab_kdam",ss_init, "log", "log")


#wr and atp 
save_2x_plots(:ω_r, :atp, wr_range, atp_range, "wr_atp",ss_init, "log", "")

#wr and kin 
save_2x_plots(:ω_r, :kin, wr_range, kin_range, "wr_kin",ss_init, "log", "")

#wr and lam 
save_2x_plots(:ω_r, :lam, wr_range, lam_range, "wr_lam",ss_init, "log", "")

#wr and kdam 
save_2x_plots(:ω_r, :kdam, wr_range, kdam_range, "wr_kdam",ss_init, "log", "log")


#atp and kin 
save_2x_plots(:atp, :kin, atp_range, kin_range, "atp_kin", ss_init,"", "")

#atp and lam 
save_2x_plots(:atp, :lam, atp_range, lam_range, "atp_lam", ss_init,"", "")

#atp and kdam 
save_2x_plots(:atp, :kdam, atp_range, kdam_range, "atp_kdam",ss_init, "", "log")


#kin and lam 
save_2x_plots(:kin, :lam, kin_range, lam_range, "kin_lam",ss_init, "", "")

#kin and kdam 
save_2x_plots(:kin, :kdam, kin_range, kdam_range, "kin_kdam",ss_init, "", "log")


#lam and kdam 
save_2x_plots(:lam, :kdam, lam_range, kdam_range, "lam_kdam", ss_init,"", "log")


sweep_paramx2_new(rtc_model, :rt, get_ssval, :atp, :lam, atp_range, lam_range)


atp_range = range(1000,stop=5000,length=101)
kdam_range = 10 .^ range(-4,stop=0,length=101)

sweep_paramx2_new(rtc_model, :rm_a, get_ssval, :kin, :kdam, kin_range, kdam_range)
sweep_paramx2_new(rtc_model, :rt, get_ssval, :kin, :kdam, kin_range, kdam_range)



sweep_paramx2_new(rtc_model, :rm_a, get_ssval, :atp, :kdam, atp_range, kdam_range)
sweep_paramx2_new(rtc_model, :rtca, get_ssval, :atp, :kdam, atp_range, kdam_range)
sweep_paramx2_new(rtc_model, :rt, get_ssval, :atp, :kdam, atp_range, kdam_range)
sweep_paramx2_new(rtc_model, :rh, get_ssval, :atp, :kdam, atp_range, kdam_range)




L_range = (range(0,stop=400,length=101))
c_range = 10 .^(range(-4,stop=0,length=101))
sweep_paramx2_new(rtc_model, :rm_a, get_ssval, :L, :c, L_range, c_range)
sweep_paramx2_new(rtc_model, :rt, get_ssval, :L, :c, L_range, c_range)


sweep_paramx2_new(rtc_model, :rt, get_ssval, :lam, :atp, lam_range, atp_range)
sweep_paramx2_new(rtc_model, :rm_a, get_ssval, :lam, :atp, lam_range, atp_range)



sweep_paramx2_new(rtc_model, :rm_r, get_ssval, :ω_r, :kdam, wr_range, kdam_range)
