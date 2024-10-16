using Parameters, LabelledArrays, StaticArrays, CSV, DataFrames, DifferentialEquations, BenchmarkTools, OrderedCollections, DataInterpolations, PlotlyJS, Statistics
include("/home/hollie_hindley/Documents/may23_rtc/functions/solving.jl"); include("/home/hollie_hindley/Documents/may23_rtc/functions/plotting.jl"); include("/home/hollie_hindley/Documents/may23_rtc/functions/sweep_params.jl"); include("/home/hollie_hindley/Documents/may23_rtc/models/rtc_orig.jl"); include("/home/hollie_hindley/Documents/may23_rtc/models/atp_lam_kin_t.jl"); include("/home/hollie_hindley/Documents/may23_rtc/analysis/t_param_setup.jl");

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
    kdam =  0.0#0.000147;#0.05; 
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
end
params_kdam = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
params_kdam = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, 0, ktag, kdeg, kin, atp, na, nb, nr, lam] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)

initial = @SVector [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, rh_0, rd_0, rt_0]




L_range = 10 .^(range(1,stop=3,length=101))
c_range = 10 .^(range(-4,stop=0,length=101))
wab_range = 10 .^range(-8,stop=-1,length=101)
wr_range = 10 .^(range(-7,stop=-2,length=101))
atp_range = range(500,stop=5000,length=101)
kin_range = range(0,stop=0.2,length=101)
lam_range = range(0.001,stop=0.04,length=101)
kdam_range = 10 .^ range(-4,stop=0,length=101)

initial = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, rh_0, rd_0, rt_0]
new_params=deepcopy(params1)
res=[]
for i in atp_range
    new_params.atp = i 
    solu = sol(rtc_model, initial, tspan, new_params)
    push!(res, get_ssval(solu, :rh))
end
res


plot(scatter(x=atp_range, y=res))


param_change = [:L, :c, :ω_ab, :ω_r, :atp, :kin, :lam, :kdam]
param_ranges = [L_range, c_range, wab_range, wr_range, atp_range, kin_range, lam_range, kdam_range]


for (i, j) in zip(param_change, param_ranges)
    save_1x_plots(j, i, "kdam = 0.0", "log", initial, get_ssval)
end

kdam_range1 = range(0, stop=1, length=10)
all_curves = []
for j in all_species
    curves = []
    for i in kdam_range1
        params[:kdam] = i
        solu = sol(rtc_model, initial, tspan, params)
        specie = get_curve(solu, j)
        push!(curves, scatter(x=solu.t, y=specie))
    end
    push!(all_curves, curves)
end

for curve in all_curves
    display(plot([i for i in curve], Layout(xaxis_type="log")))
end
# c_range = collect(0:0.01:1)
# # c_results = change_param(c_range, :c, rtc_model, initial, all_species, lam, atp, kin)
# # plot_change_param_sols(c_range, c_results, "c", "")
# # # plot_all_change_param(c_range, c_results)

# wab_range = collect(0:0.01:1)
# # wab_results = change_param(wab_range, :ω_ab, rtc_model, initial, all_species, lam, atp, kin)
# # plot_change_param_sols(wab_range, wab_results, "ω_ab", "")

# wr_range = collect(0:0.01:1)
# # wr_results = change_param(wr_range, :ω_r, rtc_model, initial, all_species, lam, atp, kin)
# # plot_change_param_sols(wr_range, wr_results, "ω_r", "")

# atp_range = collect(0:50:5000)
# # atp_results = change_param(atp_range, :atp, rtc_model, initial, all_species, lam, atp, kin)
# # plot_change_param_sols(atp_range, atp_results, "ATP", "")

# kin_range = collect(0:0.1:10)
# # kin_results = change_param(kin_range, :kin, rtc_model, initial, all_species, lam, atp, kin)
# # plot_change_param_sols(kin_range, kin_results, "kin", "")

# lam_range = collect(0.1:0.01:1.1)
# # lam_results = change_param(lam_range, :lam, rtc_model, initial, all_species, lam, atp, kin)
# # plot_change_param_sols(lam_range, lam_results, "λ", "")

# kdam_range = collect(0:0.01:1)
# kdam_results = change_param(kdam_range, :kdam, rtc_model, initial, all_species, lam, atp, kin)
# plot_change_param_sols(kdam_range, kdam_results, "kdam", "")