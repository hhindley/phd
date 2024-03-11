using Parameters, CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, Statistics
using Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra

using PlotlyJS, ProgressBars
include("/home/hollie_hindley/Documents/may23_rtc/functions/solving.jl"); include("/home/hollie_hindley/Documents/may23_rtc/functions/set_ups.jl"); include("/home/hollie_hindley/Documents/may23_rtc/functions/plotting.jl"); 
include("/home/hollie_hindley/Documents/may23_rtc/functions/sweep_params.jl"); include("/home/hollie_hindley/Documents/may23_rtc/models/rtc_orig.jl"); include("/home/hollie_hindley/Documents/may23_rtc/models/atp_lam_kin_t.jl"); 
include("/home/hollie_hindley/Documents/may23_rtc/models/single_t.jl"); include("/home/hollie_hindley/Documents/may23_rtc/models/combinations_t.jl"); 
include("/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/bf_funcs.jl");

include("/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/init_switch/funcs.jl");




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
    # atp = 4000;#2500; 
    km_a = 20; 
    km_b = 16;
    g_max = 2.0923; 
    gr_c = 0.0008856; # 0.000599; 
    kdeg = 0.001; 
    # kin = 0.054; #2.381 
    ω_ab = 4#4#0.093; #0.0828304057748932;#4; 
    ω_r = 0.0019*6#2e-7 #0.0019*6 #70.53; #0.0019*6#79.43865871861044; #0.0019*6;  
    ω_a = 4; 
    ω_b = 4;
    # kdam =  0.#0.000147;#0.05; 
    k = 2; # carrying capacity - changes depending on the data?
    # lam = 0.033;

    # rtca_0 = 0#0.00894; 
    # rtcb_0 = 0#0.0216; 
    # rh_0 = 11.29; #69.56; #69.4
    # rtcr_0 = 0# 0.0131 #0.04; # 8.67e-3; # change this based on keeping steady state level the whole time course (levels shouldn't really change)
    # rm_a_0 = 0; 
    # rm_b_0 = 0; 
    # rm_r_0 = 0#0.0131#0.04 # 0; 
    # rd_0 = 0; 
    # rt_0 = 0;
end

params1 = @LArray [L, c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 1., ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
# kdam_range = range(0.6359851,2.016996,length=500)
# logkdam_range = 10 .^ range(log10(0.6359851),log10(2.016996),length=500)
# kdam_range1 = vcat(range(0.6359851,0.643,length=200), range(0.6451,2.016996,length=10))

tspan=(0,1e9)


# lower to upper - double change 
n = 1; l = 50;
branches1 = setup_ssvals_from_bfkit(0.6359851)
ps = deepcopy(params1)
ps.kdam = 0.6359851
lower_ranges = get_all_ranges(set_ss_range_Nssval, branches1, "ss_val_off", n, l)

sss, init1s, init2s, init3s = triple_init(branches1.ss_val_off, lower_ranges, branches1.ss_val_on, "rtcr", "rtcb", "rt", ps)
# # findall(x->x==branches1.ss_val_upper[7],sss)

# p = plot_binary(sss, set_shared_range_0ton(n,l), branches1.ss_val_off[7], "rtcr", "rt", "lower to upper - kdam = 1")

# # # open("$PATHmay23_rtc/analysis/bifurcation_analysis/plots/init_switch_lowertoupper1.html", "w") do io
# # #     PlotlyBase.to_html(io, p.plot)
# # # end
# binary = get_binary(branches1.ss_val_off[7], sss)
# x, y = get_plot_vals(binary, l, n)
# plot(scatter(x=x,y=y), Layout(xaxis_range=(0,1000), yaxis_range=(0,1000)))

sss


kdam_range = range(0.6359851,2.016996,length=10)
tspan=(0,1e9)
sss_res=[]; branches=[]; xs=[]; ys=[]

for kdam_val in ProgressBar(kdam_range)
    ps = deepcopy(params1)
    ps.kdam = kdam_val
    branches1 = setup_ssvals_from_bfkit(kdam_val)

    
    n = 1; l = 50;
    lower_ranges = get_all_ranges(set_ss_range_Nssval, branches1, "ss_val_off", n, l)

    sss, init1s, init2s, init3s = triple_init(branches1.ss_val_off, lower_ranges, branches1.ss_val_on, "rtcr", "rtcb", "rt", ps)
    push!(sss_res, sss)
end

sss_res

