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

params_for_ssval_setup = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
kdeg = 0.001, kin = 0.022222222, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
kdam =  0.01, lam = 0.014)
params1 = @LArray [L, c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 1., ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)

tspan=(0,1e9)


kdam_range_offon = range(0.636,2.0175,length=250)
svals_offon = DataFrame(rm_a=[],rtca=[],rm_b=[],rtcb=[],rm_r=[],rtcr=[],rh=[],rd=[],rt=[])
for kdam_val in ProgressBar(kdam_range_offon)
    psm = deepcopy(params1)
    psm.kdam = kdam_val
    branches1 = setup_ssvals_from_bfkit(kdam_val, params_for_ssval_setup)
    # @show psm
    
    n = 10000; l = 15000;
    lower_ranges = get_all_ranges(set_ss_range_Nssval, branches1, "ss_val_off", n, l)
    # @show upper_ranges[9]
    all, init_vals = get_rh_init_switch_all_ranges(lower_ranges, branches1.ss_val_off,:rh,l,psm)
    binary = upper_or_lower(all, branches1.ss_val_off[7], l)
    inds = get_switch_ind(binary, 0)
    vals = get_switch_vals(inds, init_vals)
    push!(svals_offon.rm_a, vals[1])
    push!(svals_offon.rtca, vals[2])
    push!(svals_offon.rm_b, vals[3])
    push!(svals_offon.rtcb, vals[4])
    push!(svals_offon.rm_r, vals[5])
    push!(svals_offon.rtcr, vals[6])
    push!(svals_offon.rh, vals[7])
    push!(svals_offon.rd, vals[8])
    push!(svals_offon.rt, vals[9])
end
# switch_ind(dashed green)/(on state(green)-off state(red))
CSV.write("/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/init_switch/off_on/switch_vals.csv", svals_offon)

