using Parameters, CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, Statistics
# using Plots
using PlotlyJS
include("/home/holliehindley/phd/may23_rtc/functions/solving.jl"); include("/home/holliehindley/phd/may23_rtc/functions/set_ups.jl"); include("/home/holliehindley/phd/may23_rtc/functions/plotting.jl"); 
include("/home/holliehindley/phd/may23_rtc/functions/sweep_params.jl"); include("/home/holliehindley/phd/may23_rtc/models/rtc_orig.jl"); include("/home/holliehindley/phd/may23_rtc/models/atp_lam_kin_t.jl"); 
include("/home/holliehindley/phd/may23_rtc/models/single_t.jl"); include("/home/holliehindley/phd/may23_rtc/models/combinations_t.jl"); include("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/init_switch/funcs.jl");


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

branches1 = setup_ssvals(params1)

tspan=(0,1e9)

# upper to lower - single change - don't see anything
n1 = 1000; l1= 1000;
upper_ranges = get_all_ranges(set_ss_range_zerotossval, branches1, "ss_val_upper", n1,l1)
all1 = get_rh_init_switch_all_ranges(upper_ranges, branches1.ss_val_upper, :rh, l1, params1)
dfs1 = upper_or_lower(all1, "upper",l1)
shared_range1 = set_shared_range(n1,l1)

p = plot([scatter(x=shared_range1,y=dfs1.rm_a, name="rm_a"),scatter(x=shared_range1,y=dfs1.rtca, name="rtca"),scatter(x=shared_range1,y=dfs1.rm_b, name="rm_b"),
scatter(x=shared_range1,y=dfs1.rtcb,name="rtcb"),scatter(x=shared_range1,y=dfs1.rm_r,name="rm_r"),scatter(x=shared_range1,y=dfs1.rtcr,name="rtcr"),
scatter(x=shared_range1,y=dfs1.rh,name="rh"),scatter(x=shared_range1,y=dfs1.rd,name="rd"),scatter(x=shared_range1,y=dfs1.rt,name="rt")],
Layout(xaxis_title="% increase from initial ss val", yaxis_title="Branch",
yaxis_tickvals=[0,1], yaxis_ticktext=["lower","upper"], title="switch from upper to lower (rh)"
))




# upper to lower - double change - where we see change rtcr and rtcb
n1 = 0; l1= 100;
upper_ranges = get_all_ranges(set_ss_range_zerotossval, branches1, "ss_val_upper", n1,l1)
sss, rtcbs, rtcrs = double_init(branches1.ss_val_upper, upper_ranges, branches1.ss_val_lower, "rtcb", "rtcr", params1)
p = plot_binary(sss, range(-100,0,length=l1), branches1.ss_val_lower, "rtcb", "rtcr", "upper to lower - rh, kdam = 1")


open("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/init_switch_uppertolower.html", "w") do io
    PlotlyBase.to_html(io, p.plot)
end


params2 = @LArray [L, c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 0.64, ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)

branches2 = setup_ssvals(params2)

# upper to lower - single change - don't see anything
n1 = 1000; l1= 100;
upper_ranges = get_all_ranges(set_ss_range_zerotoNssval, branches2, "ss_val_upper", n1,l1)
all1 = get_rh_init_switch_all_ranges(upper_ranges, branches2.ss_val_upper, :rh, l1, params2)
dfs1 = upper_or_lower(all1, "upper",l1)
shared_range1 = set_shared_range(n1,l1)

p = plot([scatter(x=shared_range1,y=dfs1.rm_a, name="rm_a"),scatter(x=shared_range1,y=dfs1.rtca, name="rtca"),scatter(x=shared_range1,y=dfs1.rm_b, name="rm_b"),
scatter(x=shared_range1,y=dfs1.rtcb,name="rtcb"),scatter(x=shared_range1,y=dfs1.rm_r,name="rm_r"),scatter(x=shared_range1,y=dfs1.rtcr,name="rtcr"),
scatter(x=shared_range1,y=dfs1.rh,name="rh"),scatter(x=shared_range1,y=dfs1.rd,name="rd"),scatter(x=shared_range1,y=dfs1.rt,name="rt")],
Layout(xaxis_title="% increase from initial ss val", yaxis_title="Branch",
yaxis_tickvals=[0,1], yaxis_ticktext=["lower","upper"], title="switch from upper to lower (rh)"
))


# upper to lower - double change - where we see change rtcr and rtcb
n1 = 10000; l1= 50;
upper_ranges = get_all_ranges(set_ss_range_zerotoNssval, branches2, "ss_val_upper", n1,l1)
sss, rtcbs, rtcrs = double_init(branches2.ss_val_upper, upper_ranges, branches2.ss_val_lower, "rtcr", "rt", params2)
findall(x->x==branches2.ss_val_lower[7],sss)
sss
p = plot_binary(sss, range(-100,0,length=l1), branches2.ss_val_lower, "rt", "rtcr", "upper to lower - rh")
sss
round(sss[1];digits=6)

branches2
branches2.ss_val_upper[7]

open("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/init_switch_uppertolower.html", "w") do io
    PlotlyBase.to_html(io, p.plot)
end









params3 = @LArray [L, c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 2, ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)

branches3 = setup_ssvals(params3)

# upper to lower - single change - don't see anything
n1 = 1000; l1= 1000;
upper_ranges = get_all_ranges(set_ss_range_zerotoNssval, branches3, "ss_val_upper", n1,l1)
all1 = get_rh_init_switch_all_ranges(upper_ranges, branches3.ss_val_upper, :rh, l1, params3)
dfs1 = upper_or_lower(all1, "upper",l1)
shared_range1 = set_shared_range(n1,l1)

p = plot([scatter(x=shared_range1,y=dfs1.rm_a, name="rm_a"),scatter(x=shared_range1,y=dfs1.rtca, name="rtca"),scatter(x=shared_range1,y=dfs1.rm_b, name="rm_b"),
scatter(x=shared_range1,y=dfs1.rtcb,name="rtcb"),scatter(x=shared_range1,y=dfs1.rm_r,name="rm_r"),scatter(x=shared_range1,y=dfs1.rtcr,name="rtcr"),
scatter(x=shared_range1,y=dfs1.rh,name="rh"),scatter(x=shared_range1,y=dfs1.rd,name="rd"),scatter(x=shared_range1,y=dfs1.rt,name="rt")],
Layout(xaxis_title="% increase from initial ss val", yaxis_title="Branch",
yaxis_tickvals=[0,1], yaxis_ticktext=["lower","upper"], title="switch from upper to lower (rh) - kdam = 2"
))

open("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/init_switch_upperlower2.html", "w") do io
    PlotlyBase.to_html(io, p.plot)
end


# upper to lower - double change - where we see change rtcr and rtcb
n1 = 30; l1= 50;
upper_ranges = get_all_ranges(set_ss_range_zerotoNssval, branches3, "ss_val_upper", n1,l1)
sss, rtcbs, rtcrs = double_init(branches3.ss_val_upper, upper_ranges, branches3.ss_val_lower, "rtcr", "rt", params3)
findall(x->x==branches3.ss_val_lower[7],sss)
shared_range1 = set_shared_range(n1,l1)
p = plot_binary(sss, shared_range1, branches3.ss_val_lower, "rt", "rtcr", "upper to lower - rh, kdam = 2")
sss
round(sss[1];digits=6)

branches2
branches2.ss_val_upper[7]

open("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/init_switch_uppertolower.html", "w") do io
    PlotlyBase.to_html(io, p.plot)
end
