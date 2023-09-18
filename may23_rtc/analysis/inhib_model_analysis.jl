using Parameters, CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, Statistics
using Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra, Printf
# using Plots
using PlotlyJS, ProgressBars
include("/home/holliehindley/phd/may23_rtc/functions/solving.jl"); include("/home/holliehindley/phd/may23_rtc/functions/set_ups.jl"); include("/home/holliehindley/phd/may23_rtc/functions/plotting.jl"); 
include("/home/holliehindley/phd/may23_rtc/functions/sweep_params.jl"); include("/home/holliehindley/phd/may23_rtc/models/rtc_orig.jl"); include("/home/holliehindley/phd/may23_rtc/models/atp_lam_kin_t.jl"); 
include("/home/holliehindley/phd/may23_rtc/models/single_t.jl"); include("/home/holliehindley/phd/may23_rtc/models/combinations_t.jl"); 
include("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/bf_funcs.jl");

include("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/init_switch/funcs.jl"); include("/home/holliehindley/phd/may23_rtc/models/rtc_inhibition_model.jl");


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

tspan = (0,1e9)
k_inhib1 = 1
k_inhib2 = 0.025
inhib = 0.1
params_inhib = @LArray [L, c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 0.5, ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014, k_inhib1, k_inhib2, inhib] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam, :k_inhib1, :k_inhib2, :inhib)

initial_i = [0,0,0,0,0,0,11.29,0,0,0]


solu_rtcb = sol(rtc_inhib_model_rtcb, initial_i, tspan, params_inhib)
p1_rtcb = plotly_plot_sol(solu_rtcb, "", "", "RtcB inhib");

solu_rtca = sol(rtc_inhib_model_rtca, initial_i, tspan, params_inhib)
p1_rtca = plotly_plot_sol(solu_rtca, "", "", "RtcA inhib")

params1 = @LArray [L, c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 0.5, ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
solu_orig = sol(rtc_model, initial, tspan, params1)
p2 = plotly_plot_sol(solu_orig, "", "", "orig");


p_species_inhib = [p2 p1_rtca p1_rtcb]

p_alpha, p_fa, p_ra, p_vinit, p_tscr_a, p_tscr_b, p_tlr_a, p_tlr_b, p_tlr_r, p_rtca1, p_rtcb1, p_vrep, p_vdam, p_vtag = plot_all_vars(solu_orig)
p_alpha1, p_fa1, p_ra1, p_vinit1, p_tscr_a1, p_tscr_b1, p_tlr_a1, p_tlr_b1, p_tlr_r1, p_rtca11, p_rtcb11, p_vrep1, p_vdam1, p_vtag1 = plot_all_vars(solu_rtca)
p_alpha2, p_fa2, p_ra2, p_vinit2, p_tscr_a2, p_tscr_b2, p_tlr_a2, p_tlr_b2, p_tlr_r2, p_rtca12, p_rtcb12, p_vrep2, p_vdam2, p_vtag2 = plot_all_vars(solu_rtcb)

p_orig_all = plot([p_alpha, p_fa, p_ra, p_vinit, p_tscr_a, p_tscr_b, p_tlr_a, p_tlr_b, p_tlr_r, p_rtca1, p_rtcb1, p_vrep, p_vdam, p_vtag], Layout(title="orig",xaxis=attr(range=(0, 800))))
p_rtca_all = plot([p_alpha1, p_fa1, p_ra1, p_vinit1, p_tscr_a1, p_tscr_b1, p_tlr_a1, p_tlr_b1, p_tlr_r1, p_rtca11, p_rtcb11, p_vrep1, p_vdam1, p_vtag1], Layout(title="RtcA",xaxis=attr(range=(0, 800))))
p_rtcb_all = plot([p_alpha2, p_fa2, p_ra2, p_vinit2, p_tscr_a2, p_tscr_b2, p_tlr_a2, p_tlr_b2, p_tlr_r2, p_rtca12, p_rtcb12, p_vrep2, p_vdam2, p_vtag2], Layout(title="RtcB",xaxis=attr(range=(0, 800))))

p_all_inhib = [p_orig_all p_rtca_all p_rtcb_all]

open("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/p_species_inhib.html", "w") do io
    PlotlyBase.to_html(io, p_species_inhib.plot)
end
open("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/p_all_inhib.html", "w") do io
    PlotlyBase.to_html(io, p_all_inhib.plot)
end

solu_rtca = sol(rtc_inhib_model_rtca, initial_i, tspan, params_inhib)
p1_rtca = plotly_plot_sol(solu, "", "", "RtcA inhib")


solu_rtcb = sol(rtc_inhib_model_rtcb, initial_i, tspan, params_inhib)
p1_rtcb = plotly_plot_sol(solu, "", "", "RtcB inhib");

atp = 3578.9473684210525
rtcb = get_curve(solu, :rtcb)
rt = get_curve(solu, :rt)
rtcb_i = get_curve(solu, :rt)
rtcb1 = @. (atp*rtcb)/(atp+(km_b*rt))
Vrep = @. krep*rtcb1*rt

rd = get_curve(solu, :rd)
rtca = get_curve(solu, :rtca)
rtca1 = @. (atp*rtca)/(atp+(km_a*rd)) 
Vtag = @. ktag*rtca1*rd


params1 = @LArray [L, c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 0.5, ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
solu_orig = sol(rtc_model, initial, tspan, params1)
p2 = plotly_plot_sol(solu_orig, "", "", "orig");

rtcb_orig = get_curve(solu_orig, :rtcb)
rt_orig = get_curve(solu_orig, :rt)
rtcb1_orig = @. (atp*rtcb_orig)/(atp+(km_b*rt_orig)) 
Vrep_orig = @. krep*rtcb1_orig*rt_orig

rd_orig = get_curve(solu_orig, :rd)
rtca_orig = get_curve(solu_orig, :rtca)
rtca1_orig = @. (atp*rtca_orig)/(atp+(km_a*rd_orig)) 
Vtag_orig = @. ktag*rtca1_orig*rd_orig

[p1 p2]

plot([scatter(x=solu.t, y=rtcb_i, name="Inactive RtcB"), scatter(x=solu.t, y=rtcb, name="RtcB"), scatter(x=solu_orig.t, y=rtcb_orig, name="RtcB orig")], Layout(xaxis=attr(range=(0, 800))))

plot([scatter(x=solu.t, y=rtcb1, name="RtcB1 inhib"), scatter(x=solu_orig.t, y=rtcb1_orig, name="RtcB1 orig")], Layout(xaxis=attr(range=(0, 800))))
plot([scatter(x=solu.t, y=Vrep, name="Vrep inhib"), scatter(x=solu_orig.t, y=Vrep_orig, name="Vrep orig")], Layout(xaxis=attr(range=(0, 800))))

plot([scatter(x=solu.t,y=rt, name="Rt inhib"), scatter(x=solu_orig.t, y=rt_orig, name="Rt orig")], Layout(xaxis=attr(range=(0, 800))))

plot([scatter(x=solu.t, y=Vtag, name="Vtag inhib"), scatter(x=solu_orig.t, y=Vtag_orig, name="Vtag orig")], Layout(xaxis=attr(range=(0, 800))))



rtcb_ss=[]
rt_ss=[]
rtcb1_ss=[]
Vrep_ss=[]
rtcb_ss_orig=[]
rt_ss_orig=[]
rtcb1_ss_orig=[]
Vrep_ss_orig=[]
for i in range(0,3,length=50)
    params_inhib_new = deepcopy(params_inhib)
    params1_new = deepcopy(params1)
    params_inhib_new.kdam = i
    params1_new.kdam = i
    solu = sol(rtc_inhib_model, initial_i, tspan, params_inhib_new)
    solu_orig = sol(rtc_model, initial, tspan, params1_new)
    rtcb = get_ssval(solu, :rtcb)
    rt = get_ssval(solu, :rt)
    rtcb1 = @. (atp*rtcb)/(atp+(km_b*rt)) 
    Vrep = @. krep*rtcb1*rt
    rtcb_orig = get_ssval(solu_orig, :rtcb)
    rt_orig = get_ssval(solu_orig, :rt)
    rtcb1_orig = @. (atp*rtcb_orig)/(atp+(km_b*rt_orig)) 
    Vrep_orig = @. krep*rtcb1_orig*rt_orig
    push!(rtcb_ss, rtcb)
    push!(rt_ss, rt)
    push!(rtcb1_ss, rtcb1)
    push!(Vrep_ss, Vrep)
    push!(rtcb_ss_orig, rtcb_orig)
    push!(rt_ss_orig, rt_orig)
    push!(rtcb1_ss_orig, rtcb1_orig)
    push!(Vrep_ss_orig, Vrep_orig)
end

rtcb_ss

plot([scatter(x=range(0,3,length=50),y=rtcb_ss,name="RtcB inhib"), scatter(x=range(0,3,length=50),y=rtcb_ss_orig,name="RtcB orig")])
plot([scatter(x=range(0,3,length=50),y=rt_ss,name="Rt inhib"), scatter(x=range(0,3,length=50),y=rt_ss_orig,name="Rt orig")])
plot([scatter(x=range(0,3,length=50),y=rtcb1_ss,name="RtcB1 inhib"), scatter(x=range(0,3,length=50),y=rtcb1_ss_orig,name="RtcB1 orig")])
plot([scatter(x=range(0,3,length=50),y=Vrep_ss,name="Vrep inhib"), scatter(x=range(0,3,length=50),y=Vrep_ss_orig,name="Vrep orig")])



kinhib_range = 10 .^range(log10(0.0001),log10(1),length=5)
inhib_range = 10 .^range(log10(0.001),log10(10),length=5)

sweep_paramx2_new(rtc_inhib_model, :rtcb, get_ssval, :k_inhib, :inhib, kinhib_range, inhib_range, initial, params_inhib, "", "")

kdam_r = range(0,3,length=100)
res=[]
for i in kinhib_range
    p = deepcopy(params_inhib)
    p.k_inhib = i
    rtcb_ss=[]
    for kdam in kdam_r
        p.kdam = kdam
        solu = sol(rtc_inhib_model, initial, tspan, p)
        push!(rtcb_ss, get_ssval(solu, :rtcb))
    end
    push!(res, rtcb_ss)
end
plot([scatter(x=kdam_r, y=res[1], name="$(kinhib_range[1])"),scatter(x=kdam_r, y=res[2], name="$(kinhib_range[2])"),scatter(x=kdam_r, y=res[3], name="$(kinhib_range[3])"),scatter(x=kdam_r, y=res[4], name="$(kinhib_range[4])"),scatter(x=kdam_r, y=res[5], name="$(kinhib_range[5])")])

res2=[]
for i in inhib_range
    p = deepcopy(params_inhib)
    p.inhib = i
    rtcb_ss=[]
    for kdam in kdam_r
        p.kdam = kdam
        solu = sol(rtc_inhib_model, initial, tspan, p)
        push!(rtcb_ss, get_ssval(solu, :rtcb))
    end
    push!(res2, rtcb_ss)
end
res
plot([scatter(x=kdam_r, y=res2[1]),scatter(x=kdam_r, y=res2[2]),scatter(x=kdam_r, y=res2[3]),scatter(x=kdam_r, y=res2[4]),scatter(x=kdam_r, y=res2[5])])







rtcb_ss2=[]
for i in inhib_range
    p = deepcopy(params_inhib)
    p.inhib = i
    solu = sol(rtc_inhib_model, initial, tspan, p)
    push!(rtcb_ss2, get_ssval(solu, :rtcb))
end

plot(scatter(x=inhib_range, y=rtcb_ss2))









[p1_rtca p1_rtcb p2]



solu_rtca = sol(rtc_inhib_model_rtca, initial_i, tspan, params_inhib)
p1_rtca = plotly_plot_sol(solu, "", "", "RtcA inhib")

atp = 3578.9473684210525
rtcb = get_curve(solu, :rtcb)
rt = get_curve(solu, :rt)
rtcb_i = get_curve(solu, :rt)
rtcb1 = @. (atp*rtcb)/(atp+(km_b*rt))
Vrep = @. krep*rtcb1*rt

rd = get_curve(solu, :rd)
rtca = get_curve(solu, :rtca)
rtca1 = @. (atp*rtca)/(atp+(km_a*rd)) 
Vtag = @. ktag*rtca1*rd


params1 = @LArray [L, c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 0.5, ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
solu_orig = sol(rtc_model, initial, tspan, params1)
p2 = plotly_plot_sol(solu_orig, "", "", "");

rtcb_orig = get_curve(solu_orig, :rtcb)
rt_orig = get_curve(solu_orig, :rt)
rtcb1_orig = @. (atp*rtcb_orig)/(atp+(km_b*rt_orig)) 
Vrep_orig = @. krep*rtcb1_orig*rt_orig

rd_orig = get_curve(solu_orig, :rd)
rtca_orig = get_curve(solu_orig, :rtca)
rtca1_orig = @. (atp*rtca_orig)/(atp+(km_a*rd_orig)) 
Vtag_orig = @. ktag*rtca1_orig*rd_orig

[p1 p2]

plot([scatter(x=solu.t, y=rtcb_i, name="Inactive RtcB"), scatter(x=solu.t, y=rtcb, name="RtcB"), scatter(x=solu_orig.t, y=rtcb_orig, name="RtcB orig")], Layout(xaxis=attr(range=(0, 800))))

plot([scatter(x=solu.t, y=rtcb1, name="RtcB1 inhib"), scatter(x=solu_orig.t, y=rtcb1_orig, name="RtcB1 orig")], Layout(xaxis=attr(range=(0, 800))))
plot([scatter(x=solu.t, y=Vrep, name="Vrep inhib"), scatter(x=solu_orig.t, y=Vrep_orig, name="Vrep orig")], Layout(xaxis=attr(range=(0, 800))))

plot([scatter(x=solu.t,y=rt, name="Rt inhib"), scatter(x=solu_orig.t, y=rt_orig, name="Rt orig")], Layout(xaxis=attr(range=(0, 800))))

plot([scatter(x=solu.t, y=Vtag, name="Vtag inhib"), scatter(x=solu_orig.t, y=Vtag_orig, name="Vtag orig")], Layout(xaxis=attr(range=(0, 800))))
