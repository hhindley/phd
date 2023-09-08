using Parameters, CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, Statistics
using Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra

using PlotlyJS, ProgressBars
include("/home/holliehindley/phd/may23_rtc/functions/solving.jl"); include("/home/holliehindley/phd/may23_rtc/functions/set_ups.jl"); include("/home/holliehindley/phd/may23_rtc/functions/plotting.jl"); 
include("/home/holliehindley/phd/may23_rtc/functions/sweep_params.jl"); include("/home/holliehindley/phd/may23_rtc/models/rtc_orig.jl"); include("/home/holliehindley/phd/may23_rtc/models/atp_lam_kin_t.jl"); 
include("/home/holliehindley/phd/may23_rtc/models/single_t.jl"); include("/home/holliehindley/phd/may23_rtc/models/combinations_t.jl"); 
include("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/bf_funcs.jl");

include("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/init_switch/funcs.jl");

function rtc_inhib_model(initial, params, t) 
    L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, k_inhib, inhib = params
    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt = initial


    # dilution by growth and degradation 
    dil(species) = lam*species
    deg(species) = d*species
    
    # MWC
    alpha = rt/kr 
    fa = (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6)
    ra = fa*rtcr
    
    # transcription
    Vinit = ra*Vmax_init*atp/(Km_init+atp)
    tscr_el_a = ω_ab*atp/(θtscr+atp)
    tscr_a = Vinit*tscr_el_a
    tscr_el_b = ω_ab*atp/(θtscr+atp)
    tscr_b = Vinit*tscr_el_b
    tscr_r = ω_r*atp/(θtscr+atp)

    # translation
    tlr_el = g_max*atp/(θtlr+atp)
    tlr(rm_x, nx) = (1/nx)*rh*rm_x*tlr_el # *1/nx nx = length of RtcA, RtcB and RtcR

    # ribosomes
    rtca1 = (atp*rtca)/(atp+(km_a*rd)) 
    rtcb1 = (atp*rtcb)/(atp+(km_b*rt)) 

    # rtcb1 = k_b*atp*rtcb/(k_b*atp+krep*rt)

    Vrep = krep*rtcb1*rt
    Vdam = kdam*rh
    Vinflux = kin*tlr_el
    Vtag = ktag*rtca1*rd

    rtcb_i = inhib*k_inhib*rtcb/((k_inhib*inhib)+k_inhib)
    # ODEs
    drm_a = tscr_a - dil(rm_a) - deg(rm_a)
    drtca = tlr(rm_a, na) - dil(rtca)    
    drm_b = tscr_b - dil(rm_b) - deg(rm_b)
    drtcb = tlr(rm_b, nb) - dil(rtcb) - k_inhib*rtcb*inhib + k_inhib*rtcb_i
    drm_r = tscr_r - dil(rm_r) - deg(rm_r)
    drtcr = tlr(rm_r, nr) - dil(rtcr)
    drh = Vrep - Vdam + Vinflux - dil(rh)
    drd = Vdam - Vtag - kdeg*rd - dil(rd)
    drt = Vtag - Vrep - dil(rt)

    # @SVector [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt]
    [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt]

end


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
k_inhib = 0.001
inhib = 1
params_inhib = @LArray [L, c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 1., ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014, k_inhib, inhib] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam, :k_inhib, :inhib)

# kinhib_range = 10 .^range(log10(0.0001),log10(1),length=5)
# inhib_range = 10 .^range(log10(0.001),log10(10),length=5)

# sweep_paramx2_new(rtc_inhib_model, :rtcb, get_ssval, :k_inhib, :inhib, kinhib_range, inhib_range, initial, params_inhib, "", "")

# kdam_r = range(0,3,length=100)
# res=[]
# for i in kinhib_range
#     p = deepcopy(params_inhib)
#     p.k_inhib = i
#     rtcb_ss=[]
#     for kdam in kdam_r
#         p.kdam = kdam
#         solu = sol(rtc_inhib_model, initial, tspan, p)
#         push!(rtcb_ss, get_ssval(solu, :rtcb))
#     end
#     push!(res, rtcb_ss)
# end
# plot([scatter(x=kdam_r, y=res[1], name="$(kinhib_range[1])"),scatter(x=kdam_r, y=res[2], name="$(kinhib_range[2])"),scatter(x=kdam_r, y=res[3], name="$(kinhib_range[3])"),scatter(x=kdam_r, y=res[4], name="$(kinhib_range[4])"),scatter(x=kdam_r, y=res[5], name="$(kinhib_range[5])")])

# res2=[]
# for i in inhib_range
#     p = deepcopy(params_inhib)
#     p.inhib = i
#     rtcb_ss=[]
#     for kdam in kdam_r
#         p.kdam = kdam
#         solu = sol(rtc_inhib_model, initial, tspan, p)
#         push!(rtcb_ss, get_ssval(solu, :rtcb))
#     end
#     push!(res2, rtcb_ss)
# end
# res
# plot([scatter(x=kdam_r, y=res2[1]),scatter(x=kdam_r, y=res2[2]),scatter(x=kdam_r, y=res2[3]),scatter(x=kdam_r, y=res2[4]),scatter(x=kdam_r, y=res2[5])])







# rtcb_ss2=[]
# for i in inhib_range
#     p = deepcopy(params_inhib)
#     p.inhib = i
#     solu = sol(rtc_inhib_model, initial, tspan, p)
#     push!(rtcb_ss2, get_ssval(solu, :rtcb))
# end

# plot(scatter(x=inhib_range, y=rtcb_ss2))



k_inhib = 0.1
inhib = 1
params_inhib = @LArray [L, c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 1., ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014, k_inhib, inhib] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam, :k_inhib, :inhib)

solu = sol(rtc_inhib_model, initial, tspan, params_inhib)
p1 = plotly_plot_sol(solu, "", "", "");


rtcb = get_curve(solu, :rtcb)
rtcb_i = inhib*k_inhib*rtcb/((k_inhib*inhib)+k_inhib)

y = rtcb_i./(rtcb)
# plot([scatter(x=solu.t, y=rtcb_i, name="Inactive RtcB"), scatter(x=solu.t, y=rtcb, name="RtcB")], Layout(xaxis=attr(range=(0, 800))))

params1 = @LArray [L, c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 1., ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)

solu1 = sol(rtc_model, initial, tspan, params1)
p2 = plotly_plot_sol(solu1, "", "", "");
rtcb1 = get_curve(solu1, :rtcb)
# [p1 p2]

plot([scatter(x=solu.t, y=rtcb_i, name="Inactive RtcB"), scatter(x=solu.t, y=rtcb, name="RtcB"), scatter(x=solu1.t, y=rtcb1, name="RtcB orig")], Layout(xaxis=attr(range=(0, 800))))




function rtc_inhib_mod!(dz, z, p, t)
    @unpack L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, k_inhib, inhib = p
    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt = z


    # dilution by growth and degradation 
    dil(species) = lam*species
    deg(species) = d*species
    
    # MWC
    alpha = rt/kr 
    fa = (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6)
    ra = fa*rtcr
    
    # transcription
    Vinit = ra*Vmax_init*atp/(Km_init+atp)
    tscr_el_a = ω_ab*atp/(θtscr+atp)
    tscr_a = Vinit*tscr_el_a
    tscr_el_b = ω_ab*atp/(θtscr+atp)
    tscr_b = Vinit*tscr_el_b
    tscr_r = ω_r*atp/(θtscr+atp)

    # translation
    tlr_el = g_max*atp/(θtlr+atp)
    tlr(rm_x, nx) = (1/nx)*rh*rm_x*tlr_el # *1/nx nx = length of RtcA, RtcB and RtcR

    # ribosomes
    rtca1 = (atp*rtca)/(atp+(km_a*rd)) 
    rtcb1 = (atp*rtcb)/(atp+(km_b*rt)) 

    # rtcb1 = k_b*atp*rtcb/(k_b*atp+krep*rt)

    Vrep = krep*rtcb1*rt
    Vdam = kdam*rh
    Vinflux = kin*tlr_el
    Vtag = ktag*rtca1*rd

    rtcb_i = inhib*k_inhib*rtcb/((k_inhib*inhib)+k_inhib)


    # ODEs
    dz[1] = tscr_a - dil(rm_a) - deg(rm_a)
    dz[2] = tlr(rm_a, na) - dil(rtca)    
    dz[3] = tscr_b - dil(rm_b) - deg(rm_b)
    dz[4] = tlr(rm_b, nb) - dil(rtcb) - k_inhib*rtcb*inhib + k_inhib*rtcb_i
    dz[5] = tscr_r - dil(rm_r) - deg(rm_r)
    dz[6] = tlr(rm_r, nr) - dil(rtcr)
    dz[7] = Vrep - Vdam + Vinflux - dil(rh)
    dz[8] = Vdam - Vtag - kdeg*rd - dil(rd)
    dz[9] = Vtag - Vrep - dil(rt)
    
    dz
    
end

rtc_inhib_mod(z, p) = rtc_inhib_mod!(similar(z), z, p, 0)
