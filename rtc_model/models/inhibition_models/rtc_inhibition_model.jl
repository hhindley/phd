using Parameters, CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, Statistics
using Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra

# using PlotlyJS, ProgressBars
include("/home/holliehindley/phd/may23_rtc/functions/solving.jl"); include("/home/holliehindley/phd/may23_rtc/functions/set_ups.jl"); include("/home/holliehindley/phd/may23_rtc/functions/plotting.jl"); 
include("/home/holliehindley/phd/may23_rtc/functions/sweep_params.jl"); include("/home/holliehindley/phd/may23_rtc/models/rtc_orig.jl"); include("/home/holliehindley/phd/may23_rtc/models/atp_lam_kin_t.jl"); 
include("/home/holliehindley/phd/may23_rtc/models/single_t.jl"); include("/home/holliehindley/phd/may23_rtc/models/combinations_t.jl"); 
# include("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/bf_funcs.jl");

# include("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/init_switch/funcs.jl");

function rtc_inhib_model_rtcb(initial, params, t) 
    L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, kc, k_inhib1, k_inhib2, inhib = params
    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt, rtc_i = initial


    # dilution by growth and degradation 
    dil(species) = lam*species
    deg(species) = d*species
    
    # MWC
    alpha = rt/kr # unitless
    fa = (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6) # unitless 
    ra = fa*rtcr # uM 
    
    # transcription
    Voc = Vmax_init*atp/(Km_init+atp) # uM min-1 
    sig_o = ra*Voc/k_diss # uM

    tscr_el_a = ω_ab*atp/(θtscr+atp) # min-1
    tscr_a = sig_o*tscr_el_a # uM min-1
    tscr_el_b = ω_ab*atp/(θtscr+atp) # min-1
    tscr_b = sig_o*tscr_el_b # uM min-1
    tscr_r = ω_r*atp/(θtscr+atp) # uM min-1

    # translation
    tlr_el = g_max*atp/(θtlr+atp) # aa min-1 molec-1
    tlr(rm_x, nx) = (1/nx)*kc*rh*rm_x*tlr_el # uM min-1

    # # ribosomes
    Vrep = rtcb*rt*krep/(rt+km_b) # uM min-1 
    Vdam = kdam*rh # uM min-1
    Vinflux = kin*tlr_el # uM min-1 
    Vtag = rtca*rd*ktag/(rd+km_a) # uM min-1 

    # ODEs
    drm_a = tscr_a - dil(rm_a) - deg(rm_a)
    drtca = tlr(rm_a, na) - dil(rtca)    
    drm_b = tscr_b - dil(rm_b) - deg(rm_b)
    drtcb = tlr(rm_b, nb) - dil(rtcb) - k_inhib1*rtcb*inhib + k_inhib2*rtc_i
    drm_r = tscr_r - dil(rm_r) - deg(rm_r)
    drtcr = tlr(rm_r, nr) - dil(rtcr)
    drh = Vrep - Vdam + Vinflux - dil(rh)
    drd = Vdam - Vtag - kdeg*rd - dil(rd)
    drt = Vtag - Vrep - dil(rt)
    drtc_i = k_inhib1*rtcb*inhib - k_inhib2*rtc_i - dil(rtc_i)
    # @SVector [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt]
    [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt, drtc_i]

end


function rtc_inhib_mod_rtcb!(dz, z, p, t)
    @unpack L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, kc, k_inhib1, k_inhib2, inhib = p
    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt, rtc_i = z


    # dilution by growth and degradation 
    dil(species) = lam*species
    deg(species) = d*species
    
    # MWC
    alpha = rt/kr # unitless
    fa = (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6) # unitless 
    ra = fa*rtcr # uM 
    
    # transcription
    Voc = Vmax_init*atp/(Km_init+atp) # uM min-1 
    sig_o = ra*Voc/k_diss # uM

    tscr_el_a = ω_ab*atp/(θtscr+atp) # min-1
    tscr_a = sig_o*tscr_el_a # uM min-1
    tscr_el_b = ω_ab*atp/(θtscr+atp) # min-1
    tscr_b = sig_o*tscr_el_b # uM min-1
    tscr_r = ω_r*atp/(θtscr+atp) # uM min-1

    # translation
    tlr_el = g_max*atp/(θtlr+atp) # aa min-1 molec-1
    tlr(rm_x, nx) = (1/nx)*kc*rh*rm_x*tlr_el # uM min-1

    # # ribosomes
    Vrep = rtcb*rt*krep/(rt+km_b) # uM min-1 
    Vdam = kdam*rh # uM min-1
    Vinflux = kin*tlr_el # uM min-1 
    Vtag = rtca*rd*ktag/(rd+km_a) # uM min-1 

    # ODEs
    dz[1] = tscr_a - dil(rm_a) - deg(rm_a)
    dz[2] = tlr(rm_a, na) - dil(rtca)    
    dz[3] = tscr_b - dil(rm_b) - deg(rm_b)
    dz[4] = tlr(rm_b, nb) - dil(rtcb) - k_inhib1*rtcb*inhib + k_inhib2*rtc_i
    dz[5] = tscr_r - dil(rm_r) - deg(rm_r)
    dz[6] = tlr(rm_r, nr) - dil(rtcr)
    dz[7] = Vrep - Vdam + Vinflux - dil(rh)
    dz[8] = Vdam - Vtag - kdeg*rd - dil(rd)
    dz[9] = Vtag - Vrep - dil(rt)
    dz[10] = k_inhib1*rtcb*inhib - k_inhib2*rtc_i - dil(rtc_i)

    dz
end


rtc_inhib_mod_rtcb(z, p) = rtc_inhib_mod_rtcb!(similar(z), z, p, 0)








function rtc_inhib_model_rtca(initial, params, t) 
    L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, kc, k_inhib1, k_inhib2, inhib = params
    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt, rtc_i = initial


    # dilution by growth and degradation 
    dil(species) = lam*species
    deg(species) = d*species
    
    # MWC
    alpha = rt/kr # unitless
    fa = (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6) # unitless 
    ra = fa*rtcr # uM 
    
    # transcription
    Voc = Vmax_init*atp/(Km_init+atp) # uM min-1 
    sig_o = ra*Voc/k_diss # uM

    tscr_el_a = ω_ab*atp/(θtscr+atp) # min-1
    tscr_a = sig_o*tscr_el_a # uM min-1
    tscr_el_b = ω_ab*atp/(θtscr+atp) # min-1
    tscr_b = sig_o*tscr_el_b # uM min-1
    tscr_r = ω_r*atp/(θtscr+atp) # uM min-1

    # translation
    tlr_el = g_max*atp/(θtlr+atp) # aa min-1 molec-1
    tlr(rm_x, nx) = (1/nx)*kc*rh*rm_x*tlr_el # uM min-1

    # # ribosomes
    Vrep = rtcb*rt*krep/(rt+km_b) # uM min-1 
    Vdam = kdam*rh # uM min-1
    Vinflux = kin*tlr_el # uM min-1 
    Vtag = rtca*rd*ktag/(rd+km_a) # uM min-1 

    # ODEs
    drm_a = tscr_a - dil(rm_a) - deg(rm_a)
    drtca = tlr(rm_a, na) - dil(rtca) - k_inhib1*rtca*inhib + k_inhib2*rtc_i
    drm_b = tscr_b - dil(rm_b) - deg(rm_b)
    drtcb = tlr(rm_b, nb) - dil(rtcb) 
    drm_r = tscr_r - dil(rm_r) - deg(rm_r)
    drtcr = tlr(rm_r, nr) - dil(rtcr)
    drh = Vrep - Vdam + Vinflux - dil(rh)
    drd = Vdam - Vtag - kdeg*rd - dil(rd)
    drt = Vtag - Vrep - dil(rt)
    drtc_i = k_inhib1*rtca*inhib - k_inhib2*rtc_i - dil(rtc_i)
    # @SVector [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt]
    [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt, drtc_i]

end


function rtc_inhib_mod_rtca!(dz, z, p, t)
    @unpack L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, kc, k_inhib1, k_inhib2, inhib = p
    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt, rtc_i = z


    # dilution by growth and degradation 
    dil(species) = lam*species
    deg(species) = d*species
    
    # MWC
    alpha = rt/kr # unitless
    fa = (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6) # unitless 
    ra = fa*rtcr # uM 
    
    # transcription
    Voc = Vmax_init*atp/(Km_init+atp) # uM min-1 
    sig_o = ra*Voc/k_diss # uM

    tscr_el_a = ω_ab*atp/(θtscr+atp) # min-1
    tscr_a = sig_o*tscr_el_a # uM min-1
    tscr_el_b = ω_ab*atp/(θtscr+atp) # min-1
    tscr_b = sig_o*tscr_el_b # uM min-1
    tscr_r = ω_r*atp/(θtscr+atp) # uM min-1

    # translation
    tlr_el = g_max*atp/(θtlr+atp) # aa min-1 molec-1
    tlr(rm_x, nx) = (1/nx)*kc*rh*rm_x*tlr_el # uM min-1

    # # ribosomes
    Vrep = rtcb*rt*krep/(rt+km_b) # uM min-1 
    Vdam = kdam*rh # uM min-1
    Vinflux = kin*tlr_el # uM min-1 
    Vtag = rtca*rd*ktag/(rd+km_a) # uM min-1 

    # ODEs
    dz[1] = tscr_a - dil(rm_a) - deg(rm_a)
    dz[2] = tlr(rm_a, na) - dil(rtca) - k_inhib1*rtca*inhib + k_inhib2*rtc_i
    dz[3] = tscr_b - dil(rm_b) - deg(rm_b)
    dz[4] = tlr(rm_b, nb) - dil(rtcb) 
    dz[5] = tscr_r - dil(rm_r) - deg(rm_r)
    dz[6] = tlr(rm_r, nr) - dil(rtcr)
    dz[7] = Vrep - Vdam + Vinflux - dil(rh)
    dz[8] = Vdam - Vtag - kdeg*rd - dil(rd)
    dz[9] = Vtag - Vrep - dil(rt)
    dz[10] = k_inhib1*rtca*inhib - k_inhib2*rtc_i - dil(rtc_i)
    dz
end


rtc_inhib_mod_rtca(z, p) = rtc_inhib_mod_rtca!(similar(z), z, p, 0)





function rtc_inhib_model_rtcr(initial, params, t) 
    L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, kc, k_inhib1, k_inhib2, inhib = params
    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt, rtc_i = initial


    # dilution by growth and degradation 
    dil(species) = lam*species
    deg(species) = d*species
    
    # MWC
    alpha = rt/kr # unitless
    fa = (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6) # unitless 
    ra = fa*rtcr # uM 
    
    # transcription
    Voc = Vmax_init*atp/(Km_init+atp) # uM min-1 
    sig_o = ra*Voc/k_diss # uM

    tscr_el_a = ω_ab*atp/(θtscr+atp) # min-1
    tscr_a = sig_o*tscr_el_a # uM min-1
    tscr_el_b = ω_ab*atp/(θtscr+atp) # min-1
    tscr_b = sig_o*tscr_el_b # uM min-1
    tscr_r = ω_r*atp/(θtscr+atp) # uM min-1

    # translation
    tlr_el = g_max*atp/(θtlr+atp) # aa min-1 molec-1
    tlr(rm_x, nx) = (1/nx)*kc*rh*rm_x*tlr_el # uM min-1

    # # ribosomes
    Vrep = rtcb*rt*krep/(rt+km_b) # uM min-1 
    Vdam = kdam*rh # uM min-1
    Vinflux = kin*tlr_el # uM min-1 
    Vtag = rtca*rd*ktag/(rd+km_a) # uM min-1 

    # ODEs
    drm_a = tscr_a - dil(rm_a) - deg(rm_a)
    drtca = tlr(rm_a, na) - dil(rtca)    
    drm_b = tscr_b - dil(rm_b) - deg(rm_b)
    drtcb = tlr(rm_b, nb) - dil(rtcb) 
    drm_r = tscr_r - dil(rm_r) - deg(rm_r)
    drtcr = tlr(rm_r, nr) - dil(rtcr) - k_inhib1*rtcr*inhib + k_inhib2*rtc_i
    drh = Vrep - Vdam + Vinflux - dil(rh)
    drd = Vdam - Vtag - kdeg*rd - dil(rd)
    drt = Vtag - Vrep - dil(rt)
    drtc_i = k_inhib1*rtcr*inhib - k_inhib2*rtc_i - dil(rtc_i)
    # @SVector [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt]
    [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt, drtc_i]

end


function rtc_inhib_mod_rtcr!(dz, z, p, t)
    @unpack L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, kc, k_inhib1, k_inhib2, inhib = p
    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt, rtc_i = z


    # dilution by growth and degradation 
    dil(species) = lam*species
    deg(species) = d*species
    
    # MWC
    alpha = rt/kr # unitless
    fa = (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6) # unitless 
    ra = fa*rtcr # uM 
    
    # transcription
    Voc = Vmax_init*atp/(Km_init+atp) # uM min-1 
    sig_o = ra*Voc/k_diss # uM

    tscr_el_a = ω_ab*atp/(θtscr+atp) # min-1
    tscr_a = sig_o*tscr_el_a # uM min-1
    tscr_el_b = ω_ab*atp/(θtscr+atp) # min-1
    tscr_b = sig_o*tscr_el_b # uM min-1
    tscr_r = ω_r*atp/(θtscr+atp) # uM min-1

    # translation
    tlr_el = g_max*atp/(θtlr+atp) # aa min-1 molec-1
    tlr(rm_x, nx) = (1/nx)*kc*rh*rm_x*tlr_el # uM min-1

    # # ribosomes
    Vrep = rtcb*rt*krep/(rt+km_b) # uM min-1 
    Vdam = kdam*rh # uM min-1
    Vinflux = kin*tlr_el # uM min-1 
    Vtag = rtca*rd*ktag/(rd+km_a) # uM min-1 

    # ODEs
    dz[1] = tscr_a - dil(rm_a) - deg(rm_a)
    dz[2] = tlr(rm_a, na) - dil(rtca)    
    dz[3] = tscr_b - dil(rm_b) - deg(rm_b)
    dz[4] = tlr(rm_b, nb) - dil(rtcb)
    dz[5] = tscr_r - dil(rm_r) - deg(rm_r)
    dz[6] = tlr(rm_r, nr) - dil(rtcr) - k_inhib1*rtcr*inhib + k_inhib2*rtc_i
    dz[7] = Vrep - Vdam + Vinflux - dil(rh)
    dz[8] = Vdam - Vtag - kdeg*rd - dil(rd)
    dz[9] = Vtag - Vrep - dil(rt)
    dz[10] = k_inhib1*rtcr*inhib - k_inhib2*rtc_i - dil(rtc_i)

    dz
end


rtc_inhib_mod_rtcr(z, p) = rtc_inhib_mod_rtcr!(similar(z), z, p, 0)











function rtc_inhib_model_rt(initial, params, t) 
    L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, kc, k_inhib1, k_inhib2, inhib = params
    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt, rtc_i = initial


    # dilution by growth and degradation 
    dil(species) = lam*species
    deg(species) = d*species
    
    # MWC
    alpha = rt/kr # unitless
    fa = (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6) # unitless 
    ra = fa*rtcr # uM 
    
    # transcription
    Voc = Vmax_init*atp/(Km_init+atp) # uM min-1 
    sig_o = ra*Voc/k_diss # uM

    tscr_el_a = ω_ab*atp/(θtscr+atp) # min-1
    tscr_a = sig_o*tscr_el_a # uM min-1
    tscr_el_b = ω_ab*atp/(θtscr+atp) # min-1
    tscr_b = sig_o*tscr_el_b # uM min-1
    tscr_r = ω_r*atp/(θtscr+atp) # uM min-1

    # translation
    tlr_el = g_max*atp/(θtlr+atp) # aa min-1 molec-1
    tlr(rm_x, nx) = (1/nx)*kc*rh*rm_x*tlr_el # uM min-1

    # # ribosomes
    Vrep = rtcb*rt*krep/(rt+km_b) # uM min-1 
    Vdam = kdam*rh # uM min-1
    Vinflux = kin*tlr_el # uM min-1 
    Vtag = rtca*rd*ktag/(rd+km_a) # uM min-1 

    # ODEs
    drm_a = tscr_a - dil(rm_a) - deg(rm_a)
    drtca = tlr(rm_a, na) - dil(rtca)    
    drm_b = tscr_b - dil(rm_b) - deg(rm_b)
    drtcb = tlr(rm_b, nb) - dil(rtcb) 
    drm_r = tscr_r - dil(rm_r) - deg(rm_r)
    drtcr = tlr(rm_r, nr) - dil(rtcr) 
    drh = Vrep - Vdam + Vinflux - dil(rh)
    drd = Vdam - Vtag - kdeg*rd - dil(rd)
    drt = Vtag - Vrep - dil(rt) - k_inhib1*rt*inhib + k_inhib2*rtc_i
    drtc_i = k_inhib1*rt*inhib - k_inhib2*rtc_i - dil(rtc_i)
    # @SVector [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt]
    [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt, drtc_i]

end


function rtc_inhib_mod_rt!(dz, z, p, t)
    @unpack L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, kc, k_inhib1, k_inhib2, inhib = p
    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt, rtc_i = z


    # dilution by growth and degradation 
    dil(species) = lam*species
    deg(species) = d*species
    
    # MWC
    alpha = rt/kr # unitless
    fa = (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6) # unitless 
    ra = fa*rtcr # uM 
    
    # transcription
    Voc = Vmax_init*atp/(Km_init+atp) # uM min-1 
    sig_o = ra*Voc/k_diss # uM

    tscr_el_a = ω_ab*atp/(θtscr+atp) # min-1
    tscr_a = sig_o*tscr_el_a # uM min-1
    tscr_el_b = ω_ab*atp/(θtscr+atp) # min-1
    tscr_b = sig_o*tscr_el_b # uM min-1
    tscr_r = ω_r*atp/(θtscr+atp) # uM min-1

    # translation
    tlr_el = g_max*atp/(θtlr+atp) # aa min-1 molec-1
    tlr(rm_x, nx) = (1/nx)*kc*rh*rm_x*tlr_el # uM min-1

    # # ribosomes
    Vrep = rtcb*rt*krep/(rt+km_b) # uM min-1 
    Vdam = kdam*rh # uM min-1
    Vinflux = kin*tlr_el # uM min-1 
    Vtag = rtca*rd*ktag/(rd+km_a) # uM min-1 

    # ODEs
    dz[1] = tscr_a - dil(rm_a) - deg(rm_a)
    dz[2] = tlr(rm_a, na) - dil(rtca)    
    dz[3] = tscr_b - dil(rm_b) - deg(rm_b)
    dz[4] = tlr(rm_b, nb) - dil(rtcb)
    dz[5] = tscr_r - dil(rm_r) - deg(rm_r)
    dz[6] = tlr(rm_r, nr) - dil(rtcr)
    dz[7] = Vrep - Vdam + Vinflux - dil(rh)
    dz[8] = Vdam - Vtag - kdeg*rd - dil(rd)
    dz[9] = Vtag - Vrep - dil(rt) - k_inhib1*rt*inhib + k_inhib2*rtc_i
    dz[10] =k_inhib1*rt*inhib - k_inhib2*rtc_i - dil(rtc_i)

    dz
end


rtc_inhib_mod_rt(z, p) = rtc_inhib_mod_rt!(similar(z), z, p, 0)










function rtc_inhib_model_rtcab(initial, params, t) 
    L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, kc, k_inhib1, k_inhib2, inhib = params
    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt, rtca_i, rtcb_i = initial


    # dilution by growth and degradation 
    dil(species) = lam*species
    deg(species) = d*species
    
    # MWC
    alpha = rt/kr # unitless
    fa = (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6) # unitless 
    ra = fa*rtcr # uM 
    
    # transcription
    Voc = Vmax_init*atp/(Km_init+atp) # uM min-1 
    sig_o = ra*Voc/k_diss # uM

    tscr_el_a = ω_ab*atp/(θtscr+atp) # min-1
    tscr_a = sig_o*tscr_el_a # uM min-1
    tscr_el_b = ω_ab*atp/(θtscr+atp) # min-1
    tscr_b = sig_o*tscr_el_b # uM min-1
    tscr_r = ω_r*atp/(θtscr+atp) # uM min-1

    # translation
    tlr_el = g_max*atp/(θtlr+atp) # aa min-1 molec-1
    tlr(rm_x, nx) = (1/nx)*kc*rh*rm_x*tlr_el # uM min-1

    # # ribosomes
    Vrep = rtcb*rt*krep/(rt+km_b) # uM min-1 
    Vdam = kdam*rh # uM min-1
    Vinflux = kin*tlr_el # uM min-1 
    Vtag = rtca*rd*ktag/(rd+km_a) # uM min-1 

    # ODEs
    drm_a = tscr_a - dil(rm_a) - deg(rm_a)
    drtca = tlr(rm_a, na) - dil(rtca) - k_inhib1*rtca*inhib + k_inhib2*rtca_i
    drm_b = tscr_b - dil(rm_b) - deg(rm_b)
    drtcb = tlr(rm_b, nb) - dil(rtcb) - k_inhib1*rtcb*inhib + k_inhib2*rtcb_i
    drm_r = tscr_r - dil(rm_r) - deg(rm_r)
    drtcr = tlr(rm_r, nr) - dil(rtcr)
    drh = Vrep - Vdam + Vinflux - dil(rh)
    drd = Vdam - Vtag - kdeg*rd - dil(rd)
    drt = Vtag - Vrep - dil(rt)
    drtca_i = k_inhib1*rtca*inhib - k_inhib2*rtca_i - dil(rtca_i)
    drtcb_i = k_inhib1*rtcb*inhib - k_inhib2*rtcb_i - dil(rtcb_i)
    # @SVector [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt]
    [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt, drtca_i, drtcb_i]

end


function rtc_inhib_mod_rtcab!(dz, z, p, t)
    @unpack L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, kc, k_inhib1, k_inhib2, inhib = p
    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt, rtca_i, rtcb_i = z


    # dilution by growth and degradation 
    dil(species) = lam*species
    deg(species) = d*species
    
    # MWC
    alpha = rt/kr # unitless
    fa = (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6) # unitless 
    ra = fa*rtcr # uM 
    
    # transcription
    Voc = Vmax_init*atp/(Km_init+atp) # uM min-1 
    sig_o = ra*Voc/k_diss # uM

    tscr_el_a = ω_ab*atp/(θtscr+atp) # min-1
    tscr_a = sig_o*tscr_el_a # uM min-1
    tscr_el_b = ω_ab*atp/(θtscr+atp) # min-1
    tscr_b = sig_o*tscr_el_b # uM min-1
    tscr_r = ω_r*atp/(θtscr+atp) # uM min-1

    # translation
    tlr_el = g_max*atp/(θtlr+atp) # aa min-1 molec-1
    tlr(rm_x, nx) = (1/nx)*kc*rh*rm_x*tlr_el # uM min-1

    # # ribosomes
    Vrep = rtcb*rt*krep/(rt+km_b) # uM min-1 
    Vdam = kdam*rh # uM min-1
    Vinflux = kin*tlr_el # uM min-1 
    Vtag = rtca*rd*ktag/(rd+km_a) # uM min-1 

    # ODEs
    dz[1] = tscr_a - dil(rm_a) - deg(rm_a)
    dz[2] = tlr(rm_a, na) - dil(rtca) - k_inhib1*rtca*inhib + k_inhib2*rtca_i
    dz[3] = tscr_b - dil(rm_b) - deg(rm_b)
    dz[4] = tlr(rm_b, nb) - dil(rtcb) - k_inhib1*rtcb*inhib + k_inhib2*rtcb_i
    dz[5] = tscr_r - dil(rm_r) - deg(rm_r)
    dz[6] = tlr(rm_r, nr) - dil(rtcr)
    dz[7] = Vrep - Vdam + Vinflux - dil(rh)
    dz[8] = Vdam - Vtag - kdeg*rd - dil(rd)
    dz[9] = Vtag - Vrep - dil(rt)
    dz[10] = k_inhib1*rtca*inhib - k_inhib2*rtca_i - dil(rtca_i)
    dz[11] = k_inhib1*rtcb*inhib - k_inhib2*rtcb_i - dil(rtcb_i)

    dz
end


rtc_inhib_mod_rtcab(z, p) = rtc_inhib_mod_rtcab!(similar(z), z, p, 0)












function rtc_inhib_model_rtcbr(initial, params, t) 
    L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, kc, k_inhib1, k_inhib2, inhib = params
    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt, rtcr_i, rtcb_i = initial


    # dilution by growth and degradation 
    dil(species) = lam*species
    deg(species) = d*species
    
    # MWC
    alpha = rt/kr # unitless
    fa = (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6) # unitless 
    ra = fa*rtcr # uM 
    
    # transcription
    Voc = Vmax_init*atp/(Km_init+atp) # uM min-1 
    sig_o = ra*Voc/k_diss # uM

    tscr_el_a = ω_ab*atp/(θtscr+atp) # min-1
    tscr_a = sig_o*tscr_el_a # uM min-1
    tscr_el_b = ω_ab*atp/(θtscr+atp) # min-1
    tscr_b = sig_o*tscr_el_b # uM min-1
    tscr_r = ω_r*atp/(θtscr+atp) # uM min-1

    # translation
    tlr_el = g_max*atp/(θtlr+atp) # aa min-1 molec-1
    tlr(rm_x, nx) = (1/nx)*kc*rh*rm_x*tlr_el # uM min-1

    # # ribosomes
    Vrep = rtcb*rt*krep/(rt+km_b) # uM min-1 
    Vdam = kdam*rh # uM min-1
    Vinflux = kin*tlr_el # uM min-1 
    Vtag = rtca*rd*ktag/(rd+km_a) # uM min-1 

    # ODEs
    drm_a = tscr_a - dil(rm_a) - deg(rm_a)
    drtca = tlr(rm_a, na) - dil(rtca) 
    drm_b = tscr_b - dil(rm_b) - deg(rm_b)
    drtcb = tlr(rm_b, nb) - dil(rtcb) - k_inhib1*rtcb*inhib + k_inhib2*rtcb_i
    drm_r = tscr_r - dil(rm_r) - deg(rm_r)
    drtcr = tlr(rm_r, nr) - dil(rtcr)- k_inhib1*rtcr*inhib + k_inhib2*rtcr_i
    drh = Vrep - Vdam + Vinflux - dil(rh)
    drd = Vdam - Vtag - kdeg*rd - dil(rd)
    drt = Vtag - Vrep - dil(rt)
    drtcr_i = k_inhib1*rtcr*inhib - k_inhib2*rtcr_i - dil(rtcr_i)
    drtcb_i = k_inhib1*rtcb*inhib - k_inhib2*rtcb_i - dil(rtcb_i)
    # @SVector [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt]
    [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt, drtcr_i, drtcb_i]

end


function rtc_inhib_mod_rtcbr!(dz, z, p, t)
    @unpack L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, kc, k_inhib1, k_inhib2, inhib = p
    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt, rtcr_i, rtcb_i = z


    # dilution by growth and degradation 
    dil(species) = lam*species
    deg(species) = d*species
    
    # MWC
    alpha = rt/kr # unitless
    fa = (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6) # unitless 
    ra = fa*rtcr # uM 
    
    # transcription
    Voc = Vmax_init*atp/(Km_init+atp) # uM min-1 
    sig_o = ra*Voc/k_diss # uM

    tscr_el_a = ω_ab*atp/(θtscr+atp) # min-1
    tscr_a = sig_o*tscr_el_a # uM min-1
    tscr_el_b = ω_ab*atp/(θtscr+atp) # min-1
    tscr_b = sig_o*tscr_el_b # uM min-1
    tscr_r = ω_r*atp/(θtscr+atp) # uM min-1

    # translation
    tlr_el = g_max*atp/(θtlr+atp) # aa min-1 molec-1
    tlr(rm_x, nx) = (1/nx)*kc*rh*rm_x*tlr_el # uM min-1

    # # ribosomes
    Vrep = rtcb*rt*krep/(rt+km_b) # uM min-1 
    Vdam = kdam*rh # uM min-1
    Vinflux = kin*tlr_el # uM min-1 
    Vtag = rtca*rd*ktag/(rd+km_a) # uM min-1 

    # ODEs
    dz[1] = tscr_a - dil(rm_a) - deg(rm_a)
    dz[2] = tlr(rm_a, na) - dil(rtca)
    dz[3] = tscr_b - dil(rm_b) - deg(rm_b)
    dz[4] = tlr(rm_b, nb) - dil(rtcb) - k_inhib1*rtcb*inhib + k_inhib2*rtcb_i
    dz[5] = tscr_r - dil(rm_r) - deg(rm_r)
    dz[6] = tlr(rm_r, nr) - dil(rtcr) - k_inhib1*rtcr*inhib + k_inhib2*rtcr_i
    dz[7] = Vrep - Vdam + Vinflux - dil(rh)
    dz[8] = Vdam - Vtag - kdeg*rd - dil(rd)
    dz[9] = Vtag - Vrep - dil(rt)
    dz[10] = k_inhib1*rtcr*inhib - k_inhib2*rtcr_i - dil(rtcr_i)
    dz[11] = k_inhib1*rtcb*inhib - k_inhib2*rtcb_i - dil(rtcb_i)

    dz
end


rtc_inhib_mod_rtcbr(z, p) = rtc_inhib_mod_rtcbr!(similar(z), z, p, 0)














function rtc_inhib_model_rtcar(initial, params, t) 
    L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, kc, k_inhib1, k_inhib2, inhib = params
    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt, rtcr_i, rtca_i = initial


    # dilution by growth and degradation 
    dil(species) = lam*species
    deg(species) = d*species
    
    # MWC
    alpha = rt/kr # unitless
    fa = (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6) # unitless 
    ra = fa*rtcr # uM 
    
    # transcription
    Voc = Vmax_init*atp/(Km_init+atp) # uM min-1 
    sig_o = ra*Voc/k_diss # uM

    tscr_el_a = ω_ab*atp/(θtscr+atp) # min-1
    tscr_a = sig_o*tscr_el_a # uM min-1
    tscr_el_b = ω_ab*atp/(θtscr+atp) # min-1
    tscr_b = sig_o*tscr_el_b # uM min-1
    tscr_r = ω_r*atp/(θtscr+atp) # uM min-1

    # translation
    tlr_el = g_max*atp/(θtlr+atp) # aa min-1 molec-1
    tlr(rm_x, nx) = (1/nx)*kc*rh*rm_x*tlr_el # uM min-1

    # # ribosomes
    Vrep = rtcb*rt*krep/(rt+km_b) # uM min-1 
    Vdam = kdam*rh # uM min-1
    Vinflux = kin*tlr_el # uM min-1 
    Vtag = rtca*rd*ktag/(rd+km_a) # uM min-1 

    # ODEs
    drm_a = tscr_a - dil(rm_a) - deg(rm_a)
    drtca = tlr(rm_a, na) - dil(rtca) - k_inhib1*rtca*inhib + k_inhib2*rtca_i
    drm_b = tscr_b - dil(rm_b) - deg(rm_b)
    drtcb = tlr(rm_b, nb) - dil(rtcb) 
    drm_r = tscr_r - dil(rm_r) - deg(rm_r)
    drtcr = tlr(rm_r, nr) - dil(rtcr)- k_inhib1*rtcr*inhib + k_inhib2*rtcr_i
    drh = Vrep - Vdam + Vinflux - dil(rh)
    drd = Vdam - Vtag - kdeg*rd - dil(rd)
    drt = Vtag - Vrep - dil(rt)
    drtcr_i = k_inhib1*rtcr*inhib - k_inhib2*rtcr_i - dil(rtcr_i)
    drtca_i = k_inhib1*rtca*inhib - k_inhib2*rtca_i - dil(rtca_i)
    # @SVector [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt]
    [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt, drtcr_i, drtca_i]

end


function rtc_inhib_mod_rtcar!(dz, z, p, t)
    @unpack L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, kc, k_inhib1, k_inhib2, inhib = p
    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt, rtcr_i, rtca_i = z


    # dilution by growth and degradation 
    dil(species) = lam*species
    deg(species) = d*species
    
    # MWC
    alpha = rt/kr # unitless
    fa = (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6) # unitless 
    ra = fa*rtcr # uM 
    
    # transcription
    Voc = Vmax_init*atp/(Km_init+atp) # uM min-1 
    sig_o = ra*Voc/k_diss # uM

    tscr_el_a = ω_ab*atp/(θtscr+atp) # min-1
    tscr_a = sig_o*tscr_el_a # uM min-1
    tscr_el_b = ω_ab*atp/(θtscr+atp) # min-1
    tscr_b = sig_o*tscr_el_b # uM min-1
    tscr_r = ω_r*atp/(θtscr+atp) # uM min-1

    # translation
    tlr_el = g_max*atp/(θtlr+atp) # aa min-1 molec-1
    tlr(rm_x, nx) = (1/nx)*kc*rh*rm_x*tlr_el # uM min-1

    # # ribosomes
    Vrep = rtcb*rt*krep/(rt+km_b) # uM min-1 
    Vdam = kdam*rh # uM min-1
    Vinflux = kin*tlr_el # uM min-1 
    Vtag = rtca*rd*ktag/(rd+km_a) # uM min-1 

    # ODEs
    dz[1] = tscr_a - dil(rm_a) - deg(rm_a)
    dz[2] = tlr(rm_a, na) - dil(rtca)- k_inhib1*rtca*inhib + k_inhib2*rtca_i
    dz[3] = tscr_b - dil(rm_b) - deg(rm_b)
    dz[4] = tlr(rm_b, nb) - dil(rtcb) 
    dz[5] = tscr_r - dil(rm_r) - deg(rm_r)
    dz[6] = tlr(rm_r, nr) - dil(rtcr) - k_inhib1*rtcr*inhib + k_inhib2*rtcr_i
    dz[7] = Vrep - Vdam + Vinflux - dil(rh)
    dz[8] = Vdam - Vtag - kdeg*rd - dil(rd)
    dz[9] = Vtag - Vrep - dil(rt)
    dz[10] = k_inhib1*rtcr*inhib - k_inhib2*rtcr_i - dil(rtcr_i)
    dz[11] = k_inhib1*rtca*inhib - k_inhib2*rtca_i - dil(rtca_i)

    dz
end


rtc_inhib_mod_rtcar(z, p) = rtc_inhib_mod_rtcar!(similar(z), z, p, 0)










function rtc_inhib_model_rtcat(initial, params, t) 
    L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, kc, k_inhib1, k_inhib2, inhib = params
    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt, rt_i, rtca_i = initial


    # dilution by growth and degradation 
    dil(species) = lam*species
    deg(species) = d*species
    
    # MWC
    alpha = rt/kr # unitless
    fa = (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6) # unitless 
    ra = fa*rtcr # uM 
    
    # transcription
    Voc = Vmax_init*atp/(Km_init+atp) # uM min-1 
    sig_o = ra*Voc/k_diss # uM

    tscr_el_a = ω_ab*atp/(θtscr+atp) # min-1
    tscr_a = sig_o*tscr_el_a # uM min-1
    tscr_el_b = ω_ab*atp/(θtscr+atp) # min-1
    tscr_b = sig_o*tscr_el_b # uM min-1
    tscr_r = ω_r*atp/(θtscr+atp) # uM min-1

    # translation
    tlr_el = g_max*atp/(θtlr+atp) # aa min-1 molec-1
    tlr(rm_x, nx) = (1/nx)*kc*rh*rm_x*tlr_el # uM min-1

    # # ribosomes
    Vrep = rtcb*rt*krep/(rt+km_b) # uM min-1 
    Vdam = kdam*rh # uM min-1
    Vinflux = kin*tlr_el # uM min-1 
    Vtag = rtca*rd*ktag/(rd+km_a) # uM min-1 

    # ODEs
    drm_a = tscr_a - dil(rm_a) - deg(rm_a)
    drtca = tlr(rm_a, na) - dil(rtca) - k_inhib1*rtca*inhib + k_inhib2*rtca_i
    drm_b = tscr_b - dil(rm_b) - deg(rm_b)
    drtcb = tlr(rm_b, nb) - dil(rtcb) 
    drm_r = tscr_r - dil(rm_r) - deg(rm_r)
    drtcr = tlr(rm_r, nr) - dil(rtcr)
    drh = Vrep - Vdam + Vinflux - dil(rh)
    drd = Vdam - Vtag - kdeg*rd - dil(rd)
    drt = Vtag - Vrep - dil(rt) - k_inhib1*rt*inhib + k_inhib2*rt_i
    drt_i = k_inhib1*rt*inhib - k_inhib2*rt_i - dil(rt_i)
    drtca_i = k_inhib1*rtca*inhib - k_inhib2*rtca_i - dil(rtca_i)
    # @SVector [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt]
    [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt, drt_i, drtca_i]

end


function rtc_inhib_mod_rtcat!(dz, z, p, t)
    @unpack L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, kc, k_inhib1, k_inhib2, inhib = p
    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt, rt_i, rtca_i = z


    # dilution by growth and degradation 
    dil(species) = lam*species
    deg(species) = d*species
    
    # MWC
    alpha = rt/kr # unitless
    fa = (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6) # unitless 
    ra = fa*rtcr # uM 
    
    # transcription
    Voc = Vmax_init*atp/(Km_init+atp) # uM min-1 
    sig_o = ra*Voc/k_diss # uM

    tscr_el_a = ω_ab*atp/(θtscr+atp) # min-1
    tscr_a = sig_o*tscr_el_a # uM min-1
    tscr_el_b = ω_ab*atp/(θtscr+atp) # min-1
    tscr_b = sig_o*tscr_el_b # uM min-1
    tscr_r = ω_r*atp/(θtscr+atp) # uM min-1

    # translation
    tlr_el = g_max*atp/(θtlr+atp) # aa min-1 molec-1
    tlr(rm_x, nx) = (1/nx)*kc*rh*rm_x*tlr_el # uM min-1

    # # ribosomes
    Vrep = rtcb*rt*krep/(rt+km_b) # uM min-1 
    Vdam = kdam*rh # uM min-1
    Vinflux = kin*tlr_el # uM min-1 
    Vtag = rtca*rd*ktag/(rd+km_a) # uM min-1 

    # ODEs
    dz[1] = tscr_a - dil(rm_a) - deg(rm_a)
    dz[2] = tlr(rm_a, na) - dil(rtca)- k_inhib1*rtca*inhib + k_inhib2*rtca_i
    dz[3] = tscr_b - dil(rm_b) - deg(rm_b)
    dz[4] = tlr(rm_b, nb) - dil(rtcb) 
    dz[5] = tscr_r - dil(rm_r) - deg(rm_r)
    dz[6] = tlr(rm_r, nr) - dil(rtcr) 
    dz[7] = Vrep - Vdam + Vinflux - dil(rh)
    dz[8] = Vdam - Vtag - kdeg*rd - dil(rd)
    dz[9] = Vtag - Vrep - dil(rt)- k_inhib1*rt*inhib + k_inhib2*rt_i
    dz[10] = k_inhib1*rt*inhib - k_inhib2*rt_i - dil(rt_i)
    dz[11] = k_inhib1*rtca*inhib - k_inhib2*rtca_i - dil(rtca_i)

    dz
end


rtc_inhib_mod_rtcat(z, p) = rtc_inhib_mod_rtcat!(similar(z), z, p, 0)






function rtc_inhib_model_rtcbt(initial, params, t) 
    L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, kc, k_inhib1, k_inhib2, inhib = params
    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt, rt_i, rtcb_i = initial


    # dilution by growth and degradation 
    dil(species) = lam*species
    deg(species) = d*species
    
    # MWC
    alpha = rt/kr # unitless
    fa = (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6) # unitless 
    ra = fa*rtcr # uM 
    
    # transcription
    Voc = Vmax_init*atp/(Km_init+atp) # uM min-1 
    sig_o = ra*Voc/k_diss # uM

    tscr_el_a = ω_ab*atp/(θtscr+atp) # min-1
    tscr_a = sig_o*tscr_el_a # uM min-1
    tscr_el_b = ω_ab*atp/(θtscr+atp) # min-1
    tscr_b = sig_o*tscr_el_b # uM min-1
    tscr_r = ω_r*atp/(θtscr+atp) # uM min-1

    # translation
    tlr_el = g_max*atp/(θtlr+atp) # aa min-1 molec-1
    tlr(rm_x, nx) = (1/nx)*kc*rh*rm_x*tlr_el # uM min-1

    # # ribosomes
    Vrep = rtcb*rt*krep/(rt+km_b) # uM min-1 
    Vdam = kdam*rh # uM min-1
    Vinflux = kin*tlr_el # uM min-1 
    Vtag = rtca*rd*ktag/(rd+km_a) # uM min-1 

    # ODEs
    drm_a = tscr_a - dil(rm_a) - deg(rm_a)
    drtca = tlr(rm_a, na) - dil(rtca) 
    drm_b = tscr_b - dil(rm_b) - deg(rm_b)
    drtcb = tlr(rm_b, nb) - dil(rtcb) - k_inhib1*rtcb*inhib + k_inhib2*rtcb_i
    drm_r = tscr_r - dil(rm_r) - deg(rm_r)
    drtcr = tlr(rm_r, nr) - dil(rtcr)
    drh = Vrep - Vdam + Vinflux - dil(rh)
    drd = Vdam - Vtag - kdeg*rd - dil(rd)
    drt = Vtag - Vrep - dil(rt) - k_inhib1*rt*inhib + k_inhib2*rt_i
    drt_i = k_inhib1*rt*inhib - k_inhib2*rt_i - dil(rt_i)
    drtcb_i = k_inhib1*rtcb*inhib - k_inhib2*rtcb_i - dil(rtcb_i)
    # @SVector [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt]
    [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt, drt_i, drtcb_i]

end


function rtc_inhib_mod_rtcbt!(dz, z, p, t)
    @unpack L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, kc, k_inhib1, k_inhib2, inhib = p
    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt, rt_i, rtcb_i = z


    # dilution by growth and degradation 
    dil(species) = lam*species
    deg(species) = d*species
    
    # MWC
    alpha = rt/kr # unitless
    fa = (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6) # unitless 
    ra = fa*rtcr # uM 
    
    # transcription
    Voc = Vmax_init*atp/(Km_init+atp) # uM min-1 
    sig_o = ra*Voc/k_diss # uM

    tscr_el_a = ω_ab*atp/(θtscr+atp) # min-1
    tscr_a = sig_o*tscr_el_a # uM min-1
    tscr_el_b = ω_ab*atp/(θtscr+atp) # min-1
    tscr_b = sig_o*tscr_el_b # uM min-1
    tscr_r = ω_r*atp/(θtscr+atp) # uM min-1

    # translation
    tlr_el = g_max*atp/(θtlr+atp) # aa min-1 molec-1
    tlr(rm_x, nx) = (1/nx)*kc*rh*rm_x*tlr_el # uM min-1

    # # ribosomes
    Vrep = rtcb*rt*krep/(rt+km_b) # uM min-1 
    Vdam = kdam*rh # uM min-1
    Vinflux = kin*tlr_el # uM min-1 
    Vtag = rtca*rd*ktag/(rd+km_a) # uM min-1 

    # ODEs
    dz[1] = tscr_a - dil(rm_a) - deg(rm_a)
    dz[2] = tlr(rm_a, na) - dil(rtca)
    dz[3] = tscr_b - dil(rm_b) - deg(rm_b)
    dz[4] = tlr(rm_b, nb) - dil(rtcb) - k_inhib1*rtcb*inhib + k_inhib2*rtcb_i
    dz[5] = tscr_r - dil(rm_r) - deg(rm_r)
    dz[6] = tlr(rm_r, nr) - dil(rtcr) 
    dz[7] = Vrep - Vdam + Vinflux - dil(rh)
    dz[8] = Vdam - Vtag - kdeg*rd - dil(rd)
    dz[9] = Vtag - Vrep - dil(rt)- k_inhib1*rt*inhib + k_inhib2*rt_i
    dz[10] = k_inhib1*rt*inhib - k_inhib2*rt_i - dil(rt_i)
    dz[11] = k_inhib1*rtcb*inhib - k_inhib2*rtcb_i - dil(rtcb_i)

    dz
end


rtc_inhib_mod_rtcbt(z, p) = rtc_inhib_mod_rtcbt!(similar(z), z, p, 0)










function rtc_inhib_model_rtcrt(initial, params, t) 
    L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, kc, k_inhib1, k_inhib2, inhib = params
    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt, rt_i, rtc_i = initial


    # dilution by growth and degradation 
    dil(species) = lam*species
    deg(species) = d*species
    
    # MWC
    alpha = rt/kr # unitless
    fa = (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6) # unitless 
    ra = fa*rtcr # uM 
    
    # transcription
    Voc = Vmax_init*atp/(Km_init+atp) # uM min-1 
    sig_o = ra*Voc/k_diss # uM

    tscr_el_a = ω_ab*atp/(θtscr+atp) # min-1
    tscr_a = sig_o*tscr_el_a # uM min-1
    tscr_el_b = ω_ab*atp/(θtscr+atp) # min-1
    tscr_b = sig_o*tscr_el_b # uM min-1
    tscr_r = ω_r*atp/(θtscr+atp) # uM min-1

    # translation
    tlr_el = g_max*atp/(θtlr+atp) # aa min-1 molec-1
    tlr(rm_x, nx) = (1/nx)*kc*rh*rm_x*tlr_el # uM min-1

    # # ribosomes
    Vrep = rtcb*rt*krep/(rt+km_b) # uM min-1 
    Vdam = kdam*rh # uM min-1
    Vinflux = kin*tlr_el # uM min-1 
    Vtag = rtca*rd*ktag/(rd+km_a) # uM min-1 

    # ODEs
    drm_a = tscr_a - dil(rm_a) - deg(rm_a)
    drtca = tlr(rm_a, na) - dil(rtca) 
    drm_b = tscr_b - dil(rm_b) - deg(rm_b)
    drtcb = tlr(rm_b, nb) - dil(rtcb) 
    drm_r = tscr_r - dil(rm_r) - deg(rm_r)
    drtcr = tlr(rm_r, nr) - dil(rtcr) - k_inhib1*rtcr*inhib + k_inhib2*rtc_i
    drh = Vrep - Vdam + Vinflux - dil(rh)
    drd = Vdam - Vtag - kdeg*rd - dil(rd)
    drt = Vtag - Vrep - dil(rt) - k_inhib1*rt*inhib + k_inhib2*rt_i
    drt_i = k_inhib1*rt*inhib - k_inhib2*rt_i - dil(rt_i)
    drtcr_i = k_inhib1*rtcr*inhib - k_inhib2*rtc_i - dil(rtc_i)
    # @SVector [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt]
    [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt, drt_i, drtc_i]

end


function rtc_inhib_mod_rtcrt!(dz, z, p, t)
    @unpack L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, kc, k_inhib1, k_inhib2, inhib = p
    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt, rt_i, rtc_i = z


    # dilution by growth and degradation 
    dil(species) = lam*species
    deg(species) = d*species
    
    # MWC
    alpha = rt/kr # unitless
    fa = (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6) # unitless 
    ra = fa*rtcr # uM 
    
    # transcription
    Voc = Vmax_init*atp/(Km_init+atp) # uM min-1 
    sig_o = ra*Voc/k_diss # uM

    tscr_el_a = ω_ab*atp/(θtscr+atp) # min-1
    tscr_a = sig_o*tscr_el_a # uM min-1
    tscr_el_b = ω_ab*atp/(θtscr+atp) # min-1
    tscr_b = sig_o*tscr_el_b # uM min-1
    tscr_r = ω_r*atp/(θtscr+atp) # uM min-1

    # translation
    tlr_el = g_max*atp/(θtlr+atp) # aa min-1 molec-1
    tlr(rm_x, nx) = (1/nx)*kc*rh*rm_x*tlr_el # uM min-1

    # # ribosomes
    Vrep = rtcb*rt*krep/(rt+km_b) # uM min-1 
    Vdam = kdam*rh # uM min-1
    Vinflux = kin*tlr_el # uM min-1 
    Vtag = rtca*rd*ktag/(rd+km_a) # uM min-1 

    # ODEs
    dz[1] = tscr_a - dil(rm_a) - deg(rm_a)
    dz[2] = tlr(rm_a, na) - dil(rtca)
    dz[3] = tscr_b - dil(rm_b) - deg(rm_b)
    dz[4] = tlr(rm_b, nb) - dil(rtcb) 
    dz[5] = tscr_r - dil(rm_r) - deg(rm_r)
    dz[6] = tlr(rm_r, nr) - dil(rtcr) - k_inhib1*rtcr*inhib + k_inhib2*rtc_i
    dz[7] = Vrep - Vdam + Vinflux - dil(rh)
    dz[8] = Vdam - Vtag - kdeg*rd - dil(rd)
    dz[9] = Vtag - Vrep - dil(rt)- k_inhib1*rt*inhib + k_inhib2*rt_i
    dz[10] = k_inhib1*rt*inhib - k_inhib2*rt_i - dil(rt_i)
    dz[11] = k_inhib1*rtcr*inhib - k_inhib2*rtc_i - dil(rtc_i)

    dz
end


rtc_inhib_mod_rtcrt(z, p) = rtc_inhib_mod_rtcrt!(similar(z), z, p, 0)











