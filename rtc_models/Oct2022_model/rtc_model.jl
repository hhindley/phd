
function rtc_model(initial, params, t) 
    L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr = params
    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt = initial


    # growth rate
    lam = gr_c*rh 
    @show lam
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


    # ODEs
    drm_a = tscr_a - dil(rm_a) - deg(rm_a)
    drtca = tlr(rm_a, na) - dil(rtca)    
    drm_b = tscr_b - dil(rm_b) - deg(rm_b)
    drtcb = tlr(rm_b, nb) - dil(rtcb)
    drm_r = tscr_r - dil(rm_r) - deg(rm_r)
    drtcr = tlr(rm_r, nr) - dil(rtcr)
    drh = Vrep - Vdam + Vinflux - dil(rh)
    drd = Vdam - Vtag - kdeg*rd - dil(rd)
    drt = Vtag - Vrep - dil(rt)

    @SVector [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt]
end

function rtc_model_knockout_OD(initial, params, t) 
    L, c, kr, Vmax_init, Km_init, ω_a, ω_b, ω_r, θtscr, g_max, θtlr, km_a, km_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr = params
    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt, OD = initial


    # growth rate
    lam = gr_c*rh 
    
    # dilution by growth and degradation 
    dil(species) = lam*species
    deg(species) = d*species
    
    # MWC
    alpha = rt/kr 
    fa = (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6)
    ra = fa*rtcr
    
    # transcription
    Vinit = ra*Vmax_init*atp/(Km_init+atp)
    tscr_el_a = ω_a*atp/(θtscr+atp)
    tscr_a = Vinit*tscr_el_a
    tscr_el_b = ω_b*atp/(θtscr+atp)
    tscr_b = Vinit*tscr_el_b
    tscr_r = ω_r*atp/(θtscr+atp)

    # translation
    tlr_el = g_max*atp/(θtlr+atp)
    tlr(rm_x, nx) = (1/nx)*rh*rm_x*tlr_el # *1/nx nx = length of RtcA, RtcB and RtcR

    # ribosomes
    rtca1 = (atp*rtca)/(atp+(km_a*rd)) 
    rtcb1 = (atp*rtcb)/(atp+(km_b*rt)) 

    Vrep = krep*rtcb1*rt
    Vdam = kdam*rh
    Vinflux = kin*tlr_el
    Vtag = ktag*rtca1*rd


    # ODEs
    drm_a = tscr_a - dil(rm_a) - deg(rm_a)
    drtca = tlr(rm_a, na) - dil(rtca)    
    drm_b = tscr_b - dil(rm_b) - deg(rm_b)
    drtcb = tlr(rm_b, nb) - dil(rtcb)
    drm_r = tscr_r - dil(rm_r) - deg(rm_r)
    drtcr = tlr(rm_r, nr) - dil(rtcr)
    drh = Vrep - Vdam + Vinflux - dil(rh)
    drd = Vdam - Vtag - kdeg*rd - dil(rd)
    drt = Vtag - Vrep - dil(rt)
    dOD = lam*OD*(1-(OD/k))

    @SVector [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt, dOD]
end

function rtc_model_OD(initial, params, t) 
    L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, k = params
    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt, OD = initial


    # growth rate
    lam = gr_c*rh 
    
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

    Vrep = krep*rtcb1*rt
    Vdam = kdam*rh
    Vinflux = kin*tlr_el
    Vtag = ktag*rtca1*rd


    # ODEs
    drm_a = tscr_a - dil(rm_a) - deg(rm_a)
    drtca = tlr(rm_a, na) - dil(rtca)    
    drm_b = tscr_b - dil(rm_b) - deg(rm_b)
    drtcb = tlr(rm_b, nb) - dil(rtcb)
    drm_r = tscr_r - dil(rm_r) - deg(rm_r)
    drtcr = tlr(rm_r, nr) - dil(rtcr)
    drh = Vrep - Vdam + Vinflux - dil(rh)
    drd = Vdam - Vtag - kdeg*rd - dil(rd)
    drt = Vtag - Vrep - dil(rt)
    dOD = lam*OD*(1-(OD/k))

    @SVector [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt, dOD]
end