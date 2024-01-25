function rtc_model(initial, params, t) 
    L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, kc = params
    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt = initial


    # dilution by growth and degradation 
    dil(species) = lam*species # uM min-1
    deg(species) = d*species # uM min-1 
    
    # MWC
    alpha = rt/kr # unitless
    fa = (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6) # unitless 
    ra = fa*rtcr # uM 
    
    # transcription
    Vinit = ra*Vmax_init*atp/(Km_init+atp) # uM min-1 
    tscr_el_a = ω_ab*atp/(θtscr+atp) # unitless
    tscr_a = Vinit*tscr_el_a # uM min-1
    tscr_el_b = ω_ab*atp/(θtscr+atp) # unitless
    tscr_b = Vinit*tscr_el_b # uM min-1
    tscr_r = ω_r*atp/(θtscr+atp) # uM min-1

    # translation
    tlr_el = g_max*atp/(θtlr+atp) # aa min-1 molec-1
    tlr(rm_x, nx) = (1/nx)*kc*rh*rm_x*tlr_el # uM min-1

    # ribosomes
    rtca1 = (atp*rtca)/(atp+(km_a*rd)) # uM ???
    rtcb1 = (atp*rtcb)/(atp+(km_b*rt)) # uM ???

    # rtcb1 = k_b*atp*rtcb/(k_b*atp+krep*rt)

    Vrep = krep*rtcb1*rt # uM min-1 
    Vdam = kdam*rh # uM min-1
    Vinflux = kin*tlr_el # uM min-1 
    Vtag = ktag*rtca1*rd # uM min-1 


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

    # @SVector [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt]
    [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt]

end




function rtc_mod!(dz, z, p, t)
    @unpack L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, kc = p
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
    tlr(rm_x, nx) = (1/nx)*kc*rh*rm_x*tlr_el # uM min-1

    # ribosomes
    rtca1 = (atp*rtca)/(atp+(km_a*rd)) 
    rtcb1 = (atp*rtcb)/(atp+(km_b*rt)) 

    # rtcb1 = k_b*atp*rtcb/(k_b*atp+krep*rt)

    Vrep = krep*rtcb1*rt
    Vdam = kdam*rh
    Vinflux = kin*tlr_el
    Vtag = ktag*rtca1*rd


    # ODEs
    dz[1] = tscr_a - dil(rm_a) - deg(rm_a)
    dz[2] = tlr(rm_a, na) - dil(rtca)    
    dz[3] = tscr_b - dil(rm_b) - deg(rm_b)
    dz[4] = tlr(rm_b, nb) - dil(rtcb) 
    dz[5] = tscr_r - dil(rm_r) - deg(rm_r)
    dz[6] = tlr(rm_r, nr) - dil(rtcr)
    dz[7] = Vrep - Vdam + Vinflux - dil(rh)
    dz[8] = Vdam - Vtag - kdeg*rd - dil(rd)
    dz[9] = Vtag - Vrep - dil(rt)
    
    dz
    
end


rtc_mod(z, p) = rtc_mod!(similar(z), z, p, 0)
