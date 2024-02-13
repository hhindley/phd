
#trna

function rtc_trna_inhib_model_rtcb(initial, params, t) 
    L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, rh, thr_t, k_inhib1, k_inhib2, inhib = params
    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, trna, rd, rt, rtc_i = initial


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
    tlr_el = (g_max*atp/(θtlr+atp)) * trna/(thr_t+trna) 
    tlr(rm_x, nx) = (1/nx)*rh*rm_x*tlr_el # *1/nx nx = length of RtcA, RtcB and RtcR


    # rtcb_i = inhib*k_inhib*rtcb/((k_inhib*inhib)+k_inhib)
    # rtcb_a = rtcb - rtcb_i

    # ribosomes
    rtca1 = (atp*rtca)/(atp+(km_a*rd)) 
    rtcb1 = (atp*rtcb)/(atp+(km_b*rt)) 
    # rtcb1 = k_b*atp*rtcb/(k_b*atp+krep*rt)

    Vrep = krep*rtcb1*rt
    Vdam = kdam*trna
    Vinflux = kin*tlr_el
    Vtag = ktag*rtca1*rd

    # ODEs
    drm_a = tscr_a - dil(rm_a) - deg(rm_a)
    drtca = tlr(rm_a, na) - dil(rtca)    
    drm_b = tscr_b - dil(rm_b) - deg(rm_b)
    drtcb = tlr(rm_b, nb) - dil(rtcb) - k_inhib1*rtcb*inhib + k_inhib2*rtc_i
    drm_r = tscr_r - dil(rm_r) - deg(rm_r)
    drtcr = tlr(rm_r, nr) - dil(rtcr)
    dtrna = Vrep - Vdam + Vinflux - dil(trna)
    drd = Vdam - Vtag - kdeg*rd - dil(rd)
    drt = Vtag - Vrep - dil(rt)
    drtc_i = k_inhib1*rtcb*inhib - k_inhib2*rtc_i - dil(rtc_i)
    # @SVector [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt]
    [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, dtrna, drd, drt, drtc_i]

end


function rtc_trna_inhib_mod_rtcb!(dz, z, p, t)
    @unpack L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, rh, thr_t, k_inhib1, k_inhib2, inhib = p
    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, trna, rd, rt, rtc_i = z


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
    tlr_el = (g_max*atp/(θtlr+atp)) * trna/(thr_t+trna) 
    tlr(rm_x, nx) = (1/nx)*rh*rm_x*tlr_el # *1/nx nx = length of RtcA, RtcB and RtcR


    # rtcb_i = inhib*k_inhib*rtcb/((k_inhib*inhib)+k_inhib)
    # rtcb_a = rtcb - rtcb_i

    # ribosomes
    rtca1 = (atp*rtca)/(atp+(km_a*rd)) 
    rtcb1 = (atp*rtcb)/(atp+(km_b*rt)) 
    # rtcb1 = k_b*atp*rtcb/(k_b*atp+krep*rt)

    Vrep = krep*rtcb1*rt
    Vdam = kdam*trna
    Vinflux = kin*tlr_el
    Vtag = ktag*rtca1*rd

    # ODEs
    dz[1] = tscr_a - dil(rm_a) - deg(rm_a)
    dz[2] = tlr(rm_a, na) - dil(rtca)    
    dz[3] = tscr_b - dil(rm_b) - deg(rm_b)
    dz[4] = tlr(rm_b, nb) - dil(rtcb) - k_inhib1*rtcb*inhib + k_inhib2*rtc_i
    dz[5] = tscr_r - dil(rm_r) - deg(rm_r)
    dz[6] = tlr(rm_r, nr) - dil(rtcr)
    dz[7] = Vrep - Vdam + Vinflux - dil(trna)
    dz[8] = Vdam - Vtag - kdeg*rd - dil(rd)
    dz[9] = Vtag - Vrep - dil(rt)
    dz[10] = k_inhib1*rtcb*inhib - k_inhib2*rtc_i - dil(rtc_i)

    dz
end


rtc_trna_inhib_mod_rtcb(z, p) = rtc_trna_inhib_mod_rtcb!(similar(z), z, p, 0)








function rtc_trna_inhib_model_rtca(initial, params, t) 
    L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, rh, thr_t, k_inhib1, k_inhib2, inhib = params
    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, trna, rd, rt, rtc_i = initial


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
    tlr_el = (g_max*atp/(θtlr+atp)) * trna/(thr_t+trna) 
    tlr(rm_x, nx) = (1/nx)*rh*rm_x*tlr_el # *1/nx nx = length of RtcA, RtcB and RtcR


    # rtcb_i = inhib*k_inhib*rtcb/((k_inhib*inhib)+k_inhib)
    # rtcb_a = rtcb - rtcb_i

    # ribosomes
    rtca1 = (atp*rtca)/(atp+(km_a*rd)) 
    rtcb1 = (atp*rtcb)/(atp+(km_b*rt)) 
    # rtcb1 = k_b*atp*rtcb/(k_b*atp+krep*rt)

    Vrep = krep*rtcb1*rt
    Vdam = kdam*trna
    Vinflux = kin*tlr_el
    Vtag = ktag*rtca1*rd

    # ODEs
    drm_a = tscr_a - dil(rm_a) - deg(rm_a)
    drtca = tlr(rm_a, na) - dil(rtca) - k_inhib1*rtca*inhib + k_inhib2*rtc_i
    drm_b = tscr_b - dil(rm_b) - deg(rm_b)
    drtcb = tlr(rm_b, nb) - dil(rtcb) 
    drm_r = tscr_r - dil(rm_r) - deg(rm_r)
    drtcr = tlr(rm_r, nr) - dil(rtcr)
    dtrna = Vrep - Vdam + Vinflux - dil(trna)
    drd = Vdam - Vtag - kdeg*rd - dil(rd)
    drt = Vtag - Vrep - dil(rt)
    drtc_i = k_inhib1*rtca*inhib - k_inhib2*rtc_i - dil(rtc_i)
    # @SVector [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt]
    [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, dtrna, drd, drt, drtc_i]

end


function rtc_trna_inhib_mod_rtca!(dz, z, p, t)
    @unpack L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, rh, thr_t, k_inhib1, k_inhib2, inhib = p
    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, trna, rd, rt, rtc_i = z


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
    tlr_el = (g_max*atp/(θtlr+atp)) * trna/(thr_t+trna) 
    tlr(rm_x, nx) = (1/nx)*rh*rm_x*tlr_el # *1/nx nx = length of RtcA, RtcB and RtcR


    # rtcb_i = inhib*k_inhib*rtcb/((k_inhib*inhib)+k_inhib)
    # rtcb_a = rtcb - rtcb_i

    # ribosomes
    rtca1 = (atp*rtca)/(atp+(km_a*rd)) 
    rtcb1 = (atp*rtcb)/(atp+(km_b*rt)) 
    # rtcb1 = k_b*atp*rtcb/(k_b*atp+krep*rt)

    Vrep = krep*rtcb1*rt
    Vdam = kdam*trna
    Vinflux = kin*tlr_el
    Vtag = ktag*rtca1*rd

    # ODEs
    dz[1] = tscr_a - dil(rm_a) - deg(rm_a)
    dz[2] = tlr(rm_a, na) - dil(rtca) - k_inhib1*rtca*inhib + k_inhib2*rtc_i
    dz[3] = tscr_b - dil(rm_b) - deg(rm_b)
    dz[4] = tlr(rm_b, nb) - dil(rtcb) 
    dz[5] = tscr_r - dil(rm_r) - deg(rm_r)
    dz[6] = tlr(rm_r, nr) - dil(rtcr)
    dz[7] = Vrep - Vdam + Vinflux - dil(trna)
    dz[8] = Vdam - Vtag - kdeg*rd - dil(rd)
    dz[9] = Vtag - Vrep - dil(rt)
    dz[10] = k_inhib1*rtca*inhib - k_inhib2*rtc_i - dil(rtc_i)
    dz
end


rtc_trna_inhib_mod_rtca(z, p) = rtc_trna_inhib_mod_rtca!(similar(z), z, p, 0)






function rtc_trna_inhib_model_rtcr(initial, params, t) 
    L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, rh, thr_t, k_inhib1, k_inhib2, inhib = params
    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, trna, rd, rt, rtc_i = initial


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
    tlr_el = (g_max*atp/(θtlr+atp)) * trna/(thr_t+trna) 
    tlr(rm_x, nx) = (1/nx)*rh*rm_x*tlr_el # *1/nx nx = length of RtcA, RtcB and RtcR


    # rtcb_i = inhib*k_inhib*rtcb/((k_inhib*inhib)+k_inhib)
    # rtcb_a = rtcb - rtcb_i

    # ribosomes
    rtca1 = (atp*rtca)/(atp+(km_a*rd)) 
    rtcb1 = (atp*rtcb)/(atp+(km_b*rt)) 
    # rtcb1 = k_b*atp*rtcb/(k_b*atp+krep*rt)

    Vrep = krep*rtcb1*rt
    Vdam = kdam*trna
    Vinflux = kin*tlr_el
    Vtag = ktag*rtca1*rd

    # ODEs
    drm_a = tscr_a - dil(rm_a) - deg(rm_a)
    drtca = tlr(rm_a, na) - dil(rtca) 
    drm_b = tscr_b - dil(rm_b) - deg(rm_b)
    drtcb = tlr(rm_b, nb) - dil(rtcb) 
    drm_r = tscr_r - dil(rm_r) - deg(rm_r)
    drtcr = tlr(rm_r, nr) - dil(rtcr) -  k_inhib1*rtcr*inhib + k_inhib2*rtc_i
    dtrna = Vrep - Vdam + Vinflux - dil(trna)
    drd = Vdam - Vtag - kdeg*rd - dil(rd)
    drt = Vtag - Vrep - dil(rt)
    drtc_i = k_inhib1*rtcr*inhib - k_inhib2*rtc_i - dil(rtc_i)
    # @SVector [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt]
    [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, dtrna, drd, drt, drtc_i]

end


function rtc_trna_inhib_mod_rtcr!(dz, z, p, t)
    @unpack L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, rh, thr_t, k_inhib1, k_inhib2, inhib = p
    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, trna, rd, rt, rtc_i = z


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
    tlr_el = (g_max*atp/(θtlr+atp)) * trna/(thr_t+trna) 
    tlr(rm_x, nx) = (1/nx)*rh*rm_x*tlr_el # *1/nx nx = length of RtcA, RtcB and RtcR


    # rtcb_i = inhib*k_inhib*rtcb/((k_inhib*inhib)+k_inhib)
    # rtcb_a = rtcb - rtcb_i

    # ribosomes
    rtca1 = (atp*rtca)/(atp+(km_a*rd)) 
    rtcb1 = (atp*rtcb)/(atp+(km_b*rt)) 
    # rtcb1 = k_b*atp*rtcb/(k_b*atp+krep*rt)

    Vrep = krep*rtcb1*rt
    Vdam = kdam*trna
    Vinflux = kin*tlr_el
    Vtag = ktag*rtca1*rd

    # ODEs
    dz[1] = tscr_a - dil(rm_a) - deg(rm_a)
    dz[2] = tlr(rm_a, na) - dil(rtca) 
    dz[3] = tscr_b - dil(rm_b) - deg(rm_b)
    dz[4] = tlr(rm_b, nb) - dil(rtcb) 
    dz[5] = tscr_r - dil(rm_r) - deg(rm_r)
    dz[6] = tlr(rm_r, nr) - dil(rtcr) - k_inhib1*rtcr*inhib + k_inhib2*rtc_i
    dz[7] = Vrep - Vdam + Vinflux - dil(trna)
    dz[8] = Vdam - Vtag - kdeg*rd - dil(rd)
    dz[9] = Vtag - Vrep - dil(rt)
    dz[10] = k_inhib1*rtcr*inhib - k_inhib2*rtc_i - dil(rtc_i)
    dz
end


rtc_trna_inhib_mod_rtcr(z, p) = rtc_trna_inhib_mod_rtcr!(similar(z), z, p, 0)
