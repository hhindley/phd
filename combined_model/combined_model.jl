function combined_model(initial, params, t)
    d, kr, L, c, w_rh, w_t, w_m, w_q, w_BA, w_R, θ_rh, θ_nr, Kq, hq, Vmax_init, Km_init, nrh, nx, nR, nA, nB, gmax, Kgamma, M, kb, ku, abx, kon, koff, vt, vm, s0, Kt, Km, km_a, km_b, krep, kdam_p, ktag, kdeg, k_diss, ns = params
    m_rh, m_t, m_m, m_q, m_R, m_A, m_B, c_rh, c_t, c_m, c_q, c_R, c_A, c_B, et, em, q, R, A, B, rh, rt, rd, z_rh, z_t, z_m, z_q, z_R, z_A, z_B, si, a = initial
    
    # MWC
    alpha = rt/kr # unitless 
    fa = (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6) # unitless
    ra = fa*R # uM 

    # transcription 
    Voc = Vmax_init*a/(Km_init+a) # uM min-1 
    sig_o = ra*Voc/k_diss
    ω_rtcBA(θ_x) = sig_o * w_BA*a/(θ_x + a) # transcription rate for RtcBA # uM min-1

    ω_p(w_x, θ_x) = w_x*a/(θ_x + a) # transcription rate for rh, et, em, rtcR # uM min-1
    ω_q(θ_x) = ((w_q*a/(θ_x + a))/(1+(q/Kq)^hq)) # transcription rate for q # uM min-1

    # translation 
    gamma = gmax*a/(Kgamma+a) # aa min-1 uM-1
    v_x(c_x, nx) = c_x*gamma/nx # min-1
    ttrate = (c_q+c_rh+c_t+c_m+c_R+c_A+c_B)*gamma # aa min-1

    # growth rate 
    λ = ttrate/M # min-1

    # dilution by growth and degradation 
    dil(x) = λ*x # uM min-1 
    deg(x) = d*x # uM min-1 

    # ribosome binding/unbinding
    rh_bind(m_x) = kb*rh*m_x # uM min-1 
    rh_unbind(c_x) = ku*c_x # uM min-1 

    # zombie complex formation 
    zm(c_x) = c_x*abx*kon # uM min-1
    zm_diss(z_x) = koff*z_x # uM min-1 

    # metabolism and import 
    vimp = (et*vt*s0/(Kt+s0)) # uM min-1 
    vcat = (em*vm*si/(Km+si)) # uM min-1 

    # rtc
    Vrep = B*rt*krep/(rt+km_b) # uM min-1
    Vdam = (zm_diss(z_rh)+zm_diss(z_t)+zm_diss(z_m)+zm_diss(z_q)+zm_diss(z_R)+zm_diss(z_A)+zm_diss(z_B))*kdam_p*ku #(z_rh + z_t + z_m + z_q + z_R + z_A + z_B)*koff*kdam # uM min-1
    Vtag = A*rd*ktag/(rd+km_a) # uM min-1

    # ODEs
    # mRNA
    dm_rh = ω_p(w_rh, θ_rh) - dil(m_rh) - deg(m_rh) + v_x(c_rh, nrh) - rh_bind(m_rh) + rh_unbind(c_rh) + zm_diss(z_rh)*ku*kdam_p
    dm_t = ω_p(w_t, θ_nr) - dil(m_t) - deg(m_t) + v_x(c_t, nx) - rh_bind(m_t) + rh_unbind(c_t) + zm_diss(z_t)*ku*kdam_p
    dm_m = ω_p(w_m, θ_nr) - dil(m_m) - deg(m_m) + v_x(c_m, nx) - rh_bind(m_m) + rh_unbind(c_m) + zm_diss(z_m)*ku*kdam_p
    dm_q = ω_q(θ_nr) - dil(m_q) - deg(m_q) + v_x(c_q, nx) - rh_bind(m_q) + rh_unbind(c_q) + zm_diss(z_q)*ku*kdam_p
    dm_R = ω_p(w_R, θ_nr) - dil(m_R) - deg(m_R) + v_x(c_R, nR) - rh_bind(m_R) + rh_unbind(c_R) + zm_diss(z_R)*ku*kdam_p
    dm_A = ω_rtcBA(θ_nr) - dil(m_A) - deg(m_A) + v_x(c_A, nA) - rh_bind(m_A) + rh_unbind(c_A) + zm_diss(z_A)*ku*kdam_p
    dm_B = ω_rtcBA(θ_nr) - dil(m_B) - deg(m_B) + v_x(c_B, nB) - rh_bind(m_B) + rh_unbind(c_B) + zm_diss(z_B)*ku*kdam_p

    # mRNA:ribosome complexes 
    dc_rh = rh_bind(m_rh) - dil(c_rh) - rh_unbind(c_rh) - v_x(c_rh, nrh) - zm(c_rh) + zm_diss(z_rh)*(1-kdam_p)
    dc_t = rh_bind(m_t) - dil(c_t) - rh_unbind(c_t) - v_x(c_t, nx) - zm(c_t) + zm_diss(z_t)*(1-kdam_p)
    dc_m = rh_bind(m_m) - dil(c_m) - rh_unbind(c_m) - v_x(c_m, nx) - zm(c_m) + zm_diss(z_m)*(1-kdam_p)
    dc_q = rh_bind(m_q) - dil(c_q) - rh_unbind(c_q) - v_x(c_q, nx) - zm(c_q) + zm_diss(z_q)*(1-kdam_p)
    dc_R = rh_bind(m_R) - dil(c_R) - rh_unbind(c_R) - v_x(c_R, nR) - zm(c_R) + zm_diss(z_R)*(1-kdam_p)
    dc_A = rh_bind(m_A) - dil(c_A) - rh_unbind(c_A) - v_x(c_A, nA) - zm(c_A) + zm_diss(z_A)*(1-kdam_p)
    dc_B = rh_bind(m_B) - dil(c_B) - rh_unbind(c_B) - v_x(c_B, nB) - zm(c_B) + zm_diss(z_B)*(1-kdam_p)

    # proteins 
    det = v_x(c_t, nx) - dil(et)
    dem = v_x(c_m, nx) - dil(em)
    dq = v_x(c_q, nx) - dil(q)
    dR = v_x(c_R, nR) - dil(R)
    dA = v_x(c_A, nA) - dil(A)
    dB = v_x(c_B, nB) - dil(B)

    # ribosomes
    drh = v_x(c_rh, nrh) - dil(rh) - (rh_bind(m_rh) + rh_bind(m_t) + rh_bind(m_m) + rh_bind(m_q) + rh_bind(m_R) + rh_bind(m_A) + rh_bind(m_B)) + (rh_unbind(c_rh) + rh_unbind(c_t) + rh_unbind(c_m) + rh_unbind(c_q) + rh_unbind(c_R) + rh_unbind(c_A) + rh_unbind(c_B)) + (v_x(c_rh, nrh) + v_x(c_t, nx) + v_x(c_m, nx) + v_x(c_q, nx) + v_x(c_R, nR) + v_x(c_A, nA) + v_x(c_B, nB)) + Vrep
    drt = Vtag - Vrep - dil(rt)
    drd = Vdam - Vtag - rd*kdeg - dil(rd)

    # zombie complexes 
    dz_rh = zm(c_rh) - zm_diss(z_rh)*(1-kdam_p) - zm_diss(z_rh)*kdam_p - dil(z_rh)
    dz_t = zm(c_t) - zm_diss(z_t)*(1-kdam_p) - zm_diss(z_t)*kdam_p - dil(z_t)
    dz_m = zm(c_m) - zm_diss(z_m)*(1-kdam_p) - zm_diss(z_m)*kdam_p - dil(z_m)
    dz_q = zm(c_q) - zm_diss(z_q)*(1-kdam_p) - zm_diss(z_q)*kdam_p - dil(z_q)
    dz_R = zm(c_R) - zm_diss(z_R)*(1-kdam_p) - zm_diss(z_R)*kdam_p - dil(z_R)
    dz_A = zm(c_A) - zm_diss(z_A)*(1-kdam_p) - zm_diss(z_A)*kdam_p - dil(z_A)
    dz_B = zm(c_B) - zm_diss(z_B)*(1-kdam_p) - zm_diss(z_B)*kdam_p - dil(z_B)

    # nutrient
    dsi = vimp - vcat - dil(si)

    # energy
    da = ns*vcat - ttrate - dil(a)

    [dm_rh, dm_t, dm_m, dm_q, dm_R, dm_A, dm_B, dc_rh, dc_t, dc_m, dc_q, dc_R, dc_A, dc_B, det, dem, dq, dR, dA, dB, drh, drt, drd, dz_rh, dz_t, dz_m, dz_q, dz_R, dz_A, dz_B, dsi, da]

end



function combined_model_v2(initial, params, t)
    d, kr, L, c, w_rh, w_t, w_m, w_q, w_BA, w_R, θ_rh, θ_nr, Kq, hq, Vmax_init, Km_init, nrh, nx, nR, nA, nB, gmax, Kgamma, M, kb, ku, abx, kon, koff, vt, vm, s0, Kt, Km, km_a, km_b, krep, kdam, ktag, kdeg, kdiss, ns = params
    m_rh, m_t, m_m, m_q, m_R, m_A, m_B, c_rh, c_t, c_m, c_q, c_R, c_A, c_B, et, em, q, R, A, B, rh, rt, rd, z_rh, z_t, z_m, z_q, z_R, z_A, z_B, c_rhD, c_tD, c_mD, c_qD, c_RD, c_AD, c_BD, si, a = initial
    # m_rh, m_t, m_m, m_q, m_R, m_A, m_B, c_rh, c_t, c_m, c_q, c_R, c_A, c_B, et, em, q, R, A, B, rh, rt, rd, z_rh, z_t, z_m, z_q, z_R, z_A, z_B, si, a = initial
    
    # MWC
    alpha = rt/kr # unitless 
    fa = (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6) # unitless
    ra = fa*R # uM 

    # transcription 
    Voc = Vmax_init*a/(Km_init+a) # uM min-1 
    sig_o = ra*Voc/kdiss
    ω_rtcBA(θ_x) = (w_BA*a/(θ_x + a)) * sig_o # transcription rate for RtcBA # uM min-1

    ω_p(w_x, θ_x) = w_x*a/(θ_x + a) # transcription rate for rh, et, em, rtcR # uM min-1
    ω_q(θ_x) = ((w_q*a/(θ_x + a))/(1+(q/Kq)^hq)) # transcription rate for q # uM min-1

    # translation 
    gamma = gmax*a/(Kgamma+a) # aa min-1 uM-1
    v_x(c_x, nx) = c_x*gamma/nx # min-1
    ttrate = (c_q+c_rh+c_t+c_m+c_R+c_A+c_B)*gamma # aa min-1

    # growth rate 
    λ = ttrate/M # min-1

    # dilution by growth and degradation 
    dil(x) = λ*x # uM min-1 
    deg(x) = d*x # uM min-1 

    # ribosome binding/unbinding
    rh_bind(m_x) = kb*rh*m_x # uM min-1 
    rh_unbind(c_x) = ku*c_x # uM min-1 

    # zombie complex formation 
    zm(c_x) = c_x*abx*kon # uM min-1
    zm_diss(z_x) = koff*z_x # uM min-1 

    # metabolism and import 
    vimp = (et*vt*s0/(Kt+s0)) # uM min-1 
    vcat = (em*vm*si/(Km+si)) # uM min-1 

    # rtc
    Vrep = B*rt*krep/(rt+km_b) # uM min-1
    Vdam = (c_rhD + c_tD + c_mD + c_qD + c_RD + c_AD + c_BD) * ku 
    # Vdam = (zm_diss(z_rh)+zm_diss(z_t)+zm_diss(z_m)+zm_diss(z_q)+zm_diss(z_R)+zm_diss(z_A)+zm_diss(z_B))*kdam*ku #(z_rh + z_t + z_m + z_q + z_R + z_A + z_B)*koff*kdam # uM min-1
    Vtag = A*rd*ktag/(rd+km_a) # uM min-1

    # ODEs
    # mRNA
    dm_rh = ω_p(w_rh, θ_rh) - dil(m_rh) - deg(m_rh) + v_x(c_rh, nrh) - rh_bind(m_rh) + rh_unbind(c_rh) + rh_unbind(c_rhD)
    dm_t = ω_p(w_t, θ_nr) - dil(m_t) - deg(m_t) + v_x(c_t, nx) - rh_bind(m_t) + rh_unbind(c_t) + rh_unbind(c_tD)
    dm_m = ω_p(w_m, θ_nr) - dil(m_m) - deg(m_m) + v_x(c_m, nx) - rh_bind(m_m) + rh_unbind(c_m) + rh_unbind(c_mD)
    dm_q = ω_q(θ_nr) - dil(m_q) - deg(m_q) + v_x(c_q, nx) - rh_bind(m_q) + rh_unbind(c_q) + rh_unbind(c_qD)
    dm_R = ω_p(w_R, θ_nr) - dil(m_R) - deg(m_R) + v_x(c_R, nR) - rh_bind(m_R) + rh_unbind(c_R) + rh_unbind(c_RD)
    dm_A = ω_rtcBA(θ_nr) - dil(m_A) - deg(m_A) + v_x(c_A, nA) - rh_bind(m_A) + rh_unbind(c_A) + rh_unbind(c_AD)
    dm_B = ω_rtcBA(θ_nr) - dil(m_B) - deg(m_B) + v_x(c_B, nB) - rh_bind(m_B) + rh_unbind(c_B) + rh_unbind(c_BD)

    # mRNA:ribosome complexes 
    dc_rh = rh_bind(m_rh) - dil(c_rh) - rh_unbind(c_rh) - v_x(c_rh, nrh) - zm(c_rh) + zm_diss(z_rh)*(1-kdam)
    dc_t = rh_bind(m_t) - dil(c_t) - rh_unbind(c_t) - v_x(c_t, nx) - zm(c_t) + zm_diss(z_t)*(1-kdam)
    dc_m = rh_bind(m_m) - dil(c_m) - rh_unbind(c_m) - v_x(c_m, nx) - zm(c_m) + zm_diss(z_m)*(1-kdam)
    dc_q = rh_bind(m_q) - dil(c_q) - rh_unbind(c_q) - v_x(c_q, nx) - zm(c_q) + zm_diss(z_q)*(1-kdam)
    dc_R = rh_bind(m_R) - dil(c_R) - rh_unbind(c_R) - v_x(c_R, nR) - zm(c_R) + zm_diss(z_R)*(1-kdam)
    dc_A = rh_bind(m_A) - dil(c_A) - rh_unbind(c_A) - v_x(c_A, nA) - zm(c_A) + zm_diss(z_A)*(1-kdam)
    dc_B = rh_bind(m_B) - dil(c_B) - rh_unbind(c_B) - v_x(c_B, nB) - zm(c_B) + zm_diss(z_B)*(1-kdam)

    # proteins 
    det = v_x(c_t, nx) - dil(et)
    dem = v_x(c_m, nx) - dil(em)
    dq = v_x(c_q, nx) - dil(q)
    dR = v_x(c_R, nR) - dil(R)
    dA = v_x(c_A, nA) - dil(A)
    dB = v_x(c_B, nB) - dil(B)

    # ribosomes
    drh = v_x(c_rh, nrh) - dil(rh) - (rh_bind(m_rh) + rh_bind(m_t) + rh_bind(m_m) + rh_bind(m_q)) + (rh_unbind(c_rh) + rh_unbind(c_t) + rh_unbind(c_m) + rh_unbind(c_q)) + (v_x(c_rh, nrh) + v_x(c_t, nx) + v_x(c_m, nx) + v_x(c_q, nx)) + Vrep
    drt = Vtag - Vrep - dil(rt)
    drd = Vdam - Vtag - rd*kdeg - dil(rd)

    # zombie complexes 
    dz_rh = zm(c_rh) - zm_diss(z_rh)*(1-kdam) - zm_diss(z_rh)*kdam - dil(z_rh)
    dz_t = zm(c_t) - zm_diss(z_t)*(1-kdam) - zm_diss(z_t)*kdam - dil(z_t)
    dz_m = zm(c_m) - zm_diss(z_m)*(1-kdam) - zm_diss(z_m)*kdam - dil(z_m)
    dz_q = zm(c_q) - zm_diss(z_q)*(1-kdam) - zm_diss(z_q)*kdam - dil(z_q)
    dz_R = zm(c_R) - zm_diss(z_R)*(1-kdam) - zm_diss(z_R)*kdam - dil(z_R)
    dz_A = zm(c_A) - zm_diss(z_A)*(1-kdam) - zm_diss(z_A)*kdam - dil(z_A)
    dz_B = zm(c_B) - zm_diss(z_B)*(1-kdam) - zm_diss(z_B)*kdam - dil(z_B)

    dc_rhD = zm_diss(z_rh)*kdam - rh_unbind(c_rhD) - dil(c_rhD)
    dc_tD = zm_diss(z_t)*kdam - rh_unbind(c_tD) - dil(c_tD)
    dc_mD = zm_diss(z_m)*kdam - rh_unbind(c_mD) - dil(c_mD)
    dc_qD = zm_diss(z_q)*kdam - rh_unbind(c_qD) - dil(c_qD)
    dc_RD = zm_diss(z_R)*kdam - rh_unbind(c_RD) - dil(c_RD)
    dc_AD = zm_diss(z_A)*kdam - rh_unbind(c_AD) - dil(c_AD)
    dc_BD = zm_diss(z_B)*kdam - rh_unbind(c_BD) - dil(c_BD)

    # nutrient
    dsi = vimp - vcat - dil(si)

    # energy
    da = ns*vcat - ttrate - dil(a)

    [dm_rh, dm_t, dm_m, dm_q, dm_R, dm_A, dm_B, dc_rh, dc_t, dc_m, dc_q, dc_R, dc_A, dc_B, det, dem, dq, dR, dA, dB, drh, drt, drd, dz_rh, dz_t, dz_m, dz_q, dz_R, dz_A, dz_B, dc_rhD, dc_tD, dc_mD, dc_qD, dc_RD, dc_AD, dc_BD, dsi, da]
    # [dm_rh, dm_t, dm_m, dm_q, dm_R, dm_A, dm_B, dc_rh, dc_t, dc_m, dc_q, dc_R, dc_A, dc_B, det, dem, dq, dR, dA, dB, drh, drt, drd, dz_rh, dz_t, dz_m, dz_q, dz_R, dz_A, dz_B, dsi, da]

end