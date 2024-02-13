function rtc_prototype!(dydt, initial, params, t) # USE THIS MODEL 
    kr, L, c, sigma, atp, w_rtc, theta_rtc, max, thr, d, k1, k2, k3, k4, k5, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, d1 = params
    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, ribo_h, ribo_d, ribo_t = initial 
    drm_adt, drtcadt, drm_bdt, drtcbdt, drm_rdt, drtcrdt, dribo_hdt, dribo_ddt, dribo_tdt = zeros(length(dydt))

    rdrtca = k1_a*ribo_d*rtca/(k1_a*ribo_d+k2_a+k3_a*atp)
    rtrtcb = ka_b*ribo_t*rtcb/(ka_b*ribo_t+kb_b+kc_b*atp)

    gr = 0.01*ribo_h
    influx = 0.4
    alpha = ribo_t/kr
    fa = (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6)
    ra = fa*rtcr
    v = ra*k1*sigma*k3*atp*k5/(k1*sigma*k3*atp+k1*sigma*k4+k1*sigma*k5+k2*k4+k2*k5+k3*atp*k5)
    tscr_a = v*(w_rtc*atp/(theta_rtc+atp))
    tscr_b = v*(w_rtc*atp/(theta_rtc+atp))
    tscr_r = w_rtc*atp/(theta_rtc+atp)    
    tlr_a = ribo_h*rm_a*(max*atp/(thr+atp))
    tlr_b = ribo_h*rm_b*(max*atp/(thr+atp))
    tlr_r = ribo_h*rm_r*(max*atp/(thr+atp))

    dydt[1] = tscr_a - gr*rm_a - d*rm_a
    dydt[2] = tlr_a - gr*rtca 
    dydt[3] = tscr_b - gr*rm_b - d*rm_b
    dydt[4] = tlr_b - gr*rtcb 
    dydt[5] = tscr_r - gr*rm_r - d*rm_r
    dydt[6] = tlr_r - gr*rtcr

    dydt[7] = kc_b*(rtrtcb)*atp - k*ribo_h - gr*ribo_h + influx*(max*atp/(thr+atp))
    dydt[8] = k*ribo_h - k3_a*(rdrtca)*atp - gr*ribo_d - d1*ribo_d   
    dydt[9] = k3_a*(rdrtca)*atp - kc_b*(rtrtcb)*atp - gr*ribo_t

end