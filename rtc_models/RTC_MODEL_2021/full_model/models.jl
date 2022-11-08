function rtc_system_ribos!(dydt, initial, params, t) # 1 - all ribosomes as ODEs, ints algebraic 

    kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, rtca_tot, rtcb_tot = params

    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, ribo_d, ribo_h, ribo_t = initial 
    drm_adt, drtcadt, drm_bdt, drtcbdt, drm_rdt, drtcrdt, dribo_ddt, dribo_hdt, dribo_tdt = zeros(length(dydt))

    n = 6
    gr = 0.01*ribo_h # NEEDS CHANGING BACK TO 0.01

    rdrtca = k1_a*ribo_d*rtca_tot/(k1_a*ribo_d+k2_a+k3_a*atp)
    rtrtcb = ka_b*ribo_t*rtcb_tot/(ka_b*ribo_t+kb_b+kc_b*atp)

    alpha = ribo_t/kr
    fa = (1+alpha)^n/(L*((1+c*alpha)^n)+(1+alpha)^n) # fraction of active RtcR
    ra = fa*rtcr # amount of active RtcR
    v = ra*k1*sigma*k3*atp*k5/(k1*sigma*k3*atp+k1*sigma*(k4+k5)+k2*(k4+k5)+k3*atp*k5) # rate of open complex formation
    
    # rtca
    tscr = v*(w_rtc*atp/(theta_rtc+atp)) # transcription rate 
    tlel = max*atp/(thr+atp) # translation elongation rate
    tlr = ribo_h*rm_a*tlel # translation rate so formation of Rtc 

    # rtcb
    tscr_b = v*(w_rtc*atp/(theta_rtc+atp)) # transcription rate 
    tlel_b = max*atp/(thr+atp) # translation elongation rate
    tlr_b = ribo_h*rm_b*tlel_b # translation rate so formation of Rtc 

    # rtcr 
    rtcr_tscr = w_rtc*atp/(theta_rtc+atp)
    rtcr_tlel = max*atp/(thr+atp)
    rtcr_tlr = ribo_h*rm_r*rtcr_tlel

    dydt[1] = tscr - gr*rm_a - d*rm_a
    dydt[2] = tlr - k1_a*rtca*ribo_d + k2_a*rdrtca + k3_a*atp*rdrtca  - gr*rtca 
    dydt[3] = tscr_b - gr*rm_b - d*rm_b 
    dydt[4] = tlr_b - ka_b*rtcb*ribo_t + kb_b*rtrtcb + kc_b*atp*rtrtcb - gr*rtcb
    dydt[5] = rtcr_tscr - gr*rm_r - d*rm_r
    dydt[6] = rtcr_tlr - gr*rtcr

    dydt[7] = k*ribo_h - k1_a*rtca*ribo_d + k2_a*rdrtca
    dydt[8] = kc_b*atp*rtrtcb - k*ribo_h 
    dydt[9] = k3_a*atp*rdrtca - ka_b*rtcb*ribo_t + kb_b*rtrtcb


end

function rtc_system_ribos_h!(dydt, initial, params, t) # 2 - rh and ints as algebraic 

    kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, rtca_tot, rtcb_tot, ribo_tot = params

    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, ribo_d, ribo_t = initial 
    drm_adt, drtcadt, drm_bdt, drtcbdt, drm_rdt, drtcrdt, dribo_ddt, dribo_tdt = zeros(length(dydt))

    n = 6
    rdrtca = k1_a*ribo_d*rtca_tot/(k1_a*ribo_d+k2_a+k3_a*atp)
    rtrtcb = ka_b*ribo_t*rtcb_tot/(ka_b*ribo_t+kb_b+kc_b*atp)
   
    ribo_h = ribo_tot - ribo_d - ribo_t - rdrtca - rtrtcb
    gr = 0.01*ribo_h

    alpha = ribo_t/kr
    fa = (1+alpha)^n/(L*((1+c*alpha)^n)+(1+alpha)^n) # fraction of active RtcR
    ra = fa*rtcr # amount of active RtcR
    v = ra*k1*sigma*k3*atp*k5/(k1*sigma*k3*atp+k1*sigma*(k4+k5)+k2*(k4+k5)+k3*atp*k5) # rate of open complex formation
    
    # rtca
    tscr = v*(w_rtc*atp/(theta_rtc+atp)) # transcription rate 
    tlel = max*atp/(thr+atp) # translation elongation rate
    tlr = ribo_h*rm_a*tlel # translation rate so formation of Rtc 

    # rtcb
    tscr_b = v*(w_rtc*atp/(theta_rtc+atp)) # transcription rate 
    tlel_b = max*atp/(thr+atp) # translation elongation rate
    tlr_b = ribo_h*rm_b*tlel_b # translation rate so formation of Rtc 

    # rtcr 
    rtcr_tscr = w_rtc*atp/(theta_rtc+atp)
    rtcr_tlel = max*atp/(thr+atp)
    rtcr_tlr = ribo_h*rm_r*rtcr_tlel

    dydt[1] = tscr - gr*rm_a - d*rm_a
    dydt[2] = tlr - k1_a*rtca*ribo_d + k2_a*rdrtca + k3_a*atp*rdrtca  - gr*rtca 
    dydt[3] = tscr_b - gr*rm_b - d*rm_b 
    dydt[4] = tlr_b - ka_b*rtcb*ribo_t + kb_b*rtrtcb + kc_b*atp*rtrtcb - gr*rtcb
    dydt[5] = rtcr_tscr - gr*rm_r - d*rm_r
    dydt[6] = rtcr_tlr - gr*rtcr

    dydt[7] = k*ribo_h - k1_a*rtca*ribo_d + k2_a*rdrtca
    dydt[8] = k3_a*atp*rdrtca - ka_b*rtcb*ribo_t + kb_b*rtrtcb


end

function rtc_system_ribos_ints!(dydt, initial, params, t) # 3 - all ribosomes as ODEs and ints as odes

    kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k = params

    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rdrtca, rtrtcb, ribo_d, ribo_h, ribo_t = initial 
    drm_adt, drtcadt, drm_bdt, drtcbdt, drm_rdt, drtcrdt, drdrtcadt, drtrtcbdt, dribo_ddt, dribo_hdt, dribo_tdt = zeros(length(dydt))

    n = 6
    gr = 0.01*ribo_h # NEEDS CHANGING BACK TO 0.01
    
    alpha = ribo_t/kr
    fa = (1+alpha)^n/(L*((1+c*alpha)^n)+(1+alpha)^n) # fraction of active RtcR
    ra = fa*rtcr # amount of active RtcR
    v = ra*k1*sigma*k3*atp*k5/(k1*sigma*k3*atp+k1*sigma*(k4+k5)+k2*(k4+k5)+k3*atp*k5) # rate of open complex formation
    
    # rtca
    tscr = v*(w_rtc*atp/(theta_rtc+atp)) # transcription rate 
    tlel = max*atp/(thr+atp) # translation elongation rate
    tlr = ribo_h*rm_a*tlel # translation rate so formation of Rtc 

    # rtcb
    tscr_b = v*(w_rtc*atp/(theta_rtc+atp)) # transcription rate 
    tlel_b = max*atp/(thr+atp) # translation elongation rate
    tlr_b = ribo_h*rm_b*tlel_b # translation rate so formation of Rtc 

    # rtcr 
    rtcr_tscr = w_rtc*atp/(theta_rtc+atp)
    rtcr_tlel = max*atp/(thr+atp)
    rtcr_tlr = ribo_h*rm_r*rtcr_tlel



    dydt[1] = tscr - gr*rm_a - d*rm_a
    dydt[2] = tlr - k1_a*rtca*ribo_d + k2_a*rdrtca + k3_a*atp*rdrtca  - gr*rtca 
    dydt[3] = tscr_b - gr*rm_b - d*rm_b 
    dydt[4] = tlr_b - ka_b*rtcb*ribo_t + kb_b*rtrtcb + kc_b*atp*rtrtcb - gr*rtcb
    dydt[5] = rtcr_tscr - gr*rm_r - d*rm_r
    dydt[6] = rtcr_tlr - gr*rtcr

    dydt[7] = k1_a*rtca*ribo_d - k2_a*rdrtca - k3_a*atp*rdrtca
    dydt[8] = ka_b*rtcb*ribo_t - kb_b*rtrtcb - kc_b*atp*rtrtcb

    dydt[9] = k*ribo_h - k1_a*rtca*ribo_d + k2_a*rdrtca
    dydt[10] = kc_b*atp*rtrtcb - k*ribo_h 
    dydt[11] = k3_a*atp*rdrtca - ka_b*rtcb*ribo_t + kb_b*rtrtcb


end

function rtc_system_fullodes!(dydt, initial, params, t) # 4 - ints as odes and rh algebraic

    kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, ribo_tot = params

    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rdrtca, rtrtcb, ribo_d, ribo_t = initial 
    drm_adt, drtcadt, drm_bdt, drtcbdt, drm_rdt, drtcrdt, drdrtcadt, drtrtcbdt, dribo_ddt, dribo_tdt = zeros(length(dydt))

    n = 6
    ribo_h = ribo_tot - ribo_d - ribo_t - rdrtca - rtrtcb
    gr = 0.01*ribo_h

    alpha = ribo_t/kr
    fa = (1+alpha)^n/(L*((1+c*alpha)^n)+(1+alpha)^n) # fraction of active RtcR
    ra = fa*rtcr # amount of active RtcR
    v = ra*k1*sigma*k3*atp*k5/(k1*sigma*k3*atp+k1*sigma*(k4+k5)+k2*(k4+k5)+k3*atp*k5) # rate of open complex formation
    
    # rtca
    tscr = v*(w_rtc*atp/(theta_rtc+atp)) # transcription rate 
    tlel = max*atp/(thr+atp) # translation elongation rate
    tlr = ribo_h*rm_a*tlel # translation rate so formation of Rtc 

    # rtcb
    tscr_b = v*(w_rtc*atp/(theta_rtc+atp)) # transcription rate 
    tlel_b = max*atp/(thr+atp) # translation elongation rate
    tlr_b = ribo_h*rm_b*tlel_b # translation rate so formation of Rtc 

    # rtcr 
    rtcr_tscr = w_rtc*atp/(theta_rtc+atp)
    rtcr_tlel = max*atp/(thr+atp)
    rtcr_tlr = ribo_h*rm_r*rtcr_tlel


    dydt[1] = tscr - gr*rm_a - d*rm_a
    dydt[2] = tlr - k1_a*rtca*ribo_d + k2_a*rdrtca + k3_a*atp*rdrtca  - gr*rtca 
    dydt[3] = tscr_b - gr*rm_b - d*rm_b 
    dydt[4] = tlr_b - ka_b*rtcb*ribo_t + kb_b*rtrtcb + kc_b*atp*rtrtcb - gr*rtcb
    dydt[5] = rtcr_tscr - gr*rm_r - d*rm_r
    dydt[6] = rtcr_tlr - gr*rtcr

    dydt[7] = k1_a*rtca*ribo_d - k2_a*rdrtca - k3_a*atp*rdrtca
    dydt[8] = ka_b*rtcb*ribo_t - kb_b*rtrtcb - kc_b*atp*rtrtcb 
    
    dydt[9] = k*ribo_h - k1_a*rtca*ribo_d + k2_a*rdrtca
    dydt[10] = k3_a*atp*rdrtca - ka_b*rtcb*ribo_t + kb_b*rtrtcb

end