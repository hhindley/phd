import numpy as np

def model(t, y, params, rates): # not sure about having to include rates here

    # LHS
    r = y[0] # free ribosomes
    et = y[1] # transporter enzyme 
    em = y[2] # metabolism enzyme 
    q = y[3] # house-keeping proteins
    si = y[4] # internal nutrient 
    a = y[5] # energy
    mr = y[6] # free mRNA for ribosome
    mt = y[7] # free mRNA for transporter enzyme
    mm = y[8] # free mRNA for metabolism enzyme
    mq = y[9] # free mRNA for house-keeping proteins
    cr = y[10] # ribosome bound mRNA for free ribosome
    ct = y[11] # ribosome bound mRNA for transporter enzyme
    cm = y[12] # ribosome bound mRNA for metabolism enzyme
    cq = y[13] # ribosome bound mRNA for house-keeping proteins

    # params
    s = params[0] # external nutrient 
    dm = params[1] # mRNA-degradation rate 
    ns = params[2] # nutrient quality 
    nr = params[3] # ribosome length 
    nx = params[4] # length of nonribosomal proteins
    gmax = params[5] # max. translational elongation rate 
    Kp = params[6] # protein threshold - where did this come from? 
    vt = params[7] # max. nutrient import rate - should maybe be in rates?
    Kt = params[8] # nutrient import threshold for transporter enzyme
    vm = params[9] # max. enzymatic rate - should maybe be in rates? 
    Km = params[10] # enzymatic threshold for metabolic enzyme 
    wr = params[11] # max. ribosome transcription rate
    we = params[12] # max. enzyme transcription rate
    wq = params[13] # max. q-transcription rate 
    thetar = params[14] # ribosome transcription threshold 
    thetax = params[15] # non-ribosomal transcription threshold 
    Kq = params[16] # q-autoinhibition threshold 
    hq = params[17] # q-autoinhbition Hill coeff.
    kb = params[18] # mRNA-ribosome binding rate 
    ku = params[19] # mRNA-ribosome unbinding rate 
    M = params[20] # total cell mass 
    kcm = params[21] # chloramphenicol binding rate 
    
    
    # rates
    lam = rates[0] # dilution
    dm = rates[2] # degradation
    kb = rates[3] # ribosome binding 
    ku = rates[4] # ribosome unbinding 
    wx = rates[5] # transcription
    vx = rates[6] # translation

    # equations
    Kgamma = gmax/Kp # threshold for half maximal elongation - not sure about value or equation for Kp
    gamma = (gmax*a)/(Kgamma+a) # rate of translational elongation 
    vimp = et((vt*s)/(Kt*s)) # transporter enzyme kinetics 
    nucat = em((vm*si)/(Km+si)) # metabolic enzyme kinetics

    # RHS (dydt)
    dydt = np.zeros(14)
    # first reaction - dilution
    dydt[0] = dydt[0] - (r*lam)
    dydt[1] = dydt[1] - (et*lam)
    dydt[2] = dydt[2] - (em*lam)
    dydt[3] = dydt[3] - (q*lam)
    dydt[4] = dydt[4] - (si*lam)
    dydt[5] = dydt[5] - (a*lam)
    
    # second reaction - transcription
    dydt[4] = dydt[4] + (s*(vimp)) 
    dydt[6] = dydt[6] + ((wr*a)/thetax+a)
    dydt[7] = dydt[7] + ((we*a)/thetax+a)
    dydt[8] = dydt[8] + ((we*a)/thetax+a)
    dydt[9] = dydt[9] + ((wq*a)/thetax+a)*(1/1+((q/Kq)^hq)) 
    
    # third reaction - dilution/degradation
    dydt[6] = dydt[6] - (mr*(lam+dm))
    dydt[7] = dydt[7] - (mt*(lam+dm))
    dydt[8] = dydt[8] - (mm*(lam+dm))
    dydt[9] = dydt[9] - (mq*(lam+dm))
    dydt[4] = dydt[4] - (si*nucat) 
        
    # fourth reaction - ribosome binding 
    dydt[0] = dydt[0] - (r*mr*kb) - (r*mt*kb) - (r*mm*kb) - (r*mq*kb)
    dydt[6] = dydt[6] - (r*mr*kb)
    dydt[7] = dydt[7] - (r*mt*kb)
    dydt[8] = dydt[8] - (r*mm*kb)
    dydt[9] = dydt[9] - (r*mq*kb)
    dydt[10] = dydt[10] + (r*mr*kb)
    dydt[11] = dydt[11] + (r*mt*kb)
    dydt[12] = dydt[12] + (r*mm*kb)
    dydt[13] = dydt[13] + (r*mq*kb)
    
    # fifth reaction - ribosome unbinding 
    dydt[0] = dydt[0] + (cr*ku) + (ct*ku) + (cm*ku) + (cq*ku)
    dydt[6] = dydt[6] + (cr*ku)
    dydt[7] = dydt[7] + (ct*ku)
    dydt[8] = dydt[8] + (cm*ku)
    dydt[9] = dydt[9] + (cq*ku)
    dydt[10] = dydt[10] - (cr*ku)
    dydt[11] = dydt[11] - (ct*ku)
    dydt[12] = dydt[12] - (cm*ku)
    dydt[13] = dydt[13] - (cq*ku)
    
    # sixth reaction - dilution
    dydt[10] = dydt[10] - (cr*lam)
    dydt[11] = dydt[11] - (ct*lam)
    dydt[12] = dydt[12] - (cm*lam)
    dydt[13] = dydt[13] - (cq*lam)
    
    # seventh reaction - translation 
    dydt[0] = dydt[0] + (cr*(gamma/nr)) + (ct*(gamma/nx)) + (cm*(gamma/nx)) + (cq*(gamma/nx)) + (cr*(gamma/nr))
    dydt[6] = dydt[6] + (cr*(gamma/nr))
    dydt[7] = dydt[7] + (ct*(gamma/nx))
    dydt[8] = dydt[8] + (cm*(gamma/nx))
    dydt[9] = dydt[9] + (cq*(gamma/nx))
    dydt[1] = dydt[1] + (ct*(gamma/nx))
    dydt[2] = dydt[2] + (cm*(gamma/nx))
    dydt[3] = dydt[3] + (cq*(gamma/nx))
    dydt[10] = dydt[10] - (cr*(gamma/nr))
    dydt[11] = dydt[11] - (ct*(gamma/nx))
    dydt[12] = dydt[12] - (cm*(gamma/nx))
    dydt[13] = dydt[13] - (cq*(gamma/nx))
    dydt[5] = dydt[5] - (cr*(gamma/nr)) - (ct*(gamma/nx)) - (cm*(gamma/nx)) - (cq*(gamma/nx))

    return dydt


def rmodel(t, variables, s, d_m, n_s, n_r, n_t, n_m, n_q, gamma_max, K_gamma, v_t, K_t, v_m, K_m, w_r, w_e, w_t, w_m, \
    w_q, theta_r, theta_nr, K_q, h_q, k_b, k_u, M, k_cm):

    s_i, a, r, e_t, e_m, q, m_r, m_t, m_m, m_q, c_r, c_t, c_m, c_q =  variables

    # helper function
    def v_imp(e_t, s):
        return e_t * v_t * s / (K_t + s)


    def v_cat(e_m, s_i):
        return e_m * v_m * s_i / (K_m + s_i)


    def gamma(a):
        return gamma_max * a / (K_gamma + a)


    def v_x(x, c_x, a):
        if x == 'r':
            return c_x * gamma(a) / n_r
        elif x == 't':
            return c_x * gamma(a) / n_t
        elif x == 'm':
            return c_x * gamma(a) / n_m
        else:
            return c_x * gamma(a) / n_q

    def omega(x, a):
        if x == 'r':
            return w_r * a / (theta_r + a)
        elif x == 't':
            return w_t * a / (theta_nr + a)
        else:
            return w_m * a / (theta_nr + a)


    def iota(q):
        return 1 / (1 + (q + K_q) ** h_q)


    def omega_q(q, a):
        return w_q * a * iota(q) / (theta_nr + a)

    R_t = c_r + c_t + c_m + c_q
    growth_rate = gamma(a) * R_t / M

    dsidt = v_imp(e_t, s) - v_cat(e_m, s_i) - growth_rate * s_i
    dadt = n_s * v_cat(e_m, s_i) - n_r * v_x('r', c_r, a) - n_t * v_x('t', c_t, a) - n_m * v_x('m', c_m, a) - n_q * v_x('q', c_m, a) - growth_rate * a
    drdt = v_x('r', c_r, a) - growth_rate * r + v_x('r', c_r, a) - k_b * r * m_r + k_u * c_r + v_x('t', c_t, a) - k_b * r * m_t + k_u * c_t + v_x('m', c_m, a) - k_b * r * m_m + k_u * c_m + v_x('q', c_q, a) - k_b * r * m_q + k_u * c_q
    detdt = v_x('t', c_t, a) - growth_rate * e_t
    demdt = v_x('m', c_m, a) - growth_rate * e_m
    dqdt = v_x('q', c_q, a) - growth_rate * q
    dmrdt = omega('r', a) - (growth_rate + d_m) * m_r - k_b * r * m_r + k_u * c_r + v_x('r', c_r, a)
    dmtdt = omega('t', a) - (growth_rate + d_m) * m_t - k_b * r * m_t + k_u * c_t + v_x('t', c_t, a)
    dmmdt = omega('m', a) - (growth_rate + d_m) * m_m - k_b * r * m_m + k_u * c_m + v_x('m', c_t, a)
    dmqdt = omega_q(q, a) - (growth_rate + d_m) * m_q - k_b * r * m_q + k_u * c_q + v_x('q', c_t, a)
    dcrdt = k_b * m_r * r - k_u * c_r - v_x('r', c_r, a) - growth_rate * c_r
    dctdt = k_b * m_t * r - k_u * c_t - v_x('t', c_t, a) - growth_rate * c_t
    dcmdt = k_b * m_m * r - k_u * c_m - v_x('m', c_m, a) - growth_rate * c_m
    dcqdt = k_b * m_q * r - k_u * c_q - v_x('q', c_q, a) - growth_rate * c_q

    return [dsidt, dadt, drdt, detdt, demdt, dqdt, dmrdt, dmtdt, dmmdt, dmqdt, dcrdt, dctdt, dcmdt, dcqdt]