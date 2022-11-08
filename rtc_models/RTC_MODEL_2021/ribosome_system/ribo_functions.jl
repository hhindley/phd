function ribo_odes!(dydt, initial, params, t)

    k1, k2, k3, ka, kb, kc, k, atp = params 

    rtca, rtcb, rdrtca, rtrtcb, ribo_d, ribo_t, ribo_h = initial 

    drtcadt, drtcbdt, drdrtcadt, drtrtcbdt, dribo_ddt, dribo_tdt, dribo_hdt = zeros(length(dydt))
    
    dydt[1] = -k1*rtca*ribo_d + k2*rdrtca + k3*atp*rdrtca
    dydt[2] = -ka*rtcb*ribo_t + kb*rtrtcb + kc*atp*rtrtcb
    dydt[3] = k1*rtca*ribo_d - k2*rdrtca - k3*atp*rdrtca
    dydt[4] = ka*rtcb*ribo_t - kb*rtrtcb - kc*atp*rtrtcb
    dydt[5] = k*ribo_h - k1*rtca*ribo_d + k2*rdrtca
    dydt[6] = k3*atp*rdrtca - ka*rtcb*ribo_t + kb*rtrtcb
    dydt[7] = kc*atp*rtrtcb - k*ribo_h 
    
end

function ribo_odes_alg!(dydt, initial, params, t)  # using mass conservation to get rid of the intermediate equations - doesn't work in full model because mass conservation no longer stands 

    k1, k2, k3, ka, kb, kc, k, atp, ribo_drtca_0, rtrtcb_0, rtot = params 

    rtca, rtcb, ribo_d, ribo_h = initial 

    drtcadt, drtcbdt, dribo_ddt, dribo_hdt = zeros(length(dydt))

    rdrtca = rtca_0 + ribo_drtca_0 - rtca
    rtrtcb = rtcb_0 + rtrtcb_0 - rtcb
    
    ribo_t = rtot - ribo_d - ribo_h - rdrtca - rtrtcb

    dydt[1] = -k1*rtca*ribo_d + k2*rdrtca + k3*atp*rdrtca
    dydt[2] = -ka*rtcb*ribo_t + kb*rtrtcb + kc*atp*rtrtcb
    dydt[3] = k*ribo_h - k1*rtca*ribo_d + k2*rdrtca
    dydt[4] = kc*atp*rtrtcb - k*ribo_h 
    
end

function ribo_all_alg!(dydt, initial, params, t)  

    k1, k2, k3, ka, kb, kc, k, atp, rtca_tot, rtcb_tot, ribo_tot = params 

    rtca, rtcb, ribo_d, ribo_t = initial 

    drtcadt, drtcbdt, dribo_ddt, dribo_tdt = zeros(length(dydt))

    rdrtca = k1*ribo_d*rtca_tot/(k1*ribo_d+k2+k3*atp)
    rtrtcb = ka*ribo_t*rtcb_tot/(ka*ribo_t+kb+kc*atp)
    
    ribo_h = ribo_tot - ribo_d - ribo_t - rdrtca - rtrtcb
    
    dydt[1] = -k1*rtca*ribo_d + k2*rdrtca + k3*atp*rdrtca
    dydt[2] = -ka*rtcb*ribo_t + kb*rtrtcb + kc*atp*rtrtcb
    dydt[3] = k*ribo_h - k1*rtca*ribo_d + k2*rdrtca
    dydt[4] = k3*atp*rdrtca - ka*rtcb*ribo_t + kb*rtrtcb
end

function ribo_mix1!(dydt, initial, params, t)  

    k1, k2, k3, ka, kb, kc, k, atp, rtca_tot, rtcb_tot = params 

    rtca, rtcb, ribo_d, ribo_t, ribo_h = initial 

    drtcadt, drtcbdt, dribo_ddt, dribo_tdt, dribo_hdt = zeros(length(dydt))

    rdrtca = k1*ribo_d*rtca_tot/(k1*ribo_d+k2+k3*atp)
    rtrtcb = ka*ribo_t*rtcb_tot/(ka*ribo_t+kb+kc*atp)
        
    dydt[1] = -k1*rtca*ribo_d + k2*rdrtca + k3*atp*rdrtca
    dydt[2] = -ka*rtcb*ribo_t + kb*rtrtcb + kc*atp*rtrtcb
    dydt[3] = k*ribo_h - k1*rtca*ribo_d + k2*rdrtca
    dydt[4] = k3*atp*rdrtca - ka*rtcb*ribo_t + kb*rtrtcb
    dydt[5] = kc*atp*rtrtcb - k*ribo_h 
end

function ribo_mix2!(dydt, initial, params, t)  

    k1, k2, k3, ka, kb, kc, k, atp, ribo_tot = params 

    rtca, rtcb, ribo_d, ribo_t, rdrtca, rtrtcb = initial 

    drtcadt, drtcbdt, dribo_ddt, dribo_tdt, drdrtcadt, drdrtcbdt = zeros(length(dydt))

    ribo_h = ribo_tot - ribo_d - ribo_t - rdrtca - rtrtcb
        
    dydt[1] = -k1*rtca*ribo_d + k2*rdrtca + k3*atp*rdrtca
    dydt[2] = -ka*rtcb*ribo_t + kb*rtrtcb + kc*atp*rtrtcb
    dydt[3] = k*ribo_h - k1*rtca*ribo_d + k2*rdrtca
    dydt[4] = k3*atp*rdrtca - ka*rtcb*ribo_t + kb*rtrtcb
    dydt[5] = k1*rtca*ribo_d - k2*rdrtca - k3*atp*rdrtca
    dydt[6] = ka*rtcb*ribo_t - kb*rtrtcb - kc*atp*rtrtcb    
end

function split_all_species(species, solution)
    solDF = DataFrame([[j[i] for j in solution.u] for i=1:length(solution.u[1])], species);
    rtca = solDF[:, :rtca]
    rtcb = solDF[:, :rtcb]
    rdrtca = solDF[:, :rdrtca]
    rtrtcb = solDF[:, :rtrtcb]
    ribo_d = solDF[:, :ribo_d]
    ribo_t = solDF[:, :ribo_t]
    ribo_h = solDF[:, :ribo_h]
    return rtca, rtcb, rdrtca, rtrtcb, ribo_d, ribo_t, ribo_h
end

function split_fewer_species(species, solution)
    solDF = DataFrame([[j[i] for j in solution.u] for i=1:length(solution.u[1])], species);
    rtca = solDF[:, :rtca]
    rtcb = solDF[:, :rtcb]
    ribo_d = solDF[:, :ribo_d]
    ribo_h = solDF[:, :ribo_h]
    ribo_t = rtot .- ribo_d .-ribo_h 
    rdrtca = rtca_0 + rdrtca_0 .- rtca
    rtrtcb = rtcb_0 + rtrtcb_0 .- rtcb
    return rtca, rtcb, rdrtca, rtrtcb, ribo_d, ribo_t, ribo_h
end

function split_all_alg(species, solution)
    solDF = DataFrame([[j[i] for j in solution.u] for i=1:length(solution.u[1])], species);
    rtca = solDF[:, :rtca]
    rtcb = solDF[:, :rtcb]
    ribo_d = solDF[:, :ribo_d]
    ribo_t = solDF[:, :ribo_t]

    rdrtca = k1.*ribo_d.*rtca_tot./(k1.*ribo_d.+k2.+k3.*atp)
    rtrtcb = ka.*ribo_t.*rtcb_tot./(ka.*ribo_t.+kb.+kc.*atp)
    
    ribo_h = ribo_tot .- ribo_d .- ribo_t .- rdrtca .- rtrtcb
    return rtca, rtcb, rdrtca, rtrtcb, ribo_d, ribo_t, ribo_h
end 

function split_mix1(species, solution)
    solDF = DataFrame([[j[i] for j in solution.u] for i=1:length(solution.u[1])], species);
    rtca = solDF[:, :rtca]
    rtcb = solDF[:, :rtcb]
    ribo_d = solDF[:, :ribo_d]
    ribo_t = solDF[:, :ribo_t]
    ribo_h = solDF[:, :ribo_h]

    rdrtca = k1.*ribo_d.*rtca_tot./(k1.*ribo_d.+k2.+k3.*atp)
    rtrtcb = ka.*ribo_t.*rtcb_tot./(ka.*ribo_t.+kb.+kc.*atp)

    return rtca, rtcb, rdrtca, rtrtcb, ribo_d, ribo_t, ribo_h
end 

function split_mix2(species, solution)
    solDF = DataFrame([[j[i] for j in solution.u] for i=1:length(solution.u[1])], species);
    rtca = solDF[:, :rtca]
    rtcb = solDF[:, :rtcb]
    ribo_d = solDF[:, :ribo_d]
    ribo_t = solDF[:, :ribo_t]
    rdrtca = solDF[:, :rdrtca]
    rtrtcb = solDF[:, :rtrtcb]

    ribo_h = ribo_tot .- ribo_d .- ribo_t .- rdrtca .- rtrtcb

    return rtca, rtcb, rdrtca, rtrtcb, ribo_d, ribo_t, ribo_h
end 

function split(species, solution)
    if length(species) == 7
        rtca, rtcb, ribo_drtca, rtrtcb, ribo_d, rt, ribo_h = split_all_species(species, solution)
    else
        rtca, rtcb, ribo_drtca, rtrtcb, ribo_d, rt, ribo_h = split_fewer_species(species, solution)
    end
end