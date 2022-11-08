function simple_solve!(model, init, tspan, params)
    prob = ODEProblem(model, init, tspan, params);
    sol = solve(prob, Rodas4())#, abstol=1e-15, reltol=1e-12);
    return sol
end



function rtc_system!(dydt, initial, params, t)

    kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, rdrtca_0, rtrtcb_0, ribo_t_0 #=ribo_tot=# = params

    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, ribo_d, ribo_h = initial 
    drm_adt, drtcadt, drm_bdt, drtcbdt, drm_rdt, drtcrdt, dribo_ddt, dribo_hdt = zeros(length(dydt))

    rtca_tot = rtca + rdrtca_0
    rdrtca = k1_a*ribo_d*rtca_tot/(k1_a*ribo_d+k2_a+k3_a*atp)
    rtcb_tot = rtcb + rtrtcb_0
    rtrtcb = ka_b*ribo_t*rtcb_tot/(ka_b*ribo_t+kb_b+kc_b*atp)
    n = 6
    ribo_tot = ribo_t_0 + ribo_d_0 + ribo_h_0 - rdrtca - rtrtcb #10
    gr = 0.01*ribo_h
    ribo_t = ribo_tot - ribo_d - ribo_h
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


end



function rtc_system_ribos_rd!(dydt, initial, params, t)

    kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, gr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, rdrtca_0, rtrtcb_0 = params

    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, ribo_h, ribo_t = initial 
    drm_adt, drtcadt, drm_bdt, drtcbdt, drm_rdt, drtcrdt, dribo_hdt, dribo_tdt = zeros(length(dydt))

    n = 6
    # ribo_tot = ribo_t_0 + ribo_d_0 + ribo_h_0 #10
    # gr = 0.01*ribo_tot
    ribo_d = ribo_tot - ribo_t - ribo_h
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

    rtca_tot = rtca + rdrtca_0
    rdrtca = k1_a*ribo_d*rtca_tot/(k1_a*ribo_d+k2_a+k3_a*atp)
    rtcb_tot = rtcb + rtrtcb_0
    rtrtcb = ka_b*ribo_t*rtcb_tot/(ka_b*ribo_t+kb_b+kc_b*atp)

    dydt[1] = tscr - gr*rm_a - d*rm_a
    dydt[2] = tlr - k1_a*rtca*ribo_d + k2_a*rdrtca + k3_a*atp*rdrtca  - gr*rtca 
    dydt[3] = tscr_b - gr*rm_b - d*rm_b 
    dydt[4] = tlr_b - ka_b*rtcb*ribo_t + kb_b*rtrtcb + kc_b*atp*rtrtcb - gr*rtcb
    dydt[5] = rtcr_tscr - gr*rm_r - d*rm_r
    dydt[6] = rtcr_tlr - gr*rtcr

    dydt[7] = kc_b*atp*rtrtcb - k*ribo_h 
    dydt[8] = k3_a*atp*rdrtca - ka_b*rtcb*ribo_t + kb_b*rtrtcb


end




function split_8(species, solution)

    solDF = DataFrame([[j[i] for j in solution.u] for i=1:length(solution.u[1])], species)

    rm_a = solDF[:, :rm_a]
    rtca = solDF[:, :rtca]
    rm_b = solDF[:, :rm_b]
    rtcb = solDF[:, :rtcb]
    rm_r = solDF[:, :rm_r]
    rtcr = solDF[:, :rtcr]
    ribo_d = solDF[:, :ribo_d]
    ribo_h = solDF[:, :ribo_h]

    ribo_tot = ribo_t_0 + ribo_d_0 + ribo_h_0 #10
    ribo_t = ribo_tot .- ribo_d .- ribo_h 
    alpha = ribo_t./kr
    fa_1 = (alpha.+1).^6 # fraction of active RtcR
    fa_2 = L*(c.*alpha.+1).^6 +(alpha.+1).^6
    fa = fa_1./fa_2
    ra = fa.*rtcr # amount of active RtcR
    v = ra*k1*sigma*k3*atp*k5/(k1*sigma*k3*atp+k1*sigma*(k4+k5)+k2*(k4+k5)+k3*atp*k5) # rate of open complex formation
    tscr = v*(w_rtc*atp/(theta_rtc+atp)) # transcription rate 
    tlel = max*atp/(thr+atp) # translation elongation rate
    tlr = ribo_h.*rm_a.*tlel # translation rate so formation of Rtc 

    tscr_b = v*(w_rtc*atp/(theta_rtc+atp)) # transcription rate 
    tlel_b = max*atp/(thr+atp) # translation elongation rate
    tlr_b = ribo_h.*rm_b.*tlel_b

    # rtcr 
    rtcr_tscr = w_rtc*atp/(theta_rtc+atp)
    rtcr_tlel = max*atp/(thr+atp)
    rtcr_tlr = ribo_h.*rm_r.*rtcr_tlel

    rtca_tot = rtca .+ rdrtca_0
    rdrtca = k1_a.*ribo_d.*rtca_tot./(k1_a.*ribo_d.+k2_a.+k3_a.*atp)
    rtcb_tot = rtcb .+ rtrtcb_0
    rtrtcb = ka_b.*ribo_t.*rtcb_tot./(ka_b.*ribo_t.+kb_b.+kc_b.*atp)

    return rm_a, rtca, rm_b, rtcb, rm_r, rtcr, ribo_d, ribo_h, ribo_t, alpha, fa, ra, v, tscr, tlr, tscr_b, tlr_b, rtcr_tscr, rtcr_tlr, rdrtca, rtrtcb
end

function split_rh(species, solution) # 2

    solDF = DataFrame([[j[i] for j in solution.u] for i=1:length(solution.u[1])], species)

    rm_a = solDF[:, :rm_a]
    rtca = solDF[:, :rtca]
    rm_b = solDF[:, :rm_b]
    rtcb = solDF[:, :rtcb]
    rm_r = solDF[:, :rm_r]
    rtcr = solDF[:, :rtcr]
    ribo_d = solDF[:, :ribo_d]
    ribo_t = solDF[:, :ribo_t]

    n = 6
    rdrtca = k1_a.*ribo_d.*rtca./(k1_a.*ribo_d.+k2_a.+k3_a.*atp)
    rtrtcb = ka_b.*ribo_t.*rtcb./(ka_b.*ribo_t.+kb_b.+kc_b.*atp)
   
    ribo_tot = ribo_t_0 + ribo_d_0 + ribo_h_0 + rdrtca_0 + rtrtcb_0 #10
    ribo_h = ribo_tot .- ribo_d .- ribo_t .- rdrtca .- rtrtcb

    alpha = ribo_t./kr
    fa = ([1].+alpha).^n./(L.*(([1].+c.*alpha).^n)+([1].+alpha).^n) # fraction of active RtcR
    ra = fa.*rtcr # amount of active RtcR
    v = ra.*k1.*sigma.*k3.*atp.*k5./(k1.*sigma.*k3.*atp.+k1.*sigma.*(k4.+k5).+k2.*(k4.+k5).+k3.*atp.*k5) # rate of open complex formation
    
    # rtca
    tscr = v.*(w_rtc*atp/(theta_rtc+atp)) # transcription rate 
    tlel = max*atp/(thr+atp) # translation elongation rate
    tlr = ribo_h.*rm_a.*tlel # translation rate so formation of Rtc 

    # rtcb
    tscr_b = v.*(w_rtc*atp/(theta_rtc+atp)) # transcription rate 
    tlel_b = max*atp/(thr+atp) # translation elongation rate
    tlr_b = ribo_h.*rm_b.*tlel_b # translation rate so formation of Rtc 

    # rtcr 
    rtcr_tscr = w_rtc*atp/(theta_rtc+atp)
    rtcr_tlel = max*atp/(thr+atp)
    rtcr_tlr = ribo_h.*rm_r.*rtcr_tlel



    return rm_a, rtca, rm_b, rtcb, rm_r, rtcr, ribo_d, ribo_h, ribo_t, alpha, fa, ra, v, tscr, tlr, tscr_b, tlr_b, rtcr_tscr, rtcr_tlr, rdrtca, rtrtcb
end

function split_9(species, solution) # 1

    solDF = DataFrame([[j[i] for j in solution.u] for i=1:length(solution.u[1])], species)

    rm_a = solDF[:, :rm_a]
    rtca = solDF[:, :rtca]
    rm_b = solDF[:, :rm_b]
    rtcb = solDF[:, :rtcb]
    rm_r = solDF[:, :rm_r]
    rtcr = solDF[:, :rtcr]
    ribo_d = solDF[:, :ribo_d]
    ribo_h = solDF[:, :ribo_h]
    ribo_t = solDF[:, :ribo_t]

    n = 6

    rtca_tot = rtca_0 + rdrtca_0
    rdrtca = k1_a.*ribo_d.*rtca_tot./(k1_a.*ribo_d.+k2_a.+k3_a.*atp)
    rtcb_tot = rtcb_0 + rtrtcb_0
    rtrtcb = ka_b.*ribo_t.*rtcb_tot./(ka_b.*ribo_t.+kb_b.+kc_b.*atp)

    alpha = ribo_t./kr
    fa = ([1].+alpha).^n./(L.*(([1].+c.*alpha).^n).+([1].+alpha).^n) # fraction of active RtcR
    ra = fa.*rtcr # amount of active RtcR
    v = ra.*k1.*sigma.*k3.*atp.*k5./(k1.*sigma.*k3.*atp.+k1.*sigma.*(k4.+k5).+k2.*(k4.+k5).+k3.*atp.*k5) # rate of open complex formation
    
    # rtca
    tscr = v.*(w_rtc*atp/(theta_rtc+atp)) # transcription rate 
    tlel = max*atp/(thr+atp) # translation elongation rate
    tlr = ribo_h.*rm_a.*tlel # translation rate so formation of Rtc 

    # rtcb
    tscr_b = v.*(w_rtc*atp/(theta_rtc+atp)) # transcription rate 
    tlel_b = max*atp/(thr+atp) # translation elongation rate
    tlr_b = ribo_h.*rm_b.*tlel_b # translation rate so formation of Rtc 

    # rtcr 
    rtcr_tscr = w_rtc*atp/(theta_rtc+atp)
    rtcr_tlel = max*atp/(thr+atp)
    rtcr_tlr = ribo_h.*rm_r.*rtcr_tlel

    return rm_a, rtca, rm_b, rtcb, rm_r, rtcr, ribo_d, ribo_h, ribo_t, alpha, fa, ra, v, tscr, tlr, tscr_b, tlr_b, rtcr_tscr, rtcr_tlr, rdrtca, rtrtcb
end

function split_all(species, solution) # 4

    solDF = DataFrame([[j[i] for j in solution.u] for i=1:length(solution.u[1])], species)

    rm_a = solDF[:, :rm_a]
    rtca = solDF[:, :rtca]
    rm_b = solDF[:, :rm_b]
    rtcb = solDF[:, :rtcb]
    rm_r = solDF[:, :rm_r]
    rtcr = solDF[:, :rtcr]
    ribo_d = solDF[:, :ribo_d]
    ribo_t = solDF[:, :ribo_t]
    rdrtca = solDF[:, :rdrtca]
    rtrtcb = solDF[:, :rtrtcb]

    n = 6
    ribo_tot = ribo_t_0+ribo_d_0+ribo_h_0+rdrtca_0+rtrtcb_0
    ribo_h = ribo_tot .- ribo_d .- ribo_t .- rdrtca .- rtrtcb

    alpha = ribo_t./kr
    fa = ([1].+alpha).^n./(L.*(([1].+c.*alpha).^n).+([1].+alpha).^n) # fraction of active RtcR
    ra = fa.*rtcr # amount of active RtcR
    v = ra.*k1.*sigma.*k3.*atp.*k5./(k1.*sigma.*k3.*atp.+k1.*sigma.*(k4.+k5).+k2.*(k4.+k5).+k3.*atp.*k5) # rate of open complex formation
    
    # rtca
    tscr = v.*(w_rtc*atp/(theta_rtc+atp)) # transcription rate 
    tlel = max*atp/(thr+atp) # translation elongation rate
    tlr = ribo_h.*rm_a.*tlel # translation rate so formation of Rtc 

    # rtcb
    tscr_b = v.*(w_rtc*atp/(theta_rtc+atp)) # transcription rate 
    tlel_b = max*atp/(thr+atp) # translation elongation rate
    tlr_b = ribo_h.*rm_b.*tlel_b # translation rate so formation of Rtc 

    # rtcr 
    rtcr_tscr = w_rtc*atp/(theta_rtc+atp)
    rtcr_tlel = max*atp/(thr+atp)
    rtcr_tlr = ribo_h.*rm_r.*rtcr_tlel
    return rm_a, rtca, rm_b, rtcb, rm_r, rtcr, ribo_d, ribo_h, ribo_t, alpha, fa, ra, v, tscr, tlr, tscr_b, tlr_b, rtcr_tscr, rtcr_tlr, rdrtca, rtrtcb
end




function split(species, solution)
    if length(species) == 8
        rm_a, rtca, rm_b, rtcb, rm_r, rtcr, ribo_d, ribo_h, ribo_t, alpha, fa, ra, v, tscr, tlr, tscr_b, tlr_b, rtcr_tscr, rtcr_tlr, rdrtca, rtrtcb = split_8(species, solution)
    elseif  length(species) == 9
        rm_a, rtca, rm_b, rtcb, rm_r, rtcr, ribo_d, ribo_h, ribo_t, alpha, fa, ra, v, tscr, tlr, tscr_b, tlr_b, rtcr_tscr, rtcr_tlr, rdrtca, rtrtcb = split_9(species, solution)
    else
        rm_a, rtca, rm_b, rtcb, rm_r, rtcr, ribo_d, ribo_h, ribo_t, alpha, fa, ra, v, tscr, tlr, tscr_b, tlr_b, rtcr_tscr, rtcr_tlr, rdrtca, rtrtcb = split_all(species, solution)
    end
end


function vary_parameter_for_full_model(varying, param, range, model, init, tspan, params, species)
    rm_a_res = []
    rtca_res = []
    rm_b_res = []
    rtcb_res = []
    rm_r_res = []
    rtcr_res = []
    ribo_d_res = []
    ribo_h_res = []
    ribo_t_res = []
    rdrtca_res = []
    rtrtcb_res = []
    solutions = [];
    solutions_t = [];
    for x in range
        varying[param] = x
        sol = simple_solve!(model, init, tspan, params)
        rm_a, rtca, rm_b, rtcb, rm_r, rtcr, ribo_d, ribo_h, ribo_t, alpha, fa, ra, v, tscr, tlr, tscr_b, tlr_b, rtcr_tscr, rtcr_tlr, rdrtca, rtrtcb = split(species, sol)
        push!(rm_a_res, rm_a)
        push!(rtca_res, rtca)
        push!(rm_b_res, rm_b)
        push!(rtcb_res, rtcb)
        push!(rm_r_res, rm_r)
        push!(rtcr_res, rtcr)
        push!(ribo_d_res, ribo_d)
        push!(ribo_h_res, ribo_h)
        push!(ribo_t_res, ribo_t)
        push!(rdrtca_res, rdrtca)
        push!(rtrtcb_res, rtrtcb)
        push!(solutions, sol)
        push!(solutions_t, sol.t)
    end
    sol1 = solutions[1]
    sol2 = solutions[2]
    sol3 = solutions[3]
    sol4 = solutions[4]
    sol5 = solutions[5]
    sol6 = solutions[6]
    return rm_a_res, rtca_res, rm_b_res, rtcb_res, rm_r_res, rtcr_res, ribo_d_res, ribo_h_res, ribo_t_res, rdrtca_res, rtrtcb_res, solutions, solutions_t, sol1, sol2, sol3, sol4, sol5, sol6
end

function plots_for_param(v_param, varying, param, range, model, init, tspan, params, species) # varying should be init or param
    rm_a_res, rtca_res, rm_b_res, rtcb_res, rm_r_res, rtcr_res, ribo_d_res, ribo_h_res, ribo_t_res, rdrtca_res, rtrtcb_res, solutions, solutions_t, sol1, sol2, sol3, sol4, sol5, sol6 = vary_parameter_for_full_model(varying, param, range, model, init, tspan, params, species);

    p1 = plot(sol1, labels=labels, palette=:tab10, title="$v_param = $(collect(range)[1])")
    p2 = plot(sol2, labels=labels, palette=:tab10, title="$v_param = $(collect(range)[2])")
    p3 = plot(sol3, labels=labels, palette=:tab10, title="$v_param = $(collect(range)[3])")
    p4 = plot(sol4, labels=labels, palette=:tab10, title="$v_param = $(collect(range)[4])")
    p5 = plot(sol5, labels=labels, palette=:tab10, title="$v_param = $(collect(range)[5])")
    p6 = plot(sol6, labels=labels, palette=:tab10, title="$v_param = $(collect(range)[6])")

    p_main = plot(p1, p2, p3, p4, p5, p6, size=(1000,800))

    labels1 = ["$(collect(range)[1])" "$(collect(range)[2])" "$(collect(range)[3])" "$(collect(range)[4])" "$(collect(range)[5])" "$(collect(range)[6])"]
        
    p1_a = plot(solutions_t, rtcr_res, title="RtcR", xlabel="t", palette=:tab10, labels=labels1)
    p2_a = plot(solutions_t, rtca_res, title="RtcA", xlabel="t", palette=:tab10, labels=labels1)
    p3_a = plot(solutions_t, rtcb_res, title="RtcB", xlabel="t", palette=:tab10, labels=labels1)
    p4_a = plot(solutions_t, ribo_d_res, title="Rd", xlabel="t", palette=:tab10, labels=labels1)
    p5_a = plot(solutions_t, ribo_t_res, title="Rt", xlabel="t", palette=:tab10, labels=labels1)
    p6_a = plot(solutions_t, ribo_h_res, title="Rh", xlabel="t", palette=:tab10, labels=labels1)

    p_a_main = plot(p1_a, p2_a, p3_a, p4_a, p5_a, p6_a, size=(1000,800), plot_title="Varying $v_param")
    return p_main, p_a_main
end

function logplots_for_param(v_param, varying, param, range, model, init, tspan, params, species, labels) # varying should be init or param
    rm_a_res, rtca_res, rm_b_res, rtcb_res, rm_r_res, rtcr_res, ribo_d_res, ribo_h_res, ribo_t_res, rdrtca_res, rtrtcb_res, solutions, solutions_t, sol1, sol2, sol3, sol4, sol5, sol6 = vary_parameter_for_full_model(varying, param, range, model, init, tspan, params, species);

    p1 = plot(sol1[2:end], labels=labels, palette=:tab10, title="$v_param = $(collect(range)[1])",  xaxis=(:log10, (1,Inf)), yaxis=(:log10, (1,Inf)))
    p2 = plot(sol2[2:end], labels=labels, palette=:tab10, title="$v_param = $(collect(range)[2])",  xaxis=(:log10, (1,Inf)), yaxis=(:log10, (1,Inf)))
    p3 = plot(sol3[2:end], labels=labels, palette=:tab10, title="$v_param = $(collect(range)[3])", xaxis=(:log10, (1,Inf)), yaxis=(:log10, (1,Inf)))
    p4 = plot(sol4[2:end], labels=labels, palette=:tab10, title="$v_param = $(collect(range)[4])", xaxis=(:log10, (1,Inf)), yaxis=(:log10, (1,Inf)))
    p5 = plot(sol5[2:end], labels=labels, palette=:tab10, title="$v_param = $(collect(range)[5])", xaxis=(:log10, (1,Inf)), yaxis=(:log10, (1,Inf)))
    p6 = plot(sol6[2:end], labels=labels, palette=:tab10, title="$v_param = $(collect(range)[6])", xaxis=(:log10, (1,Inf)), yaxis=(:log10, (1,Inf)))

    p_main = plot(p1, p2, p3, p4, p5, p6, size=(1000,800))

    labels1 = ["$(collect(range)[1])" "$(collect(range)[2])" "$(collect(range)[3])" "$(collect(range)[4])" "$(collect(range)[5])" "$(collect(range)[6])"]
        
    p1_a = plot(solutions_t, rtcr_res, title="RtcR", xlabel="t", palette=:tab10, labels=labels1, xaxis=(:log10, (1,Inf)))
    p2_a = plot(solutions_t, rtca_res, title="RtcA", xlabel="t", palette=:tab10, labels=labels1, xaxis=(:log10, (1,Inf)))
    p3_a = plot(solutions_t, rtcb_res, title="RtcB", xlabel="t", palette=:tab10, labels=labels1, xaxis=(:log10, (1,Inf)))
    p4_a = plot(solutions_t, ribo_d_res, title="Rd", xlabel="t", palette=:tab10, labels=labels1, xaxis=(:log10, (1,Inf)))
    p5_a = plot(solutions_t, ribo_t_res, title="Rt", xlabel="t", palette=:tab10, labels=labels1, xaxis=(:log10, (1,Inf)))
    p6_a = plot(solutions_t, ribo_h_res, title="Rh", xlabel="t", palette=:tab10, labels=labels1, xaxis=(:log10, (1,Inf)))

    p_a_main = plot(p1_a, p2_a, p3_a, p4_a, p5_a, p6_a, size=(1000,800), plot_title="Varying $v_param")
    return p_main, p_a_main
end

function vary_param_xaxis(long_range, varying, param, species)
    rm_a_res = []
    rtca_res = []
    rm_b_res = []
    rtcb_res = []
    rm_r_res = []
    rtcr_res = []
    ribo_d_res = []
    ribo_h_res = []
    ribo_t_res = []
    rdrtca_res = []
    rtrtcb_res = []
    for i in long_range
        varying[param] = i
        sol1 = simple_solve!(rtc_system!, init, tspan, params)
        solDF = DataFrame([[j[i] for j in sol1.u] for i=1:length(sol1.u[1])], species)
        rm_a = solDF[end, :rm_a]
        rtca = solDF[end, :rtca]
        rm_b = solDF[end, :rm_b]
        rtcb = solDF[end, :rtcb]
        rm_r = solDF[end, :rm_r]
        rtcr = solDF[end, :rtcr]
        ribo_d = solDF[end, :ribo_d]
        ribo_h = solDF[end, :ribo_h]
        # rdrtca = solDF[end, :rdrtca]
        # rtrtcb = solDF[end, :rtrtcb]
        ribo_t = ribo_tot .- ribo_d .- ribo_h 
        push!(rm_a_res, rm_a)
        push!(rtca_res, rtca)
        push!(rm_b_res, rm_b)
        push!(rtcb_res, rtcb)
        push!(rm_r_res, rm_r)
        push!(rtcr_res, rtcr)
        push!(ribo_d_res, ribo_d)
        push!(ribo_h_res, ribo_h)
        push!(ribo_t_res, ribo_t)
        # push!(rdrtca_res, rdrtca)
        # push!(rtrtcb_res, rtrtcb)
    end
    return rm_a_res, rtca_res, rm_b_res, rtcb_res, rm_r_res, rtcr_res, ribo_d_res, ribo_h_res, ribo_t_res
end

function checking_difference_with_mm_and_odes()
    all = [rm_a, rtca, rm_b, rtcb, rm_r, rtcr, ribo_d, ribo_h, ribo_t, alpha, fa, ra, v, tscr, tlr, tscr_b, tlr_b, rtcr_tscr, rtcr_tlr, rdrtca, rtrtcb]

    all_max = []
    for i in all
        push!(all_max, (maximum(i)))
    end

    println(typeof(all_max)) 

    all_ode = [rm_a_ode, rtca_ode, rm_b_ode, rtcb_ode, rm_r_ode, rtcr_ode, ribo_d_ode, ribo_h_ode, ribo_t_ode, alpha_ode, fa_ode, ra_ode, v_ode, tscr_ode, tlr_ode, tscr_b_ode, tlr_b_ode, rtcr_tscr_ode, rtcr_tlr_ode, rdrtca_ode, rtrtcb_ode]
    all_max_ode = []
    for j in all_ode
        push!(all_max_ode, (maximum(j)))
    end

    println(all_max_ode)

    dif = all_max_ode.-all_max

    names = ["rm_a", "rtca", "rm_b", "rtcb", "rm_r", "rtcr", "ribo_d", "ribo_h", "ribo_t", "alpha", "fa", "ra", "v", "tscr", "tlr", "tscr_b", "tlr_b", "rtcr_tscr", "rtcr_tlr", "rdrtca", "rtrtcb"]

    df = DataFrame(species=names, all_max=all_max, all_max_ode=all_max_ode, difference=dif)
end










function vary_parameter_for_full_model_rh(varying, param, range, model, init, tspan, params, species)
    rm_a_res = []
    rtca_res = []
    rm_b_res = []
    rtcb_res = []
    rm_r_res = []
    rtcr_res = []
    ribo_d_res = []
    ribo_h_res = []
    ribo_t_res = []
    rdrtca_res = []
    rtrtcb_res = []
    solutions = [];
    solutions_t = [];
    for x in range
        varying[param] = x
        sol = simple_solve!(model, init, tspan, params)
        rm_a, rtca, rm_b, rtcb, rm_r, rtcr, ribo_d, ribo_h, ribo_t, alpha, fa, ra, v, tscr, tlr, tscr_b, tlr_b, rtcr_tscr, rtcr_tlr, rdrtca, rtrtcb = split_rh(species, sol)
        push!(rm_a_res, rm_a)
        push!(rtca_res, rtca)
        push!(rm_b_res, rm_b)
        push!(rtcb_res, rtcb)
        push!(rm_r_res, rm_r)
        push!(rtcr_res, rtcr)
        push!(ribo_d_res, ribo_d)
        push!(ribo_h_res, ribo_h)
        push!(ribo_t_res, ribo_t)
        push!(rdrtca_res, rdrtca)
        push!(rtrtcb_res, rtrtcb)
        push!(solutions, sol)
        push!(solutions_t, sol.t)
    end
    sol1 = solutions[1]
    sol2 = solutions[2]
    sol3 = solutions[3]
    sol4 = solutions[4]
    sol5 = solutions[5]
    sol6 = solutions[6]
    return rm_a_res, rtca_res, rm_b_res, rtcb_res, rm_r_res, rtcr_res, ribo_d_res, ribo_h_res, ribo_t_res, rdrtca_res, rtrtcb_res, solutions, solutions_t, sol1, sol2, sol3, sol4, sol5, sol6
end

function plots_for_param_rh(v_param, varying, param, range, model, init, tspan, params, species) # varying should be init or param
    rm_a_res, rtca_res, rm_b_res, rtcb_res, rm_r_res, rtcr_res, ribo_d_res, ribo_h_res, ribo_t_res, rdrtca_res, rtrtcb_res, solutions, solutions_t, sol1, sol2, sol3, sol4, sol5, sol6 = vary_parameter_for_full_model_rh(varying, param, range, model, init, tspan, params, species);

    p1 = plot(sol1, labels=labels, palette=:tab10, title="$v_param = $(collect(range)[1])")
    p2 = plot(sol2, labels=labels, palette=:tab10, title="$v_param = $(collect(range)[2])")
    p3 = plot(sol3, labels=labels, palette=:tab10, title="$v_param = $(collect(range)[3])")
    p4 = plot(sol4, labels=labels, palette=:tab10, title="$v_param = $(collect(range)[4])")
    p5 = plot(sol5, labels=labels, palette=:tab10, title="$v_param = $(collect(range)[5])")
    p6 = plot(sol6, labels=labels, palette=:tab10, title="$v_param = $(collect(range)[6])")

    p_main = plot(p1, p2, p3, p4, p5, p6, size=(1000,800))

    labels1 = ["$(collect(range)[1])" "$(collect(range)[2])" "$(collect(range)[3])" "$(collect(range)[4])" "$(collect(range)[5])" "$(collect(range)[6])"]
        
    p1_a = plot(solutions_t, rtcr_res, title="RtcR", xlabel="t", palette=:tab10, labels=labels1)
    p2_a = plot(solutions_t, rtca_res, title="RtcA", xlabel="t", palette=:tab10, labels=labels1)
    p3_a = plot(solutions_t, rtcb_res, title="RtcB", xlabel="t", palette=:tab10, labels=labels1)
    p4_a = plot(solutions_t, ribo_d_res, title="Rd", xlabel="t", palette=:tab10, labels=labels1)
    p5_a = plot(solutions_t, ribo_t_res, title="Rt", xlabel="t", palette=:tab10, labels=labels1)
    p6_a = plot(solutions_t, ribo_h_res, title="Rh", xlabel="t", palette=:tab10, labels=labels1)

    p_a_main = plot(p1_a, p2_a, p3_a, p4_a, p5_a, p6_a, size=(1000,800), plot_title="Varying $v_param")
    return p_main, p_a_main
end


function split_rd(species, solution)

    solDF = DataFrame([[j[i] for j in solution.u] for i=1:length(solution.u[1])], species)

    rm_a = solDF[:, :rm_a]
    rtca = solDF[:, :rtca]
    rm_b = solDF[:, :rm_b]
    rtcb = solDF[:, :rtcb]
    rm_r = solDF[:, :rm_r]
    rtcr = solDF[:, :rtcr]
    ribo_h = solDF[:, :ribo_h]
    ribo_t = solDF[:, :ribo_t]

    ribo_tot = ribo_t_0 + ribo_d_0 + ribo_h_0 #10
    ribo_d = ribo_tot .- ribo_h .- ribo_t 
    alpha = ribo_t./kr
    fa_1 = (alpha.+1).^6 # fraction of active RtcR
    fa_2 = L*(c.*alpha.+1).^6 +(alpha.+1).^6
    fa = fa_1./fa_2
    ra = fa.*rtcr # amount of active RtcR
    v = ra*k1*sigma*k3*atp*k5/(k1*sigma*k3*atp+k1*sigma*(k4+k5)+k2*(k4+k5)+k3*atp*k5) # rate of open complex formation
    tscr = v*(w_rtc*atp/(theta_rtc+atp)) # transcription rate 
    tlel = max*atp/(thr+atp) # translation elongation rate
    tlr = ribo_h.*rm_a.*tlel # translation rate so formation of Rtc 

    tscr_b = v*(w_rtc*atp/(theta_rtc+atp)) # transcription rate 
    tlel_b = max*atp/(thr+atp) # translation elongation rate
    tlr_b = ribo_h.*rm_b.*tlel_b

    # rtcr 
    rtcr_tscr = w_rtc*atp/(theta_rtc+atp)
    rtcr_tlel = max*atp/(thr+atp)
    rtcr_tlr = ribo_h.*rm_r.*rtcr_tlel

    rtca_tot = rtca .+ rdrtca_0
    rdrtca = k1_a.*ribo_d.*rtca_tot./(k1_a.*ribo_d.+k2_a.+k3_a.*atp)
    rtcb_tot = rtcb .+ rtrtcb_0
    rtrtcb = ka_b.*ribo_t.*rtcb_tot./(ka_b.*ribo_t.+kb_b.+kc_b.*atp)

    return rm_a, rtca, rm_b, rtcb, rm_r, rtcr, ribo_d, ribo_h, ribo_t, alpha, fa, ra, v, tscr, tlr, tscr_b, tlr_b, rtcr_tscr, rtcr_tlr, rdrtca, rtrtcb
end

function vary_parameter_for_full_model_rd(varying, param, range, model, init, tspan, params, species)
    rm_a_res = []
    rtca_res = []
    rm_b_res = []
    rtcb_res = []
    rm_r_res = []
    rtcr_res = []
    ribo_d_res = []
    ribo_h_res = []
    ribo_t_res = []
    rdrtca_res = []
    rtrtcb_res = []
    solutions = [];
    solutions_t = [];
    for x in range
        varying[param] = x
        sol = simple_solve!(model, init, tspan, params)
        rm_a, rtca, rm_b, rtcb, rm_r, rtcr, ribo_d, ribo_h, ribo_t, alpha, fa, ra, v, tscr, tlr, tscr_b, tlr_b, rtcr_tscr, rtcr_tlr, rdrtca, rtrtcb = split_rd(species, sol)
        push!(rm_a_res, rm_a)
        push!(rtca_res, rtca)
        push!(rm_b_res, rm_b)
        push!(rtcb_res, rtcb)
        push!(rm_r_res, rm_r)
        push!(rtcr_res, rtcr)
        push!(ribo_d_res, ribo_d)
        push!(ribo_h_res, ribo_h)
        push!(ribo_t_res, ribo_t)
        push!(rdrtca_res, rdrtca)
        push!(rtrtcb_res, rtrtcb)
        push!(solutions, sol)
        push!(solutions_t, sol.t)
    end
    sol1 = solutions[1]
    sol2 = solutions[2]
    sol3 = solutions[3]
    sol4 = solutions[4]
    sol5 = solutions[5]
    sol6 = solutions[6]
    return rm_a_res, rtca_res, rm_b_res, rtcb_res, rm_r_res, rtcr_res, ribo_d_res, ribo_h_res, ribo_t_res, rdrtca_res, rtrtcb_res, solutions, solutions_t, sol1, sol2, sol3, sol4, sol5, sol6
end

function plots_for_param_rd(v_param, varying, param, range, model, init, tspan, params, species) # varying should be init or param
    rm_a_res, rtca_res, rm_b_res, rtcb_res, rm_r_res, rtcr_res, ribo_d_res, ribo_h_res, ribo_t_res, rdrtca_res, rtrtcb_res, solutions, solutions_t, sol1, sol2, sol3, sol4, sol5, sol6 = vary_parameter_for_full_model_rd(varying, param, range, model, init, tspan, params, species);

    p1 = plot(sol1, labels=labels, palette=:tab10, title="$v_param = $(collect(range)[1])")
    p2 = plot(sol2, labels=labels, palette=:tab10, title="$v_param = $(collect(range)[2])")
    p3 = plot(sol3, labels=labels, palette=:tab10, title="$v_param = $(collect(range)[3])")
    p4 = plot(sol4, labels=labels, palette=:tab10, title="$v_param = $(collect(range)[4])")
    p5 = plot(sol5, labels=labels, palette=:tab10, title="$v_param = $(collect(range)[5])")
    p6 = plot(sol6, labels=labels, palette=:tab10, title="$v_param = $(collect(range)[6])")

    p_main = plot(p1, p2, p3, p4, p5, p6, size=(1000,800))

    labels1 = ["$(collect(range)[1])" "$(collect(range)[2])" "$(collect(range)[3])" "$(collect(range)[4])" "$(collect(range)[5])" "$(collect(range)[6])"]
        
    p1_a = plot(solutions_t, rtcr_res, title="RtcR", xlabel="t", palette=:tab10, labels=labels1)
    p2_a = plot(solutions_t, rtca_res, title="RtcA", xlabel="t", palette=:tab10, labels=labels1)
    p3_a = plot(solutions_t, rtcb_res, title="RtcB", xlabel="t", palette=:tab10, labels=labels1)
    p4_a = plot(solutions_t, ribo_d_res, title="Rd", xlabel="t", palette=:tab10, labels=labels1)
    p5_a = plot(solutions_t, ribo_t_res, title="Rt", xlabel="t", palette=:tab10, labels=labels1)
    p6_a = plot(solutions_t, ribo_h_res, title="Rh", xlabel="t", palette=:tab10, labels=labels1)

    p_a_main = plot(p1_a, p2_a, p3_a, p4_a, p5_a, p6_a, size=(1000,800), plot_title="Varying $v_param")
    return p_main, p_a_main
end









# new 
function split_rtc!(solution, species, ribo_tot)
    solDF = DataFrame([[j[i] for j in solution.u] for i=1:length(solution.u[1])], species)
    rdrtca = k1_a.*solDF[:, :ribo_d].*solDF[:, :rtca]./(k1_a.*solDF[:, :ribo_d].+k2_a.+k3_a.*atp)
    rtrtcb = ka_b.*solDF[:, :ribo_t].*solDF[:, :rtcb]./(ka_b.*solDF[:, :ribo_t].+kb_b.+kc_b.*atp)
    ribo_h = ribo_tot .- solDF[:, :ribo_d] .- solDF[:, :ribo_t] .- rdrtca .- rtrtcb
    return solDF[:, :rm_a], solDF[:, :rtca], solDF[:, :rm_b], solDF[:, :rtcb], solDF[:, :rm_r], solDF[:, :rtcr], rdrtca, rtrtcb, ribo_h, solDF[:, :ribo_d], solDF[:, :ribo_t]
end

function empty_res()
    sols=[]
    sols_t=[]
    rm_a=[]
    rtca=[]
    rm_b=[]
    rtcb=[]
    rm_r=[]
    rtcr=[]
    ribo_h=[]
    ribo_d=[]
    ribo_t=[]
    rdrtca=[]
    rtrtcb=[]
    ribo_tot=[]
    return sols, sols_t, rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rdrtca, rtrtcb, ribo_h, ribo_d, ribo_t, ribo_tot
end

function fill_empty_arrays(sols, sols_t, rm_a_sols, rtca_sols, rm_b_sols, rtcb_sols, rm_r_sols, rtcr_sols, rdrtca_sols, rtrtcb_sols, ribo_h_sols, ribo_d_sols, ribo_t_sols, ribo_tot_sols, sol, rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rdrtca, rtrtcb, ribo_h, ribo_d, ribo_t, ribo_tot)
    push!(sols, sol)
    push!(sols_t, sol.t)
    push!(rm_a_sols, rm_a)
    push!(rtca_sols, rtca)
    push!(rm_b_sols, rm_b)
    push!(rtcb_sols, rtcb)
    push!(rm_r_sols, rm_r)
    push!(rtcr_sols, rtcr)
    push!(ribo_h_sols, ribo_h)
    push!(ribo_d_sols, ribo_d)
    push!(ribo_t_sols, ribo_t)
    push!(rdrtca_sols, rdrtca)
    push!(rtrtcb_sols, rtrtcb)
    push!(ribo_tot_sols, ribo_tot)

    return sols, sols_t, rm_a_sols, rtca_sols, rm_b_sols, rtcb_sols, rm_r_sols, rtcr_sols, rdrtca_sols, rtrtcb_sols, ribo_h_sols, ribo_d_sols, ribo_t_sols, ribo_tot_sols
end


function vary_param(range, param_or_init, x, model, init, params, species)
    sols, sols_t, rm_a_sols, rtca_sols, rm_b_sols, rtcb_sols, rm_r_sols, rtcr_sols, rdrtca_sols, rtrtcb_sols, ribo_h_sols, ribo_d_sols, ribo_t_sols, ribo_tot_sols = empty_res()
    for i in range
        param_or_init[x] = i
        ribo_tot = init_new[11] + init_new[10] + init_new[9] + rdrtca_0 + rtrtcb_0 
        sol = simple_solve!(model, init, tspan, params)
        rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rdrtca, rtrtcb, ribo_h, ribo_d, ribo_t = split_full_model(species, sol)
        sols, sols_t, rm_a_sols, rtca_sols, rm_b_sols, rtcb_sols, rm_r_sols, rtcr_sols, rdrtca_sols, rtrtcb_sols, ribo_h_sols, ribo_d_sols, ribo_t_sols, ribo_tot_sols = fill_empty_arrays(sols, sols_t, rm_a_sols, rtca_sols, rm_b_sols, rtcb_sols, rm_r_sols, rtcr_sols, rdrtca_sols, rtrtcb_sols, ribo_h_sols, ribo_d_sols, ribo_t_sols, ribo_tot_sols, sol, rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rdrtca, rtrtcb, ribo_h, ribo_d, ribo_t, ribo_tot)
    end
    return sols, sols_t, rm_a_sols, rtca_sols, rm_b_sols, rtcb_sols, rm_r_sols, rtcr_sols, rdrtca_sols, rtrtcb_sols, ribo_h_sols, ribo_d_sols, ribo_t_sols, ribo_tot_sols
end

function logplots(range, param_or_init, x, model, init, params, species, labels, param_name) # varying should be init or param
    sols, sols_t, rm_a_sols, rtca_sols, rm_b_sols, rtcb_sols, rm_r_sols, rtcr_sols, rdrtca_sols, rtrtcb_sols, ribo_h_sols, ribo_d_sols, ribo_t_sols = vary_param(range, param_or_init, x, model, init, params, species)

    p1 = plot(sols[1][2:end], palette=:tab10, title="$param_name = $(collect(range)[1])", ylabel="amount of species", xaxis=(:log10, (1,Inf)), yaxis=(:log10, (1,Inf)), legend=false)
    p2 = plot(sols[2][2:end], palette=:tab10, title="$param_name = $(collect(range)[2])", ylabel="amount of species", xaxis=(:log10, (1,Inf)), yaxis=(:log10, (1,Inf)), legend=false)
    p3 = plot(sols[3][2:end], palette=:tab10, title="$param_name = $(collect(range)[3])", ylabel="amount of species", xaxis=(:log10, (1,Inf)), yaxis=(:log10, (1,Inf)), legend=false)
    p4 = plot(sols[4][2:end], palette=:tab10, title="$param_name = $(collect(range)[4])", ylabel="amount of species", xaxis=(:log10, (1,Inf)), yaxis=(:log10, (1,Inf)), legend=false)
    p5 = plot(sols[5][2:end], palette=:tab10, title="$param_name = $(collect(range)[5])", ylabel="amount of species", xaxis=(:log10, (1,Inf)), yaxis=(:log10, (1,Inf)), legend=false)
    p6 = plot(sols[6][2:end], palette=:tab10, title="$param_name = $(collect(range)[6])", ylabel="amount of species", xaxis=(:log10, (1,Inf)), yaxis=(:log10, (1,Inf)), legend=false)
    legend = plot([0 0 0 0 0 0 0 0 0 0 0], palette=:tab10, axis=([], false), foreground_color_legend = nothing, showaxislabels=false, showaxis=false, grid=false, label=labels, legend=:inside)
    plot_1 = plot(p1, p2, p3, p4, p5, p6, legend, size=(1000,800), layout=@layout([[A B C; D E F] G{0.2w}]))
    # plot_1 = plot(p1, p3, p5, legend, size=(800,400), layout=@layout([[A B C] G{0.2w}]), margin=2mm)
    # plot_1 = plot(p1, legend, size=(1000,800))

    p7 = plot(sols_t, rtcr_sols, title="RtcR", xlabel="t", ylabel="amount of species", palette=:tab10, yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)))
    p8 = plot(sols_t, rtca_sols, title="RtcA", xlabel="t", ylabel="amount of species", palette=:tab10, yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)))
    # p8 = plot(sols_t, rtca_sols, title="RtcA/B", xlabel="t", ylabel="amount of species", palette=:tab10, yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)))

    p9 = plot(sols_t, rtcb_sols, title="RtcB", xlabel="t", ylabel="amount of species", palette=:tab10, yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)))
    p10 = plot(sols_t, ribo_d_sols, title="Rd", xlabel="t", ylabel="amount of species", palette=:tab10, xaxis=(:log10, (1,Inf)))
    p11 = plot(sols_t, ribo_t_sols, title="Rt", xlabel="t", ylabel="amount of species", palette=:tab10, xaxis=(:log10, (1,Inf)))
    p12 = plot(sols_t, ribo_h_sols, title="Rh", xlabel="t", ylabel="amount of species", palette=:tab10, xaxis=(:log10, (1,Inf)))

    plot_2 = plot(p7, p8, p9, p10, p11, p12, labels=["$(collect(range)[1])" "$(collect(range)[2])" "$(collect(range)[3])" "$(collect(range)[4])" "$(collect(range)[5])" "$(collect(range)[6])"], size=(1000,800), plot_title="Varying $param_name")
    # plot_2 = plot(p7, p8, labels=["$(collect(range)[1])" "$(collect(range)[2])" "$(collect(range)[3])" "$(collect(range)[4])" "$(collect(range)[5])" "$(collect(range)[6])"], size=(600,400), plot_title="Varying $param_name")

    return plot_1, plot_2
end


function split_rh(species, solution) # 2

    solDF = DataFrame([[j[i] for j in solution.u] for i=1:length(solution.u[1])], species)

    rm_a = solDF[:, :rm_a]
    rtca = solDF[:, :rtca]
    rm_b = solDF[:, :rm_b]
    rtcb = solDF[:, :rtcb]
    rm_r = solDF[:, :rm_r]
    rtcr = solDF[:, :rtcr]
    ribo_d = solDF[:, :ribo_d]
    ribo_t = solDF[:, :ribo_t]

    n = 6
    rdrtca = k1_a.*solDF[:, :ribo_d].*solDF[:, :rtca]./(k1_a.*solDF[:, :ribo_d].+k2_a.+k3_a.*atp)
    rtrtcb = ka_b.*solDF[:, :ribo_t].*solDF[:, :rtcb]./(ka_b.*solDF[:, :ribo_t].+kb_b.+kc_b.*atp)
    ribo_tot = ribo_t_0 + ribo_d_0 + ribo_h_0 + rdrtca_0 + rtrtcb_0
    ribo_h = ribo_tot .- solDF[:, :ribo_d] .- solDF[:, :ribo_t] .- rdrtca .- rtrtcb

    alpha = solDF[:, :ribo_t]./kr
    fa = ([1].+alpha).^n./(L.*(([1].+c.*alpha).^n)+([1].+alpha).^n) # fraction of active RtcR
    ra = fa.*solDF[:, :rtcr] # amount of active RtcR
    v = ra.*k1.*sigma.*k3.*atp.*k5./(k1.*sigma.*k3.*atp.+k1.*sigma.*(k4.+k5).+k2.*(k4.+k5).+k3.*atp.*k5) # rate of open complex formation
    
    # rtca
    tscr = v.*(w_rtc*atp/(theta_rtc+atp)) # transcription rate 
    tlel = max*atp/(thr+atp) # translation elongation rate
    tlr = ribo_h.*solDF[:, :rm_a].*tlel # translation rate so formation of Rtc 

    # rtcb
    tscr_b = v.*(w_rtc*atp/(theta_rtc+atp)) # transcription rate 
    tlel_b = max*atp/(thr+atp) # translation elongation rate
    tlr_b = ribo_h.*solDF[:, :rm_b].*tlel_b # translation rate so formation of Rtc 

    # rtcr 
    rtcr_tscr = w_rtc*atp/(theta_rtc+atp)
    rtcr_tlel = max*atp/(thr+atp)
    rtcr_tlr = ribo_h.*solDF[:, :rm_r].*rtcr_tlel



    return rm_a, rtca, rm_b, rtcb, rm_r, rtcr, ribo_d, ribo_h, ribo_t, alpha, fa, ra, v, tscr, tlr, tscr_b, tlr_b, rtcr_tscr, rtcr_tlr, rdrtca, rtrtcb
end

function split_riboint(species, solution) # 3

    solDF = DataFrame([[j[i] for j in solution.u] for i=1:length(solution.u[1])], species)

    rm_a = solDF[:, :rm_a]
    rtca = solDF[:, :rtca]
    rm_b = solDF[:, :rm_b]
    rtcb = solDF[:, :rtcb]
    rm_r = solDF[:, :rm_r]
    rtcr = solDF[:, :rtcr]
    ribo_d = solDF[:, :ribo_d]
    ribo_h = solDF[:, :ribo_h]
    ribo_t = solDF[:, :ribo_t]
    rdrtca = solDF[:, :rdrtca]
    rtrtcb = solDF[:, :rtrtcb]

    n = 6

    alpha = solDF[:, :ribo_t]./kr
    fa = ([1].+alpha).^n./(L.*(([1].+c.*alpha).^n)+([1].+alpha).^n) # fraction of active RtcR
    ra = fa.*solDF[:, :rtcr] # amount of active RtcR
    v = ra.*k1.*sigma.*k3.*atp.*k5./(k1.*sigma.*k3.*atp.+k1.*sigma.*(k4.+k5).+k2.*(k4.+k5).+k3.*atp.*k5) # rate of open complex formation
    
    # rtca
    tscr = v.*(w_rtc*atp/(theta_rtc+atp)) # transcription rate 
    tlel = max*atp/(thr+atp) # translation elongation rate
    tlr = solDF[:, :ribo_h].*solDF[:, :rm_a].*tlel # translation rate so formation of Rtc 

    # rtcb
    tscr_b = v.*(w_rtc*atp/(theta_rtc+atp)) # transcription rate 
    tlel_b = max*atp/(thr+atp) # translation elongation rate
    tlr_b = solDF[:, :ribo_h].*solDF[:, :rm_b].*tlel_b # translation rate so formation of Rtc 

    # rtcr 
    rtcr_tscr = w_rtc*atp/(theta_rtc+atp)
    rtcr_tlel = max*atp/(thr+atp)
    rtcr_tlr = solDF[:, :ribo_h].*solDF[:, :rm_r].*rtcr_tlel


    return rm_a, rtca, rm_b, rtcb, rm_r, rtcr, ribo_d, ribo_h, ribo_t, alpha, fa, ra, v, tscr, tlr, tscr_b, tlr_b, rtcr_tscr, rtcr_tlr, rdrtca, rtrtcb
end


function split_full_model(species, solution)
    solDF = DataFrame([[j[i] for j in solution.u] for i=1:length(solution.u[1])], species)
    return solDF[:, :rm_a], solDF[:, :rtca], solDF[:, :rm_b], solDF[:, :rtcb], solDF[:, :rm_r], solDF[:, :rtcr], solDF[:, :rdrtca], solDF[:, :rtrtcb], solDF[:, :ribo_h], solDF[:, :ribo_d], solDF[:, :ribo_t]
end


function split_reduced_model(species, solution)
    solDF = DataFrame([[j[i] for j in solution.u] for i=1:length(solution.u[1])], species)
    rdrtca = k1_a.*solDF[:, :ribo_d].*solDF[:, :rtca]./(k1_a.*solDF[:, :ribo_d].+k2_a.+k3_a.*atp)
    rtrtcb = ka_b.*solDF[:, :ribo_t].*solDF[:, :rtcb]./(ka_b.*solDF[:, :ribo_t].+kb_b.+kc_b.*atp)
    return solDF[:, :rm_a], solDF[:, :rtca], solDF[:, :rm_b], solDF[:, :rtcb], solDF[:, :rm_r], solDF[:, :rtcr], rdrtca, rtrtcb, solDF[:, :ribo_h], solDF[:, :ribo_d], solDF[:, :ribo_t]
end



function split_reduced_model_rh(species, solution)
    solDF = DataFrame([[j[i] for j in solution.u] for i=1:length(solution.u[1])], species)
    rdrtca = k1_a.*solDF[:, :ribo_d].*solDF[:, :rtca]./(k1_a.*solDF[:, :ribo_d].+k2_a.+k3_a.*atp)
    rtrtcb = ka_b.*solDF[:, :ribo_t].*solDF[:, :rtcb]./(ka_b.*solDF[:, :ribo_t].+kb_b.+kc_b.*atp)
    ribo_h = ribo_tot .- solDF[:, :ribo_d] .- solDF[:, :ribo_t] .- rdrtca .- rtrtcb 
    # rdrtca = solDF[:, :rtca].*k1.*solDF[:, :ribo_d].*(k4.+k3_a)./(k2.*k4.+k2.*k3_a.+k3.*atp.*k3_a.+k1.*solDF[:, :ribo_d].*(k4.+k3_a.+k3.*atp))
    # rtrtcb = solDF[:, :rtcb].*k1.*solDF[:, :ribo_t].*(k4.+k3_a)./(k2.*k4.+k2.*k3_a.+k3.*atp.*k3_a.+k1.*solDF[:, :ribo_t].*(k4.+k3_a.+k3.*atp))
    return solDF[:, :rm_a], solDF[:, :rtca], solDF[:, :rm_b], solDF[:, :rtcb], solDF[:, :rm_r], solDF[:, :rtcr], rdrtca, rtrtcb, ribo_h, solDF[:, :ribo_d], solDF[:, :ribo_t]
end





function vary_param_reduced(range, param_or_init, x, model, init, params, species)
    sols, sols_t, rm_a_sols, rtca_sols, rm_b_sols, rtcb_sols, rm_r_sols, rtcr_sols, rdrtca_sols, rtrtcb_sols, ribo_h_sols, ribo_d_sols, ribo_t_sols, ribo_tot_sols = empty_res()
    for i in range
        param_or_init[x] = i
        ribo_tot = init_new[7] + init_new[8] + init_new[9] + rdrtca_0 + rtrtcb_0 
        sol = simple_solve!(model, init, tspan, params)
        rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rdrtca, rtrtcb, ribo_h, ribo_d, ribo_t = split_reduced_model(species, sol)
        sols, sols_t, rm_a_sols, rtca_sols, rm_b_sols, rtcb_sols, rm_r_sols, rtcr_sols, rdrtca_sols, rtrtcb_sols, ribo_h_sols, ribo_d_sols, ribo_t_sols, ribo_tot_sols = fill_empty_arrays(sols, sols_t, rm_a_sols, rtca_sols, rm_b_sols, rtcb_sols, rm_r_sols, rtcr_sols, rdrtca_sols, rtrtcb_sols, ribo_h_sols, ribo_d_sols, ribo_t_sols, ribo_tot_sols, sol, rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rdrtca, rtrtcb, ribo_h, ribo_d, ribo_t, ribo_tot)
    end
    return sols, sols_t, rm_a_sols, rtca_sols, rm_b_sols, rtcb_sols, rm_r_sols, rtcr_sols, rdrtca_sols, rtrtcb_sols, ribo_h_sols, ribo_d_sols, ribo_t_sols, ribo_tot_sols
end

function logplots_reduced(range, param_or_init, x, model, init, params, species, labels, param_name) # varying should be init or param
    sols, sols_t, rm_a_sols, rtca_sols, rm_b_sols, rtcb_sols, rm_r_sols, rtcr_sols, rdrtca_sols, rtrtcb_sols, ribo_h_sols, ribo_d_sols, ribo_t_sols = vary_param_reduced(range, param_or_init, x, model, init, params, species)

    p1 = plot(sols[1][2:end], palette=:tab10, title="$param_name = $(collect(range)[1])", ylabel="amount of species", xaxis=(:log10, (1,Inf)), yaxis=(:log10, (1,Inf)), legend=false)
    p2 = plot(sols[2][2:end], palette=:tab10, title="$param_name = $(collect(range)[2])", ylabel="amount of species", xaxis=(:log10, (1,Inf)), yaxis=(:log10, (1,Inf)), legend=false)
    p3 = plot(sols[3][2:end], palette=:tab10, title="$param_name = $(collect(range)[3])", ylabel="amount of species", xaxis=(:log10, (1,Inf)), yaxis=(:log10, (1,Inf)), legend=false)
    p4 = plot(sols[4][2:end], palette=:tab10, title="$param_name = $(collect(range)[4])", ylabel="amount of species", xaxis=(:log10, (1,Inf)), yaxis=(:log10, (1,Inf)), legend=false)
    p5 = plot(sols[5][2:end], palette=:tab10, title="$param_name = $(collect(range)[5])", ylabel="amount of species", xaxis=(:log10, (1,Inf)), yaxis=(:log10, (1,Inf)), legend=false)
    p6 = plot(sols[6][2:end], palette=:tab10, title="$param_name = $(collect(range)[6])", ylabel="amount of species", xaxis=(:log10, (1,Inf)), yaxis=(:log10, (1,Inf)), legend=false)
    legend = plot([0 0 0 0 0 0 0 0 0], palette=:tab10, axis=([], false), foreground_color_legend = nothing, showaxislabels=false, showaxis=false, grid=false, label=labels, legend=:inside)
    plot_1 = plot(p1, p2, p3, p4, p5, p6, legend, size=(1000,800), layout=@layout([[A B C; D E F] G{0.2w}]))
    # plot_1 = plot(p1, p3, p5, legend, size=(800,400), layout=@layout([[A B C] G{0.2w}]), margin=2mm)
    # plot_1 = plot(p1, legend, size=(1000,800))
    
    p7 = plot(sols_t, rtcr_sols, title="RtcR", xlabel="t", ylabel="amount of species", palette=:tab10, yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)))
    p8 = plot(sols_t, rtca_sols, title="RtcA", xlabel="t", ylabel="amount of species", palette=:tab10, yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)))
    # p8 = plot(sols_t, rtca_sols, title="RtcA/B", xlabel="t", ylabel="amount of species", palette=:tab10, yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)))

    p9 = plot(sols_t, rtcb_sols, title="RtcB", xlabel="t", ylabel="amount of species", palette=:tab10, yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)))
    p10 = plot(sols_t, ribo_d_sols, title="Rd", xlabel="t", ylabel="amount of species", palette=:tab10, xaxis=(:log10, (1,Inf)))
    p11 = plot(sols_t, ribo_t_sols, title="Rt", xlabel="t", ylabel="amount of species", palette=:tab10, xaxis=(:log10, (1,Inf)))
    p12 = plot(sols_t, ribo_h_sols, title="Rh", xlabel="t", ylabel="amount of species", palette=:tab10, xaxis=(:log10, (1,Inf)))

    plot_2 = plot(p7, p8, p9, p10, p11, p12, labels=["$(collect(range)[1])" "$(collect(range)[2])" "$(collect(range)[3])" "$(collect(range)[4])" "$(collect(range)[5])" "$(collect(range)[6])"], size=(1000,800), plot_title="Varying $param_name")
    # plot_2 = plot(p7, p8, labels=["$(collect(range)[1])" "$(collect(range)[2])" "$(collect(range)[3])" "$(collect(range)[4])" "$(collect(range)[5])" "$(collect(range)[6])"], size=(600,400), plot_title="Varying $param_name")

    return plot_1, plot_2
end

function split_reduced_model_steady_state(species, solution)
    solDF = DataFrame([[j[i] for j in solution.u] for i=1:length(solution.u[1])], species)
    rm_a = solDF[end, :rm_a]
    rtca = solDF[end, :rtca]
    rm_b = solDF[end, :rm_b]
    rtcb = solDF[end, :rtcb]
    rm_r = solDF[end, :rm_r]
    rtcr = solDF[end, :rtcr]
    ribo_h = solDF[end, :ribo_h]
    ribo_d = solDF[end, :ribo_d]
    ribo_t = solDF[end, :ribo_t]
    rdrtca = @. k1_a*solDF[end, :ribo_d]*solDF[end, :rtca]/(k1_a*solDF[end, :ribo_d]+k2_a+k3_a*atp)
    rtrtcb = @. ka_b*solDF[end, :ribo_t]*solDF[end, :rtcb]/(ka_b*solDF[end, :ribo_t]+kb_b+kc_b*atp)
    return rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rdrtca, rtrtcb, ribo_h, ribo_d, ribo_t
end

function split_full_model_steady_state(species, solution)
    solDF = DataFrame([[j[i] for j in solution.u] for i=1:length(solution.u[1])], species)
    return solDF[end, :rm_a], solDF[end, :rtca], solDF[end, :rm_b], solDF[end, :rtcb], solDF[end, :rm_r], solDF[end, :rtcr], solDF[end, :rdrtca], solDF[end, :rtrtcb], solDF[end, :ribo_h], solDF[end, :ribo_d], solDF[end, :ribo_t]
end

function param_vs_ss(long_range, param_or_init, x, model, init, params, species)
    sols, sols_t, rm_a_sols, rtca_sols, rm_b_sols, rtcb_sols, rm_r_sols, rtcr_sols, rdrtca_sols, rtrtcb_sols, ribo_h_sols, ribo_d_sols, ribo_t_sols, ribo_tot_sols = empty_res()
    for i in ProgressBar(long_range)
        param_or_init[x] = i
        sol = simple_solve!(model, init, tspan, params)
        rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rdrtca, rtrtcb, ribo_h, ribo_d, ribo_t = split_reduced_model_steady_state(species, sol)
        sols, sols_t, rm_a_sols, rtca_sols, rm_b_sols, rtcb_sols, rm_r_sols, rtcr_sols, rdrtca_sols, rtrtcb_sols, ribo_h_sols, ribo_d_sols, ribo_t_sols, ribo_tot_sols = fill_empty_arrays(sols, sols_t, rm_a_sols, rtca_sols, rm_b_sols, rtcb_sols, rm_r_sols, rtcr_sols, rdrtca_sols, rtrtcb_sols, ribo_h_sols, ribo_d_sols, ribo_t_sols, ribo_tot_sols, sol, rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rdrtca, rtrtcb, ribo_h, ribo_d, ribo_t, ribo_tot)
    end
    return rm_a_sols, rtca_sols, rm_b_sols, rtcb_sols, rm_r_sols, rtcr_sols, rdrtca_sols, rtrtcb_sols, ribo_h_sols, ribo_d_sols, ribo_t_sols
end

function param_vs_ss_full_model(long_range, param_or_init, x, model, init, params, species)
    sols, sols_t, rm_a_sols, rtca_sols, rm_b_sols, rtcb_sols, rm_r_sols, rtcr_sols, rdrtca_sols, rtrtcb_sols, ribo_h_sols, ribo_d_sols, ribo_t_sols, ribo_tot_sols = empty_res()
    for i in long_range
        param_or_init[x] = i
        sol = simple_solve!(model, init, tspan, params)
        rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rdrtca, rtrtcb, ribo_h, ribo_d, ribo_t = split_full_model_steady_state(species, sol)
        sols, sols_t, rm_a_sols, rtca_sols, rm_b_sols, rtcb_sols, rm_r_sols, rtcr_sols, rdrtca_sols, rtrtcb_sols, ribo_h_sols, ribo_d_sols, ribo_t_sols, ribo_tot_sols = fill_empty_arrays(sols, sols_t, rm_a_sols, rtca_sols, rm_b_sols, rtcb_sols, rm_r_sols, rtcr_sols, rdrtca_sols, rtrtcb_sols, ribo_h_sols, ribo_d_sols, ribo_t_sols, ribo_tot_sols, sol, rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rdrtca, rtrtcb, ribo_h, ribo_d, ribo_t, ribo_tot)
    end
    return rm_a_sols, rtca_sols, rm_b_sols, rtcb_sols, rm_r_sols, rtcr_sols, rdrtca_sols, rtrtcb_sols, ribo_h_sols, ribo_d_sols, ribo_t_sols
end


function plot_long_range_results(long_range, param_or_init, x, model, init, params, species, param_name)
    rm_a_sols, rtca_sols, rm_b_sols, rtcb_sols, rm_r_sols, rtcr_sols, rdrtca_sols, rtrtcb_sols, ribo_h_sols, ribo_d_sols, ribo_t_sols = param_vs_ss(long_range, param_or_init, x, model, init, params, species);
    # p = plot(long_range, [rm_a_sols, rtca_sols, rm_b_sols, rtcb_sols, rm_r_sols, rtcr_sols, ribo_h_sols, ribo_d_sols, ribo_t_sols], yaxis=(:log10, (1,Inf)), labels=labels_rtc, title="Reduced model - $param_name")
    p1 = plot(long_range, [ribo_h_sols], labels="Rh", xaxis="$param_name", yaxis="amount of species", labelfontsize=8)
    p2 = plot(long_range, [ribo_d_sols], labels="Rd", xaxis="$param_name", yaxis="amount of species", labelfontsize=8)
    p3 = plot(long_range, [ribo_t_sols], labels="Rt", xaxis="$param_name", yaxis="amount of species", labelfontsize=8)
    p4 = plot(long_range, [rtca_sols, rtcb_sols], labels=["RtcA" "RtcB"], xaxis="$param_name", yaxis="amount of species", labelfontsize=8)
    p5 = plot(long_range, [rtcr_sols], labels="RtcR", xaxis="$param_name", yaxis="amount of species", labelfontsize=8)
    p6 = plot(long_range, [rdrtca_sols, rtrtcb_sols], labels=["RdRtcA" "RtRtcB"], xaxis="$param_name", yaxis="amount of species", labelfontsize=8)
    p7 = plot(long_range, [rm_a_sols, rm_b_sols], labels=["mRNA-RtcA" "mRNA-RtcB"], xaxis="$param_name", yaxis="amount of species", labelfontsize=8)
    p8 = plot(long_range, [rm_r_sols], labels="mRNA-RtcR", xaxis="$param_name", yaxis="amount of species", labelfontsize=8)
    plot(p1, p2, p3, p4, p5, p6, p7, p8, plot_title="$param_name effect", layout=(2,4), size=(1000,600))
end

function full_plot_long_range_results(long_range, param_or_init, x, model, init, params, species, param_name)
    rm_a_sols, rtca_sols, rm_b_sols, rtcb_sols, rm_r_sols, rtcr_sols, rdrtca_sols, rtrtcb_sols, ribo_h_sols, ribo_d_sols, ribo_t_sols = param_vs_ss(long_range, param_or_init, x, model, init, params, species);
    # p = plot(long_range, [rm_a_sols, rtca_sols, rm_b_sols, rtcb_sols, rm_r_sols, rtcr_sols, ribo_h_sols, ribo_d_sols, ribo_t_sols], yaxis=(:log10, (1,Inf)), labels=labels_rtc, title="Reduced model - $param_name")
    p1 = plot(long_range, [ribo_h_sols], labels="Rh", xaxis="$param_name", yaxis="amount of species", labelfontsize=8)
    p2 = plot(long_range, [ribo_d_sols], labels="Rd", xaxis="$param_name", yaxis="amount of species", labelfontsize=8)
    p3 = plot(long_range, [ribo_t_sols], labels="Rt", xaxis="$param_name", yaxis="amount of species", labelfontsize=8)
    p4 = plot(long_range, [rtca_sols, rtcb_sols], labels=["RtcA" "RtcB"], xaxis="$param_name", yaxis="amount of species", labelfontsize=8)
    p5 = plot(long_range, [rtcr_sols], labels="RtcR", xaxis="$param_name", yaxis="amount of species", labelfontsize=8)
    p6 = plot(long_range, [rdrtca_sols, rtrtcb_sols], labels=["RdRtcA" "RtRtcB"], xaxis="$param_name", yaxis="amount of species", labelfontsize=8)
    p7 = plot(long_range, [rm_a_sols, rm_b_sols], labels=["mRNA-RtcA" "mRNA-RtcB"], xaxis="$param_name", yaxis="amount of species", labelfontsize=8)
    p8 = plot(long_range, [rm_r_sols], labels="mRNA-RtcR", xaxis="$param_name", yaxis="amount of species", labelfontsize=8)
    plot(p1, p2, p3, p4, p5, p6, p7, p8, plot_title="$param_name effect - full", layout=(2,4), size=(1000,600))
end















# for long ranges for varying parameter this code is quicker than above as it doesn't have the unnessecary bits

function empty_res1()
    rm_a=[]
    rtca=[]
    rm_b=[]
    rtcb=[]
    rm_r=[]
    rtcr=[]
    ribo_h=[]
    ribo_d=[]
    ribo_t=[]
    rdrtca=[]
    rtrtcb=[]
    return rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rdrtca, rtrtcb, ribo_h, ribo_d, ribo_t
end

function fill_empty_arrays1(rm_a_sols, rtca_sols, rm_b_sols, rtcb_sols, rm_r_sols, rtcr_sols, rdrtca_sols, rtrtcb_sols, ribo_h_sols, ribo_d_sols, ribo_t_sols, rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rdrtca, rtrtcb, ribo_h, ribo_d, ribo_t)
    push!(rm_a_sols, rm_a)
    push!(rtca_sols, rtca)
    push!(rm_b_sols, rm_b)
    push!(rtcb_sols, rtcb)
    push!(rm_r_sols, rm_r)
    push!(rtcr_sols, rtcr)
    push!(ribo_h_sols, ribo_h)
    push!(ribo_d_sols, ribo_d)
    push!(ribo_t_sols, ribo_t)
    push!(rdrtca_sols, rdrtca)
    push!(rtrtcb_sols, rtrtcb)

    return rm_a_sols, rtca_sols, rm_b_sols, rtcb_sols, rm_r_sols, rtcr_sols, rdrtca_sols, rtrtcb_sols, ribo_h_sols, ribo_d_sols, ribo_t_sols
end

function split_reduced_model_steady_state(species, solution)
    solDF = DataFrame([[j[i] for j in solution.u] for i=1:length(solution.u[1])], species)
    rm_a = solDF[end, :rm_a]
    rtca = solDF[end, :rtca]
    rm_b = solDF[end, :rm_b]
    rtcb = solDF[end, :rtcb]
    rm_r = solDF[end, :rm_r]
    rtcr = solDF[end, :rtcr]
    ribo_h = solDF[end, :ribo_h]
    ribo_d = solDF[end, :ribo_d]
    ribo_t = solDF[end, :ribo_t]
    rdrtca = @. k1_a*solDF[end, :ribo_d]*solDF[end, :rtca]/(k1_a*solDF[end, :ribo_d]+k2_a+k3_a*atp)
    rtrtcb = @. ka_b*solDF[end, :ribo_t]*solDF[end, :rtcb]/(ka_b*solDF[end, :ribo_t]+kb_b+kc_b*atp)
    return rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rdrtca, rtrtcb, ribo_h, ribo_d, ribo_t
end


function param_vs_ss(long_range, param_or_init, x, model, init, params, species)
    rm_a_sols, rtca_sols, rm_b_sols, rtcb_sols, rm_r_sols, rtcr_sols, rdrtca_sols, rtrtcb_sols, ribo_h_sols, ribo_d_sols, ribo_t_sols = empty_res1()
    for i in ProgressBar(long_range)
        param_or_init[x] = i
        sol = simple_solve!(model, init, tspan, params)
        rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rdrtca, rtrtcb, ribo_h, ribo_d, ribo_t = split_reduced_model_steady_state(species, sol)
        rm_a_sols, rtca_sols, rm_b_sols, rtcb_sols, rm_r_sols, rtcr_sols, rdrtca_sols, rtrtcb_sols, ribo_h_sols, ribo_d_sols, ribo_t_sols = fill_empty_arrays1(rm_a_sols, rtca_sols, rm_b_sols, rtcb_sols, rm_r_sols, rtcr_sols, rdrtca_sols, rtrtcb_sols, ribo_h_sols, ribo_d_sols, ribo_t_sols, rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rdrtca, rtrtcb, ribo_h, ribo_d, ribo_t)
    end
    return rm_a_sols, rtca_sols, rm_b_sols, rtcb_sols, rm_r_sols, rtcr_sols, rdrtca_sols, rtrtcb_sols, ribo_h_sols, ribo_d_sols, ribo_t_sols
end
