# transcription analysis 
function vary_all(rt, kr, n, L, c, rtcr, k1, k2, k3, k4, k5, sigma, atp)
    alpha = rt/kr;
    fa = (1+alpha)^n/(L*((1+c*alpha)^n)+(1+alpha)^n);
    ra = fa*rtcr; 
    v = ra*k1*sigma*k3*atp*k5/(k1*sigma*k3*atp+k1*sigma*(k4+k5)+k2*(k4+k5)+k3*atp*k5);
    return v
end

function vary_all_transcription(rt, kr, n, L, c, rtcr, k1, k2, k3, k4, k5, sigma, atp, w_rtc, theta_rtc)
    alpha = rt/kr;
    fa = (1+alpha)^n/(L*((1+c*alpha)^n)+(1+alpha)^n);
    ra = fa*rtcr; 
    v = ra*k1*sigma*k3*atp*k5/(k1*sigma*k3*atp+k1*sigma*(k4+k5)+k2*(k4+k5)+k3*atp*k5)
    tscr = v*(w_rtc*atp/(theta_rtc+atp))
    return tscr
end

function vary_all_fraction_active_rtcr(rt, kr, L, c)
    alpha = rt/kr;
    n = 6
    fa = (1+alpha)^n/(L*((1+c*alpha)^n)+(1+alpha)^n);
    return fa
end 

# michaelis menten - gives same results as partition analysis (slightly different numbers, same trend)
function mm_rt(rt_range)
    oc = [];
    for i in rt_range
        alpha = i/kr;
        fa = (1+alpha)^n/(L*((1+c*alpha)^n)+(1+alpha)^n);
        ra = fa*rtcr; 
        # v = ra*k1*sigma*k3*atp*k5/(k1*sigma*k3*atp+k1*sigma*(k4+k5)+k2*(k4+k5)+k3*atp*k5)
        v = k3*atp*k1*sigma*ra/(k1*sigma+k2+k3*atp)
        push!(oc, v)
    end 
    return oc
end

function mm_rtcr(rtcr_range)
    oc = [];
    for i in rtcr_range
        alpha = rt/kr;
        fa = (1+alpha)^n/(L*((1+c*alpha)^n)+(1+alpha)^n);
        ra = fa*i; 
        v = k3*atp*k1*sigma*ra/(k1*sigma+k2+k3*atp)
        push!(oc, v)
    end 
    return oc
end


# plotting functions
function plot_rtc_comp1(func, rtc_range, title, xlabel, ylabel)
    return plot(rtcr_range, func(rtc_range), title=title, xlabel=xlabel, ylabel=ylabel, labels=false);
end


function plotting1(rtc_range, y, title, xlabel, ylabel)
    return plot(rtc_range, y, title=title, xlabel=xlabel, ylabel=ylabel, labels=false)
end

function plotting2(rtc_range, y, title, xlabel, ylabel, labels)
    y = reshape(y, (100,5));
    y = [d[:] for d in eachcol(y)];
    return plot(rtc_range, y, title=title, xlabel=xlabel, ylabel=ylabel, labels=labels)
end 

# model functions - translation analysis 
function rtc_expression!(dydt, initial, params, t)

    rt, kr, L, c, rtcr, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, r, gr, d = params

    rm, rtc = initial 
    drmdt, drtcdt = zeros(length(dydt))

    n = 6
    alpha = rt/kr;
    fa = (1+alpha)^n/(L*((1+c*alpha)^n)+(1+alpha)^n); # fraction of active RtcR
    ra = fa*rtcr; # amount of active RtcR
    v = ra*k1*sigma*k3*atp*k5/(k1*sigma*k3*atp+k1*sigma*(k4+k5)+k2*(k4+k5)+k3*atp*k5); # rate of open complex formation
    tscr = v*(w_rtc*atp/(theta_rtc+atp)); # transcription rate 
    tlel = max*atp/(thr+atp) # translation elongation rate
    tlr = r*rm*tlel # translation rate so formation of Rtc 

    dydt[1] = tscr - gr*rm - d*rm # rm = rtc mRNA 
    dydt[2] = tlr - gr*rtc # - 0.1*rtc # rtc = rtcA and rtcB 

end

function rtc_expression_with_rtcr!(dydt, initial, params, t)

    rt, kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, r, gr, d = params

    rm, rtc, rm_r, rtcr = initial 
    drmdt, drtcdt, drm_rdt, drtcrdt = zeros(length(dydt))

    n = 6
    alpha = rt/kr;
    fa = (1+alpha)^n/(L*((1+c*alpha)^n)+(1+alpha)^n); # fraction of active RtcR
    ra = fa*rtcr; # amount of active RtcR
    v = ra*k1*sigma*k3*atp*k5/(k1*sigma*k3*atp+k1*sigma*(k4+k5)+k2*(k4+k5)+k3*atp*k5); # rate of open complex formation
    tscr = v*(w_rtc*atp/(theta_rtc+atp)); # transcription rate 
    tlel = max*atp/(thr+atp) # translation elongation rate
    tlr = r*rm*tlel # translation rate so formation of Rtc 

    # rtcr 
    rtcr_tscr = w_rtc*atp/(theta_rtc+atp)
    rtcr_tlel = max*atp/(thr+atp)
    rtcr_tlr = r*rm_r*rtcr_tlel

    dydt[1] = tscr - gr*rm - d*rm # rm = rtc mRNA 
    dydt[2] = tlr - gr*rtc # - 0.1*rtc # rtc = rtcA and rtcB 
    dydt[3] = rtcr_tscr - gr*rm_r - d*rm_r
    dydt[4] = rtcr_tlr - gr*rtcr #- ra*rtcr

end

function simple_solve!(model, init, tspan, params)
    prob = ODEProblem(model, init, tspan, params);
    sol = solve(prob, alg=Rodas4())#, abstol=1e-9, reltol=1e-6);
    return sol
end

# vary parameters for model 
function vary_param_for_rtcab_model_plot!(x_range, param, model, init, tspan, params, xlabel, ylabel)
    res_rtc = [];
    res_rm = [];
    for x in x_range
        params[param] = x
        sol = simple_solve!(model, init, tspan, params)
        solDF = DataFrame([[j[i] for j in sol.u] for i=1:length(sol.u[1])], species);
        rm = solDF[end, :rm];
        rtc = solDF[end, :rtc];
        push!(res_rtc, rtc)
        push!(res_rm, rm)
    end
    # plot(x_range, res_rtc, xlabel=xlabel, ylabel=ylabel, labels=false, color=:blue)
    # plot(x_range, res_rm, xlabel=xlabel, ylabel=ylabel, labels=false, color=:orange)
    plot(x_range, [res_rtc, res_rm], xlabel=xlabel, ylabel=ylabel, labels=["rtc" "rm"])
end

function vary_param_for_rtcRab_plot!(x_range, param, model, init, tspan, params, species, xlabel, ylabel)
    res_rtc = [];
    res_rm = [];
    res_rm_r = [];
    res_rtcr = [];
    for x in x_range
        params[param] = x
        sol = simple_solve!(model, init, tspan, params)
        solDF = DataFrame([[j[i] for j in sol.u] for i=1:length(sol.u[1])], species);
        rm = solDF[end, :rm];
        rtc = solDF[end, :rtc];
        rm_r = solDF[end, :rm_r];
        rtcr = solDF[end, :rtcr];
        push!(res_rtc, rtc)
        push!(res_rm, rm)
        push!(res_rm_r, rm_r)
        push!(res_rtcr, rtcr)
    end
    plot(x_range, [res_rtc, res_rm, res_rm_r, res_rtcr], xlabel=xlabel, ylabel=ylabel, labels=["RtcAB" "RtcAB-mRNA" "RtcR-mRNA" "RtcR"])
end

function vary_param_for_rtcRab_model!(long_range, param, model, init, tspan, params, species)
    res_rtc = [];
    res_rm = [];
    res_rm_r = [];
    res_rtcr = [];
    for x in long_range
        params[param] = x
        sol = simple_solve!(model, init, tspan, params)
        solDF = DataFrame([[j[i] for j in sol.u] for i=1:length(sol.u[1])], species);
        rm = solDF[end, :rm];
        rtc = solDF[end, :rtc];
        rm_r = solDF[end, :rm_r];
        rtcr = solDF[end, :rtcr];
        push!(res_rtc, rtc)
        push!(res_rm, rm)
        push!(res_rm_r, rm_r)
        push!(res_rtcr, rtcr)
    end
    return res_rm, res_rtc, res_rm_r, res_rtcr
end

function vary_parameter_total_expression!(range, param, model, init, tspan, params, species)
    res_rtc = [];
    res_rm = [];
    res_rm_r = [];
    res_rtcr = [];
    solutions = [];
    solutions_t = [];
    for x in range
        params[param] = x
        sol = simple_solve!(model, init, tspan, params)
        solDF = DataFrame([[j[i] for j in sol.u] for i=1:length(sol.u[1])], species);
        rm = solDF[:, :rm];
        rtc = solDF[:, :rtc];
        rm_r = solDF[:, :rm_r];
        rtcr = solDF[:, :rtcr];
        push!(res_rtc, rtc)
        push!(res_rm, rm)
        push!(res_rm_r, rm_r)
        push!(res_rtcr, rtcr)
        push!(solutions_t, sol.t)
        push!(solutions, sol)
    end
    atp1 = solutions[1]
    atp2 = solutions[2]
    atp3 = solutions[3]
    atp4 = solutions[4]
    atp5 = solutions[5]
    return res_rm, res_rtc, res_rm_r, res_rtcr, atp1, atp2, atp3, atp4, atp5, solutions_t
end

function all_species!(long_range, param, model, init, tspan, params, species)
    res_rtc = [];
    res_rm = [];
    res_rm_r = [];
    res_rtcr = [];
    res_fa = [];
    res_ra = [];
    for x in long_range
        params[param] = x
        sol = simple_solve!(model, init, tspan, params)
        solDF = DataFrame([[j[i] for j in sol.u] for i=1:length(sol.u[1])], species);
        rm = solDF[end, :rm];
        rtc = solDF[end, :rtc];
        rm_r = solDF[end, :rm_r];
        rtcr = solDF[end, :rtcr];
        push!(res_rtc, rtc)
        push!(res_rm, rm)
        push!(res_rm_r, rm_r)
        push!(res_rtcr, rtcr)
        n=6
        alpha = rt/kr;
        fa = (1+alpha)^n/(L*((1+c*alpha)^n)+(1+alpha)^n)
        ra = fa*rtcr
        push!(res_fa, fa)
        push!(res_ra, ra)
    end
    return res_rm, res_rtc, res_rm_r, res_rtcr, res_fa, res_ra
end

function find_intercept!(long_range, x)
    res_rm_ss, res_rtc_ss, res_rm_r_ss, res_rtcr_ss = vary_param_for_rtcRab_model!(long_range, x, rtc_expression_with_rtcr!, init_rtcr, tspan, params_rtcr, species_rtcr);
    spl = Spline1D(long_range, res_rtcr_ss.-res_rtc_ss)
    x0 = roots(spl)
    y0 = Spline1D(long_range, res_rtc_ss)(x0)
    print(x0)
    print(y0)
    return res_rtc_ss, res_rtcr_ss, x0, y0 
end

function log_intercept_plot!(long_range, x, param, param_range)
    res_rtc_ss, res_rtcr_ss, x0, y0 = find_intercept!(long_range, x);
    plot(long_range, [res_rtc_ss, res_rtcr_ss], yaxis=(:log10, (1,Inf)), ylabel="log(conc)", xlabel=("$param"), labels=["RtcAB" "RtcR"], legend=:best, title=("Varying $param from $param_range"))#, annotate=(x0.+2,y0.+10,text("($x0, $y0)", 8)))#, ylim=(0,0.011), xlim=(0,0.0001))
    x0 = round.(x0; sigdigits=3)
    y0 = round.(y0; sigdigits=3)
    scatter!([x0],[y0], markersize=5, color="black", labels="$x0, $y0")
end

function intercept_plot!(long_range, x, param, param_range)
    res_rtc_ss, res_rtcr_ss, x0, y0 = find_intercept!(long_range, x);
    plot(range, [res_rtc_ss, res_rtcr_ss], ylabel="conc", xlabel=("$param"), labels=["RtcAB" "RtcR"], legend=:best, title=("Varying $param from $param_range"))#, annotate=(x0.+2,y0.+10,text("($x0, $y0)", 8)))#, ylim=(0,0.011), xlim=(0,0.0001))
    x0 = round.(x0; sigdigits=3)
    y0 = round.(y0; sigdigits=3)
    scatter!([x0],[y0], markersize=5, color="black", labels="$x0, $y0")
end


function param_explore_plots!(range, x, v_param, v_param_range, long_range)

    mrna_ab, rtc_ab, mrna_r, rtc_r, param1, param2, param3, param4, param5, solutions_t = vary_parameter_total_expression!(range, x, rtc_expression_with_rtcr!, init_rtcr, tspan, params_rtcr, species_rtcr);

    p1 = plot(param1, title="$v_param = $(collect(range)[1])");
    p2 = plot(param2, title="$v_param = $(collect(range)[2])");
    p3 = plot(param3, title="$v_param = $(collect(range)[3])");
    p4 = plot(param4, title="$v_param = $(collect(range)[4])");
    p5 = plot(param5, title="$v_param = $(collect(range)[5])");
    p_vary_params_solutions = plot(p1, p2, p3, p4, p5, layout=(5,1), size=(800, 1000), labels=["RtcAB-mRNA" "RtcAB" "RtcR-mRNA" "RtcR"], margin=2mm);

    labels1 = ["$(collect(range)[1])" "$(collect(range)[2])" "$(collect(range)[3])" "$(collect(range)[4])" "$(collect(range)[5])"]
    mrna_r_1 = plot(solutions_t, [mrna_r[1], mrna_r[2], mrna_r[3], mrna_r[4], mrna_r[5]], title="RtcR-mRNA", xlabel="t", ylabel="conc", labels=labels1);
    mrna_ab_2 = plot(solutions_t, [mrna_ab[1], mrna_ab[2], mrna_ab[3], mrna_ab[4], mrna_ab[5]], title="RtcAB-mRNA", xlabel="t", ylabel="conc", labels=labels1);
    rtc_r_3 = plot(solutions_t, [rtc_r[1], rtc_r[2], rtc_r[3], rtc_r[4], rtc_r[5]], title="RtcR", xlabel="t", ylabel="conc", labels=labels1);
    rtc_ab_4 = plot(solutions_t, [rtc_ab[1], rtc_ab[2], rtc_ab[3], rtc_ab[4], rtc_ab[5]], title="RtcAB", xlabel="t", ylabel="conc", labels=labels1);
    p_vary_params_species = plot(mrna_r_1, mrna_ab_2, rtc_r_3, rtc_ab_4, size=(800,600), plot_title="Different values of $v_param");

    res_rm, res_rtc, res_rm_r, res_rtcr = vary_param_for_rtcRab_model!(long_range, x, rtc_expression_with_rtcr!, init_rtcr, tspan, params_rtcr, species_rtcr);
    p_vary_param_long_sol = plot(long_range, [res_rm, res_rtc, res_rm_r, res_rtcr], xlabel=("$v_param"), ylabel=("conc"), title=("Varying $v_param from $v_param_range"), labels=species_rtcr_labels, legend=:best)
    p_vary_param_long_sol_log = plot(long_range, [res_rm, res_rtc, res_rm_r, res_rtcr], yaxis=(:log10, (1,Inf)), xlabel=("$v_param"), ylabel=("log(conc)"), title=("Varying $v_param from $v_param_range"), labels=species_rtcr_labels, legend=:best)
    
    return p_vary_params_solutions, p_vary_params_species, p_vary_param_long_sol, p_vary_param_long_sol_log
end


