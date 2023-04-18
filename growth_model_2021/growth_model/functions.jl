function simple_solve!(model, init, tspan, params)
    prob = ODEProblem(model, init, tspan, params);
    solu = solve(prob, Rodas4(), isoutofdomain=(y,p,t)->any(x->x<0,y))#, abstol=1e-15, reltol=1e-12);
    return solu
end

function get_ss_init(solu)
    species = [:cr, :em, :cp, :cq, :ct, :et, :cm, :mt, :mm, :q, :p, :si, :mq, :mp, :mr, :r, :a]
    initsolDF = DataFrame([[j[i] for j in solu.u] for i=1:length(solu.u[1])], species)
    sscr = initsolDF[end,:][:cr]
    ssem = initsolDF[end,:][:em]
    sscp = initsolDF[end,:][:cp]
    sscq = initsolDF[end,:][:cq]
    ssct = initsolDF[end,:][:ct]
    sscm = initsolDF[end,:][:cm]
    ssmt = initsolDF[end,:][:mt]
    ssmm = initsolDF[end,:][:mm]
    ssq = initsolDF[end,:][:q]
    ssp = initsolDF[end,:][:p]
    sssi = initsolDF[end,:][:si]
    ssmq = initsolDF[end,:][:mq]
    ssmp = initsolDF[end,:][:mp]
    ssmr = initsolDF[end,:][:mr]
    ssr = initsolDF[end,:][:r]
    ssa = initsolDF[end,:][:a]
    # sset = initsolDF[end,:][:et]
    # sss = initsolDF[end,:][:s]
    # ssN = initsolDF[end,:][:N]

    sset = 0.00001#0.1
    return @LArray [sscr, ssem, sscp, sscq, ssct, sset, sscm, ssmt, ssmm, ssq, ssp, sssi, ssmq, ssmp, ssmr, ssr, ssa, s_0, N_0] (:sscr, :ssem, :sscp, :sscq, :ssct, :sset, :sscm, :ssmt, :ssmm, :ssq, :ssp, :sssi, :ssmq, :ssmp, :ssmr, :ssr, :ssa, :s_0, :N_0)
end

function SSinitVals()
    include("initial.jl")
    initsol = initial_model(odemodelfull!, initfull, tspan, params)
    species = [:cr, :em, :cp, :cq, :ct, :et, :cm, :mt, :mm, :q, :p, :si, :mq, :mp, :mr, :r, :a, :zmr, :zmp, :zmq, :zmt, :zmm]
    initsolDF = DataFrame([[j[i] for j in initsol.u] for i=1:length(initsol.u[1])], species)
    sscr = initsolDF[end,:][:cr]
    ssem = initsolDF[end,:][:em]
    sscp = initsolDF[end,:][:cp]
    sscq = initsolDF[end,:][:cq]
    ssct = initsolDF[end,:][:ct]
    sset = initsolDF[end,:][:et]
    sscm = initsolDF[end,:][:cm]
    ssmt = initsolDF[end,:][:mt]
    ssmm = initsolDF[end,:][:mm]
    ssq = initsolDF[end,:][:q]
    ssp = initsolDF[end,:][:p]
    sssi = initsolDF[end,:][:si]
    ssmq = initsolDF[end,:][:mq]
    ssmp = initsolDF[end,:][:mp]
    ssmr = initsolDF[end,:][:mr]
    ssr = initsolDF[end,:][:r]
    ssa = initsolDF[end,:][:a]
    sszmr = initsolDF[end,:][:zmr]
    sszmp = initsolDF[end,:][:zmp]
    sszmq = initsolDF[end,:][:zmq]
    sszmt = initsolDF[end,:][:zmt]
    sszmm = initsolDF[end,:][:zmm]
    ssinit = [sscr, ssem, sscp, sscq, ssct, sset, sscm, ssmt, ssmm, ssq, ssp, sssi, ssmq, ssmp, ssmr, ssr, ssa, sszmr, sszmp, sszmq, sszmt, sszmm]
    return ssinit
end

function calcGrowthrate(systemState)
    Kgamma = gmax/Kp
    a = systemState[:a]
    gamma = gmax*a/(Kgamma+a)
    ttrate = sum(systemState[[:cq, :cr, :cp, :ct, :cm]])*gamma
    return ttrate/M
end

function calcRMF(systemState)
	rmf = nr*sum(systemState[[:r, :cr, :cp, :ct, :cm, :cq, :zmr, :zmp, :zmt, :zmm, :zmq]]) / nr*sum(systemState[[:r, :cr, :cp, :ct, :cm, :cq, :zmr, :zmp, :zmt, :zmm, :zmq]]) + nx*sum(systemState[[:p, :q, :et, :em]])
    return rmf
end

function get_curve(sol, species)
    df = DataFrame(sol)
    if length(sol[1]) == 17
        rename!(df, [:time, :cr, :em, :cp, :cq, :ct, :et, :cm, :mt, :mm, :q, :p, :si, :mq, :mp, :mr, :r, :a])
    else
        rename!(df, [:time, :cr, :em, :cp, :cq, :ct, :et, :cm, :mt, :mm, :q, :p, :si, :mq, :mp, :mr, :r, :a, :s, :N])
    end 
    species = df[:, species]
    return species
end

function plotly_plot_sol(sol, log, log1)
    cr = get_curve(sol, :cr); em = get_curve(sol, :em); cp = get_curve(sol, :cp); cq = get_curve(sol, :cq); 
    ct = get_curve(sol, :ct); et = get_curve(sol, :et); cm = get_curve(sol, :cm); mt = get_curve(sol, :mt); 
    mm = get_curve(sol, :mm); q = get_curve(sol, :q); p = get_curve(sol, :p); si = get_curve(sol, :si);
    mq = get_curve(sol, :mq); mp = get_curve(sol, :mp); mr = get_curve(sol, :mr); r = get_curve(sol, :r); 
    a = get_curve(sol, :a)

    cr_curve = scatter(x=sol.t, y=cr, name="cr")
    em_curve = scatter(x=sol.t, y=em, name="em")
    cp_curve = scatter(x=sol.t, y=cp, name="cp")
    cq_curve = scatter(x=sol.t, y=cq, name="cq")
    ct_curve = scatter(x=sol.t, y=ct, name="ct")
    et_curve = scatter(x=sol.t, y=et, name="et")
    cm_curve = scatter(x=sol.t, y=cm, name="cm")
    mt_curve = scatter(x=sol.t, y=mt, name="mt")
    mm_curve = scatter(x=sol.t, y=mm, name="mm")
    q_curve = scatter(x=sol.t, y=q, name="q")
    p_curve = scatter(x=sol.t, y=p, name="p")
    si_curve = scatter(x=sol.t, y=si, name="si")
    mq_curve = scatter(x=sol.t, y=mq, name="mq")
    mp_curve = scatter(x=sol.t, y=mp, name="mp")
    mr_curve = scatter(x=sol.t, y=mr, name="mr")
    r_curve = scatter(x=sol.t, y=r, name="r")
    a_curve = scatter(x=sol.t, y=a, name="a")
    

    return (plot([cr_curve, em_curve, cp_curve, cq_curve, ct_curve, et_curve, cm_curve, mt_curve, mm_curve, q_curve, p_curve, si_curve, mq_curve, mp_curve, mr_curve, r_curve, a_curve] ,Layout(xaxis_type=log, yaxis_type=log1)))
end

function plotly_plot_sol_pop(sol, log, log1)
    cr = get_curve(sol, :cr); em = get_curve(sol, :em); cp = get_curve(sol, :cp); cq = get_curve(sol, :cq); 
    ct = get_curve(sol, :ct); et = get_curve(sol, :et); cm = get_curve(sol, :cm); mt = get_curve(sol, :mt); 
    mm = get_curve(sol, :mm); q = get_curve(sol, :q); p = get_curve(sol, :p); si = get_curve(sol, :si);
    mq = get_curve(sol, :mq); mp = get_curve(sol, :mp); mr = get_curve(sol, :mr); r = get_curve(sol, :r); 
    a = get_curve(sol, :a); s = get_curve(sol, :s); N = get_curve(sol, :N);

    cr_curve = scatter(x=sol.t, y=cr, name="cr")
    em_curve = scatter(x=sol.t, y=em, name="em")
    cp_curve = scatter(x=sol.t, y=cp, name="cp")
    cq_curve = scatter(x=sol.t, y=cq, name="cq")
    ct_curve = scatter(x=sol.t, y=ct, name="ct")
    et_curve = scatter(x=sol.t, y=et, name="et")
    cm_curve = scatter(x=sol.t, y=cm, name="cm")
    mt_curve = scatter(x=sol.t, y=mt, name="mt")
    mm_curve = scatter(x=sol.t, y=mm, name="mm")
    q_curve = scatter(x=sol.t, y=q, name="q")
    p_curve = scatter(x=sol.t, y=p, name="p")
    si_curve = scatter(x=sol.t, y=si, name="si")
    mq_curve = scatter(x=sol.t, y=mq, name="mq")
    mp_curve = scatter(x=sol.t, y=mp, name="mp")
    mr_curve = scatter(x=sol.t, y=mr, name="mr")
    r_curve = scatter(x=sol.t, y=r, name="r")
    a_curve = scatter(x=sol.t, y=a, name="a")
    s_curve = scatter(x=sol.t, y=s, name="s")
    N_curve = scatter(x=sol.t, y=N, name="N")


    return (plot([cr_curve, em_curve, cp_curve, cq_curve, ct_curve, et_curve, cm_curve, mt_curve, mm_curve, q_curve, p_curve, si_curve, mq_curve, mp_curve, mr_curve, r_curve, a_curve, s_curve, N_curve] ,Layout(xaxis_type=log, yaxis_type=log1)))
end