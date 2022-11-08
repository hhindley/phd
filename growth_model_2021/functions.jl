function SSinitVals(final_time)
    initsol = initial_model(odemodelfull!, initfull, params_init, final_time)
    species = [:cr, :em, :cp, :cq, :ct, :et, :cm, :zmr, :zmp, :zmq, :zmt, :zmm, :mt, :mm, :q, :p, :si, :mq, :mp, :mr, :r, :a]
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
    ssinit = [sscr, ssem, sscp, sscq, ssct, sset, sscm, sszmr, sszmp, sszmq, sszmt, sszmm, ssmt, ssmm, ssq, ssp, sssi, ssmq, ssmp, ssmr, ssr, ssa]
    return ssinit
end

function calcGrowthrate(systemState, Kgamma)
    #Kgamma = gmax/Kp
    a = systemState[:a]
    gamma = gmax*a/(Kgamma+a)
    ttrate = sum(systemState[[:cq, :cr, :cp, :ct, :cm]])*gamma
    return ttrate/M
end

function calcRMF(systemState)
	rmf = nr*sum(systemState[[:r, :cr, :cp, :ct, :cm, :cq, :zmr, :zmp, :zmt, :zmm, :zmq]]) / (nr*sum(systemState[[:r, :cr, :cp, :ct, :cm, :cq, :zmr, :zmp, :zmt, :zmm, :zmq]]) + nx*sum(systemState[[:p, :q, :et, :em]]))
    return rmf
end

function calcNucat(systemState)
    em = systemState[:em]
    si = systemState[:si]
    nucat = (em*si*vm/(Km+si))
    return nucat
end

function calcVimp(systemState)
    et = systemState[:et]
    vimp = (et*vt*s0/(Kt+s0))
    return vimp
end

function initial_model(model, initvals, params, final_time)
    tspan = (0., final_time)

    prob = ODEProblem(model, initvals, tspan, params)
    sol = solve(prob, alg=Rosenbrock23())
end
