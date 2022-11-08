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
SSinitVals
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
