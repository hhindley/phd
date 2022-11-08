function minimalRibofinal!(dydt, init, parameters, t)
    rho, Kp, thetar, kq, thetax, wr0, wq0, we0 = parameters
    a, em, et, mm, mq, mr, mt, eq, er, rm1m, rm1q, rm1r, rm1t, rm2m, rm2q, rm2r, rm2t, si, zmm, zmq, zmr, zmt = init
    
    #define change arrays
    da, dem, det, dmm, dmq, dmr, dmt, deq, der, drm1m, drm1q, drm1r, drm1t, drm2m, drm2q, drm2r, drm2t, dsi, dzmm, dzmq, dzmr, dzmt = zeros(length(dydt))
    
    
    # import
    # s0 -> si; following a Michaelis-Menten dynamic: et*vt*s0/(Kt+s0)
    imp = (et*vt*s0/(Kt + s0))
    
    dsi = dsi + imp
    
    # metabolism
    # si -> a; following a Michaelis-Menten dynamic: em*vm*si/(Km + si)
    nucat= em*vm*si/(Km + si)
    
    da = da +nucat * ns
    dsi = dsi -nucat
        
    # transcription 
    # 0 -> rx; following a Michaelis-Menten dynamic: x*a/(theta+a)
    function transcriptionMM(wx,theta) wx*a/(theta+a) end
    
    dmm = dmm +transcriptionMM(we0,thetax)
    dmr = dmr +transcriptionMM(wr0,thetar)
    dmt = dmt +transcriptionMM(we0,thetax)
    
    # q-compartment inhibition
    # In the case of compartment q, included autoinhibition 1+(q/kq^nq)
    qMM = 1+(eq/kq)^nq

    dmq = dmq + transcriptionMM(wq0,thetax)/qMM   
    
    # translation
    
    ## binding
    # rx + r <-> rm1x; at forward rate: kb*rx*r
    #                  at backward rate: ku*rm1x

    function bind(mx) mx*er*kb end
    function unbind(rm1x) rm1x*ku end

    dmm = dmm + unbind(rm1m) -bind(mm)
    dmq = dmq + unbind(rm1q) -bind(mq)
    dmr = dmr + unbind(rm1r) -bind(mr)
    dmt = dmt + unbind(rm1t) -bind(mt)
    der = der + unbind(rm1r+rm1t+rm1m+rm1q) -bind(mr+mt+mm+mq)
    drm1m = drm1m -unbind(rm1m) +bind(mm)
    drm1q = drm1q -unbind(rm1q) +bind(mq)
    drm1r = drm1r -unbind(rm1r) +bind(mr)
    drm1t = drm1t -unbind(rm1t) +bind(mt)
    
    ## translation initiation
    # rm1x -> rm2x + mx; at forward rate: rho*rm1x
    function transInit(riborna) rho*riborna end

    dmm = dmm +transInit(rm1m)
    dmq = dmq +transInit(rm1q)
    dmr = dmr +transInit(rm1r)
    dmt = dmt +transInit(rm1t)
    drm1m = drm1m -transInit(rm1m)
    drm1q = drm1q -transInit(rm1q)
    drm1r = drm1r -transInit(rm1r)
    drm1t = drm1t -transInit(rm1t)
    drm2m = drm2m +transInit(rm1m)
    drm2q = drm2q +transInit(rm1q)
    drm2r = drm2r +transInit(rm1r)
    drm2t = drm2t +transInit(rm1t)
    
    ## main part and release
    # rm2x + nx * a ->  x + r; at forward rate: gmax*a/(Kg + a)
    
    Kg= gmax/Kp
    gamma= gmax*a/(Kg + a)

    function transRate(transRibo,protLen) transRibo*gamma/protLen end
    
    da = da -gamma*(rm2q + rm2r + rm2t + rm2m)
    dem = dem +  transRate(rm2m,nx)
    deq = deq +  transRate(rm2q,nx)
    der = der +2*transRate(rm2r,nr) +transRate(rm2t+rm2m+rm2q,nx)
    det = det +  transRate(rm2t,nx)
    drm2m = drm2m -transRate(rm2m,nx)
    drm2q = drm2q -transRate(rm2q,nx)
    drm2r = drm2r -transRate(rm2r,nr)
    drm2t = drm2t -transRate(rm2t,nx)

    # chloramphenicol
    # rmx <-> zmx ; at a forward rate cl*k_cm and backward b=0
    f = cl*k_cm
    
    function inhibit(rmx) f*rmx end
    function uninhibit(zmx) b*zmx end
    
    drm2m = drm2m -inhibit(rm2m) +uninhibit(zmm)
    drm2q = drm2q -inhibit(rm2q) +uninhibit(zmq)
    drm2r = drm2r -inhibit(rm2r) +uninhibit(zmr)
    drm2t = drm2t -inhibit(rm2t) +uninhibit(zmt)
    dzmm = dzmm +inhibit(rm2m) -uninhibit(zmm)
    dzmq = dzmq +inhibit(rm2q) -uninhibit(zmq)
    dzmr = dzmr +inhibit(rm2r) -uninhibit(zmr)
    dzmt = dzmt +inhibit(rm2t) -uninhibit(zmt)
    
    # degradation
    # mx -> 0; at a rate dm*mx
    function degrade(mx) dm*mx end
    dmm = dmm -degrade(mm)
    dmq = dmq -degrade(mq)
    dmr = dmr -degrade(mr)
    dmt = dmt -degrade(mt)
    
    # growth rate
    ttrate= (rm2q + rm2r + rm2t + rm2m)*gamma
    lam= ttrate/aatot

    #dilution
    # x -> 0; at a rate lam*x
    function dilute(x) lam*x end

    da = da -dilute(a)
    dem = dem -dilute(em)
    det = det -dilute(et)
    dmm = dmm -dilute(mm)
    dmq = dmq -dilute(mq)
    dmr = dmr -dilute(mr)
    dmt = dmt -dilute(mt)
    deq = deq -dilute(eq)
    der = der -dilute(er)
    drm1m = drm1m -dilute(rm1m)
    drm1q = drm1q -dilute(rm1q)
    drm1r = drm1r -dilute(rm1r)
    drm1t = drm1t -dilute(rm1t)
    drm2m = drm2m -dilute(rm2m)
    drm2q = drm2q -dilute(rm2q)
    drm2r = drm2r -dilute(rm2r)
    drm2t = drm2t -dilute(rm2t)
    dsi = dsi -dilute(si)
    dzmm = dzmm - dilute(dzmm)
    dzmq = dzmq - dilute(dzmq)
    dzmr = dzmr - dilute(dzmr)
    dzmt = dzmt - dilute(dzmt)
    
    dydt[1] = da
    dydt[2] = dem
    dydt[3] = det
    dydt[4] = dmm
    dydt[5] = dmq
    dydt[6] = dmr
    dydt[7] = dmt
    dydt[8] = deq
    dydt[9] = der
    dydt[10] = drm1m
    dydt[11] = drm1q
    dydt[12] = drm1r
    dydt[13] = drm1t
    dydt[14] = drm2m
    dydt[15] = drm2q
    dydt[16] = drm2r
    dydt[17] = drm2t
    dydt[18] = dsi
    dydt[19] = zmm
    dydt[20] = zmq
    dydt[21] = zmr
    dydt[22] = zmt
end 

const dm= 0.1
const kb= 1.
const ku= 1.0
rho= 0.5

Kp= 7.0 #parameterised
thetar= 426.87 # parameterised
kq= 152219. # parameterised
thetax= 4.38 # parameterised
wr0= 93. # parameterised
wq0= 949. # parameterised
we0= 4.38 # parameterised
const gmax= 1260.0
const aatot= 1.0e8
const vt= 726.0
const Kt= 1.0e3
const s0= 1.0e4
const vm= 5800.0
const Km= 1.0e3
global ns= 0.5
const nq= 4.
const nr= 7549.0
const nx= 300.0
const b = 0.
const cl= 0.
const k_cm = 0.00554752/662.435565 # parameterised

rates = [dm, kb, ku, rho]
parameters=  [Kp, thetar, kq, thetax, wr0, wq0, we0, gmax, aatot, vt, Kt, s0, vm, Km, ns, nq, nr, nx]
minimalParameters = [rho, Kp, thetar, kq, thetax, wr0, wq0, we0]

# define initial conditions
a_0= 1000.0
em_0= 0.
et_0= 0.
mm_0= 0.
mq_0= 0.
mr_0= 0.
mt_0= 0.
eq_0= 0.
er_0= 10.0
rm1m_0= 0.
rm1t_0= 0.
rm1q_0= 0.
rm1r_0= 0.
rm2m_0= 0.
rm2t_0= 0.
rm2q_0= 0.
rm2r_0= 0.
si_0= 0.
zmm_0= 0.
zmq_0= 0.
zmr_0= 0.
zmt_0= 0.

init= [a_0, em_0, et_0, mm_0, mq_0, mr_0, mt_0, eq_0, er_0, rm1m_0, rm1q_0, rm1r_0, rm1t_0, rm2m_0, rm2q_0, rm2r_0, rm2t_0, si_0, zmm_0, zmq_0, zmr_0, zmt_0]

using DifferentialEquations, DataFrames, Plots, DiffEqCallbacks, Random, BlackBoxOptim, Sundials, Dates

include("model.jl")
include("parameters.jl")
include("functions.jl")

species = [:a, :em, :et, :mm, :mq, :mr, :mt, :eq, :er, :rm1m, :rm1q, :rm1r, :rm1t, :rm2m, :rm2q, :rm2r, :rm2t, :si, :zmm, :zmq, :zmr, :zmt]
tspan = (0.,1e5)

cb = TerminateSteadyState(1e-6, 1e-4)
prob = ODEProblem(ribofinal!,init,tspan, (rates,parameters))
using BenchmarkTools
# Seems to be considerably slower (2x), so I'll manually solve using Rosenbrock
#@btime sol = solve(prob, alg=DynamicSS(Rosenbrock23(), abstol=1e-6, reltol=1e-4, tspan=Inf), save_start=false, save_everystep=false)
@btime sol = solve(prob, alg=Rosenbrock23(), abstol=1e-4, save_start=false, save_everystep=false; callback=cb)
