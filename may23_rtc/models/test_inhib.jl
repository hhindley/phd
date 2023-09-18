using Parameters, CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, Statistics
using Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra

using PlotlyJS, ProgressBars

include("/home/holliehindley/phd/may23_rtc/functions/solving.jl"); include("/home/holliehindley/phd/may23_rtc/functions/set_ups.jl"); include("/home/holliehindley/phd/may23_rtc/functions/plotting.jl"); 
include("/home/holliehindley/phd/may23_rtc/functions/sweep_params.jl"); include("/home/holliehindley/phd/may23_rtc/models/rtc_orig.jl"); include("/home/holliehindley/phd/may23_rtc/models/atp_lam_kin_t.jl"); 
include("/home/holliehindley/phd/may23_rtc/models/single_t.jl"); include("/home/holliehindley/phd/may23_rtc/models/combinations_t.jl"); 
include("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/bf_funcs.jl");

include("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/init_switch/funcs.jl");


# RtcB + i <=> rtcb_i

function test_mod(initial, params, t) 
    k = params
    rtcb, i, rtcb_i = initial

    # ODEs
    drtcb = -rtcb*i*k + rtcb_i*k
    di = -rtcb*i*k + rtcb_i*k
    drtcb_i = rtcb*i*k - rtcb_i*k

    [drtcb, di, drtcb_i]

end

init1 = [10, 1, 0]
p1 = 1
tspan = (0,1)
solu = sol(test_mod, init1, tspan, p1)
sols = DataFrame(solu)
rename!(sols, [:time, :rtcb, :i, :rtcb_i])

plot([scatter(x=solu.t, y=sols.rtcb), scatter(x=solu.t, y=sols.rtcb_i)])


plot([scatter(x=solu.t, y=sols.rtcb), scatter(x=solu.t, y=sols.i), scatter(x=solu.t, y=sols.rtcb_i)])

rtcb_tot = @. sols.rtcb+sols.rtcb_i

y = @. sols.rtcb_i/rtcb_tot

plot(scatter(x=solu.t, y=y))
plot(scatter(x=sols.i, y=y))


function test_mod_2(initial, params, t) 
    k, i = params
    rtcb = initial

    # rtcb_i = i*k*rtcb/((k*i)+k)


    # ODEs
    drtcb = -rtcb*i*k + rtcb_i*k
    # di = -rtcb*i*k + rtcb_i*k
    # drtcb_i = rtcb*i*k - rtcb_i*k

    drtcb

end

init2 = 10
k=1;i=1;
p2 = [k,i]
tspan = (0,1)
solu2 = sol(test_mod_2, init2, tspan, p2)
sols2 = DataFrame(solu2)
rename!(sols2, [:time, :rtcb])

rtcb = sols2.rtcb
rtcb_i = i*k*rtcb/((k*i)+k)
rtcb_i = rtcb*i
plot([scatter(x=solu2.t, y=sols2.rtcb), scatter(x=solu2.t, y=rtcb_i)])









function test_mod_2(initial, params, t) 
    k, i = params
    rtcb, rtcb_i = initial

    # rtcb_i = i*k*rtcb/((k*i)+k)


    # ODEs
    drtcb = -rtcb*i*k + rtcb_i*k
    # di = -rtcb*i*k + rtcb_i*k
    drtcb_i = rtcb*i*k - rtcb_i*k

    [drtcb, drtcb_i]

end

init2 = [10,1]
k=1;i=1;
p2 = [k,i]
tspan = (0,1)
solu2 = sol(test_mod_2, init2, tspan, p2)
sols2 = DataFrame(solu2)
rename!(sols2, [:time, :rtcb, :rtcb_i])

rtcb = sols2.rtcb
rtcb_i = i*k*rtcb/((k*i)+k)
rtcb_i = rtcb*i
plot([scatter(x=solu2.t, y=sols2.rtcb), scatter(x=solu2.t, y=sols.rtcb_i)])

r_tot = sols.rtcb+sols.rtcb_i

y = @. sols.rtcb_i/r_tot

plot(scatter(x=sols.rtcb, y=y))