using DifferentialEquations, PlotlyJS, DataFrames, Measures

include("model.jl")
include("parameters.jl")
include("functions.jl")

tspan = (0,1e9)

solu = simple_solve!(odemodel!, init, tspan, params)

plotly_plot_sol(solu, "", "")

species = [:cr, :em, :cp, :cq, :ct, :et, :cm, :mt, :mm, :q, :p, :si, :mq, :mp, :mr, :r, :a]
initsolDF = DataFrame([[j[i] for j in solu.u] for i=1:length(solu.u[1])], species)
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

ssinit = [sscr, ssem, sscp, sscq, ssct, sset, sscm, ssmt, ssmm, ssq, ssp, sssi, ssmq, ssmp, ssmr, ssr, ssa]

solu_ss = simple_solve!(odemodel!, ssinit, tspan, params)

plotly_plot_sol(solu_ss, "", "")