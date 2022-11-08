using DifferentialEquations, Plots, DataFrames, Measures

include("model.jl")
include("parameters.jl")
include("functions.jl")

tspan = (0,1e9)

sol = simple_solve!(odemodel!, init, tspan, params)

labels = ["cr" "em" "cp" "cq" "ct" "et" "cm" "mt" "mm" "q" "p" "si" "mq" "mp" "mr" "r" "a"]
plot(sol, labels=labels, margin=5mm)

species = [:cr, :em, :cp, :cq, :ct, :et, :cm, :mt, :mm, :q, :p, :si, :mq, :mp, :mr, :r, :a]
solDF = DataFrame([[j[i] for j in sol.u] for i=1:length(sol.u[1])], species)
a = solDF[!, :a]
cq = solDF[!, :cq]
cr = solDF[!, :cr]
cp = solDF[!, :cp]
ct = solDF[!, :ct]
cm = solDF[!, :cm]
si = solDF[!, :si]
em = solDF[!, :em]
et = solDF[!, :et]
mt = solDF[!, :mt]
mm = solDF[!, :mm]
mq = solDF[!, :mq]
mp = solDF[!, :mp]
mr = solDF[!, :mr]
r = solDF[!, :r]
p = solDF[!, :p]
q = solDF[!, :q]

plot(sol.t, r)

plot(sol.t, cp)