using DifferentialEquations, Plots, DataFrames  

include("model.jl")
include("parameters.jl")

function initial_model(model, initvals, params)
    tspan = (0., 1e9)

    prob = ODEProblem(model, initvals, tspan, params)
    sol = solve(prob, alg_hints=[:stiff])
end
solution = initial_model(odemodel!, init, params_init)  
plot(solution)

# species = [:cr, :em, :cp, :cq, :ct, :et, :cm, :mt, :mm, :q, :p, :si, :mq, :mp, :mr, :r, :a]
# solDF = DataFrame([[j[i] for j in solution.u] for i=1:length(solution.u[1])], species)
# a = solDF[!, :a]
# cq = solDF[!, :cq]
# cr = solDF[!, :cr]
# cp = solDF[!, :cp]
# ct = solDF[!, :ct]
# cm = solDF[!, :cm]
# si = solDF[!, :si]
# em = solDF[!, :em]
# et = solDF[!, :et]
# mt = solDF[!, :mt]
# mm = solDF[!, :mm]
# mq = solDF[!, :mq]
# mp = solDF[!, :mp]
# mr = solDF[!, :mr]
# r = solDF[!, :r]
# p = solDF[!, :p]
# q = solDF[!, :q]
# plot(solution.t, mt)
