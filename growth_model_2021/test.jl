using DifferentialEquations, Plots, DataFrames, BenchmarkTools

include("model.jl")
include("parameters.jl")
include("functions.jl")


tspan = (0.0, 1e9)


species = [:cr, :em, :cp, :cq, :ct, :et, :cm, :mt, :mm, :q, :p, :si, :mq, :mp, :mr, :r, :a, :zmr, :zmp, :zmt, :zmm, :zmq]
prob = ODEProblem(odemodelfull!,initfull,tspan,params_init)
sol = solve(prob, alg=Rosenbrock23())
solDF = DataFrame([[j[i] for j in sol.u] for i=1:length(sol.u[1])], species)
# a = solDF[!, :zmp]
plot(sol)



# speciesn = [:cr, :em, :cp, :cq, :ct, :et, :cm, :mt, :mm, :q, :p, :si, :mq, :mp, :mr, :r, :a]

# probn = ODEProblem(odemodel!,init,tspan,params)
# soln = solve(probn, alg=Rosenbrock23())
# solDF = DataFrame([[j[i] for j in soln.u] for i=1:length(soln.u[1])], speciesn)
# # a = solDF[!, :zmp]
# plot(soln)