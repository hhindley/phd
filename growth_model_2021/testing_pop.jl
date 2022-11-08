using DifferentialEquations, Plots, DataFrames, BenchmarkTools

include("model.jl")
include("parameters.jl")
include("functions.jl")

ssinit = SSinitVals(1e9) # change this to run the initial values for a longer time period, standard is 1e9 

global ns = 100 # 100
pop_params[25] = 100

ssinit[6] = 1

push!(ssinit, 1)
push!(ssinit, 1e10)
# pushfirst!(ssinit, 1e10)
# pushfirst!(ssinit, 1.)
# insert!(ssinit, 8, 1)
# insert!(ssinit, 9, 1e10)

tspan1 = (0.0, 600.0)

species = [:cr, :em, :cp, :cq, :ct, :et, :cm, :zmr, :zmp, :zmt, :zmm, :zmq, :mt, :mm, :q, :p, :si, :mq, :mp, :mr, :r, :a, :N, :s0]

prob = ODEProblem(popodemodelfull!, ssinit, tspan1, pop_params)
sol = solve(prob, alg_hints=[:stiff], abstol=1e-6, reltol=1e-6)#Rosenbrock23()) # alg=Rodas4(), abstol=1e-12, reltol=1e-12)
solDF = DataFrame([[j[i] for j in sol.u] for i=1:length(sol.u[1])], species)
N = solDF[!,:N]
s0 = solDF[!,:s0]
a = solDF[!,:a]
plot(sol.t, N, legend= false, title="N", xlabel="Time (minutes)", ylabel="Number of cells")
plot(sol.t, s0, title="s0")
plot(sol.t, a, title="a")

# lam = calcGrowthrate(solDF[end,:], 3e9)
# vimp = calcVimp(solDF[end,:])

# pop_init = [0, 1e10]
# pparam = [vimp, lam]
# species1 = [:N, :s0]

# prob1 = ODEProblem(pop_odes!, pop_init, tspan1, pparam)
# sol1 = solve(prob1, alg=Rodas4(), abstol=1e-12, reltol=1e-12)
# solDF1 = DataFrame([[j[i] for j in sol1.u] for i=1:length(sol1.u[1])], species1)

# N1 = solDF1[!,:N]
# plot(sol.t, N1)

