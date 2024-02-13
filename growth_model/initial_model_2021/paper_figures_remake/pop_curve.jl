using DifferentialEquations, Plots, DataFrames, BenchmarkTools

include("model.jl")
include("parameters.jl")
include("functions.jl")
# include("initial.jl")

global ns = 50 # 100
params[26] = 50
global s0 = 1e10 # 1e10
params[9] = 1e10

tspan = collect(1:600)
ssinit = SSinitVals(10) # change this to run the initial values for a longer time period 

ssinit[6] = 1
# initfull[6] = 1

function pop(tspan)
    N = []
    species = [:cr, :em, :cp, :cq, :ct, :et, :cm, :zmr, :zmp, :zmt, :zmm, :zmq, :mt, :mm, :q, :p, :si, :mq, :mp, :mr, :r, :a]
    for i in tspan
        prob = ODEProblem(odemodelfull!,ssinit,i,params) # kinda works better if use initial values set to 0 so next need to try running the ss values for a shorter time span to see if that makes a difference 
        sol = solve(prob, alg=Rodas4(), abstol=1e-12, reltol=1e-12) #Rodas4() works
        solDF = DataFrame([[j[i] for j in sol.u] for i=1:length(sol.u[1])], species)
        gr = calcGrowthrate(solDF[end,:], 3e8) # 3e8
        nucat = calcNucat(solDF[end,:])
        push!(N, gr*nucat)
    end
    return N    
end

@time population = pop(tspan)

plot(tspan, population)

