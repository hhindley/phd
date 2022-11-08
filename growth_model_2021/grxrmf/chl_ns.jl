using DifferentialEquations, Plots, DataFrames, BenchmarkTools

include("model.jl")
include("parameters.jl")
include("functions.jl")

tspan = (0.0, 1e9)
nutrientQ = 10 .^ range(log10(0.08), stop=log10(0.5), length=6)
chloram = [12, 8, 4, 2, 0]

# species = [:cr, :em, :cp, :cq, :ct, :et, :cm, :mt, :mm, :q, :p, :si, :mq, :mp, :mr, :r, :a, :zmr, :zmp, :zmt, :zmm, :zmq]
# solDF = DataFrame([[j[i] for j in sol.u] for i=1:length(sol.u[1])], species)
# a = solDF[!, :zmp]
# prob = ODEProblem(odemodelfull!,initfull,tspan,params)
# sol = solve(prob, alg=Rosenbrock23())

# plot(sol)

function grrmfCurve(nutrientQ, chloram, tspan)
    gr = []
    rmf = []
    species = [:cr, :em, :cp, :cq, :ct, :et, :cm, :zmr, :zmp, :zmt, :zmm, :zmq, :mt, :mm, :q, :p, :si, :mq, :mp, :mr, :r, :a]

    for i in nutrientQ
        for j in chloram
            global ns = i
            params[26] = i
            global cl = j 
            params[5] = j 
            global f = j*k_cm
            params[7] = j*k_cm
            println(f)
            prob = ODEProblem(odemodelfull!,initfull,tspan,params)
            sol = solve(prob, alg=Rosenbrock23())
            solDF = DataFrame([[j[i] for j in sol.u] for i=1:length(sol.u[1])], species)
            # display(solDF[end,:][:mr])
            push!(gr, calcGrowthrate(solDF[end,:]))
            push!(rmf, calcRMF(solDF[end,:]))
        end
    end
    return gr, rmf
end 

@time grVal, rmfVal = grrmfCurve(nutrientQ, chloram, tspan)

# println(length(grVal))

newgr = reshape(grVal, (5,6))
newrmf = reshape(rmfVal, (5,6))

# println(newrmf)
# show(size(newgr))

# print(size(newgr[:,1]))

# show(typeof(newgr))

# show(typeof(newrmf))
# plot(newgr[:,4], newrmf[:,4])
# plot(newgr, newrmf)

p = plot()
for (c, r) in zip(eachcol(newgr), eachcol(newrmf))
    plot!(c, r)
end
display(p)

# show((eachcol(newgr)))
# show(size(newrmf[:]))