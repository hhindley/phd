using DifferentialEquations, Plots, DataFrames, BenchmarkTools

# include("initial.jl")
# initial_model

include("model.jl")
include("parameters.jl")
include("functions.jl")

tspan = (0.0, 100.0)
nutrient = collect(1:5e4)

function monodCurve(nutrient, tspan)
    gr = []
    species = [:cr, :em, :cp, :cq, :ct, :et, :cm, :mt, :mm, :q, :p, :si, :mq, :mp, :mr, :r, :a]

    for i in nutrient
        # global s0 = i
        params[8] = i
        
        prob = ODEProblem(odemodel!,init,tspan,params)
        sol = solve(prob, alg=Rosenbrock23())
        solDF = DataFrame([[j[i] for j in sol.u] for i=1:length(sol.u[1])], species)
        # display(solDF[end,:][:mr])
        push!(gr, calcGrowthrate(solDF[end,:]))

    end
    return gr
end 

@time growth = monodCurve(nutrient, tspan)

plot(nutrient, growth)
