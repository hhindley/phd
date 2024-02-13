using DifferentialEquations, Plots, DataFrames, BenchmarkTools

include("model.jl")
include("parameters.jl")
include("functions.jl")

tspan = (0.0, 1e9)
nutrientQ = 10 .^ range(log10(0.08), stop=log10(0.5), length=6)
chloram = [12, 8, 4, 2, 0]

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

newgr = reshape(grVal, (5,6))
newrmf = reshape(rmfVal, (5,6))

p = plot(legend=false)
for (c, r) in zip(eachcol(newgr), eachcol(newrmf))
    plot!(c, r)
end
display(p)
xlabel!("growth rate [1/min]")
ylabel!("ribosomal mass fraction")
plot!([0.0018,0.02],[0.1,0.18], arrow=true, color=:black)

#= to plot using Gadfly
using Cairo, Fontconfig, Gadfly
img = plot(Theme(background_color = "white"),
layer(x=newgr[:,1], y=newrmf[:,1], Geom.point, Geom.line, Theme(default_color=colorant"green")),
layer(x=newgr[:,2], y=newrmf[:,2], Geom.point, Geom.line, Theme(default_color=colorant"blue")),  
layer(x=newgr[:,3], y=newrmf[:,3], Geom.point, Geom.line, Theme(default_color=colorant"pink")),
layer(x=newgr[:,4], y=newrmf[:,4], Geom.point, Geom.line, Theme(default_color=colorant"yellow")),
layer(x=newgr[:,5], y=newrmf[:,5], Geom.point, Geom.line, Theme(default_color=colorant"red")),
layer(x=newgr[:,6], y=newrmf[:,6], Geom.point, Geom.line, Theme(default_color=colorant"purple"))
)

draw(PNG("chl_ns.png", 15cm, 9cm), img) =#

