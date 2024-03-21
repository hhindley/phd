using DifferentialEquations
using Plots


function determine_partitioning(u, prop, thresh, FixDetReact)
    # reactions with a propensity below the threshold are modeled stochastically
    isStochReact = prop(u) .< thresh

    if !isempty(FixDetReact)
        isStochReact[FixDetReact] .= 0
    end

    return isStochReact
end

function odefunction!(du, u, isStochReact, t)
    # println([t u])

    # Enforce non-negativity
    u[u .< 0] .= 0
    
    # isStochReact = getpartitioning(u)
    a = prop(u)
    # println(u)

    # Sum of stochastic propensities
    a0_s = sum(a[isStochReact])
    
    # Rates (=propensities) of deterministic reactions
    a_d = a .* .~isStochReact
    
    dfposdt = []
    
    if length(u) > I.nrOfItems
        fpos = u[(I.nrOfItems + 1):end]  # Fork position
        dfposdt = ones(length(fpos))
    end

    dySpecies = S * a_d        # ODEs with deterministic reactions
    
    dySpecies[end] = a0_s       # Modify total propensities
    du .= vcat(dySpecies, dfposdt)'   # Additional replication forks

    return du   
end

function conditionStoch(u, t, integrator) #Stochastic reaction 

    if  u[sidx("TotProp")] >= xi[1]
        global stochReact = true 
    else 
        global stochReact = false 
    end 

    # u[I.nrOfItems] - 1
    stochReact
    # u[I.nrOfItems] - xi[1]
end


function conditionSampling(u, t, integrator) #Stochastic reaction 

    if  t >= t_save
        global savedat = true 
    else 
        global savedat = false 
    end 

    savedat
end

    
function affect_term!(integrator) 
    terminate!(integrator)
end
   



