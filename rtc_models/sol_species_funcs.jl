# get curve of species that we want to compare data/model to 

species = [:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt]

function get_curve(sol, species)
    df = DataFrame(sol)
    rename!(df, [:time, :rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt])
    species = df[:, species]
    return species
end


# solving function
function sol(model, init, tspan, params, t)
    prob = ODEProblem(model, init, tspan, params)
    sol = solve(prob, Rodas4(), saveat=t)
    return sol
end
