# get curve of species that we want to compare data/model to 

all_species = [:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt]

function get_curve(sol, species)
    df = DataFrame(sol)
    rename!(df, [:time, :rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt])
    species = df[:, species]
    return species
end

function get_sscurve(sol, species)
    df = DataFrame(sol)
    rename!(df, [:time, :rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt])
    species = df[end, species]
    return species
end

function get_all_curves(sol, species) 
    dict_res = OrderedDict("rm_a"=>[], "rtca"=>[], "rm_b"=>[], "rtcb"=>[], "rm_r"=>[], "rtcr"=>[], "rh"=>[], "rd"=>[], "rt"=>[])
    for (i,j) in zip(values(dict_res), species)
        push!(i, get_curve(sol, j))
    end
    return dict_res
end


# solving functions
function sol(model, init, tspan, params)
    prob = ODEProblem(model, init, tspan, params)
    sol = solve(prob, Rodas4())
    return sol
end

function sol_with_t(model, init, tspan, params, t)
    prob = ODEProblem(model, init, tspan, params)
    sol = solve(prob, Rodas4(), saveat=t)
    return sol
end

param_dict = OrderedDict("L"=>L, "c"=>c, "kr"=>kr, "Vmax_init"=>Vmax_init, "Km_init"=>Km_init, "ω_ab"=>ω_ab, "ω_r"=>ω_r, "θtscr"=>θtscr, "g_max"=>g_max, "θtlr"=>θtlr, "km"=>km, "k_b"=>k_b, "gr_c"=>gr_c, "d"=>d, "krep"=>krep, "kdam"=>kdam, "ktag"=>ktag, "kdeg"=>kdeg, "kin"=>kin, "atp"=>atp)
function change_param(param_range, parameter)
    dict_res = OrderedDict("rm_a"=>[], "rtca"=>[], "rm_b"=>[], "rtcb"=>[], "rm_r"=>[], "rtcr"=>[], "rh"=>[], "rd"=>[], "rt"=>[])
    for val in param_range  
        param_dict[parameter] = val
        params = values(param_dict) # try get this to be an Svector at some point
        solu = sol(rtc_model, init, tspan, params)
        for (i,j) in zip(values(dict_res), all_species)
            push!(i, get_sscurve(solu, j))
        end
    end
    return results
end

