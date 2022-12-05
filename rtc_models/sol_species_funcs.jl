# get curve of species that we want to compare data/model to 

all_species = [:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt]
all_species_OD = [:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt, :OD]

function get_curve(sol, species)
    df = DataFrame(sol)
    if length(sol[1]) == 9
        rename!(df, [:time, :rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt])
    else
        rename!(df, [:time, :rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt, :OD])
    end
    species = df[:, species]
    return species
end

function get_ssval(sol, species)
    df = DataFrame(sol)
    if length(sol[1]) == 9
        rename!(df, [:time, :rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt])
    else
        rename!(df, [:time, :rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt, :OD])
    end
    species = df[end, species]
    return species
end

function get_all_curves(sol, species) 
    dict_res = OrderedDict(name => [] for name in all_species)
    for (i,j) in zip(values(dict_res), species)
        push!(i, get_curve(sol, j))
    end
    return dict_res
end


# solving functions
function choose_param_vector(model)
    if model == rtc_model
        params = param_dict
        init = initial
    elseif model == rtc_model_OD
        params =  param_dict_OD
        init = init_OD
    else 
        params = param_dict_ko
        init = init_OD
    end
    return (@SVector [values(params)])[1], init
end



function sol(model, init, tspan)
    params, init = choose_param_vector(model)
    prob = ODEProblem(model, init, tspan, params)
    sol = solve(prob, Rodas4())
    return sol
end

function sol_with_t(model, init, params, tspan, t)
    # params, init = choose_param_vector(model)
    prob = ODEProblem(model, init, tspan, params)
    sol = solve(prob, Rodas4(), saveat=t)
    return sol
end

function change_param(param_range, parameter, model, init, species, param_dict)
    # param_dict = OrderedDict("L"=>L, "c"=>c, "kr"=>kr, "Vmax_init"=>Vmax_init, "Km_init"=>Km_init, "ω_ab"=>ω_ab, "ω_r"=>ω_r, "θtscr"=>θtscr, "g_max"=>g_max, "θtlr"=>θtlr, "km_a"=>km_a, "km_b"=>km_b, "gr_c"=>gr_c, "d"=>d, "krep"=>krep, "kdam"=>kdam, "ktag"=>ktag, "kdeg"=>kdeg, "kin"=>kin, "atp"=>atp, "na"=>na, "nb"=>nb, "nr"=>nr, "k"=>k)
    dict_res = OrderedDict(name => [] for name in species)
    for val in param_range  
        param_dict[parameter] = val
        params = values(param_dict) # try get this to be an Svector at some point
        solu = sol(model, init, tspan)
        for (i,j) in zip(values(dict_res), species)
            push!(i, get_ssval(solu, j))
        end
    end
    return dict_res
end

