function sweep_param(model, ss_init, params, param_range, param, species, title)
    new_params = deepcopy(params)
    new_species = deepcopy(species)
    ssvals=[]
    for i in ProgressBar(param_range)
        new_params[param] = i
        ss = steady_states(model, ss_init, new_params)
        ssvals_dict = Dict([i => j for (i,j) in zip(species, ss)])
        lam = calc_lam(new_params, ssvals_dict)
        ss_lam = collect(ss)
        rmf = calc_rmf(new_params, ssvals_dict)
        push!(ss_lam, lam)
        push!(ss_lam, rmf)
        push!(ssvals, ss_lam)
    end
    # return ssvals
    df_ssvals = DataFrame(vcat(transpose(ssvals)...), :auto)
    push!(new_species, :lam)
    push!(new_species, :rmf)
    rename!(df_ssvals, new_species)
    return plot([scatter(x=param_range, y=col, name="$(names(df_ssvals)[i])") for (col,i) in zip(eachcol(df_ssvals), range(1,length(names(df_ssvals))))],
    Layout(xaxis_title="$param", title=title))#, xaxis_type="log"))

end


function sweep_paramx2(model, ss_init, parameters, species, param1, param2, param_range1, param_range2)
    all_res = []
    params = deepcopy(parameters)
    for i in ProgressBar(param_range1)
        params[param1] = i
        res1 = []
        for val in param_range2 
            params[param2] = val
            ss = steady_states(model, ss_init, params)
            ss_lam = ([ss...])
            # ss_lam = collect(ss1)
            ssvals_dict = Dict([i => j for (i,j) in zip(species, ss)])
            push!(ss_lam, calc_lam(params, ssvals_dict))
            push!(ss_lam, calc_rmf(params, ssvals_dict))
            push!(res1, ss_lam)
        end
        push!(all_res, res1)
    end
    # @show size(all_res)
    new_species = deepcopy(species)
    push!(new_species, :lam)
    push!(new_species, :rmf)
    df = DataFrame([name => [] for name in new_species])
    for i in range(1, length(all_res))
        for j in all_res[i]
            push!(df, j)
        end
    end
    return df
end

function plot_contour(res, specie, param_range1, param_range2, param1, param2)
    rh_res = reshape(res[!,specie], length(param_range1),length(param_range2))
    return plot(contour(x=param_range1, y=param_range2, z=rh_res, colorbar=attr(title="$specie", titleside="right")), Layout(xaxis_title=param1, yaxis_title=param2, title="$specie"))
end


# function sweep_paramx2(model, parameters, all_species, species, func, param1, param2, param_range1, param_range2)
#     all_res = []
#     # params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, atp, na, nb, nr, lam] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
#     params = deepcopy(parameters)
#     for i in param_range1
#         params[param1] = i
#         res1 = []
#         for val in param_range2 
#             params[param2] = val
#             solu = sol(model, initial, tspan, params)
#             push!(res1, func(solu, species, all_species))
#         end
#         push!(all_res, res1)

#     end
#     # @show size(all_res)
#     vec = []
#     for i in (1:length(param_range1))
#         append!(vec, values(all_res[i]))
#     end
#     # @show (vec)
#     vec = reshape(vec, (length(param_range1),length(param_range1)))
#     @show (vec)
#     return plot(contour(x=param_range1, y=param_range2, z=vec, colorbar=attr(title="$species", titleside="right")), Layout(xaxis_title="$param1", yaxis_title="$param2", title="$func of $species"))
# end

function sweep_paramx3(model, lam, species, func, param1, param2, param3, param_range)
    all_res = []
    params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
    for v in param_range
        params[param3] = v
        res1 = []
        for i in param_range
            params[param1] = i
            res2 = []
            for val in param_range  
                params[param2] = val
                solu = sol(model, init, tspan, params)
                push!(res2, func(solu, species))
            end
            push!(res1, res2)
        end
        push!(all_res, res1)
    end
    vec = []
    for i in (1:length(param_range))
        append!(vec, values(all_res[i]))
    end
    vec1 = []
    for i in (1:length(vec))
        append!(vec1, values(vec[i]))
    end
    
    vec2 = reshape(vec1, (length(param_range),length(param_range),length(param_range)))
    vec2 = convert(Array{Float64}, vec2)

    x,y,z = mgrid(param_range, param_range, param_range)

    return plot(volume(x=x[:], y=y[:], z=z[:], value=vec2[:], opacity=0.1, surface_count=21, slices_z=attr(show=true, locations=[0.4], fill=1),
    colorbar=attr(title="$species", titleside="right")), Layout(scene=attr(xaxis_title="$param1", yaxis_title="$param2", zaxis_title="$param3", title="$func of $species")))

end

function change_param_atp(param_range, parameter, model, init, species)#, params)
    # param_dict = OrderedDict("L"=>L, "c"=>c, "kr"=>kr, "Vmax_init"=>Vmax_init, "Km_init"=>Km_init, "ω_ab"=>ω_ab, "ω_r"=>ω_r, "θtscr"=>θtscr, "g_max"=>g_max, "θtlr"=>θtlr, "km_a"=>km_a, "km_b"=>km_b, "gr_c"=>gr_c, "d"=>d, "krep"=>krep, "kdam"=>kdam, "ktag"=>ktag, "kdeg"=>kdeg, "kin"=>kin, "atp"=>atp, "na"=>na, "nb"=>nb, "nr"=>nr, "k"=>k)
    dict_res = OrderedDict(name => [] for name in species)
    params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, na, nb, nr, d_atp] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :na, :nb, :nr, :d_atp)
    for val in param_range  
        params[parameter] = val
        param = values(params)
        # @show params[:kdam]
        solu = sol(model, init, tspan, param)
        # solu = solcb(model, init, tspan, param, cb)
        for (i,j) in zip(values(dict_res), species)
            push!(i, get_ssval(solu, j))
        end
    end
    return dict_res
end

function change_param(param_range, parameter, model, init, func, species, params)
    # param_dict = OrderedDict("L"=>L, "c"=>c, "kr"=>kr, "Vmax_init"=>Vmax_init, "Km_init"=>Km_init, "ω_ab"=>ω_ab, "ω_r"=>ω_r, "θtscr"=>θtscr, "g_max"=>g_max, "θtlr"=>θtlr, "km_a"=>km_a, "km_b"=>km_b, "gr_c"=>gr_c, "d"=>d, "krep"=>krep, "kdam"=>kdam, "ktag"=>ktag, "kdeg"=>kdeg, "kin"=>kin, "atp"=>atp, "na"=>na, "nb"=>nb, "nr"=>nr, "k"=>k)
    dict_res = OrderedDict(name => [] for name in species)
    # params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
    # println(params)
    # params = @LArray [L, c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 1, ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)

    new_params = deepcopy(params)
    for val in param_range
        new_params[parameter] = val
        # @show params[parameter]
        @show new_params[parameter]
        # param = values(params)
        # @show params
        solu = sol(model, init, tspan, new_params)
        for (i,j) in zip(values(dict_res), species)
            push!(i, func(solu, j, species))
        end
    end
    return dict_res
end

function sweep_paramx2_new(model, species, func, param1, param2, param_range1, param_range2, initial, params, xlog, ylog)
    all_res = []
    # params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
    new_params = deepcopy(params)
    for i in param_range1
        new_params[param1] = i
        # @show params[:kdam]

        res1 = []
        for val in param_range2 
            new_params[param2] = val
            # @show params[:kin]
            solu = sol(model, initial, tspan, new_params)
            push!(res1, func(solu, species))
        end
        push!(all_res, res1)

    end
    # @show size(all_res)
    vec = []
    for i in (1:length(param_range1))
        append!(vec, values(all_res[i]))
    end
    # @show (vec)
    vec = reshape(vec, (length(param_range1),length(param_range1)))
    return plot(contour(x=param_range1, y=param_range2, z=(vec), colorbar=attr(title="$species", titleside="right", exponentformat="power")), Layout(xaxis_title="$param1", yaxis_title="$param2", title="$func of $species", xaxis_type=xlog, yaxis_type=ylog))
end

function sweep_paramx2_few_t(model, tspan, t, lam, atp, kin, species, func, param1, param2, param_range1, param_range2)
    all_res = []
    params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
    for i in param_range1
        params[param1] = i
        res1 = []
        for val in param_range2 
            params[param2] = val
            solu = sol_with_t(model, initial, params, tspan, t)
            push!(res1, func(solu, species))
        end
        push!(all_res, res1)

    end
    # @show size(all_res)
    vec = []
    for i in (1:length(param_range1))
        append!(vec, values(all_res[i]))
    end
    # @show (vec)
    vec = reshape(vec, (length(param_range1),length(param_range1)))
    @show (vec)
    return plot(contour(x=param_range1, y=param_range2, z=vec, colorbar=attr(title="$species", titleside="right")), Layout(xaxis_title="$param1", yaxis_title="$param2", title="$func of $species"))
end



function save_1x_plots(range, param, title, log, initial, func)
    # print(params)
    params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
    results = change_param(range, param, rtc_model, initial, func, all_species, params)
    p = display(plot_change_param_sols(range, results, "$param", title, log))
    # return p 
    # open("/home/holliehindley/phd/may23_rtc/analysis/results/1x_param_sweep/without_damage/$param.html", "w") do io
    #     PlotlyBase.to_html(io, p.plot)
    # end
end

function save_2x_plots(param1, param2, param1_range, param2_range, folder, ss_init, xlog, ylog)
    for i in all_species
        p = display(sweep_paramx2_new(rtc_model, i, get_ssval, param1, param2, param1_range, param2_range, ss_init, xlog, ylog))
        # open("/home/holliehindley/phd/may23_rtc/analysis/results/2x_param_sweep/with_damage/$folder/$i.html", "w") do io
        #     PlotlyBase.to_html(io, p.plot)
        # end
    end
end



function solve_plot_param_comparison(rtc_model, initial, tspan, params, species, param, param_list, other_param, log)
    curves=[]
    for i in param_list
        params[param] = i
        solu = sol(rtc_model, initial, tspan, params)
        specie = get_curve(solu, species)
        push!(curves, scatter(x=solu.t, y=specie, name="$i"))
    end
    return (plot([i for i in curves], Layout(xaxis_type=log, yaxis_title="$species", xaxis_title="time", title="changing $param, $other_param")))
end

function solve_plot_param_comparison2(rtc_model, initial, tspan, params, species, param, param_list, log)
    curves=[]
    for (i,j) in zip(param_list[1], param_list[2])
        params[param[1]] = i
        params[param[2]] = j
        solu = sol(rtc_model, initial, tspan, params)
        specie = get_curve(solu, species)
        push!(curves, scatter(x=solu.t, y=specie, name="$i,$j"))
    end
    return (plot([i for i in curves], Layout(xaxis_type=log, yaxis_title="$species", xaxis_title="time", title="changing $param")))
end

function param_comparison(rtc_model, initial, tspan, params, species, param, param_list, other_param, log)
    if length(param_list) == 2
        solve_plot_param_comparison2(rtc_model, initial, tspan, params, species, param, param_list, log)
    else
        solve_plot_param_comparison(rtc_model, initial, tspan, params, species, param, param_list, other_param, log)
    end
end


function param2x_plot_same_species(param_range, param, param2_range, param2, specie)
    curves = []
    for i in param_range
        results = []
        params[param] = i
        results = change_param(param2_range, param2, rtc_model, initial, all_species, params)
        push!(curves, scatter(x=param2_range, y=results[specie], name="$param = $i"))
    end

    return plot([i for i in curves], Layout(title="changing $param and $param2", xaxis_title="$param2", yaxis_title="$specie (μM)"))
end
function save_plot(p, title, folder)
    open("/home/holliehindley/phd/may23_rtc/analysis/results/rabbit_holes/$folder/$title.html", "w") do io
        PlotlyBase.to_html(io, p.plot)
    end
end

function param2x_same_species(param_range, param, param2_range, param2, specie)
    results = []
    for i in param_range
        # results = []
        params[param] = i
        result = change_param(param2_range, param2, rtc_model, initial, all_species, params)
        push!(results, result[specie])
    end

    return results
end




function change_param_timevars(param_range, parameter, model, init, species, atp, lam, kin)#, params)
    # param_dict = OrderedDict("L"=>L, "c"=>c, "kr"=>kr, "Vmax_init"=>Vmax_init, "Km_init"=>Km_init, "ω_ab"=>ω_ab, "ω_r"=>ω_r, "θtscr"=>θtscr, "g_max"=>g_max, "θtlr"=>θtlr, "km_a"=>km_a, "km_b"=>km_b, "gr_c"=>gr_c, "d"=>d, "krep"=>krep, "kdam"=>kdam, "ktag"=>ktag, "kdeg"=>kdeg, "kin"=>kin, "atp"=>atp, "na"=>na, "nb"=>nb, "nr"=>nr, "k"=>k)
    dict_res = OrderedDict(name => [] for name in species)
    # params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
    # println(params)
    params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
    new_params = deepcopy(params)
    for val in param_range
        new_params[parameter] = val
        # @show params[parameter]
        @show new_params[:kdam]
        # param = values(params)
        # @show params
        solu = sol(model, init, tspan, new_params)
        for (i,j) in zip(values(dict_res), species)
            push!(i, get_ssval(solu, j))
        end
    end
    return dict_res
end


