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

function check_get_ssval(sol, species, n=3)
    df = DataFrame(sol)
    if length(sol[1]) == 9
        rename!(df, [:time, :rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt])
    else
        rename!(df, [:time, :rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt, :OD])
    end
    # print(df[length(sol.t)-n:length(sol.t), :rtca])
    return std(df[length(sol.t)-n:length(sol.t), species])
end

function get_all_curves(sol, species) 
    dict_res = OrderedDict(name => [] for name in species)
    for (i,j) in zip(values(dict_res), species)
        push!(i, get_curve(sol, j))
    end
    return ((dict_res))
end

function get_all_curves_df(sol, species) 
    df_res = DataFrame(name => [] for name in species)
    for (i,j) in zip(values(df_res), species)
        push!(i, get_curve(sol, j))
    end
    return df_res
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



function sol(model, init, tspan, params)
    # params, init = choose_param_vector(model)
    prob = ODEProblem(model, init, tspan, params)
    solu = solve(prob, Rodas4())
    return solu
end

function sol_with_t(model, init, params, tspan, t)
    # params, init = choose_param_vector(model)
    prob = ODEProblem(model, init, tspan, params)
    solu = solve(prob, Rodas4(), saveat=t)
    return solu
end

function change_param(param_range, parameter, model, init, species, params)
    # param_dict = OrderedDict("L"=>L, "c"=>c, "kr"=>kr, "Vmax_init"=>Vmax_init, "Km_init"=>Km_init, "ω_ab"=>ω_ab, "ω_r"=>ω_r, "θtscr"=>θtscr, "g_max"=>g_max, "θtlr"=>θtlr, "km_a"=>km_a, "km_b"=>km_b, "gr_c"=>gr_c, "d"=>d, "krep"=>krep, "kdam"=>kdam, "ktag"=>ktag, "kdeg"=>kdeg, "kin"=>kin, "atp"=>atp, "na"=>na, "nb"=>nb, "nr"=>nr, "k"=>k)
    dict_res = OrderedDict(name => [] for name in species)
    for val in param_range  
        params[parameter] = val
        param = values(params)
        solu = sol(model, init, tspan, param)
        for (i,j) in zip(values(dict_res), species)
            push!(i, get_ssval(solu, j))
        end
    end
    return dict_res
end

function plotly_plot_sol(sol)
    rm_a = get_curve(sol, :rm_a); rm_b = get_curve(sol, :rm_b); rm_r = get_curve(sol, :rm_r); rtca = get_curve(sol, :rtca); rtcb = get_curve(sol, :rtcb); rtcr = get_curve(sol, :rtcr); rh = get_curve(sol, :rh); rt = get_curve(sol, :rt); rd = get_curve(sol, :rd);

    rma_curve = scatter(x=sol.t, y=rm_a, name="mRNA RtcA")
    rmb_curve = scatter(x=sol.t, y=rm_b, name="mRNA RtcB")
    rmr_curve = scatter(x=sol.t, y=rm_r, name="mRNA RtcR")
    rtca_curve = scatter(x=sol.t, y=rtca, name="RtcA")
    rtcb_curve = scatter(x=sol.t, y=rtcb, name="RtcB")
    rtcr_curve = scatter(x=sol.t, y=rtcr, name="RtcR")
    rh_curve = scatter(x=sol.t, y=rh, name="Rh")
    rt_curve = scatter(x=sol.t, y=rt, name="Rt")
    rd_curve = scatter(x=sol.t, y=rd, name="Rd")
    return (plot([rma_curve, rmb_curve, rmr_curve, rtca_curve, rtcb_curve, rtcr_curve, rh_curve, rt_curve, rd_curve] ,Layout(xaxis_type="log")))
end

function plotly_plot_sol_withdata(sol)
    rm_a = get_curve(sol, :rm_a); rm_b = get_curve(sol, :rm_b); rm_r = get_curve(sol, :rm_r); rtca = get_curve(sol, :rtca); rtcb = get_curve(sol, :rtcb); rtcr = get_curve(sol, :rtcr); rh = get_curve(sol, :rh); rt = get_curve(sol, :rt); rd = get_curve(sol, :rd);

    rma_curve = scatter(x=sol.t, y=rm_a, name="mRNA RtcA")
    rmb_curve = scatter(x=sol.t, y=rm_b, name="mRNA RtcB")
    rmr_curve = scatter(x=sol.t, y=rm_r, name="mRNA RtcR")
    rtca_curve = scatter(x=sol.t, y=rtca, name="RtcA")
    rtcb_curve = scatter(x=sol.t, y=rtcb, name="RtcB")
    rtcr_curve = scatter(x=sol.t, y=rtcr, name="RtcR")
    rh_curve = scatter(x=sol.t, y=rh, name="Rh")
    rt_curve = scatter(x=sol.t, y=rt, name="Rt")
    rd_curve = scatter(x=sol.t, y=rd, name="Rd")
    return display(plot([rma_curve, rmb_curve, rmr_curve, rtca_curve, rtcb_curve, rtcr_curve, rh_curve, rt_curve, rd_curve, data_plot]))
end

function plotly_plot_sol_OD(sol)
    rm_a = get_curve(sol, :rm_a); rm_b = get_curve(sol, :rm_b); rm_r = get_curve(sol, :rm_r); rtca = get_curve(sol, :rtca); rtcb = get_curve(sol, :rtcb); rtcr = get_curve(sol, :rtcr); rh = get_curve(sol, :rh); rt = get_curve(sol, :rt); rd = get_curve(sol, :rd); OD = get_curve(sol, :OD)

    rma_curve = scatter(x=sol.t, y=rm_a, name="mRNA RtcA")
    rmb_curve = scatter(x=sol.t, y=rm_b, name="mRNA RtcB")
    rmr_curve = scatter(x=sol.t, y=rm_r, name="mRNA RtcR")
    rtca_curve = scatter(x=sol.t, y=rtca, name="RtcA")
    rtcb_curve = scatter(x=sol.t, y=rtcb, name="RtcB")
    rtcr_curve = scatter(x=sol.t, y=rtcr, name="RtcR")
    rh_curve = scatter(x=sol.t, y=rh, name="Rh")
    rt_curve = scatter(x=sol.t, y=rt, name="Rt")
    rd_curve = scatter(x=sol.t, y=rd, name="Rd")
    OD_curve = scatter(x=sol.t, y=OD, name="OD")
    return display(plot([rma_curve, rmb_curve, rmr_curve, rtca_curve, rtcb_curve, rtcr_curve, rh_curve, rt_curve, rd_curve, OD_curve]))
end

function plot_from_dict(dict, sol)
    rma_curve = scatter(x=sol.t, y=dict[:rm_a], name="mRNA RtcA")
    rmb_curve = scatter(x=sol.t, y=dict[:rm_b], name="mRNA RtcB")
    rmr_curve = scatter(x=sol.t, y=dict[:rm_r], name="mRNA RtcR")
    rtca_curve = scatter(x=sol.t, y=dict[:rtca], name="RtcA")
    rtcb_curve = scatter(x=sol.t, y=dict[:rtcb], name="RtcB")
    rtcr_curve = scatter(x=sol.t, y=dict[:rtcr], name="RtcR")
    rh_curve = scatter(x=sol.t, y=dict[:rh], name="Rh")
    rt_curve = scatter(x=sol.t, y=dict[:rt], name="Rt")
    rd_curve = scatter(x=sol.t, y=dict[:rd], name="Rd")
    return display(plot([rma_curve, rmb_curve, rmr_curve, rtca_curve, rtcb_curve, rtcr_curve, rh_curve, rt_curve, rd_curve])    )
end


function sweep_paramx2(model, params, species, func, param1, param2, param_range)
    all_res = []
    for i in param_range
        params[param1] = i
        res1 = []
        for val in param_range  
            params[param2] = val
            solu = sol(model, initial, tspan, params)
            push!(res1, func(solu, species))
        end
        push!(all_res, res1)

    end
    vec = []
    for i in (1:length(param_range))
        append!(vec, values(all_res[i]))
    end
    vec = reshape(vec, (length(param_range),length(param_range)))
    return plot(contour(x=param_range, y=param_range, z=vec, colorbar=attr(title="$species")), Layout(xaxis_title="$param2", yaxis_title="$param1", title="$func of $species"))

end

function sweep_paramx3(model, params, species, func, param1, param2, param3, param_range)
    all_res = []
    for v in param_range
        params[param3] = v
        res1 = []
        for i in param_range
            params[param1] = i
            res2 = []
            for val in param_range  
                params[param2] = val
                solu = sol(model, initial, tspan, params)
                push!(res2, func(solu, species))
            end
            push!(res1, res2)
        end
        push!(all_res, res1)
    end

    all_res = reshape(all_res, length(param_range), length(param_range))
    @show size(all_res)

    # rtcas = []
    # for i in (1:length(param_range))
    #     for i in (1:length(param_range))
    #         push!(rtcas, all_res[i][species])
    #     end
    # end
    # @show size(rtcas)

    # values(rtcas)
    # vec = []
    # for i in (1:length(param_range))
    #     append!(vec, values(rtcas[i]))
    # end
    # @show size(vec)
    # # vec = reshape(vec, (length(param_range),length(param_range)))
    # @show size(all_res)

    return all_res
    # return plot(contour(x=ω_ab_range, y=ω_ab_range, z=vec, colorbar=attr(title="$species")), Layout(xaxis_title="ω_ab", yaxis_title="ω_r", title="ss val of $species"))

end

function get_lambda(solu, lam)
    lambda = []
    for t in solu.t
        push!(lambda, lam(t))
    end
    return lambda
end

function plot_sol_and_lam(solu, lam)
    lambda = get_lambda(solu, lam)
    p15 = plot(scatter(x=solu.t, y=lambda), Layout(title="λ"))
    p = plotly_plot_sol(solu)
    return [p p15]
end

function plot_dilution(solu, lam)
    res = get_all_curves(solu, all_species)
    dil_dict = OrderedDict(name => [] for name in all_species)
    for (dict, species) in zip(values(dil_dict), all_species)
        for (t, i) in zip(solu.t, collect(1:length(solu.t)))
            push!(dict, lam(t)*res[species][1][i])
        end
    end
    return plot_from_dict(dil_dict, solu)
end


function plot_degradation(solu)
    res = get_all_curves(solu, all_species)
    deg_dict = OrderedDict(name => [] for name in all_species)
    for (dict, species) in zip(values(deg_dict), all_species)
        for i in collect(1:length(solu.t))
            push!(dict, d*res[species][1][i])
        end
    end
    return plot_from_dict(deg_dict, solu)
end

function plot_all_variables(solu, lam)
    res = get_all_curves(solu, all_species)
    lambda = get_lambda(solu, lam)
    alpha = res[:rt][1]/kr
    fa = @. (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6)
    ra = @. fa*res[:rtcr][1]
    Vinit = @. ra*Vmax_init*atp/(Km_init+atp)
    tscr_el_a = ω_ab*atp/(θtscr+atp)
    tscr_a = Vinit*tscr_el_a
    tscr_el_b = ω_ab*atp/(θtscr+atp)
    tscr_b = Vinit*tscr_el_b
    tscr_r = ω_r*atp/(θtscr+atp)
    tlr_el = g_max*atp/(θtlr+atp)
    tlr(rm_x, nx) = @. (1/nx)*res[:rh][1]*rm_x*tlr_el
    tlr_r = tlr(res[:rm_r][1], nr); tlr_a = tlr(res[:rm_a][1], na); tlr_b = tlr(res[:rm_b][1], nb);
    rtca1 = @. (atp*res[:rtca][1])/(atp+(km_a*res[:rd][1])) 
    rtcb1 = @. (atp*res[:rtcb][1])/(atp+(km_b*res[:rt][1])) 
    Vrep = @. krep*rtcb1*res[:rt][1]
    Vdam = @. kdam*res[:rh][1]
    Vinflux = kin*tlr_el
    Vtag = @. ktag*rtca1*res[:rd][1]


    p1 = plot(scatter(x=solu.t, y=alpha), Layout(title="alpha"));
    p2 = plot(scatter(x=solu.t, y=fa), Layout(title="fa"));
    p3 = plot(scatter(x=solu.t, y=ra), Layout(title="ra"));
    p4 = plot(scatter(x=solu.t, y=Vinit), Layout(title="Vinit"));
    p5 = plot(scatter(x=solu.t, y=tscr_a), Layout(title="tscr_a"));
    p6 = plot(scatter(x=solu.t, y=tscr_b), Layout(title="tscr_b"));
    p7 = plot(scatter(x=solu.t, y=tlr_r), Layout(title="tlr_r"));
    p8 = plot(scatter(x=solu.t, y=tlr_a), Layout(title="tlr_a"));
    p9 = plot(scatter(x=solu.t, y=tlr_b), Layout(title="tlr_b"));
    p10 = plot(scatter(x=solu.t, y=rtca1), Layout(title="rtca1"));
    p11 = plot(scatter(x=solu.t, y=rtcb1), Layout(title="rtcb1"));
    p12 = plot(scatter(x=solu.t, y=Vrep), Layout(title="Vrep"));
    p13 = plot(scatter(x=solu.t, y=Vdam), Layout(title="Vdam"));
    p14 = plot(scatter(x=solu.t, y=Vtag), Layout(title="Vtag"));
    p15 = plot(scatter(x=solu.t, y=lambda), Layout(title="λ"))

   return [p1 p2 p3; p4 p5 p6; p7 p8 p9; p10 p11 p12; p13 p14 p15]
end


function scale_lam(csv, species)
    rtcas = []
    lam_range = collect(1:1:10)
    for i in lam_range
        lam = QuadraticInterpolation((csv."gr")*i,csv."t")
        params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, atp, na, nb, nr, lam] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
        solu = sol(rtc_model1!, initial, tspan, params)
        push!(rtcas, check_get_ssval(solu, species))
    end
return plot(scatter(x=lam_range, y=rtcas), Layout(xaxis_title="λ", yaxis_title="std last vals $species"))
end