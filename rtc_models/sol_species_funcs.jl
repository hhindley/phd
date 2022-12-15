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
    return (plot([rma_curve, rmb_curve, rmr_curve, rtca_curve, rtcb_curve, rtcr_curve, rh_curve, rt_curve, rd_curve])    )
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