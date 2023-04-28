# get curve of species that we want to compare data/model to 

all_species = [:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt]
all_species_OD = [:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt, :OD]
all_species_atp = [:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt, :atp]

function get_curve(sol, species)
    df = DataFrame(sol)
    if length(sol[1]) == 9
        rename!(df, [:time, :rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt])
    else
        rename!(df, [:time, :rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt, :atp])
    end
    species = df[:, species]
    return species
end

function get_ssval(sol, species)
    df = DataFrame(sol)
    if length(sol[1]) == 9
        rename!(df, [:time, :rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt])
    else
        rename!(df, [:time, :rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt, :atp])
    end
    species = df[end, species]
    return species
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

function sweep_paramx2_new(model, lam, atp, kin, species, func, param1, param2, param_range1, param_range2)
    all_res = []
    params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
    for i in param_range1
        params[param1] = i
        res1 = []
        for val in param_range2 
            params[param2] = val
            solu = sol(model, initial, tspan, params)
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

function check_get_ssval(sol, species, n=3)
    df = DataFrame(sol)
    if length(sol[1]) == 9
        rename!(df, [:time, :rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt])
    else
        rename!(df, [:time, :rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt, :OD])
    end

    ss_vals = Dict(:rm_a=>[], :rtca=>[], :rm_b=>[], :rtcb=>[], :rm_r=>[], :rtcr=>[], :rh=>[], :rd=>[], :rt=>[])
    for i in collect(0:2)
        alpha = df[!,:rt][end-i]/kr 
        fa = (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6)
        ra = fa*df[!,:rtcr][end-i]
        
        Vinit = ra*Vmax_init*atp/(Km_init+atp)
        tscr_el_a = ω_ab*atp/(θtscr+atp)
        tscr_a = Vinit*tscr_el_a
        tscr_el_b = ω_ab*atp/(θtscr+atp)
        tscr_b = Vinit*tscr_el_b
        tscr_r = ω_r*atp/(θtscr+atp)

        tlr_el = g_max*atp/(θtlr+atp)

        rtca1 = (atp*df[!,:rtca][end-i])/(atp+(km_a*df[!,:rd][end-i])) 
        rtcb1 = (atp*df[!,:rtcb][end-i])/(atp+(km_b*df[!,:rt][end-i])) 

        Vrep = krep*rtcb1*df[!,:rt][end-i]
        Vdam = kdam*df[!,:rh][end-i]
        Vinflux = kin*tlr_el 
        Vtag = ktag*rtca1*df[!,:rd][end-i]

        push!(ss_vals[:rm_a], tscr_a - lam(lam[end-i])*df[!,:rm_a][end-i] - d*df[!,:rm_a][end-i])
        push!(ss_vals[:rtca], (1/na)*df[!,:rh][end-i]*df[!,:rm_a][end-i]*tlr_el - lam(lam[end-i])*df[!,:rtca][end-i])  
        push!(ss_vals[:rm_b], tscr_b - lam(lam[end-i])*df[!,:rm_b][end-i] - d*df[!,:rm_b][end-i])
        push!(ss_vals[:rtcb], (1/nb)*df[!,:rh][end-i]*df[!,:rm_b][end-i]*tlr_el - lam(lam[end-i])*df[!,:rtcb][end-i])
        push!(ss_vals[:rm_r], tscr_r - lam(lam[end-i])*df[!,:rm_r][end-i] - d*df[!,:rm_r][end-i])
        push!(ss_vals[:rtcr], (1/nr)*df[!,:rh][end-i]*df[!,:rm_r][end-i]*tlr_el - lam(lam[end-i])*df[!,:rtcr][end-i])
        push!(ss_vals[:rh], Vrep - Vdam + Vinflux - lam(lam[end-i])*df[!,:rh][end-i])
        push!(ss_vals[:rd], Vdam - Vtag - kdeg*df[!,:rd][end-i] - lam(lam[end-i])*df[!,:rd][end-i])
        push!(ss_vals[:rt], Vtag - Vrep - lam(lam[end-i])*df[!,:rt][end-i])
    end
    return std(ss_vals[species])
    # return (std(df[length(sol.t)-n:length(sol.t), species])/mean(df[length(sol.t)-n:length(sol.t), species]))*100
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
    # solu = solve(prob, Rodas5(), isoutofdomain=(y,p,t)->any(x->x<0,y), abstol=1e-15, reltol=1e-12);

    # solu = solve(prob, alg_hints=[:auto], isoutofdomain=(y,p,t)->any(x->x<0,y))#, abstol=1e-15, reltol=1e-12);

    return solu
end

function solcb(model, init, tspan, params, callback)
    # params, init = choose_param_vector(model)
    prob = ODEProblem(model, init, tspan, params, callback=callback)
    solu = solve(prob, Rodas4())
    return solu
end

function sol_with_t(model, init, params, tspan, t)
    # params, init = choose_param_vector(model)
    prob = ODEProblem(model, init, tspan, params)
    solu = solve(prob, Rodas4(), saveat=t)
    return solu
end

function change_param_OD(param_range, parameter, model, species, lam, data)#, params)
    # param_dict = OrderedDict("L"=>L, "c"=>c, "kr"=>kr, "Vmax_init"=>Vmax_init, "Km_init"=>Km_init, "ω_ab"=>ω_ab, "ω_r"=>ω_r, "θtscr"=>θtscr, "g_max"=>g_max, "θtlr"=>θtlr, "km_a"=>km_a, "km_b"=>km_b, "gr_c"=>gr_c, "d"=>d, "krep"=>krep, "kdam"=>kdam, "ktag"=>ktag, "kdeg"=>kdeg, "kin"=>kin, "atp"=>atp, "na"=>na, "nb"=>nb, "nr"=>nr, "k"=>k)
    dict_res = OrderedDict(name => [] for name in species)
    k = set_k(data)
    OD_0 = set_OD0(data)
    init = @SVector [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, rh_0, rd_0, rt_0, OD_0];
    params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, k, lam] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :k, :lam)
    for val in param_range  
        params[parameter] = val
        param = values(params)
        # @show params
        solu = sol(model, init, tspan, param)
        for (i,j) in zip(values(dict_res), species)
            push!(i, get_ssval(solu, j))
        end
    end
    return dict_res
end

function change_param(param_range, parameter, model, init, species, lam)#, params)
    # param_dict = OrderedDict("L"=>L, "c"=>c, "kr"=>kr, "Vmax_init"=>Vmax_init, "Km_init"=>Km_init, "ω_ab"=>ω_ab, "ω_r"=>ω_r, "θtscr"=>θtscr, "g_max"=>g_max, "θtlr"=>θtlr, "km_a"=>km_a, "km_b"=>km_b, "gr_c"=>gr_c, "d"=>d, "krep"=>krep, "kdam"=>kdam, "ktag"=>ktag, "kdeg"=>kdeg, "kin"=>kin, "atp"=>atp, "na"=>na, "nb"=>nb, "nr"=>nr, "k"=>k)
    dict_res = OrderedDict(name => [] for name in species)
    params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
    for val in param_range  
        params[parameter] = val
        param = values(params)
        # @show params[:kdam]
        solu = sol(model, init, tspan, param)
        for (i,j) in zip(values(dict_res), species)
            push!(i, get_ssval(solu, j))
        end
    end
    return dict_res
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

function plotly_plot_sol(sol, log, log1, title)
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
    rtot = scatter(x=sol.t, y=rh+rd+rt, name="rtot")
    alpha = @. rt/kr 
    fa = @. (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6)
    ra = scatter(x=sol.t, y=fa.*rtcr, name="A_RtcR")
    return (plot([rma_curve, rmb_curve, rmr_curve, rtca_curve, rtcb_curve, rtcr_curve, rh_curve, rt_curve, rd_curve, rtot] ,Layout(xaxis_type=log, yaxis_type=log1, title=title)))#, xaxis_range=(0,1320))))
end

function plotly_plot_sol_atp(sol, log, log1, title, show_leg, atp_end)
    rm_a = get_curve(sol, :rm_a); rm_b = get_curve(sol, :rm_b); rm_r = get_curve(sol, :rm_r); rtca = get_curve(sol, :rtca); rtcb = get_curve(sol, :rtcb); rtcr = get_curve(sol, :rtcr); rh = get_curve(sol, :rh); rt = get_curve(sol, :rt); rd = get_curve(sol, :rd); atp = get_curve(sol, :atp);

    alpha = @. rt/kr
    fa = @. (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6)
    ra = @. fa*rtcr
    Vinit = @. ra*Vmax_init*atp/(Km_init+atp)
    tscr_el_a = @. ω_ab*atp/(θtscr+atp)
    tscr_a = @. Vinit*tscr_el_a
    tscr_el_b = @. ω_ab*atp/(θtscr+atp)
    tscr_b = @. Vinit*tscr_el_b
    tscr_r = @. ω_r*atp/(θtscr+atp)
    tlr_el = @. g_max*atp/(θtlr+atp)
    tlr(rm_x, nx) = @. (1/nx)*rh*rm_x*tlr_el
    tlr_r = tlr(rm_r, nr); tlr_a = tlr(rm_a, na); tlr_b = tlr(rm_b, nb);
    rtca1 = @. (atp*rtca)/(atp+(km_a*rd)) 
    rtcb1 = @. (atp*rtcb)/(atp+(km_b*rt)) 
    Vrep = @. krep*rtcb1*rt
    Vdam = @. kdam*rh
    Vinflux = @. kin*tlr_el
    Vtag = @. ktag*rtca1*rd

    rma_curve = scatter(x=sol.t, y=rm_a, name="mRNA RtcA", showlegend=show_leg)
    rmb_curve = scatter(x=sol.t, y=rm_b, name="mRNA RtcB", showlegend=show_leg)
    rmr_curve = scatter(x=sol.t, y=rm_r, name="mRNA RtcR", showlegend=show_leg)
    rtca_curve = scatter(x=sol.t, y=rtca, name="RtcA", showlegend=show_leg)
    rtcb_curve = scatter(x=sol.t, y=rtcb, name="RtcB", showlegend=show_leg)
    rtcr_curve = scatter(x=sol.t, y=rtcr, name="RtcR", showlegend=show_leg)
    rh_curve = scatter(x=sol.t, y=rh, name="Rh", showlegend=show_leg)
    rt_curve = scatter(x=sol.t, y=rt, name="Rt", showlegend=show_leg)
    rd_curve = scatter(x=sol.t, y=rd, name="Rd", showlegend=show_leg)
    atp_curve = scatter(x=sol.t, y=atp, name="ATP", showlegend=show_leg)
    atp_end_curve = scatter(x=fill(atp_end, length(range(0,atp_ss))) , y=range(0,atp_ss))
    alpha_curve = scatter(x=sol.t, y=alpha, name="alpha")
    fa_curve = scatter(x=sol.t, y=fa, name="fa")
    ra_curve = scatter(x=sol.t, y=ra, name="ra")
    Vinit_curve = scatter(x=sol.t, y=Vinit, name="Vinit")
    tscr_el_a_curve = scatter(x=sol.t, y=tscr_el_a, name="tscr_el_a")
    tscr_a_curve = scatter(x=sol.t, y=tscr_a, name="tscr_a")
    tscr_r_curve = scatter(x=sol.t, y=tscr_r, name="tscr_r")
    tlr_el_curve = scatter(x=sol.t, y=tlr_el, name="tlr_el")
    tlr_a_curve = scatter(x=sol.t, y=tlr_a, name="tlr_a")
    tlr_b_curve = scatter(x=sol.t, y=tlr_b, name="tlr_b")
    tlr_r_curve = scatter(x=sol.t, y=tlr_r, name="tlr_r")
    rtca1_curve = scatter(x=sol.t, y=rtca1, name="rtca1")
    rtcb1_curve = scatter(x=sol.t, y=rtcb1, name="rtcb1")
    Vrep_curve = scatter(x=sol.t, y=Vrep, name="Vrep")
    Vdam_curve = scatter(x=sol.t, y=Vdam, name="Vdam")
    Vinflux_curve = scatter(x=sol.t, y=Vinflux, name="Vinflux")
    Vtag_curve = scatter(x=sol.t, y=Vtag, name="Vtag")

    # return (plot([rma_curve, rmb_curve, rmr_curve, rtca_curve, rtcb_curve, rtcr_curve, rh_curve, rt_curve, rd_curve, atp_curve, atp_end_curve,
    # alpha_curve, fa_curve, ra_curve, Vinit_curve, tscr_el_a_curve, tscr_a_curve, tscr_r_curve, tlr_el_curve, tlr_a_curve, tlr_b_curve, tlr_r_curve,
    # rtca1_curve, rtcb1_curve, Vrep_curve, Vdam_curve, Vinflux_curve,Vtag_curve], Layout(xaxis_type=log, yaxis_type=log1, title=title)))
    return (plot([rma_curve, rmb_curve, rmr_curve, rtca_curve, rtcb_curve, rtcr_curve, rh_curve, rt_curve, rd_curve, atp_curve], Layout(xaxis_type=log, yaxis_type=log1, title=title)))
end

function plotly_plot_sol_timepoints(sol)
    rm_a = get_curve(sol, :rm_a); rm_b = get_curve(sol, :rm_b); rm_r = get_curve(sol, :rm_r); rtca = get_curve(sol, :rtca); rtcb = get_curve(sol, :rtcb); rtcr = get_curve(sol, :rtcr); rh = get_curve(sol, :rh); rt = get_curve(sol, :rt); rd = get_curve(sol, :rd);

    rma_curve = scatter(x=collect(0: length(sol.t)), y=rm_a, name="mRNA RtcA")
    rmb_curve = scatter(x=collect(0: length(sol.t)), y=rm_b, name="mRNA RtcB")
    rmr_curve = scatter(x=collect(0: length(sol.t)), y=rm_r, name="mRNA RtcR")
    rtca_curve = scatter(x=collect(0: length(sol.t)), y=rtca, name="RtcA")
    rtcb_curve = scatter(x=collect(0: length(sol.t)), y=rtcb, name="RtcB")
    rtcr_curve = scatter(x=collect(0: length(sol.t)), y=rtcr, name="RtcR")
    rh_curve = scatter(x=collect(0: length(sol.t)), y=rh, name="Rh")
    rt_curve = scatter(x=collect(0: length(sol.t)), y=rt, name="Rt")
    rd_curve = scatter(x=collect(0: length(sol.t)), y=rd, name="Rd")
    return (plot([rma_curve, rmb_curve, rmr_curve, rtca_curve, rtcb_curve, rtcr_curve, rh_curve, rt_curve, rd_curve]))# ,Layout(xaxis_type="log")))
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

function plotly_plot_sol_OD(sol, log)
    rm_a = get_curve(sol, :rm_a); rm_b = get_curve(sol, :rm_b); rm_r = get_curve(sol, :rm_r); rtca = get_curve(sol, :rtca); rtcb = get_curve(sol, :rtcb); rtcr = get_curve(sol, :rtcr); rh = get_curve(sol, :rh); rt = get_curve(sol, :rt); rd = get_curve(sol, :rd); OD = get_curve(sol, :OD)

    rma_curve = scatter(x=sol.t, y=rm_a, name="mRNA RtcA")
    rmb_curve = scatter(x=sol.t, y=rm_b, name="mRNA RtcB")
    rmr_curve = scatter(x=sol.t, y=rm_r, name="mRNA RtcR")
    rtca_curve = scsatter(x=sol.t, y=rtca, name="RtcA")
    rtcb_curve = scatter(x=sol.t, y=rtcb, name="RtcB")
    rtcr_curve = scatter(x=sol.t, y=rtcr, name="RtcR")
    rh_curve = scatter(x=sol.t, y=rh, name="Rh")
    rt_curve = scatter(x=sol.t, y=rt, name="Rt")
    rd_curve = scatter(x=sol.t, y=rd, name="Rd")
    OD_curve = scatter(x=sol.t, y=OD, name="OD")
    return display(plot([rma_curve, rmb_curve, rmr_curve, rtca_curve, rtcb_curve, rtcr_curve, rh_curve, rt_curve, rd_curve, OD_curve], Layout(xaxis_type=log)))
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


function sweep_paramx2(model, lam, species, func, param1, param2, param_range1, param_range2)
    all_res = []
    params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, atp, na, nb, nr, lam] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
    for i in param_range1
        params[param1] = i
        res1 = []
        for val in param_range2 
            params[param2] = val
            solu = sol(model, initial, tspan, params)
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
# contours_start=0, contours_end=1000000,


function sweep_paramx3(model, lam, species, func, param1, param2, param3, param_range)
    all_res = []
    params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, atp, na, nb, nr, lam] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
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

    return plot(volume(x=x[:], y=y[:], z=z[:], value=vec2[:], opacity=0.1, surface_count=21, #slices_z=attr(show=true, locations=[0.4], fill=1),
    colorbar=attr(title="$species", titleside="right")), Layout(scene=attr(xaxis_title="$param1", yaxis_title="$param2", zaxis_title="$param3", title="$func of $species")))

end

function get_lambda(solu, lam)
    lambda = []
    for t in solu.t
        push!(lambda, lam(t))
    end
    return lambda
end

function get_atp(solu, atp)
    atp_list = []
    for t in solu.t
        push!(atp_list, atp(t))
    end
    return atp
end

function plot_sol_and_lam(solu, lam)
    lambda = get_lambda(solu, lam)
    p15 = plot(scatter(x=solu.t, y=lambda), Layout(title="λ"))
    p = plotly_plot_sol(solu, "log")
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

function extend_gr_curve(csv)
    mean_gr = mean((csv."gr"[Int64((length(csv."t")*2/3)+1):end]))
    df = DataFrame(t=Float64[], gr=Float64[])
    for t in collect(csv."t"[end]+10:5000:1e9)
        push!(df, [t, mean_gr])
    end    
    new_df = vcat(csv, df)
    lam = QuadraticInterpolation(new_df."gr",new_df."t")
    return lam, new_df
end

function extend_gr_curve_low(csv)
    mean_gr = 0.0001
    df = DataFrame(t=Float64[], gr=Float64[])
    for t in collect(csv."t"[end]+10:5000:1e9)
        push!(df, [t, mean_gr])
    end    
    new_df = vcat(csv, df)
    lam = QuadraticInterpolation(new_df."gr",new_df."t")
    return lam, new_df
end  

function extend_atp_curve(csv)
    mean_gr = mean((csv."atp"[Int64((length(csv."t")*2/3)+1):end]))
    df = DataFrame(t=Float64[], atp=Float64[])
    for t in collect(csv."t"[end]+10:5000:1e9)
        push!(df, [t, mean_gr])
    end    
    new_df = vcat(csv, df)
    lam = QuadraticInterpolation(new_df."atp",new_df."t")
    return lam, new_df
end  

function plot_all_change_param(range, res)
    # res = get_all_curves(solu, all_species)
    # lambda = get_lambda(solu, lam)
    alpha = res[:rt]/kr
    fa = @. (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6)
    ra = @. fa*res[:rtcr]
    Vinit = @. ra*Vmax_init*atp/(Km_init+atp)
    tscr_el_a = ω_ab*atp/(θtscr+atp)
    tscr_a = Vinit*tscr_el_a
    tscr_el_b = ω_ab*atp/(θtscr+atp)
    tscr_b = Vinit*tscr_el_b
    tscr_r = ω_r*atp/(θtscr+atp)
    tlr_el = g_max*atp/(θtlr+atp)
    tlr(rm_x, nx) = @. (1/nx)*res[:rh]*rm_x*tlr_el
    tlr_r = tlr(res[:rm_r], nr); tlr_a = tlr(res[:rm_a], na); tlr_b = tlr(res[:rm_b], nb);
    rtca1 = @. (atp*res[:rtca])/(atp+(km_a*res[:rd])) 
    rtcb1 = @. (atp*res[:rtcb])/(atp+(km_b*res[:rt])) 
    Vrep = @. krep*rtcb1*res[:rt]
    Vdam = @. kdam*res[:rh]
    Vinflux = kin*tlr_el
    Vtag = @. ktag*rtca1*res[:rd]


    p1 = plot(scatter(x=range, y=alpha), Layout(title="alpha"));
    p2 = plot(scatter(x=range, y=fa), Layout(title="fa"));
    p3 = plot(scatter(x=range, y=ra), Layout(title="ra"));
    p4 = plot(scatter(x=range, y=Vinit), Layout(title="Vinit"));
    p5 = plot(scatter(x=range, y=tscr_a), Layout(title="tscr_a"));
    p6 = plot(scatter(x=range, y=tscr_b), Layout(title="tscr_b"));
    p7 = plot(scatter(x=range, y=tlr_r), Layout(title="tlr_r"));
    p8 = plot(scatter(x=range, y=tlr_a), Layout(title="tlr_a"));
    p9 = plot(scatter(x=range, y=tlr_b), Layout(title="tlr_b"));
    p10 = plot(scatter(x=range, y=rtca1), Layout(title="rtca1"));
    p11 = plot(scatter(x=range, y=rtcb1), Layout(title="rtcb1"));
    p12 = plot(scatter(x=range, y=Vrep), Layout(title="Vrep"));
    p13 = plot(scatter(x=range, y=Vdam), Layout(title="Vdam"));
    p14 = plot(scatter(x=range, y=Vtag), Layout(title="Vtag"));
    # p15 = plot(scatter(x=range, y=lambda), Layout(title="λ"))

   return [p1 p2 p3; p4 p5 p6; p7 p8 p9; p10 p11 p12; p13 p14 ]
end

function plot_change_param_sols(range, results, param)

    rma = scatter(x=range, y=results[:rm_a], name="rm_a");
    rmb = scatter(x=range, y=results[:rm_b], name="rm_b");
    rmr = scatter(x=range, y=results[:rm_r], name="rm_r");
    rtca = scatter(x=range, y=results[:rtca], name="rtca");
    rtcb = scatter(x=range, y=results[:rtcb], name="rtcb");
    rtcr = scatter(x=range, y=results[:rtcr], name="rtcr");
    rh = scatter(x=range, y=results[:rh], name="rh");
    rd = scatter(x=range, y=results[:rd], name="rd");
    rt = scatter(x=range, y=results[:rt], name="rt");

    return plot([rma, rmb, rmr, rtca, rtcb, rtcr, rh, rd, rt], Layout(xaxis_title="$param", yaxis_title="species (μM)"))
end

function plot_change_param_sols_atp(range, results, param, log, log1)

    rma = scatter(x=range, y=results[:rm_a], name="rm_a");
    rmb = scatter(x=range, y=results[:rm_b], name="rm_b");
    rmr = scatter(x=range, y=results[:rm_r], name="rm_r");
    rtca = scatter(x=range, y=results[:rtca], name="rtca");
    rtcb = scatter(x=range, y=results[:rtcb], name="rtcb");
    rtcr = scatter(x=range, y=results[:rtcr], name="rtcr");
    rh = scatter(x=range, y=results[:rh], name="rh");
    rd = scatter(x=range, y=results[:rd], name="rd");
    rt = scatter(x=range, y=results[:rt], name="rt");
    atp = scatter(x=range, y=results[:atp], name="ATP");

    return plot([rma, rmb, rmr, rtca, rtcb, rtcr, rh, rd, rt, atp], Layout(xaxis_title="$param", yaxis_title="species (μM)", xaxis_type=log, yaxis_type=log1))
end

function plot_change_param_sols_OD(range, results, param)

    rma = scatter(x=range, y=results[:rm_a], name="rm_a");
    rmb = scatter(x=range, y=results[:rm_b], name="rm_b");
    rmr = scatter(x=range, y=results[:rm_r], name="rm_r");
    rtca = scatter(x=range, y=results[:rtca], name="rtca");
    rtcb = scatter(x=range, y=results[:rtcb], name="rtcb");
    rtcr = scatter(x=range, y=results[:rtcr], name="rtcr");
    rh = scatter(x=range, y=results[:rh], name="rh");
    rd = scatter(x=range, y=results[:rd], name="rd");
    rt = scatter(x=range, y=results[:rt], name="rt");
    OD = scatter(x=range, y=results[:OD], name="OD");

    return plot([rma, rmb, rmr, rtca, rtcb, rtcr, rh, rd, rt, OD], Layout(xaxis_title="$param", yaxis_title="species (μM)"))
end