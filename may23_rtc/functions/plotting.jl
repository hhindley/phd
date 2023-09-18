using LabelledArrays

# function solvePlot_time(rtc_model, init, params, tspan, title, log, logy) 
#     # params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
#     @show params
#     solu = sol(rtc_model, init, tspan, params)
#     print(solu.retcode)
#     return plotly_plot_sol(solu, log, logy, "$title")
# end

function solvePlot_time(rtc_model, lam, atp, kin, kdam, title, log) 
    params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
    solu = sol(rtc_model, initial, tspan, params)
    print(solu.retcode)
    return plotly_plot_sol(solu, "", log, "$title")
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
    # alpha = @. rt/kr 
    # fa = @. (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6)
    # ra = scatter(x=sol.t, y=fa.*rtcr, name="A_RtcR")
    return (plot([rma_curve, rmb_curve, rmr_curve, rtca_curve, rtcb_curve, rtcr_curve, rh_curve, rt_curve, rd_curve, rtot] ,Layout(xaxis_type=log, yaxis_type=log1, title=title, xaxis_range=(0,1320), xaxis_title="Time (minutes)", yaxis_title="Concentration (μM)")))
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

function plot_time_vars(lam, atp, kin)
    if atp == "atp_t"
        p = make_subplots(rows=3, cols=1, shared_xaxes=true, vertical_spacing=0.08, subplot_titles=["λ" "ATP" "kin"])
        add_trace!(p, (scatter(x=csv_atp."t", y=lam)), row=1, col=1)
        add_trace!(p, (scatter(x=csv_atp."t", y=atp)), row=2, col=1)
        add_trace!(p, (scatter(x=csv_atp."t", y=kin)), row=3, col=1)
        relayout!(p, showlegend=false, xaxis_type="log", xaxis2_type="log", xaxis3_type="log")
        return p
    else
        p = make_subplots(rows=3, cols=1, shared_xaxes=true, vertical_spacing=0.08, subplot_titles=["λ" "ATP from gr data" "kin"])
        add_trace!(p, (scatter(x=csv_atp."t", y=lam)), row=1, col=1)
        add_trace!(p, (scatter(x=csv_atp."t", y=atp)), row=2, col=1)
        add_trace!(p, (scatter(x=csv_atp."t", y=kin)), row=3, col=1)
        relayout!(p, showlegend=false, xaxis_type="log", xaxis2_type="log", xaxis3_type="log")
        return p
    end

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





function plot_all_vars(solu)
    rm_a = get_curve(solu, :rm_a)
    rtca = get_curve(solu, :rtca)
    rm_b = get_curve(solu, :rm_b)
    rtcb = get_curve(solu, :rtcb)
    rm_r = get_curve(solu, :rm_r)
    rtcr = get_curve(solu, :rtcr)
    rh = get_curve(solu, :rh)
    rd = get_curve(solu, :rd)
    rt = get_curve(solu, :rt)
    # rtc_i = get_curve(solu, :rtca_i)

    atp = 3578.9473684210525
    kdam = 0.5
    kin = 0.022222222

    alpha = @. rt/kr
    fa = @. (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6)
    ra = @. fa*rtcr
    Vinit = @. ra*Vmax_init*atp/(Km_init+atp)
    tscr_el_a = ω_ab*atp/(θtscr+atp)
    tscr_a = Vinit*tscr_el_a
    tscr_el_b = ω_ab*atp/(θtscr+atp)
    tscr_b = Vinit*tscr_el_b
    tscr_r = ω_r*atp/(θtscr+atp)
    tlr_el = g_max*atp/(θtlr+atp)
    tlr(rm_x, nx) = @. (1/nx)*rh*rm_x*tlr_el
    tlr_r = tlr(rm_r, nr); tlr_a = tlr(rm_a, na); tlr_b = tlr(rm_b, nb);
    rtca1 = @. (atp*rtca)/(atp+(km_a*rd)) 
    rtcb1 = @. (atp*rtcb)/(atp+(km_b*rt)) 
    Vrep = @. krep*rtcb1*rt
    Vdam = @. kdam*rh
    Vinflux = kin*tlr_el
    Vtag = @. ktag*rtca1*rd


    p_alpha = scatter(x=solu.t, y=alpha, name="alpha")
    p_fa = scatter(x=solu.t, y=fa, name="fa")
    p_ra = scatter(x=solu.t, y=ra, name="ra")
    p_vinit = scatter(x=solu.t, y=Vinit, name="Vinit")
    p_tscr_a = scatter(x=solu.t, y=tscr_a, name="tscr_a")
    p_tscr_b = scatter(x=solu.t, y=tscr_b, name="tscr_b")
    p_tlr_r = scatter(x=solu.t, y=tlr_r, name="tlr_r")
    p_tlr_a = scatter(x=solu.t, y=tlr_a, name="tlr_a")
    p_tlr_b = scatter(x=solu.t, y=tlr_b, name="tlr_b")
    p_rtca1 = scatter(x=solu.t, y=rtca1, name="rtca1")
    p_rtcb1 = scatter(x=solu.t, y=rtcb1, name="rtcb1")
    p_vrep = scatter(x=solu.t, y=Vrep, name="Vrep")
    p_vdam = scatter(x=solu.t, y=Vdam, name="Vdam")
    p_vtag = scatter(x=solu.t, y=Vtag, name="Vtag")

   return p_alpha, p_fa, p_ra, p_vinit, p_tscr_a, p_tscr_b, p_tlr_a, p_tlr_b, p_tlr_r, p_rtca1, p_rtcb1, p_vrep, p_vdam, p_vtag
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

function plot_change_param_sols(range, results, param, title, log)

    rma = scatter(x=range, y=results[:rm_a], name="rm_a");
    rmb = scatter(x=range, y=results[:rm_b], name="rm_b");
    rmr = scatter(x=range, y=results[:rm_r], name="rm_r");
    rtca = scatter(x=range, y=results[:rtca], name="rtca");
    rtcb = scatter(x=range, y=results[:rtcb], name="rtcb");
    rtcr = scatter(x=range, y=results[:rtcr], name="rtcr");
    rh = scatter(x=range, y=results[:rh], name="rh");
    rd = scatter(x=range, y=results[:rd], name="rd");
    rt = scatter(x=range, y=results[:rt], name="rt");

    return plot([rma, rmb, rmr, rtca, rtcb, rtcr, rh, rd, rt], Layout(xaxis_title="$param", yaxis_title="species (μM)", xaxis_type=log, title=title))
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