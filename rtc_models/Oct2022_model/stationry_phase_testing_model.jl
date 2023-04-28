using CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, PlotlyJS, Statistics
include("/home/holliehindley/phd/rtc_models/Oct2022_model/rtc_model.jl")
include("/home/holliehindley/phd/rtc_models/sol_species_funcs.jl")
include("/home/holliehindley/phd/rtc_models/params_init_tspan.jl")
include("/home/holliehindley/phd/Param_inf/inf_setup.jl")


# load csv for growth rate and atp curves 
# csv_gr = DataFrame(CSV.File("/home/holliehindley/phd/data/results_colD_grfit.csv")) # read csv to a datafram
# csv_gr = select!(csv_gr, Not(["log(OD)", "log(OD) error", "gr error", "od"]))
csv_atp = DataFrame(CSV.File("/home/holliehindley/phd/data/atp_for_rtcmodel.csv"))

csv_atp.atp = csv_atp.atp/5


# set atp and lam curves from files 
# lam_t, new_df = extend_gr_curve(csv_atp)
# lam_t[lam_t.< 0] .= 0 #zero(eltype(lam_colD))
lam_t = QuadraticInterpolation(csv_atp."gr",csv_atp."t")
atp_t = QuadraticInterpolation(csv_atp."atp",csv_atp."t")

p_atp = plot(scatter(x=csv_atp."t", y=atp_t), Layout(xaxis_type="", xaxis_range=(0,1500)))
p_lam = plot(scatter(x=csv_atp."t", y=lam_t), Layout(xaxis_type="", xaxis_range=(0,1500)));


# # set atp curve based on growth rate data 
# atp_t1 = csv_gr."gr"*13213#136000                               # 13213 is the maximum I can get it to go before it starts giving unstable results with more than one variable over time, 136000 is the maximum it can go before it gives unstable results with just ATP over time 
# df1 = DataFrame(t=csv_gr."t", atp=atp_t1)
# df = DataFrame(t=Float64[], atp=Float64[])
# for t in collect(csv_gr."t"[end]+10:5000:1e9)
#     push!(df, [t, atp_t1[end]])
# end    
# new_df = vcat(df1, df)
# atp_t_data = QuadraticInterpolation(new_df."atp",new_df."t")
# atp_t_data[atp_t_data.< 0] .= 0 #zero(eltype(lam_colD))

# p_atp_data = plot(scatter(x=new_df."t", y=atp_t_data), Layout(xaxis_type="log"));


# set kin curve  
rh = 11.29
kin_model = @. lam_t*rh/g_max
kin_t = QuadraticInterpolation(kin_model, csv_atp."t") 
p_kin = plot(scatter(x=csv_atp."t", y=kin_t), Layout(xaxis_type="", xaxis_range=(0,1500)));


# plot all variables over time 
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

p = plot_time_vars(lam_t, atp_t, kin_t)
# p2 = plot_time_vars(lam_t, atp_t_data, kin_t)

p = make_subplots(rows=3, cols=1, shared_xaxes=true, vertical_spacing=0.08, subplot_titles=["λ" "ATP" "kin"])
add_trace!(p, (scatter(x=csv_atp."t", y=lam_t)), row=1, col=1)
add_trace!(p, (scatter(x=csv_atp."t", y=atp_t)), row=2, col=1)
add_trace!(p, (scatter(x=csv_atp."t", y=kin_t)), row=3, col=1)
relayout!(p, showlegend=false, xaxis_range=(0,1500))
p

open("./rtc_model_figs/variables_atp_growth_model.html", "w") do io
    PlotlyBase.to_html(io, p.plot)
end
# open("./rtc_model_figs/variables_atp_from_data.html", "w") do io
    # PlotlyBase.to_html(io, p2.plot)
# end


tspan = (0, lam_t[end])
# lam = 0.033; atp = 4000; kin = 0.054;
function tot_ribo(solu)
    rh = get_curve(solu, :rh); rd = get_curve(solu, :rd); rt = get_curve(solu, :rt)
    rtot = @. rh+rd+rt
    return plot(scatter(x=solu.t,y=rtot), Layout(xaxis_type="log"))
end

function solvePlot_time(rtc_model, lam, atp, kin, title, log) 
    params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
    solu = sol(rtc_model, initial, tspan, params)
    print(solu.retcode)
    return plotly_plot_sol(solu, log, "", "$title")
end


# nothing varied over time 
p_none = solvePlot_time(rtc_model, 0.033, 4000, 0.054, "nothing varied over time", "")

open("./rtc_model_figs/nothing_varied.html", "w") do io
    PlotlyBase.to_html(io, p_none.plot)
end


# just lam varied over time
p_justlam = solvePlot_time(rtc_model1!, lam_t, 4000, 0.054, "just λ varied over time", "")

open("./rtc_model_figs/just_lam.html", "w") do io
    PlotlyBase.to_html(io, p_justlam.plot)
end


# just vary kin over time 
p_justkin = solvePlot_time(rtc_model_kin, 0.033, 4000, kin_t, "just kin varied over time", "")

open("./rtc_model_figs/justkin.html", "w") do io
    PlotlyBase.to_html(io, p_justkin.plot)
end


# just atp varied over time 
p_justatp = solvePlot_time(atp_t!, 0.033, atp_t, 0.054, "just ATP varied over time", "")

open("./rtc_model_figs/justatp.html", "w") do io
    PlotlyBase.to_html(io, p_justatp.plot)
end


# atp and lam varied over time 
p_atplam = solvePlot_time(atp_gr_t!, lam_t, atp_t, 0.054, "λ and ATP varied over time", "")

open("./rtc_model_figs/atplam.html", "w") do io
    PlotlyBase.to_html(io, p_atplam.plot)
end


# atp and kin varied over time 
p_atpkin = solvePlot_time(atp_kin_t!, 0.033, atp_t, kin_t, "ATP and kin varied over time", "")

open("./rtc_model_figs/atpkin.html", "w") do io
    PlotlyBase.to_html(io, p_atpkin.plot)
end


# kin and lam varied over time 
p_lamkin = solvePlot_time(lam_kin_t, lam_t, 4000, kin_t, "kin and λ varied over time", "")

open("./rtc_model_figs/kinlam.html", "w") do io
    PlotlyBase.to_html(io, p_lamkin.plot)
end


#all 
p_all = solvePlot_time(rtc_all_t!, lam_t, atp_t, kin_t, "λ, ATP and kin varied over time", "log")

open("./rtc_model_figs/all.html", "w") do io
    PlotlyBase.to_html(io, p_all.plot)
end


sweep_paramx2_new(rtc_all_t!, lam_t, atp_t, kin_t, :rm_a, get_ssval, :ω_ab, :ω_r, collect(0:10:100), collect(0:10:100))
sweep_paramx2_new(rtc_all_t!, lam_t, atp_t, kin_t, :rm_r, get_ssval, :ω_ab, :ω_r, collect(0:10:100), collect(0:10:100))

# lam_model = get_lambda(solu_atplam, lam_t)
# atp_model = get_atp(solu_atplam, atp_t)
# plot(scatter(x=solu_atp_lam.t, y=lam_model), Layout(xaxis_type="log"))
# plot(scatter(x=solu_atp_lam.t, y=atp_model), Layout(xaxis_type="log"))
