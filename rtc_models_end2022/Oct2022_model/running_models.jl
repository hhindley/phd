using CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, PlotlyJS, Statistics
include("/home/holliehindley/phd/rtc_models/Oct2022_model/rtc_model.jl")
include("/home/holliehindley/phd/rtc_models/sol_species_funcs.jl")
include("/home/holliehindley/phd/rtc_models/params_init_tspan.jl")
include("/home/holliehindley/phd/Param_inf/inf_setup.jl")

csv = DataFrame(CSV.File("/home/holliehindley/phd/data/results_colD_grfit.csv")) # read csv to a datafram
csv = select!(csv, Not(["log(OD)", "log(OD) error", "gr error", "od"]))
csv_wt = DataFrame(CSV.File("/home/holliehindley/phd/data/results_rtcOFF_grfit.csv"))
csv_wt = select!(csv_wt, Not(["log(OD)", "log(OD) error", "gr error", "od"]))
csv_new = DataFrame(CSV.File("/home/holliehindley/phd/data/results_od1303_grfit.csv"))
csv_new = select!(csv_new, Not(["log(OD)", "log(OD) error", "gr error", "od"]))
csv_new2 = DataFrame(CSV.File("/home/holliehindley/phd/data/results_od1303_grfit2.csv"))
csv_new2 = select!(csv_new2, Not(["log(OD)", "log(OD) error", "gr error", "od"]))

lam_colD, new_df = extend_gr_curve(csv)
lam_wt, new_df_wt = extend_gr_curve(csv_wt)
lam_colD[lam_colD.< 0] .= 0 #zero(eltype(lam_colD))
lam_colD

lam, new_df1 = extend_gr_curve_low(csv_new)
lam[lam.< 0] .= 0 #zero(eltype(lam_colD))

lam2, new_df2 = extend_gr_curve_low(csv_new2)
lam2[lam2.< 0] .= 0 #zero(eltype(lam_colD))

tspan_wt = (0, lam_wt[end])
params_wt = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam_wt] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
solu_wt = sol(rtc_model1!, initial, tspan_wt, params_wt)
p = plotly_plot_sol(solu_wt, "log", "")

tspan_colD = (0, lam_colD[end])
params_colD = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam_colD] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
solu_colD = sol(rtc_model1!, initial, tspan_colD, params_colD)
p = plotly_plot_sol(solu_colD, "log", "", "colD data")

rm_a = get_curve(solu_colD, :rm_a)

p = scatter(x=solu_colD.t, y=rm_a)
p1 = scatter(x=(300,300), y=(0,10200))
stat = plot([p,p1], Layout(xaxis_range=[100,1000], yaxis_range=[0,20]))
savefig(stat,"stationary_phase_vis.svg")




tspan_new = (0, lam[end])
params_new = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
solu_new = sol(rtc_model1!, initial, tspan_new, params_new)
p = plotly_plot_sol(solu_new, "log", "", "unsmoothed gr curve")
lam = get_lambda(solu_new, lam)
plot(scatter(x=solu_new.t, y=lam), Layout(xaxis_type="log"))

tspan_new2 = (0, lam2[end])
params_new2 = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam2] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
solu_new2 = sol(rtc_model1!, initial, tspan_new2, params_new2)
p = plotly_plot_sol(solu_new2, "log", "", "smoothed gr curve")
lam2 = get_lambda(solu_new2, lam2)
plot(scatter(x=solu_new2.t, y=lam2), Layout(xaxis_type="log"))

rm_a1 = get_curve(solu_new2, :rm_a)
p = scatter(x=solu_new2.t, y=rm_a1)
p1 = scatter(x=(600,600), y=(0,23500))
stat = plot([p,p1], Layout(xaxis_range=[log(1),log(50)], yaxis_range=[0,200], xaxis_type="log"))

savefig(stat,"stationary_phase_vis.svg")



params_kin_kdam = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, ktag, kdeg, atp, na, nb, nr, lam_colD] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :ktag, :kdeg, :atp, :na, :nb, :nr, :lam)


solu_kin_kdam = sol(rtc_model_kin_kdam_t!, initial, tspan_colD, params_kin_kdam)
p = plotly_plot_sol(solu_kin_kdam, "log", "log")
res = get_all_curves(solu_kin_kdam, all_species)
lam = get_lambda(solu_wt, lam_wt)
kdam = @. (kin/res[:rh][1]) - lam
kin1 = @. lam*res[:rh][1]/g_max
vinf = @. kin1*(g_max*atp/(θtlr+atp))
rtcb1 = @. (atp*res[:rtcb][1])/(atp+(km_b*res[:rt][1])) 
vrep = @. krep*rtcb1*res[:rt][1]
vdam = @. kdam*res[:rh][1]
p1 = plot(scatter(x=solu_kin_kdam.t, y=kdam), Layout(title="/g_max", xaxis_type="log", yaxis_title="kdam", xaxis_title="time"))
p2 = plot(scatter(x=solu_kin_kdam.t, y=kin1), Layout(title="/g_max", xaxis_type="log", yaxis_title="kin", xaxis_title="time"))
p3 = plot(scatter(x=solu_wt.t, y=lam), Layout(title="/g_max", xaxis_type="log", yaxis_title="λ", xaxis_title="time"))
p4 = plot(scatter(x=solu_kin_kdam.t, y=res[:rh][1]), Layout(title="/g_max", xaxis_type="log", yaxis_title="rh", xaxis_title="time"))
p5 = plot(scatter(x=solu_kin_kdam.t, y=vinf), Layout(title="/g_max", xaxis_type="log", yaxis_title="V_influx", xaxis_title="time"))
p_rep = plot(scatter(x=solu_kin_kdam.t, y=vrep), Layout(title="/g_max", xaxis_type="log", yaxis_title="V_rep", xaxis_title="time"))
p_dam = plot(scatter(x=solu_kin_kdam.t, y=vdam), Layout(title="/g_max", xaxis_type="log", yaxis_title="V_dam", xaxis_title="time"))

solu_kin_kdam_tlr = sol(rtc_model_kin_kdam_tlr!, initial, tspan_colD, params_kin_kdam)
p_tlr = plotly_plot_sol(solu_kin_kdam_tlr, "log", "log")
res_tlr = get_all_curves(solu_kin_kdam_tlr, all_species)
lam_tlr = get_lambda(solu_kin_kdam_tlr, lam_colD)
kdam_tlr = @. (kin/res_tlr[:rh][1]) - lam_tlr
kin_tlr = @. lam_tlr*res_tlr[:rh][1]/((g_max*atp/(θtlr+atp)))
vinf_tlr = @. kin_tlr*(g_max*atp/(θtlr+atp))
rtcb1_tlr = @. (atp*res_tlr[:rtcb][1])/(atp+(km_b*res_tlr[:rt][1])) 
vrep_tlr = @. krep*rtcb1_tlr*res_tlr[:rt][1]
vdam_tlr = @. kdam_tlr*res_tlr[:rh][1]
p6 = plot(scatter(x=solu_kin_kdam_tlr.t, y=kdam_tlr), Layout(xaxis_type="log", yaxis_title="kdam", xaxis_title="time", title="/tlr_el"))
p7 = plot(scatter(x=solu_kin_kdam_tlr.t, y=kin_tlr), Layout(xaxis_type="log", yaxis_title="kin", xaxis_title="time", title="/tlr_el"))
p8 = plot(scatter(x=solu_kin_kdam_tlr.t, y=lam_tlr), Layout(xaxis_type="log", yaxis_title="λ", xaxis_title="time", title="/tlr_el"))
p9 = plot(scatter(x=solu_kin_kdam_tlr.t, y=res_tlr[:rh][1]), Layout(xaxis_type="log", yaxis_title="rh", xaxis_title="time", title="tlr_el"))
p10 = plot(scatter(x=solu_kin_kdam_tlr.t, y=vinf_tlr), Layout(title="/tlr_el", xaxis_type="log", yaxis_title="V_influx", xaxis_title="time"))
p_rep1 = plot(scatter(x=solu_kin_kdam_tlr.t, y=vrep_tlr), Layout(title="/tlr_el", xaxis_type="log", yaxis_title="V_rep", xaxis_title="time"))
p_dam1 = plot(scatter(x=solu_kin_kdam_tlr.t, y=vdam_tlr), Layout(title="/tlr_el", xaxis_type="log", yaxis_title="V_dam", xaxis_title="time"))


[p1;p6]
[p2;p7]
[p3;p8]
[p4;p9]
[p5;p10]
[p_rep;p_rep1]
[p_dam;p_dam1]

solu_V = sol(rtc_model_Vt, initial, tspan_colD, params_colD)
p = plotly_plot_sol(solu_V, "log", "log")

res = get_all_curves(solu_V, all_species)
lam = get_lambda(solu_V, lam_colD)
vinf = @. kin*lam
vdam = @. kdam*res[:rh][1]*lam
p3 = plot(scatter(x=solu_V.t, y=lam), Layout(title="/g_max", xaxis_type="log", yaxis_title="λ", xaxis_title="time"))
p4 = plot(scatter(x=solu_V.t, y=res[:rh][1]), Layout(title="/g_max", xaxis_type="log", yaxis_title="rh", xaxis_title="time"))
p5 = plot(scatter(x=solu_V.t, y=vinf), Layout(title="/g_max", xaxis_type="log", yaxis_title="V_influx", xaxis_title="time"))
p_dam = plot(scatter(x=solu_V.t, y=vdam), Layout(title="/g_max", xaxis_type="log", yaxis_title="V_dam", xaxis_title="time"))



solu_colD = sol(rtc_model1!, initial, tspan_colD, params_colD)
p = plotly_plot_sol(solu_colD, "log", "")

res_cd = get_all_curves(solu_colD, all_species)
lam_cd = get_lambda(solu_colD, lam_colD)
# vinf_cd = @. kin*lam_cd
dam_cd = @. kdam*res_cd[:rh][1]
p3_cd = plot(scatter(x=solu_colD.t, y=lam_cd), Layout(title="/g_max", xaxis_type="log", yaxis_title="λ", xaxis_title="time"))
p4_cd = plot(scatter(x=solu_colD.t, y=res_cd[:rh][1]), Layout(title="/g_max", xaxis_type="log", yaxis_title="rh", xaxis_title="time"))
p_dam_cd = plot(scatter(x=solu_colD.t, y=dam_cd), Layout(title="/g_max", xaxis_type="log", yaxis_title="V_dam", xaxis_title="time"))

[p3;p3_cd]
[p4;p4_cd]

[p_dam;p_dam_cd]

[p;p3_cd]


dilu = @. lam_cd*res_cd[:rm_a][1]
p_dilu = plot(scatter(x=solu_colD.t, y=dilu), Layout(xaxis_type="log"))

alpha = @. res_cd[:rt][1]/kr 
fa = @. (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6)
ra = @. fa*res_cd[:rtcr][1]

# transcription
Vinit = @. ra*Vmax_init*atp/(Km_init+atp)
tscr_el_a = @. ω_ab*atp/(θtscr+atp)
tscr_a = @. Vinit*tscr_el_a

# ODEs
drm_a = @. tscr_a - d*res_cd[:rm_a][1]

rest = plot(scatter(x=solu_colD.t, y=drm_a), Layout(xaxis_type="log"))

[p_dilu;rest;p]