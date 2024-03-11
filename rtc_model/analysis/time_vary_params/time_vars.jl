using Parameters, CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, Statistics
using Plots
using PlotlyJS
include("$PATHmay23_rtc/functions/solving.jl"); include("$PATHmay23_rtc/functions/set_ups.jl"); include("$PATHmay23_rtc/functions/plotting.jl"); 
include("$PATHmay23_rtc/functions/sweep_params.jl"); include("$PATHmay23_rtc/models/rtc_orig.jl"); include("$PATHmay23_rtc/models/atp_lam_kin_t.jl"); 
include("$PATHmay23_rtc/models/single_t.jl"); include("$PATHmay23_rtc/models/combinations_t.jl");



t, atp_t, lam_t, kin_t = set_time_vars("$PATHdata/atp_for_rtcmodel.csv")

p = make_subplots(rows=3, cols=1, shared_xaxes=true, vertical_spacing=0.08, subplot_titles=["λ" "ATP" "kin"])
add_trace!(p, (scatter(x=t, y=lam_t)), row=1, col=1)
add_trace!(p, (scatter(x=t, y=atp_t)), row=2, col=1)
add_trace!(p, (scatter(x=t, y=kin_t)), row=3, col=1)
relayout!(p, showlegend=false, xaxis_range=(0,1500))
p

tspan = (0, lam_t[end])
# lam = 0.033; atp = 4000; kin = 0.054;

@consts begin   
    # L = 10; #10 
    # c = 0.001; 
    # kr = 0.125; 
    # Vmax_init = 39.51; 
    # Km_init = 250; 
    θtscr = 160.01;  
    θtlr = 255.73; 
    # k_b = 17.7; 
    na = 338; 
    nb = 408; 
    nr = 532*6;
    d = 0.2; 
    krep = 137; 
    ktag = 9780;#0.1; 
    # atp = 4000;#2500; 
    km_a = 20; 
    km_b = 16;
    g_max = 2.0923; 
    gr_c = 0.0008856; # 0.000599; 
    kdeg = 0.001; 
    # kin = 0.054; #2.381 
    ω_ab = 4#4#0.093; #0.0828304057748932;#4; 
    ω_r = 0.0019*6#2e-7 #0.0019*6 #70.53; #0.0019*6#79.43865871861044; #0.0019*6;  
    ω_a = 4; 
    ω_b = 4;
    # kdam =  0.#0.000147;#0.05; 
    k = 2; # carrying capacity - changes depending on the data?
    # lam = 0.033;

    rtca_0 = 0#0.00894; 
    rtcb_0 = 0#0.0216; 
    rh_0 = 11.29; #69.56; #69.4
    rtcr_0 = 0# 0.0131 #0.04; # 8.67e-3; # change this based on keeping steady state level the whole time course (levels shouldn't really change)
    rm_a_0 = 0; 
    rm_b_0 = 0; 
    rm_r_0 = 0#0.0131#0.04 # 0; 
    rd_0 = 0; 
    rt_0 = 0;
end
initial = @SVector [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, rh_0, rd_0, rt_0]

# nothing varied over time 
p_none = solvePlot_time(rtc_model, 0.033, 4000, 0.054, 0, "nothing varied over time - without damage", "")
p_none = solvePlot_time(rtc_model, 0.033, 4000, 0.054, 1, "nothing varied over time - with damage at 1", "")
p_none = solvePlot_time(rtc_model, 0.033, 4000, 0.054, 2, "nothing varied over time - with damage at 2", "")

open("./rtc_model_figs/nothing_varied.html", "w") do io
    PlotlyBase.to_html(io, p_none.plot)
end


# just lam varied over time
p_justlam = solvePlot_time(rtc_model_lam_t!, lam_t, 4000, 0.054, "just λ varied over time", "log")

open("./rtc_model_figs/just_lam.html", "w") do io
    PlotlyBase.to_html(io, p_justlam.plot)
end


# just vary kin over time 
p_justkin = solvePlot_time(rtc_model_kin, 0.033, 4000, kin_t, "just kin varied over time", "log")

open("./rtc_model_figs/justkin.html", "w") do io
    PlotlyBase.to_html(io, p_justkin.plot)
end


# just atp varied over time 
p_justatp = solvePlot_time(atp_t!, 0.033, atp_t, 0.054, "just ATP varied over time", "log")

open("./rtc_model_figs/justatp.html", "w") do io
    PlotlyBase.to_html(io, p_justatp.plot)
end


# atp and lam varied over time 
p_atplam = solvePlot_time(atp_gr_t!, lam_t, atp_t, 0.054, "λ and ATP varied over time", "log")

open("./rtc_model_figs/atplam.html", "w") do io
    PlotlyBase.to_html(io, p_atplam.plot)
end


# atp and kin varied over time 
p_atpkin = solvePlot_time(atp_kin_t!, 0.033, atp_t, kin_t, "ATP and kin varied over time", "log")

open("./rtc_model_figs/atpkin.html", "w") do io
    PlotlyBase.to_html(io, p_atpkin.plot)
end


# kin and lam varied over time 
p_lamkin = solvePlot_time(lam_kin_t, lam_t, 4000, kin_t, "kin and λ varied over time", "")

open("./rtc_model_figs/kinlam.html", "w") do io
    PlotlyBase.to_html(io, p_lamkin.plot)
end


#all 
tspan = (0,1440)
p_all = solvePlot_time(rtc_all_t!, lam_t, atp_t, kin_t, "λ, ATP and kin varied over time", "")

savefig(p_all,"p_all.svg")


open("./rtc_model_figs/all.html", "w") do io
    PlotlyBase.to_html(io, p_all.plot)
end


sweep_paramx2_new(rtc_all_t!, lam_t, atp_t, kin_t, :rm_a, get_ssval, :ω_ab, :ω_r, collect(0:10:100), collect(0:10:100))
sweep_paramx2_new(rtc_all_t!, lam_t, atp_t, kin_t, :rm_r, get_ssval, :ω_ab, :ω_r, collect(0:10:100), collect(0:10:100))

# lam_model = get_lambda(solu_atplam, lam_t)
# atp_model = get_atp(solu_atplam, atp_t)
# plot(scatter(x=solu_atp_lam.t, y=lam_model), Layout(xaxis_type="log"))
# plot(scatter(x=solu_atp_lam.t, y=atp_model), Layout(xaxis_type="log"))

kdam_range = range(0, 10,length=50)
# res = change_param_timevars(kdam_range, :kdam, rtc_all_t!, initial, all_species, atp_t, lam_t, kin_t)


tspan = (0,1e9)
params = @LArray [L, c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 1, ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)


res = change_param(kdam_range, :kdam, rtc_model, initial, get_ssval, all_species, params)
#plots plotting
p3 = plot(kdam_range, res[:rm_a], label="mRNA RtcBA")
p3 = plot!(kdam_range, res[:rtca], label="RtcA")
p3 = plot!(kdam_range, res[:rtcb], label="RtcB")
# plot!(kdam_range, res[:rm_b], label="RtcA")
p3 = plot!(kdam_range, res[:rtcr], label="RtcR")
# p3 = plot!([2.01774465,2.01774465],[0,0.061])
# p3 = plot!([0.63571146,0.63571146],[0,0.061])

p4 = plot(kdam_range, res[:rh], label="Rh")
p4 = plot!(kdam_range, res[:rt], label="Rt")
p4 = plot!(kdam_range, res[:rd], label="Rd")
# p4 = plot!([2.01774465,2.01774465],[0,3.1])
# p4 = plot!([0.63571146,0.63571146],[0,3.1])

#plotly plotting
trace1 = scatter(x=kdam_range, y=res[:rtca], name="RtcA")
trace2 = scatter(x=kdam_range, y=res[:rtcb], name="RtcB")
trace3 = scatter(x=kdam_range, y=res[:rm_a], name="mRNA RtcA")
trace4 = scatter(x=kdam_range, y=res[:rm_b], name="mRNA RtcB")
trace5 = scatter(x=kdam_range, y=res[:rtcr], name="RtcR")
trace6 = scatter(x=kdam_range, y=res[:rm_r], name="mRNA RtcR")
trace7 = scatter(x=kdam_range, y=res[:rt], name="Rt")
trace8 = scatter(x=kdam_range, y=res[:rd], name="Rd")
trace9 = scatter(x=kdam_range, y=res[:rh], name="Rh")

plot([trace1,trace2,trace3,trace4,trace5,trace6,trace7,trace8,trace9], Layout(title="Steady-state", xaxis_title="kdam", yaxis_title="concentration"))


function get_max_val(sol, species)
    df = DataFrame(sol)
    rename!(df, [:time, :rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt])
    max_, ind = findmax(df[:, species])
    timepoint_max = df[ind, :time]
    return max_, timepoint_max
end


res = change_param(kdam_range, :kdam, rtc_model, initial, get_max_val, all_species, 4000, 0.033, 0.054, 1)

max_vals = OrderedDict(name => [] for name in all_species)
tps = OrderedDict(name => [] for name in all_species)
for (max_, tp, specie) in zip(values(max_vals), values(tps), all_species)
    push!(max_, ([[res[specie][i][1] for i in range(1,length(res[:rtca]))]...]...))
    push!(tp, ([[res[specie][i][2] for i in range(1,length(res[:rtca]))]...]...))
end

max1 = scatter(x=kdam_range, y=max_vals[:rtca], name="RtcA")
max2 = scatter(x=kdam_range, y=max_vals[:rtcb], name="RtcB")
max3 = scatter(x=kdam_range, y=max_vals[:rm_a], name="mRNA RtcA")
max4 = scatter(x=kdam_range, y=max_vals[:rm_b], name="mRNA RtcB")
max5 = scatter(x=kdam_range, y=max_vals[:rtcr], name="RtcR")
max6 = scatter(x=kdam_range, y=max_vals[:rm_r], name="mRNA RtcR")
max7 = scatter(x=kdam_range, y=max_vals[:rt], name="Rt")
max8 = scatter(x=kdam_range, y=max_vals[:rd], name="Rd")
max9 = scatter(x=kdam_range, y=max_vals[:rh], name="Rh")

plot([max1,max2,max3,max4,max5,max6,max7,max8,max9], Layout(title="Max vals", xaxis_title="kdam", yaxis_title="concentration"))

tp1 = scatter(x=kdam_range, y=tps[:rtca], name="RtcA")
tp2 = scatter(x=kdam_range, y=tps[:rtcb], name="RtcB")
tp3 = scatter(x=kdam_range, y=tps[:rm_a], name="mRNA RtcA")
tp4 = scatter(x=kdam_range, y=tps[:rm_b], name="mRNA RtcB")
tp5 = scatter(x=kdam_range, y=tps[:rtcr], name="RtcR")
tp6 = scatter(x=kdam_range, y=tps[:rm_r], name="mRNA RtcR")
tp7 = scatter(x=kdam_range, y=tps[:rt], name="Rt")
tp8 = scatter(x=kdam_range, y=tps[:rd], name="Rd")
tp9 = scatter(x=kdam_range, y=tps[:rh], name="Rh")

plot([tp1,tp2,tp3,tp4,tp5,tp6,tp7,tp8,tp9], Layout(title="Timepoint of max val", xaxis_title="kdam", yaxis_title="time", yaxis_type="log"))

plot([tp1,tp2,tp3,tp4,tp5,tp6,tp7,tp8,tp9], Layout(title="Timepoint of max val", xaxis_title="kdam", yaxis_title="time"))





















kdam=0
kin=0.054
lam=0.033
atp=4000
params1 = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)

prob = ODEProblem(rtc_model, initial, tspan, params1)
solu = solve(prob, Rodas4())

max_, tp_max = get_max_val(solu, :rm_a)

df = DataFrame(solu)
rename!(df, [:time, :rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt])
df
max_, ind = findmax(df[:, :rt])

df[24, :timestamp]