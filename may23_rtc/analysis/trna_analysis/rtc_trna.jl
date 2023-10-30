using Parameters, CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, Statistics
using Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra, Printf
# using Plots
using PlotlyJS, ProgressBars
include("/home/holliehindley/phd/may23_rtc/functions/solving.jl"); include("/home/holliehindley/phd/may23_rtc/functions/set_ups.jl"); include("/home/holliehindley/phd/may23_rtc/functions/plotting.jl"); 
include("/home/holliehindley/phd/may23_rtc/functions/sweep_params.jl"); include("/home/holliehindley/phd/may23_rtc/models/rtc_orig.jl"); include("/home/holliehindley/phd/may23_rtc/models/atp_lam_kin_t.jl"); 
include("/home/holliehindley/phd/may23_rtc/models/single_t.jl"); include("/home/holliehindley/phd/may23_rtc/models/combinations_t.jl"); 
include("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/bf_funcs.jl");
include("/home/holliehindley/phd/colors_plotly.jl")
include("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/init_switch/funcs.jl"); include("/home/holliehindley/phd/may23_rtc/models/rtc_inhibition_model.jl");

@consts begin
    L = 10; #10 
    c = 0.001; 
    kr = 0.125; 
    Vmax_init = 39.51; 
    Km_init = 250; 
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
    kin = 0.022222222 #0.054; #2.381 
    ω_ab = 4#4#0.093; #0.0828304057748932;#4; 
    ω_r = 0.0019*6#2e-7 #0.0019*6 #70.53; #0.0019*6#79.43865871861044; #0.0019*6;  
    ω_a = 4; 
    ω_b = 4;
    # kdam =  0.#0.000147;#0.05; 
    k = 2; # carrying capacity - changes depending on the data?
    # lam = 0.033;

    # rtca_0 = 0#0.00894; 
    # rtcb_0 = 0#0.0216; 
    # rh_0 = 11.29; #69.56; #69.4
    # rtcr_0 = 0# 0.0131 #0.04; # 8.67e-3; # change this based on keeping steady state level the whole time course (levels shouldn't really change)
    # rm_a_0 = 0; 
    # rm_b_0 = 0; 
    # rm_r_0 = 0#0.0131#0.04 # 0; 
    # rd_0 = 0; 
    # rt_0 = 0;
end

tspan = (0,1e9)
rh = 11.29 #75 # conc of ribosomes in exponential phase 
thr_t = 5#30 # was at 5 before to get saved plots # needs to be less than 30 
kin_trna = 1

params_orig = @LArray [10., c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 0.5, ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
solu_orig = sol(rtc_model, initial, tspan, params1)
df_orig = create_solu_df(solu_orig, all_species)
p2 = plotly_plot_sol(solu_orig, all_species, "", "", "damaged ribosomes")

trna_species = [:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :trna, :rd, :rt]
init_trna = [0,0,0,0,0,0,135.5,0,0] # tRNA initial conc = 135.5
params_trna = @LArray [10., c, kr*12, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 0.5, ktag, kdeg, kin_trna, 3578.9473684210525, na, nb, nr, 0.014, rh, thr_t] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam, :rh, :thr_t)
solu_trna = sol(rtc_model_trna, init_trna, tspan, params_trna)
p_trna = plotly_plot_sol(solu_trna, trna_species, "", "", "damaged tRNA")

[p2 p_trna]


rh_range = range(0,100,length=100)
res = change_param(rh_range, :rh, rtc_model_trna, init_trna, get_ssval, trna_species, params_trna)

plot([scatter(x=rh_range, y=res[:trna], name="tRNA"),scatter(x=rh_range, y=res[:rtca], name="RtcA")], Layout(xaxis_title="Rh", yaxis_title="conc", title="θ tRNA = 5"))

kin_range = range(0,2,length=100)
res_kin = change_param(kin_range, :kin, rtc_model_trna, init_trna, get_ssval, trna_species, params_trna)

plot(scatter(x=kin_range, y=res_kin[:trna]))

thrt_range = range(0,100,length=100)
res_t = change_param(thrt_range, :thr_t, rtc_model_trna, init_trna, get_ssval, trna_species, params_trna)

plot(scatter(x=thrt_range, y=res_t[:trna]), Layout(xaxis_title="θ tRNA", yaxis_title="tRNA")) # shows that thr_t needs to be less than 30 

kr_range = range(0.0125,1.25,length=1000)
res_kr = change_param(kr_range, :kr, rtc_model_trna, init_trna, get_ssval, trna_species, params_trna)
plot(scatter(x=kr_range, y=res_kr[:rtca]), Layout(xaxis_title="Kr", yaxis_title="tRNA")) # shows that thr_t needs to be less than 30 


kdam_range = range(0,4,length=1000)
res = change_param(kdam_range, :kdam, rtc_model_trna, init_trna, get_ssval, trna_species, params_trna)

trna_h_p = (scatter(x=kdam_range, y=res[:trna], name="tRNA"))
rtca_p = scatter(x=kdam_range, y=res[:rtca], name="RtcA")
rtcb_p = scatter(x=kdam_range, y=res[:rtcb], name="RtcB")
mrna_p = scatter(x=kdam_range, y=res[:rm_a], name="mRNA AB")
rtcr_p = scatter(x=kdam_range, y=res[:rtcr], name="RtcR")
trna_d_p = scatter(x=kdam_range, y=res[:rd], name="tRNA damaged")
trna_t_p = scatter(x=kdam_range, y=res[:rt], name="tRNA tagged")

trna_kdam = plot([trna_h_p, rtca_p, rtcb_p, mrna_p, rtcr_p, trna_d_p, trna_t_p])





params_trna2 = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923,
kdeg = 0.001, kin = kin_trna, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
kdam =  0.01, lam = 0.014, rh = rh, thr_t = thr_t)

include("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/bf_funcs.jl");

br = get_br(rtc_mod_trna, params_trna2, init_trna, 5.)
# bs = plot_all_curves_bistable(br, colors2, colorsr, "Original, L = 10")

df = create_br_df(br)

trna_h_p1 = plot(scatter(x=df.kdam, y=df.rh, name="tRNA"))
rtca_p1 = scatter(x=df.kdam, y=df.rtca, name="RtcA")
rtcb_p1 = scatter(x=df.kdam, y=df.rtcb, name="RtcB")
mrna_p1 = scatter(x=df.kdam, y=df.rm_a, name="mRNA AB")
rtcr_p1 = scatter(x=df.kdam, y=df.rtcr, name="RtcR")
trna_d_p1 = scatter(x=df.kdam, y=df.rd, name="tRNA damaged")
trna_t_p1 = scatter(x=df.kdam, y=df.rt, name="tRNA tagged")

plot([trna_h_p1, rtca_p1, rtcb_p1, mrna_p1, rtcr_p1, trna_d_p1, trna_t_p1])



kdam_range1 = range(0,3,length=1000)
res1 = change_param(kdam_range1, :kdam, rtc_model, initial, get_ssval, all_species, params_orig)

rh_p = scatter(x=kdam_range1, y=res1[:rh], name="Rh")
rtca_p = scatter(x=kdam_range1, y=res1[:rtca], name="RtcA")
rtcb_p = scatter(x=kdam_range1, y=res1[:rtcb], name="RtcB")
mrna_p = scatter(x=kdam_range1, y=res1[:rm_a], name="mRNA AB")
rtcr_p = scatter(x=kdam_range1, y=res1[:rtcr], name="RtcR")
trna_d_p = scatter(x=kdam_range1, y=res1[:rd], name="Rd")
trna_t_p = scatter(x=kdam_range1, y=res1[:rt], name="Rt")

plot([rh_p, rtca_p, rtcb_p, mrna_p, rtcr_p, trna_d_p, trna_t_p])








s1 = sol(rtc_model_trna, init_trna, tspan, params_trna)
s2 = sol(rtc_mod_trna!, init_trna, tspan, params_trna)

p1 = plotly_plot_sol(s1, trna_species, "", "","")
p2 = plotly_plot_sol(s2, trna_species, "", "","")

[p1 p2]





function change_param_get_tlr(param_range, parameter, model, init, species, params)
    tlr = []
    tlrel = []
    new_params = deepcopy(params)
    for val in param_range
        new_params[parameter] = val
        solu = sol(model, init, tspan, new_params)
        if length(new_params) == 23
            tlrel_orig = new_params[:g_max]*new_params[:atp]/(new_params[:θtlr]+new_params[:atp])
            # @show tlrel_orig
            rma_orig = get_ssval(solu, :rm_a, species)
            rh_orig = get_ssval(solu, :rh, species)
            tlr_orig = @. (1/new_params[:na])*rh_orig*rma_orig*tlrel_orig
        else
            trna = get_ssval(solu, :trna, species)
            tlrel_orig = @. (new_params[:g_max]*new_params[:atp]/(new_params[:θtlr]+new_params[:atp])) * trna/(new_params[:thr_t]+trna)
            # @show tlrel_orig
            rma_orig = get_ssval(solu, :rm_a, species)
            rh_orig = get_ssval(solu, :trna, species)
            tlr_orig = @. (1/new_params[:na])*rh_orig*rma_orig*tlrel_orig
        end
        push!(tlr, tlr_orig)
        push!(tlrel, tlrel_orig)
    end
    return tlr, tlrel
end

tlr_orig, tlrel_orig = change_param_get_tlr(kdam_range, :kdam, rtc_model, initial, all_species, params_orig)
tlrel_orig
plot(scatter(x=kdam_range, y=tlr_orig))


tlr_trna, tlrel_trna = change_param_get_tlr(rh_range, :rh, rtc_model_trna, init_trna, trna_species, params_trna)

plot(scatter(x=rh_range, y=tlr_trna))
plot(scatter(x=kdam_range, y=tlrel_trna))



sweep_paramx2(rtc_model_trna, params_trna, trna_species, :trna, get_ssval, :kdam, :rh, kdam_range, rh_range)
sweep_paramx2(rtc_model_trna, params_trna, trna_species, :trna, get_ssval, :kdam, :thr_t, kdam_range, thrt_range)
sweep_paramx2(rtc_model_trna, params_trna, trna_species, :trna, get_ssval, :kdam, :kin, kdam_range, kin_range)




kdam_range = range(0,400,length=1000)
kdam_range2 = range(400,0,length=1000)
kdam_range_orig = range(0,3,length=1000)
kdam_range_orig2 = range(3,0,length=1000)

res = change_param(kdam_range, :kdam, rtc_model_trna, init_trna, get_ssval, trna_species, params_trna)
trna_h_p = (scatter(x=kdam_range, y=res[:trna], name="Vary kdam", legendgroup=1, line=attr(color="#2ca02c")))

res_orig = change_param(kdam_range_orig, :kdam, rtc_model, initial, get_ssval, all_species, params_orig)
rh_p = (scatter(x=kdam_range_orig, y=res_orig[:rh], name="Vary kdam", legendgroup=1, line=attr(color="#2ca02c")))



params2 = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
kdeg = 0.001, kin = 0.022222222, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
kdam =  0., lam = 0.014)

br_orig = get_br(rtc_mod, params2, initial, 3.)
df_orig = create_br_df(br_orig)
bifurc_orig = (scatter(x=df_orig.kdam, y=df_orig.rh, name="BifurcationKit", legendgroup=2, line=attr(color="#1f77b4")))

params_trna2 = (L = 10., c = 0.001, kr = 0.125*12, Vmax_init = 39.51, Km_init = 250.,
θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923,
kdeg = 0.001, kin = kin_trna, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
kdam =  0.01, lam = 0.014, rh = rh, thr_t = thr_t)

br = get_br(rtc_mod_trna, params_trna2, init_trna, 400.) # dsmax = 0.05, 0.08, 0.09 for to finish bifurcation at 0 
br1 = get_br(rtc_mod_trna, params_trna2, init_trna, 400.) # dsmax > 0.1 to finish bifurcation at kdam = 400 (see negative values)  
Plots.plot(br, vars=(:param, :trna), linewidthstable=5)

df = create_br_df(br)
df1 = create_br_df(br1)
bf = bf_point_df(br)
bf1 = bf_point_df(br1)
kdam1 = findall(x->x==bf.kdam[1],df.kdam)[1]
kdam2 = findall(x->x==bf1.kdam[1],df1.kdam)[1]
stable_trna = df.rh[1:kdam1]
unstable_trna = df.rh[kdam1:end]
stable_trna1 = df1.rh[1:kdam2]
unstable_trna1 = df1.rh[kdam2:end]
stable1 = scatter(x=df.kdam[1:kdam1],y=stable_trna, line=attr(color="#1f77b4"), name="BifurcationKit to 0", legendgroup=5)
unstable1 = scatter(x=df.kdam[kdam1:end],y=unstable_trna, line=attr(color="#1f77b4", dash="dash"), name="BifurcationKit to 0", legendgroup=5)
stable2 = scatter(x=df1.kdam[1:kdam2],y=stable_trna1, line=attr(color="#ff7f0e"), name="BifurcationKit to 400", legendgroup=6)
unstable2 = scatter(x=df1.kdam[kdam2:end],y=unstable_trna1, line=attr(color="#ff7f0e", dash="dash"), name="BifurcationKit to 400", legendgroup=6)
bf_p = (scatter(x=[df.kdam[kdam1]], y=[df.rh[kdam1]], legendgroup=5, line=attr(color="#1f77b4"), mode="markers", showlegend=false))
bf_p1 = (scatter(x=[df1.kdam[kdam2]], y=[df1.rh[kdam2]], legendgroup=6, line=attr(color="#ff7f0e"), mode="markers", showlegend=false))
# trna_h_p1 = (scatter(x=df.kdam, y=df.rh, name="BifurcationKit", legendgroup=2))
# trna_h_p2 = (scatter(x=df1.kdam, y=df1.rh, name="BifurcationKit", legendgroup=7))
plot([stable1, unstable1], Layout(title="θ tRNA = 30"))

plot([trna_h_p1, trna_h_p2, stable1, unstable1, stable2, unstable2, bf_p, bf_p1])


res_orig1 = numerical_bistability_analysis(rtc_model, params_orig, initial, :rh, all_species, kdam_range_orig)
res_orig2 = numerical_bistability_analysis(rtc_model, params_orig, initial, :rh, all_species, kdam_range_orig2)
porig1 = scatter(x=kdam_range_orig, y=res_orig1, name="↑ kdam", legendgroup=3, line=attr(color="#e377c2"))
porig2 = scatter(x=kdam_range_orig2, y=res_orig2, name="↓ kdam", legendgroup=4, line=attr(color="#9467bd"))

orig = plot([porig1, porig2, rh_p, bifurc_orig], Layout(title="original model - Rh"))

params_trna = @LArray [10., c, kr*12, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 0.5, ktag, kdeg, kin_trna, 3578.9473684210525, na, nb, nr, 0.014, rh, thr_t] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam, :rh, :thr_t)
res_trna1 = numerical_bistability_analysis(rtc_model_trna, params_trna, init_trna, :trna, trna_species, kdam_range)
res_trna2 = numerical_bistability_analysis(rtc_model_trna, params_trna, init_trna, :trna, trna_species, kdam_range2)
ptrna1 = scatter(x=kdam_range, y=res_trna1, name="↑ kdam", legendgroup=3, line=attr(color="#e377c2"))
ptrna2 = scatter(x=kdam_range2, y=res_trna2, name="↓ kdam", legendgroup=4, line=attr(color="#9467bd"))
trna = plot([ptrna1, ptrna2])
trna2 = plot([ptrna1, ptrna2, trna_h_p, stable1, unstable1, stable2, unstable2, bf_p, bf_p1], Layout(title="tRNA model - tRNA"))

p = [orig trna2]

open("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/ss_proof.html", "w") do io
    PlotlyBase.to_html(io, p.plot)
end



kdam = 300
params_trna = @LArray [10., c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin_trna, 3578.9473684210525, na, nb, nr, 0.014, rh, thr_t] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam, :rh, :thr_t)
solu_trna = sol(rtc_model_trna, init_trna, tspan, params_trna)
ssvals = ss_init_vals(solu_trna, trna_species)
ssvals2 = ssvals*1000000 # no matter what the initial conditions of the model here the ss will always be at zero 
solu_trna2 = sol(rtc_model_trna, ssvals2, tspan, params_trna)
ss_init_vals(solu_trna2, trna_species) 


kdam2 = 250
params_trna3 = @LArray [10., c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam2, ktag, kdeg, kin_trna, 3578.9473684210525, na, nb, nr, 0.014, rh, thr_t] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam, :rh, :thr_t)
solu_trna3 = sol(rtc_model_trna, init_trna, tspan, params_trna3)
ssvals3 = ss_init_vals(solu_trna3, trna_species)
ssvalsA = ssvals3*1000000
solu_trna2 = sol(rtc_model_trna, ssvalsA, tspan, params_trna3)
ss_init_vals(solu_trna2, trna_species) 



params_trna2 = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923,
kdeg = 0.001, kin = kin_trna, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
kdam =  0.01, lam = 0.014, rh = rh, thr_t = thr_t)

br = get_br(rtc_mod_trna, params_trna2, init_trna, 400.)

df = create_br_df(br)
bf = bf_point_df(br)
kdam1 = findall(x->x==bf.kdam[1],df.kdam)[1]
unstable_firsthalf = df[!,:rh][kdam1:end][1:2094]
unstable1 = plot(scatter(x=df.kdam[kdam1:end],y=unstable_firsthalf, line=attr(color="#1f77b4", dash="dash"),showlegend=false))


kdam_range_short = range(20,350,length=10)
new_ps = deepcopy(params_trna)
first_init = DataFrame(rm_a=[],rtca=[],rm_b=[],rtcb=[],rm_r=[],rtcr=[],trna=[],rd=[],rt=[]); 
new_init = DataFrame(rm_a=[],rtca=[],rm_b=[],rtcb=[],rm_r=[],rtcr=[],trna=[],rd=[],rt=[]); 
final_ss = DataFrame(rm_a=[],rtca=[],rm_b=[],rtcb=[],rm_r=[],rtcr=[],trna=[],rd=[],rt=[]);
for i in kdam_range_short
    for (s, f, n, l) in zip(range(1,9),eachcol(first_init), eachcol(new_init), eachcol(final_ss))
        new_ps.kdam = i    
        solu = sol(rtc_model_trna, init_trna, tspan, new_ps)
        ssvals = ss_init_vals(solu, trna_species)
        add_vals=[]
        for spec in all_species
            push!(add_vals, 0.1*maximum(df[!,spec]))
        end
        ssvalsA = ssvals+add_vals
        solu2 = sol(rtc_model_trna, ssvalsA, tspan, new_ps)
        new_ssvals = ss_init_vals(solu2, trna_species)
        push!(f, ssvals[s])
        push!(n, ssvalsA[s])
        push!(l, new_ssvals[s])
    end
end

function plot_bf_species(br, specie, specie2, kdam_range2, kdam_range_short, first_init, new_init, final_ss)
    df = create_br_df(br)
    bf = bf_point_df(br)
    kdam1 = findall(x->x==bf.kdam[1],df.kdam)[1]
    stable_trna = df[!,specie][1:kdam1]
    unstable_trna = df[!,specie][kdam1:end]
    stable1 = (scatter(x=df.kdam[1:kdam1],y=stable_trna, line=attr(color="#1f77b4"),showlegend=false))
    unstable1 = (scatter(x=df.kdam[kdam1:end],y=unstable_trna, line=attr(color="#1f77b4", dash="dash"),showlegend=false))
    res_trna2 = checking_bistability(rtc_model_trna, params_trna, init_trna, specie2, trna_species, kdam_range2)
    zeross = scatter(x=kdam_range2, y=res_trna2, name="↓ kdam", legendgroup=4, line=attr(color="#9467bd"),showlegend=false)
    diffs_ss = []
    for (i,j, kdam_val) in zip(first_init[!,specie2], new_init[!,specie2], kdam_range_short)
        push!(diffs_ss,scatter(x=[kdam_val,kdam_val], y=[i,j],line=attr(color="#e377c2"),showlegend=false))
    end
    final_ss_p = scatter(x=kdam_range_short, y=final_ss[!,specie2],showlegend=false)
    plot([stable1, unstable1, zeross, diffs_ss[1], diffs_ss[2],
     diffs_ss[3], diffs_ss[4], diffs_ss[5], diffs_ss[6], diffs_ss[7], diffs_ss[8], 
    diffs_ss[9], diffs_ss[10], final_ss_p], Layout(title="$specie"))
end

prma = plot_bf_species(br, :rm_a, :rm_a, kdam_range2, kdam_range_short, first_init, new_init, final_ss)
prtca = plot_bf_species(br, :rtca, :rtca, kdam_range2, kdam_range_short, first_init, new_init, final_ss)    
prmb = plot_bf_species(br, :rm_b, :rm_b, kdam_range2, kdam_range_short, first_init, new_init, final_ss)
prtcb = plot_bf_species(br, :rtcb, :rtcb, kdam_range2, kdam_range_short, first_init, new_init, final_ss)     
prtcr = plot_bf_species(br, :rtcr, :rtcr, kdam_range2, kdam_range_short, first_init, new_init, final_ss) 
ptrna = plot_bf_species(br, :rh, :trna, kdam_range2, kdam_range_short, first_init, new_init, final_ss)  
prd = plot_bf_species(br, :rd, :rd, kdam_range2, kdam_range_short, first_init, new_init, final_ss)    
prt = plot_bf_species(br, :rt, :rt, kdam_range2, kdam_range_short, first_init, new_init, final_ss) 
    
p = [prma prtca prmb prtcb; prtcr ptrna prd prt]










params_trna2 = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923,
kdeg = 0.001, kin = kin_trna, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
kdam =  0.01, lam = 0.014, rh = rh, thr_t = thr_t)

br = get_br(rtc_mod_trna, params_trna2, init_trna, 400.)

rh_range = range(0,100,length=10)
kin_range = range(0,2,length=10)
thrt_range = range(0,100,length=10)

function bistable_search(param, param_range)
    params_trna = deepcopy(params_trna2)
    for i in param_range
        params_trna = merge(params_trna, (param=>i,))
        br = get_br(rtc_mod_trna, params_trna, init_trna, 400.)
        if length(br.specialpoint) >= 3
            display(Plots.plot(br, vars=(:param, :trna), linewidthstable=5))
            @show i
        end
    end
end
bistable_search(:rh, rh_range)
bistable_search(:kin, kin_range)
bistable_search(:thr_t, thrt_range)

