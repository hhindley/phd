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
    # kin = 0.054; #2.381 
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
thr_t = 5 # needs to be less than 30 
kin_trna = 1



params_trna2 = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923,
kdeg = 0.001, kin = kin_trna, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
kdam =  0.01, lam = 0.014, rh = rh, thr_t = thr_t)


trna_species = [:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :trna, :rd, :rt]
init_trna = [0,0,0,0,0,0,135.5,0,0] # tRNA initial conc = 135.5

rh_range = range(0,100,length=100)
kin_range = range(0,2,length=100)
thrt_range = range(0,100,length=100)

function bistable_search(param, param_range)
    params_trna = deepcopy(params_trna2)
    for i in param_range
        params_trna = merge(params_trna, (param=>i,))
        br = get_br(rtc_mod_trna, params_trna, init_trna, 5.)
        if length(br.specialpoint) == 4
            @show i 
        end
    end
end
bistable_search(:rh, rh_range)
bistable_search(:kin, kin_range)
bistable_search(:thr_t, thrt_range)

kin_trna = 0.06060606060606061
params_trna2 = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923,
kdeg = 0.001, kin = kin_trna, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
kdam =  0.01, lam = 0.014, rh = rh, thr_t = thr_t)

br = get_br(rtc_mod_trna, params_trna2, init_trna, 4.)
# bs = plot_all_curves_bistable(br, colors2, colorsr, "Original, L = 10")

df = create_br_df(br)

trna_h_p1 = (scatter(x=df.kdam, y=df.rh, name="tRNA"))
rtca_p1 = scatter(x=df.kdam, y=df.rtca, name="RtcA")
rtcb_p1 = scatter(x=df.kdam, y=df.rtcb, name="RtcB")
mrna_p1 = scatter(x=df.kdam, y=df.rm_a, name="mRNA AB")
rtcr_p1 = scatter(x=df.kdam, y=df.rtcr, name="RtcR")
trna_d_p1 = scatter(x=df.kdam, y=df.rd, name="tRNA damaged")
trna_t_p1 = scatter(x=df.kdam, y=df.rt, name="tRNA tagged")

bistab = plot([trna_h_p1, rtca_p1, rtcb_p1, mrna_p1, rtcr_p1, trna_d_p1, trna_t_p1])



params_trna = @LArray [10., c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 0.01, ktag, kdeg, kin_trna, 3578.9473684210525, na, nb, nr, 0.014, rh, thr_t] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam, :rh, :thr_t)
solu_trna = sol(rtc_model_trna, init_trna, tspan, params_trna)
p_trna = plotly_plot_sol(solu_trna, trna_species, "", "", "damaged tRNA")




kdam_range = range(0,4,length=1000)
res = change_param(kdam_range, :kdam, rtc_model_trna, init_trna, get_ssval, trna_species, params_trna)

trna_h_p = (scatter(x=kdam_range, y=res[:trna], name="tRNA"))
rtca_p = scatter(x=kdam_range, y=res[:rtca], name="RtcA")
rtcb_p = scatter(x=kdam_range, y=res[:rtcb], name="RtcB")
mrna_p = scatter(x=kdam_range, y=res[:rm_a], name="mRNA AB")
rtcr_p = scatter(x=kdam_range, y=res[:rtcr], name="RtcR")
trna_d_p = scatter(x=kdam_range, y=res[:rd], name="tRNA damaged")
trna_t_p = scatter(x=kdam_range, y=res[:rt], name="tRNA tagged")

kdam_vary = plot([trna_h_p, rtca_p, rtcb_p, mrna_p, rtcr_p, trna_d_p, trna_t_p])



[bistab kdam_vary]