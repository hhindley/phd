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

atp=3578.9473684210525

params_orig = @LArray [10., c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 0.5, ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
solu_orig = sol(rtc_model, initial, tspan, params1)
df_orig = create_solu_df(solu_orig, all_species)
tlrel_orig = (g_max*atp/(θtlr+atp))
tlr_orig = @. (1/na)*df_orig.rh*df_orig.rm_a*tlrel_orig
plot(scatter(x=df_orig.rh,y=tlr_orig))

trna_species = [:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :trna, :rd, :rt]
init_trna = [0,0,0,0,0,0,135.5,0,0] # tRNA initial conc = 135.5
params_trna = @LArray [10., c, kr*12, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 0.5, ktag, kdeg, kin_trna, 3578.9473684210525, na, nb, nr, 0.014, rh, thr_t] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam, :rh, :thr_t)
solu_trna = sol(rtc_model_trna, init_trna, tspan, params_trna)
df_trna = create_solu_df(solu_trna, trna_species)
tlrel_trna = @. (g_max*atp/(θtlr+atp)) * df_trna.trna/(thr_t+df_trna.trna) 
tlr_trna = @. (1/na)*rh*df_trna.rm_a*tlrel_trna
plot(scatter(x=df_trna.trna,y=tlrel_trna))


rh_range = range(0,100,length=100)
res=[];
for i in rh_range
    push!(res, (1/na)*i*0.001*tlrel_orig)
end
res
plot(scatter(x=rh_range,y=res))