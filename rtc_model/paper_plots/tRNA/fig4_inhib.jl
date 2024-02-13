using Parameters, CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, Statistics
using Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra, Printf
# using Plots
using PlotlyJS, ProgressBars
include("/home/holliehindley/phd/may23_rtc/functions/solving.jl"); include("/home/holliehindley/phd/may23_rtc/functions/set_ups.jl"); include("/home/holliehindley/phd/may23_rtc/functions/plotting.jl"); 
include("/home/holliehindley/phd/may23_rtc/functions/sweep_params.jl"); include("/home/holliehindley/phd/may23_rtc/models/rtc_orig.jl"); include("/home/holliehindley/phd/may23_rtc/models/atp_lam_kin_t.jl"); 
include("/home/holliehindley/phd/may23_rtc/models/single_t.jl"); include("/home/holliehindley/phd/may23_rtc/models/combinations_t.jl"); 
include("/home/holliehindley/phd/may23_rtc/functions/bf_funcs/bf_funcs.jl");
include("/home/holliehindley/phd/colors_plotly.jl"); include("/home/holliehindley/phd/may23_rtc/models/inhibition_models/rtc_inhibition_model.jl");
include("/home/holliehindley/phd/may23_rtc/models/inhibition_models/trna_inhib.jl"); include("/home/holliehindley/phd/may23_rtc/models/rtc_trna_model.jl");

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
thr_t = 5#30 # was at 5 before to get saved plots # needs to be less than 30 
kin_trna = 1

k_inhib1 = 0.3
k_inhib2 = 0.0025
inhib = 0.1
k_inhib1a = 0.5
k_inhib1b = 1


trna_species_inhib = [:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :trna, :rd, :rt, :rtc_i]
init_trna_inhib = [0,0,0,0,0,0,135.5,0,0,0] # tRNA initial conc = 135.5
init_trna = [0,0,0,0,0,0,135.5,0,0] # tRNA initial conc = 135.5
params_trna_inhib = @LArray [10., c, kr*12, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 0.5, ktag, kdeg, kin_trna, 3578.9473684210525, na, nb, nr, 0.014, rh, thr_t,k_inhib1, k_inhib2, inhib] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam, :rh, :thr_t, :k_inhib1, :k_inhib2, :inhib)
params_trna = @LArray [10., c, kr*12, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 0.5, ktag, kdeg, kin_trna, 3578.9473684210525, na, nb, nr, 0.014, rh, thr_t] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam, :rh, :thr_t)

kdam_range = range(0,400,length=1000)
kdam_range2 = range(400,0,length=1000)


initial_i = [0,0,0,0,0,0,11.29,0,0,0]

params_for_ssval_setup_inhib = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
kdeg = 0.001, kin = 0.022222222, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
kdam =  0.01, lam = 0.014, k_inhib1=k_inhib1, k_inhib2=k_inhib2, inhib=inhib)

params2 = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
kdeg = 0.001, kin = 0.022222222, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
kdam =  0.01, lam = 0.014)


function different_levels_inhibition_trna(rtc_inhib_mod, k_inhib1, k_inhib2, inhib)
        
    params_trna_inhib = (L = 10., c = 0.001, kr = 0.125*12, Vmax_init = 39.51, Km_init = 250.,
    θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
    krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
    kdeg = 0.001, kin = kin_trna, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
    kdam =  0.01, lam = 0.014, rh=rh, thr_t=thr_t, k_inhib1=k_inhib1, k_inhib2=k_inhib2, inhib=inhib)

    br = get_br(rtc_inhib_mod, params_trna_inhib, init_trna_inhib, 400.)
    bf = bf_point_df_inhib(br)
    df = create_br_df_inhib(br)
    # kdam1 = findall(x->x==bf.kdam[1],df.kdam)[1]
    # kdam2 = findall(x->x==bf.kdam[2],df.kdam)[1]
    return bf, df
end
function trna_plotting_inhib(bf, df, specie, color, i)
    # res_trna1 = numerical_bistability_analysis(rtc_model_trna, params_trna, init_trna, :trna, trna_species, kdam_range)
    res_trna2 = numerical_bistability_analysis(rtc_model_trna, params_trna, init_trna, specie, all_species, kdam_range2)
    # ptrna1 = scatter(x=kdam_range, y=res_trna1, name="Healthy tRNA", legendgroup=3, yaxis="y2", line=attr(color=:gold,linewidth=3))
    ptrna2 = scatter(x=kdam_range2, y=res_trna2, name="", legendgroup=i, showlegend=false, line=attr(color=color,width=5))
# res_rtcb1 = numerical_bistability_analysis(rtc_model_trna, params_trna, init_trna, :rtcb, trna_species, kdam_range)
# res_rtcb2 = numerical_bistability_analysis(rtc_model_trna, params_trna, init_trna, :rtcb, trna_species, kdam_range2)
# prtcb1 = scatter(x=kdam_range, y=res_rtcb1, name="RtcB", legendgroup=3, line=attr(color=:plum,linewidth=3))
# prtcb2 = scatter(x=kdam_range2, y=res_rtcb2, name="", showlegend=false, legendgroup=3, line=attr(color=:plum,linewidth=3))
    # df = create_br_df(br)
    # df1 = create_br_df(br1)
    # bf = bf_point_df(br)
    # bf1 = bf_point_df(br1)
    kdam1 = findall(x->x==bf.kdam[1],df.kdam)[1]
    # kdam2 = findall(x->x==bf1.kdam[1],df1.kdam)[1]
    stable_trna = df[!,specie][1:kdam1]
    unstable_trna = df[!,specie][kdam1:end]
    stable1 = scatter(x=df.kdam[1:kdam1],y=stable_trna, line=attr(color=color, width=5), name="$specie", legendgroup=i)
    unstable1 = scatter(x=df.kdam[kdam1:end],y=unstable_trna, line=attr(color=color, width=5, dash="dash"), name="", showlegend=false, legendgroup=i)
    return ptrna2, stable1, unstable1
end
function creating_rtc_inhib_plot_trna(rtc_inhib_mod, specie)
    params_trna2 = (L = 10., c = 0.001, kr = 0.125*12, Vmax_init = 39.51, Km_init = 250.,
    θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
    krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923,
    kdeg = 0.001, kin = kin_trna, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
    kdam =  0.01, lam = 0.014, rh = rh, thr_t = thr_t)



    br = get_br(rtc_mod_trna, params_trna2, init_trna, 400.)
    bf0 = bf_point_df(br)
    df0 = create_br_df(br)
    kdam01 = findall(x->x==bf0.kdam[1],df0.kdam)[1]
    # kdam02 = findall(x->x==bf0.kdam[2],df0.kdam)[1]

    rtcb_01, rtcb_02, rtcb_03 = trna_plotting_inhib(bf0, df0, specie, :blue, 1)


    bf, df = different_levels_inhibition_trna(rtc_inhib_mod, k_inhib1, k_inhib2, inhib)
    rtcb1, rtcb2, rtcb3 = trna_plotting_inhib(bf, df, specie, :purple, 2)

    bfa, dfa = different_levels_inhibition_trna(rtc_inhib_mod, k_inhib1a, k_inhib2, inhib)
    rtcb1a, rtcb2a, rtcb3a = trna_plotting_inhib(bfa, dfa, specie, :orange, 3)

    bfb, dfb = different_levels_inhibition_trna(rtc_inhib_mod, k_inhib1b, k_inhib2, inhib)
    rtcb1b, rtcb2b, rtcb3b = trna_plotting_inhib(bfb, dfb, specie, :red, 4)

    return [rtcb_01, rtcb_02, rtcb_03, rtcb1, rtcb2, rtcb3, rtcb1a, rtcb2a, rtcb3a, rtcb1b, rtcb2b, rtcb3b]
end



rtcb_traces = creating_rtc_inhib_plot_trna(rtc_trna_inhib_mod_rtcb, :rtcb)
p_rtcb = plot([i for i in rtcb_traces],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="RtcB (μM)",
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white", font=attr(size=20, color="black", family="sans-serif")))


rtca_traces = creating_rtc_inhib_plot_trna(rtc_trna_inhib_mod_rtca, :rtca)
p_rtca = plot([i for i in rtca_traces],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="RtcA (μM)",
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=20, color="black", family="sans-serif")))

rtcr_traces = creating_rtc_inhib_plot_trna(rtc_trna_inhib_mod_rtcr, :rtcr)
p_rtcr = plot([i for i in rtcr_traces],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="RtcR (μM)", yaxis_range=(-0,0.023),
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=20, color="black", family="sans-serif")))



rtcb_rh_traces = creating_rtc_inhib_plot_trna(rtc_trna_inhib_mod_rtcb, :rh)
p_rtcb_rh = plot([i for i in rtcb_rh_traces],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="Rh (μM)",
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=20, color="black", family="sans-serif")))

rtca_rh_traces = creating_rtc_inhib_plot_trna(rtc_trna_inhib_mod_rtca, :rh)
p_rtca_rh = plot([i for i in rtca_rh_traces],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="Rh (μM)",
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=20, color="black", family="sans-serif")))

rtcr_rh_traces = creating_rtc_inhib_plot_trna(rtc_trna_inhib_mod_rtcr, :rh)
p_rtcr_rh = plot([i for i in rtcr_rh_traces],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="Rh (μM)",
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=20, color="black", family="sans-serif")))


[p_rtcb p_rtcb_rh; p_rtca p_rtca_rh; p_rtcr p_rtcr_rh]


savefig(p_rtcb, "/home/holliehindley/phd/may23_rtc/paper_plots/tRNA/rtcb_inhib.svg")
savefig(p_rtca, "/home/holliehindley/phd/may23_rtc/paper_plots/tRNA/rtca_inhib.svg")
savefig(p_rtcr, "/home/holliehindley/phd/may23_rtc/paper_plots/tRNA/rtcr_inhib.svg")
savefig(p_rtcb_rh, "/home/holliehindley/phd/may23_rtc/paper_plots/tRNA/rtcb_rh_inhib.svg")
savefig(p_rtca_rh, "/home/holliehindley/phd/may23_rtc/paper_plots/tRNA/rtca_rh_inhib.svg")
savefig(p_rtcr_rh, "/home/holliehindley/phd/may23_rtc/paper_plots/tRNA/rtcr_rh_inhib.svg")





function bf_size(rtc_inhib_mod)
    params2 = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
    θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
    krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
    kdeg = 0.001, kin = 0.022222222, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
    kdam =  0.01, lam = 0.014)

    br = get_br(rtc_mod, params2, initial, 3.)
    bf0 = bf_point_df(br)
    df0 = create_br_df(br)
    kdam01 = findall(x->x==bf0.kdam[1],df0.kdam)[1]
    kdam02 = findall(x->x==bf0.kdam[2],df0.kdam)[1]

    bf, df, kdam1, kdam2 = different_levels_inhibition(rtc_inhib_mod, k_inhib1, k_inhib2, inhib)

    bfa, dfa, kdam1a, kdam2a = different_levels_inhibition(rtc_inhib_mod, k_inhib1a, k_inhib2, inhib)

    bfb, dfb, kdam1b, kdam2b = different_levels_inhibition(rtc_inhib_mod, k_inhib1b, k_inhib2, inhib)

    s0 = bf0.kdam[1]-bf0.kdam[2]
    s1 = bf.kdam[1]-bf.kdam[2]
    s2 = bfa.kdam[1]-bfa.kdam[2]
    s3 = bfb.kdam[1]-bfb.kdam[2]

    rtcb_sizes = [s0, s2, s3, s1]
    percentage_of_original_size = [100*(rtcb_sizes[i]/rtcb_sizes[1]) for i in range(2,4)]
    # percentage_decrease_rtcb = [100*((rtcb_sizes[i]-rtcb_sizes[1])/rtcb_sizes[1]) for i in range(2,4)]
    return percentage_of_original_size
    # return percentage_decrease_rtcb
end

percentage_size_rtcb = bf_size(rtc_inhib_mod_rtcb)
percentage_size_rtca = bf_size(rtc_inhib_mod_rtca)
percentage_size_rtcr = bf_size(rtc_inhib_mod_rtcr)


function protein_decrease(rtc_inhib_mod, specie)
    params2 = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
    θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
    krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
    kdeg = 0.001, kin = 0.022222222, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
    kdam =  0.01, lam = 0.014)

    br = get_br(rtc_mod, params2, initial, 3.)
    bf0 = bf_point_df(br)
    df0 = create_br_df(br)
    kdam01 = findall(x->x==bf0.kdam[1],df0.kdam)[1]
    kdam02 = findall(x->x==bf0.kdam[2],df0.kdam)[1]

    bf, df, kdam1, kdam2 = different_levels_inhibition(rtc_inhib_mod, k_inhib1, k_inhib2, inhib)

    bfa, dfa, kdam1a, kdam2a = different_levels_inhibition(rtc_inhib_mod, k_inhib1a, k_inhib2, inhib)

    bfb, dfb, kdam1b, kdam2b = different_levels_inhibition(rtc_inhib_mod, k_inhib1b, k_inhib2, inhib)

    df=Float64.(df)
    df0=Float64.(df0)
    dfa=Float64.(dfa)
    dfb=Float64.(dfb)

    int_orig = QuadraticInterpolation(df0[!,specie][1:kdam01], df0.kdam[1:kdam01])
    inhib_int_rtcba = QuadraticInterpolation(dfa[!,specie][1:kdam1a], dfa.kdam[1:kdam1a])
    inhib_int_rtcbb = QuadraticInterpolation(dfb[!,specie][1:kdam1b], dfb.kdam[1:kdam1b])
    inhib_int_rtcb1 = QuadraticInterpolation(df[!,specie][1:kdam1], df.kdam[1:kdam1])

    orig_rtcb = [int_orig(i) for i in range(0,df0.kdam[kdam01],length=1000)]
    inhib_rtcba = [inhib_int_rtcba(i) for i in range(0,dfa.kdam[kdam1a],length=1000)]
    inhib_rtcbb = [inhib_int_rtcbb(i) for i in range(0,dfb.kdam[kdam1b],length=1000)]
    inhib_rtcb1 = [inhib_int_rtcb1(i) for i in range(0,df.kdam[kdam1],length=1000)]

    perc_deca = [(i/o)*100 for (i,o) in zip(inhib_rtcba, orig_rtcb)]
    perc_decb = [(i/o)*100 for (i,o) in zip(inhib_rtcbb, orig_rtcb)]
    perc_dec1 = [(i/o)*100 for (i,o) in zip(inhib_rtcb1, orig_rtcb)]

    return [perc_deca, perc_decb, perc_dec1]
end

rtcb_dec = protein_decrease(rtc_inhib_mod_rtcb, :rtcb)
# rtcb_dec_rh = protein_decrease(rtc_inhib_mod_rtcb, :rh)
rtca_dec = protein_decrease(rtc_inhib_mod_rtca, :rtca)
# rtca_dec_rh = protein_decrease(rtc_inhib_mod_rtca, :rh)
rtcr_dec = protein_decrease(rtc_inhib_mod_rtcr, :rtcr)
# rtcr_dec_rh = protein_decrease(rtc_inhib_mod_rtcr, :rh)

p_rtcb_dec = plot([scatter(x=range(0,1000,length=1000),y=rtcb_dec[1], name="kinhib = 0.3"),scatter(x=range(0,1000,length=1000),y=rtcb_dec[2], name="kinhib=0.5"),scatter(x=range(0,1000,length=1000),y=rtcb_dec[3],name="kinhib = 1")])
# p_rtcb_rh_dec = plot([scatter(x=range(0,1000,length=1000),y=rtcb_dec_rh[1], name="kinhib = 0.3"),scatter(x=range(0,1000,length=1000),y=rtcb_dec_rh[2], name="kinhib=0.5"),scatter(x=range(0,1000,length=1000),y=rtcb_dec_rh[3],name="kinhib = 1")])
p_rtca_dec = plot([scatter(x=range(0,1000,length=1000),y=rtca_dec[1], name="kinhib = 0.3"),scatter(x=range(0,1000,length=1000),y=rtca_dec[2], name="kinhib=0.5"),scatter(x=range(0,1000,length=1000),y=rtca_dec[3],name="kinhib = 1")])
# p_rtca_rh_dec = plot([scatter(x=range(0,1000,length=1000),y=rtca_dec_rh[1], name="kinhib = 0.3"),scatter(x=range(0,1000,length=1000),y=rtca_dec_rh[2], name="kinhib=0.5"),scatter(x=range(0,1000,length=1000),y=rtca_dec_rh[3],name="kinhib = 1")])
p_rtcr_dec = plot([scatter(x=range(0,1000,length=1000),y=rtcr_dec[1], name="kinhib = 0.3"),scatter(x=range(0,1000,length=1000),y=rtcr_dec[2], name="kinhib=0.5"),scatter(x=range(0,1000,length=1000),y=rtcr_dec[3],name="kinhib = 1")])
# p_rtcr_rh_dec = plot([scatter(x=range(0,1000,length=1000),y=rtcr_dec_rh[1], name="kinhib = 0.3"),scatter(x=range(0,1000,length=1000),y=rtcr_dec_rh[2], name="kinhib=0.5"),scatter(x=range(0,1000,length=1000),y=rtcr_dec_rh[3],name="kinhib = 1")])

[p_rtcb_dec p_rtcb_rh_dec; p_rtca_dec p_rtca_rh_dec; p_rtcr_dec p_rtcr_rh_dec]


av_dec_rtcb = [mean(rtcb_dec[i]) for i in range(1,3)]
# min_dec_rtcb_rh = [minimum(rtcb_dec_rh[i]) for i in range(1,3)]
# max_dec_rtcb_rh = [maximum(rtcb_dec_rh[i]) for i in range(1,3)]
# av_dec_rtcb_rh = [mean(rtcb_dec_rh[i]) for i in range(1,3)]

av_dec_rtca = [mean(rtca_dec[i]) for i in range(1,3)]
# min_dec_rtca_rh = [minimum(rtca_dec_rh[i]) for i in range(1,3)]
# max_dec_rtca_rh = [maximum(rtca_dec_rh[i]) for i in range(1,3)]
# av_dec_rtca_rh = [mean(rtca_dec_rh[i]) for i in range(1,3)]

av_dec_rtcr = [mean(rtcr_dec[i]) for i in range(1,3)]
# min_dec_rtcr_rh = [minimum(rtcr_dec_rh[i]) for i in range(1,3)]
# max_dec_rtcr_rh = [maximum(rtcr_dec_rh[i]) for i in range(1,3)]
# av_dec_rtcr_rh = [mean(rtcr_dec_rh[i]) for i in range(1,3)]


# res_rtcb = DataFrame("kinhib"=>[k_inhib1a, k_inhib1b, k_inhib1], "perc_change_bs_rtcb"=>percentage_decrease_rtcb, "perc_protein_change_rtcb"=>av_dec_rtcb, "perc_rh_protein_change rtcb"=>av_dec_rtcb_rh, "min_perc_rh_protein_change_rtcb"=>min_dec_rtcb_rh, "max_perc_rh_protein_change_rtcb"=>max_dec_rtcb_rh)

# res_rtca = DataFrame("kinhib"=>[k_inhib1a, k_inhib1b, k_inhib1], "perc_change_bs_rtca"=>percentage_decrease_rtca, "perc_protein_change_rtca"=>av_dec_rtca, "perc_rh_protein_change_rtca"=>av_dec_rtca_rh, "min_perc_rh_protein_change_rtca"=>min_dec_rtca_rh, "max_perc_rh_protein_change_rtca"=>max_dec_rtca_rh)

# res_rtcr = DataFrame("kinhib"=>[k_inhib1a, k_inhib1b, k_inhib1], "perc_change_bs_rtcr"=>percentage_decrease_rtcr, "perc_protein_change_rtcr"=>av_dec_rtcr, "perc_rh_protein_change_rtcr"=>av_dec_rtcr_rh, "min_perc_rh_protein_change_rtcr"=>min_dec_rtcr_rh, "max_perc_rh_protein_change_rtcr"=>max_dec_rtcr_rh)

# perc_dec_bs = DataFrame("kinhib"=>["kinhib=$k_inhib1a", "kinhib=$k_inhib1b", "kinhib=$k_inhib1"], "rtcb"=>-percentage_decrease_rtcb, "rtca"=>-percentage_decrease_rtca, "rtcr"=>-percentage_decrease_rtcr)
# perc_dec_protein = DataFrame("kinhib"=>["kinhib=$k_inhib1a", "kinhib=$k_inhib1b", "kinhib=$k_inhib1"], "rtcb"=>-av_dec_rtcb, "rtca"=>-av_dec_rtca, "rtcr"=>-av_dec_rtcr)
# perc_dec_rh = DataFrame("kinhib"=>["kinhib=$k_inhib1a", "kinhib=$k_inhib1b", "kinhib=$k_inhib1"], "rtcb"=>-av_dec_rtcb_rh, "rtca"=>-av_dec_rtca_rh, "rtcr"=>-av_dec_rtcr_rh)
# min_dec_rh = DataFrame("kinhib"=>["kinhib=$k_inhib1a", "kinhib=$k_inhib1b", "kinhib=$k_inhib1"], "rtcb"=>-min_dec_rtcb_rh, "rtca"=>-min_dec_rtca_rh, "rtcr"=>-min_dec_rtcr_rh)
# max_dec_rh = DataFrame("kinhib"=>["kinhib=$k_inhib1a", "kinhib=$k_inhib1b", "kinhib=$k_inhib1"], "rtcb"=>-max_dec_rtcb_rh, "rtca"=>-max_dec_rtca_rh, "rtcr"=>-max_dec_rtcr_rh)

# prot_dec_tot = DataFrame("kinhib"=>["kinhib=$k_inhib1a", "kinhib=$k_inhib1b", "kinhib=$k_inhib1"],"rtcb"=>-av_dec_rtcb, "rtca"=>-av_dec_rtca, "rtcr"=>-av_dec_rtcr, "rtcb_rh"=>-av_dec_rtcb_rh, "rtca_rh"=>-av_dec_rtca_rh, "rtcr_rh"=>-av_dec_rtcr_rh)

# protbs_dec_tot = DataFrame("kinhib"=>["kinhib=$k_inhib1a", "kinhib=$k_inhib1b", "kinhib=$k_inhib1"],"rtcb"=>-av_dec_rtcb, "rtca"=>-av_dec_rtca, "rtcr"=>-av_dec_rtcr, "rtcb_rh"=>-av_dec_rtcb_rh, "rtca_rh"=>-av_dec_rtca_rh, "rtcr_rh"=>-av_dec_rtcr_rh, "rtcb_bs"=>-percentage_decrease_rtcb, "rtca_bs"=>-percentage_decrease_rtca, "rtcr_bs"=>-percentage_decrease_rtcr)


# plot([bar(perc_dec_bs, x=:kinhib, y=y, name=String(y)) for y in [:rtcb, :rtca, :rtcr]], Layout(yaxis_title="% ↓ from original bs region"))
# plot([bar(perc_dec_protein, x=:kinhib, y=y, name=String(y)) for y in [:rtcb, :rtca, :rtcr]], Layout(yaxis_title="% ↓ from original protein"))
# plot([bar(perc_dec_rh, x=:kinhib, y=y, name=String(y)) for y in [:rtcb, :rtca, :rtcr]], Layout(yaxis_title="% ↓ from original rh"))
# plot([bar(min_dec_rh, x=:kinhib, y=y, name=String(y)) for y in [:rtcb, :rtca, :rtcr]], Layout(yaxis_title="% ↓ from original"))
# plot([bar(max_dec_rh, x=:kinhib, y=y, name=String(y)) for y in [:rtcb, :rtca, :rtcr]], Layout(yaxis_title="% ↓ from original"))

# plot([bar(prot_dec_tot, x=:kinhib, y=y, name=String(y)) for y in propertynames(prot_dec_tot)[2:end]], Layout(yaxis_title="% ↓ from original", yaxis_type="log"))

# plot([bar(protbs_dec_tot, x=:kinhib, y=y, name=String(y)) for y in propertynames(protbs_dec_tot)[2:end]], Layout(yaxis_title="% ↓ from original", yaxis_type="log"))


p_rtcb = plot([rtcb_traces[1], rtcb_traces[5], rtcb_traces[9], rtcb_traces[13]],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="RtcB (μM)",
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white", font=attr(size=20, color="black", family="sans-serif")))

# function trapezoidal_rule(x, y)
#     n = length(x)
#     area = 0.0

#     for i in 1:n-1
#         dx = x[i+1] - x[i]
#         area += (y[i] + y[i+1]) * dx / 2
#     end

#     return area
# end

br = get_br(rtc_mod, params2, initial, 3.)
bf0 = bf_point_df(br)
df0 = create_br_df(br)
kdam01 = findall(x->x==bf0.kdam[1],df0.kdam)[1]
kdam02 = findall(x->x==bf0.kdam[2],df0.kdam)[1]

bf, df, kdam1, kdam2 = different_levels_inhibition(rtc_inhib_mod_rtcb, k_inhib1, k_inhib2, inhib)

bfa, dfa, kdam1a, kdam2a = different_levels_inhibition(rtc_inhib_mod_rtcb, k_inhib1a, k_inhib2, inhib)

bfb, dfb, kdam1b, kdam2b = different_levels_inhibition(rtc_inhib_mod_rtcb, k_inhib1b, k_inhib2, inhib)

# result0 = trapezoidal_rule(df0.kdam, df0.rh)
# result = trapezoidal_rule(df.kdam, df.rh)
# resulta = trapezoidal_rule(dfa.kdam, dfa.rh)
# resultb = trapezoidal_rule(dfb.kdam, dfb.rh)

using QuadGK, Interpolations

function area_under_curve_rh(df,kdam1)
    df=Float64.(df)
    x = df.kdam[1:kdam1]
    y = df.rh[1:kdam1]
    int_orig = Interpolations.LinearInterpolation(x, y)
    f(x) = int_orig(x)
    a = minimum(x)
    b = maximum(x)

    result, error = quadgk(f, a, b)
    return @LArray [x,y,f,result] (:x,:y,:f1,:result)
end

function all_area_under_curve_rh(rtc_inhib_mod)
    br = get_br(rtc_mod, params2, initial, 3.)
    bf0 = bf_point_df(br)
    df0 = create_br_df(br)
    kdam01 = findall(x->x==bf0.kdam[1],df0.kdam)[1]
    kdam02 = findall(x->x==bf0.kdam[2],df0.kdam)[1]

    bf, df, kdam1, kdam2 = different_levels_inhibition(rtc_inhib_mod, k_inhib1, k_inhib2, inhib)

    bfa, dfa, kdam1a, kdam2a = different_levels_inhibition(rtc_inhib_mod, k_inhib1a, k_inhib2, inhib)

    bfb, dfb, kdam1b, kdam2b = different_levels_inhibition(rtc_inhib_mod, k_inhib1b, k_inhib2, inhib)

    res0 = area_under_curve_rh(df0,kdam01)
    res = area_under_curve_rh(df,kdam1)
    res1 = area_under_curve_rh(dfa,kdam1a)
    res2 = area_under_curve_rh(dfb,kdam1b)

    return [res0,res1,res2,res]
end

rtcb_auc = all_area_under_curve_rh(rtc_inhib_mod_rtcb)
rtcr_auc = all_area_under_curve_rh(rtc_inhib_mod_rtcr)
rtca_auc = all_area_under_curve_rh(rtc_inhib_mod_rtca)

rtcb_fill1 = scatter(x=rtcb_auc[1].x,y=rtcb_auc[1].f1.(rtcb_auc[1].x), fill="tozeroy", showlegend=false, mode="none")
rtcb_fill2 = scatter(x=rtcb_auc[2].x,y=rtcb_auc[2].f1.(rtcb_auc[2].x), fill="tozeroy", showlegend=false, mode="none")
rtcb_fill3 = scatter(x=rtcb_auc[3].x,y=rtcb_auc[3].f1.(rtcb_auc[3].x), fill="tozeroy", showlegend=false, mode="none")
rtcb_fill4 = scatter(x=rtcb_auc[4].x,y=rtcb_auc[4].f1.(rtcb_auc[4].x), fill="tozeroy", showlegend=false, mode="none")

rtca_fill1 = scatter(x=rtca_auc[1].x,y=rtca_auc[1].f1.(rtca_auc[1].x), fill="tozeroy", showlegend=false, mode="none")
rtca_fill2 = scatter(x=rtca_auc[2].x,y=rtca_auc[2].f1.(rtca_auc[2].x), fill="tozeroy", showlegend=false, mode="none")
rtca_fill3 = scatter(x=rtca_auc[3].x,y=rtca_auc[3].f1.(rtca_auc[3].x), fill="tozeroy", showlegend=false, mode="none")
rtca_fill4 = scatter(x=rtca_auc[4].x,y=rtca_auc[4].f1.(rtca_auc[4].x), fill="tozeroy", showlegend=false, mode="none")

rtcr_fill1 = scatter(x=rtcr_auc[1].x,y=rtcr_auc[1].f1.(rtcr_auc[1].x), fill="tozeroy", showlegend=false, mode="none")
rtcr_fill2 = scatter(x=rtcr_auc[2].x,y=rtcr_auc[2].f1.(rtcr_auc[2].x), fill="tozeroy", showlegend=false, mode="none")
rtcr_fill3 = scatter(x=rtcr_auc[3].x,y=rtcr_auc[3].f1.(rtcr_auc[3].x), fill="tozeroy", showlegend=false, mode="none")
rtcr_fill4 = scatter(x=rtcr_auc[4].x,y=rtcr_auc[4].f1.(rtcr_auc[4].x), fill="tozeroy", showlegend=false, mode="none")



rtcb_rh_traces = creating_rtc_inhib_plot(rtc_inhib_mod_rtcb, :rh)
push!(rtcb_rh_traces, rtcb_fill1, rtcb_fill2, rtcb_fill3, rtcb_fill4)
p_rtcb_rh = plot([i for i in rtcb_rh_traces],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="Rh (μM)",
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=20, color="black", family="sans-serif")))

rtca_rh_traces = creating_rtc_inhib_plot(rtc_inhib_mod_rtca, :rh)
push!(rtca_rh_traces, rtca_fill1, rtca_fill2, rtca_fill3, rtca_fill4)
p_rtca_rh = plot([i for i in rtca_rh_traces],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="Rh (μM)",
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=20, color="black", family="sans-serif")))

rtcr_rh_traces = creating_rtc_inhib_plot(rtc_inhib_mod_rtcr, :rh)
push!(rtcr_rh_traces, rtcr_fill1, rtcr_fill2, rtcr_fill3, rtcr_fill4)
p_rtcr_rh = plot([i for i in rtcr_rh_traces],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="Rh (μM)",
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=20, color="black", family="sans-serif")))


savefig(p_rtcb_rh, "/home/holliehindley/phd/may23_rtc/paper_plots/rtcb_rh_plusfill.svg")
savefig(p_rtca_rh, "/home/holliehindley/phd/may23_rtc/paper_plots/rtca_rh_plusfill.svg")
savefig(p_rtcr_rh, "/home/holliehindley/phd/may23_rtc/paper_plots/rtcr_rh_plusfill.svg")





areas_rtcb = [res0,res1,res2,res]
percs_rtcb = [100*(areas_rtcb[i]/areas_rtcb[1]) for i in range(2,4)]

areas_rtcb = [rtcb_auc[i][:result] for i in range(1,4)]
areas_rtcr = [rtcr_auc[i][:result] for i in range(1,4)]
areas_rtca = [rtca_auc[i][:result] for i in range(1,4)]

percs_rtcb = [100*(areas_rtcb[i]/areas_rtcb[1]) for i in range(2,4)]
percs_rtcr = [100*(areas_rtcr[i]/areas_rtcr[1]) for i in range(2,4)]
percs_rtca = [100*(areas_rtca[i]/areas_rtca[1]) for i in range(2,4)]



# new_df = DataFrame("data"=>["bs_region","rtc_conc","rh_conc"], "rtcb"=>[-percentage_decrease_rtcb[3],-av_dec_rtcb[3],-av_dec_rtcb_rh[3]], "rtca"=>[-percentage_decrease_rtca[3],-av_dec_rtca[3],-av_dec_rtca_rh[3]],"rtcr"=>[-percentage_decrease_rtcr[3],-av_dec_rtcr[3],-av_dec_rtcr_rh[3]])
# rdf = DataFrame("data"=>["bs_region","rtc_conc","rh_conc"], "rtcb"=>[round(-percentage_decrease_rtcb[3],digits=2),round(-av_dec_rtcb[3],digits=2),round(-av_dec_rtcb_rh[3],digits=2)], "rtca"=>[round(-percentage_decrease_rtca[3],digits=2),round(-av_dec_rtca[3],digits=2),round(-av_dec_rtca_rh[3],digits=2)],"rtcr"=>[round(-percentage_decrease_rtcr[3],digits=2),round(-av_dec_rtcr[3],digits=2),round(-av_dec_rtcr_rh[3],digits=2)])

rdf = DataFrame("data"=>["Rtc conc.","Bistability region","Growth Capacity"], 
"rtcb"=>[round(av_dec_rtcb[3],digits=2),round(percentage_size_rtcb[3],digits=2),round(percs_rtcb[3],digits=2)],
"rtcr"=>[round(av_dec_rtcr[3],digits=2),round(percentage_size_rtcr[3],digits=2),round(percs_rtcr[3],digits=2)],
"rtca"=>[round(av_dec_rtca[3],digits=2),round(percentage_size_rtca[3],digits=2),round(percs_rtca[3],digits=2)])


p = plot([bar(rdf, x=:data, y=:rtcb, text=:rtcb, textposition="auto", name=String(:rtcb), marker_color=["#7e5c94ff","#7e5c94ff","#7e5c94ff"]),
bar(rdf, x=:data, y=:rtcr, text=:rtcr, textposition="auto", name=String(:rtcr), marker_color=["#28726dff","#28726dff","#28726dff"]),
bar(rdf, x=:data, y=:rtca, text=:rtca, textposition="auto", name=String(:rtca), marker_color=["#a1403fff","#a1403fff","#a1403fff"])], 
Layout(yaxis_title="% of original",yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=20, color="black", family="sans-serif")))


savefig(p,"/home/holliehindley/phd/may23_rtc/paper_plots/bar.svg")















atp_range = range(1000,stop=5000,length=3)
lam_range = range(0.01,stop=0.04,length=3)

copyparams = deepcopy(params2)

bfs=[]; dfs=[];
for i in ProgressBar(atp_range)
    copyparams = deepcopy(params_for_ssval_setup)
    params = merge(copyparams, (:atp=>i,))
    br = get_br(rtc_mod, params, initial, 3.)
    bf = bf_point_df(br)
    df = create_br_df(br)
    push!(bfs, bf)
    push!(dfs, df)
end



rtcb1_atp1, rtcb2_atp1, rtcb3_atp1, bf_rtcb_atp1 = plot_rtc_bf(bfs[1], dfs[1], findall(x->x==bfs[1].kdam[1],dfs[1].kdam)[1], findall(x->x==bfs[1].kdam[2],dfs[1].kdam)[1], :rtcb)
rtcb1_atp2, rtcb2_atp2, rtcb3_atp2, bf_rtcb_atp2 = plot_rtc_bf(bfs[2], dfs[2], findall(x->x==bfs[2].kdam[1],dfs[2].kdam)[1], findall(x->x==bfs[2].kdam[2],dfs[2].kdam)[1], :rtcb)
rtcb1_atp3, rtcb2_atp3, rtcb3_atp3, bf_rtcb_atp3 = plot_rtc_bf(bfs[3], dfs[3], findall(x->x==bfs[3].kdam[1],dfs[3].kdam)[1], findall(x->x==bfs[3].kdam[2],dfs[3].kdam)[1], :rtcb)

rtcbp = plot([rtcb1_atp3, rtcb2_atp3, rtcb3_atp3, bf_rtcb_atp3 , rtcb1_atp1, rtcb2_atp1, rtcb3_atp1, bf_rtcb_atp1, rtcb1_atp2, rtcb2_atp2, rtcb3_atp2, bf_rtcb_atp2])

savefig(rtcbp, "/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/plots/rtcb_atp.svg")


lam_range = range(0.014,stop=0.04,length=3)

lam_range = [0.01, 0.014, 0.02]
bfs=[]; dfs=[];
for i in (lam_range)
    copyparams = deepcopy(params_for_ssval_setup)
    params = merge(copyparams, (:lam=>i,))
    br = get_br(rtc_mod, params, initial, 3.)
    bf = bf_point_df(br)
    df = create_br_df(br)
    push!(bfs, bf)
    push!(dfs, df)
end

rtcb1_2, rtcb2_2, rtcb3_2, rtcr1_2, rtcr2_2, rtcr3_2, rh1_2, rh2_2, rh3_2, rt1_2, rt2_2, rt3_2, rd1_2, rd2_2, rd3_2, rtca1_2, rtca2_2, rtca3_2, rma1_2, rma2_2, rma3_2, rmr1_2, rmr2_2, rmr3_2, rmb1_2, rmb2_2, rmb3_2, bf_rtcb_2, bf_rtca_2, bf_rtcr_2, bf_rma_2, bf_rmb_2, bf_rmr_2, bf_rh_2, bf_rd_2, bf_rt_2 = plot_all_bf(bfs[2],dfs[2], :mediumorchid, lam_range[2], "λ")
rtcb1_3, rtcb2_3, rtcb3_3, rtcr1_3, rtcr2_3, rtcr3_3, rh1_3, rh2_3, rh3_3, rt1_3, rt2_3, rt3_3, rd1_3, rd2_3, rd3_3, rtca1_3, rtca2_3, rtca3_3, rma1_3, rma2_3, rma3_3, rmr1_3, rmr2_3, rmr3_3, rmb1_3, rmb2_3, rmb3_3, bf_rtcb_3, bf_rtca_3, bf_rtcr_3, bf_rma_3, bf_rmb_3, bf_rmr_3, bf_rh_3, bf_rd_3, bf_rt_3 = plot_all_bf(bfs[3],dfs[3], :purple, lam_range[3], "λ")

a = (scatter(x=dfs[1].kdam, y=dfs[1].rtcb, line=attr(width=3, color=:plum), name="λ = $(lam_range[1])"))

rtcb_lam = plot([a, rtcb1_2, rtcb2_2, rtcb3_2, bf_rtcb_2, rtcb1_3, rtcb2_3, rtcb3_3, bf_rtcb_3],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", yaxis_title="RtcB (μM)",
yaxis=attr(showline=true,linewidth=1,linecolor="black",mirror=true),xaxis=attr(showline=true,linewidth=1,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,plot_bgcolor="white",legend=attr(x=0.75,y=1)))

savefig(rtcb_lam, "/home/hollie_hindley/Documents/may23_rtc/analysis/bifurcation_analysis/plots/rtcb_lam.svg")
