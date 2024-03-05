using Parameters, CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, Statistics
using Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra, Printf, ProgressBars, LabelledArrays, DataFrames, PlotlyJS, ModelingToolkit


include("/home/holliehindley/phd/general_funcs/solving.jl")
include("/home/holliehindley/phd/rtc_model/models/rtc_orig.jl")
include("/home/holliehindley/phd/rtc_model/models/rtc_trna_model.jl")
include("/home/holliehindley/phd/rtc_model/functions/bf_funcs/bf_funcs.jl")

br = get_br(rtc_trna_model, ssvals_trna, params_trna, 20.)
bf = bf_point_df(br)
df = create_br_df(br)
kdam1 = findall(x->x==bf.kdam[1],df.kdam)[1]
kdam2 = findall(x->x==bf.kdam[2],df.kdam)[1]

rtcb1, rtcb2, rtcb3 = plot_rtc_bf(df, kdam1, kdam2, :rtcb, "1", "b693ccff", "RtcB")
rtcr1, rtcr2, rtcr3 = plot_rtc_bf(df, kdam1, kdam2, :rtcr, "2", "4ca7a2ff", "RtcR")
rtca1, rtca2, rtca3 = plot_rtc_bf(df, kdam1, kdam2, :rtca, "3", "e48080ff", "RtcA")

p = plot([rtcb1, rtcb2, rtcb3, rtcr1, rtcr2, rtcr3, rtca1, rtca2, rtca3],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="Rtc protein (μM)",
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))

rh1, rh2, rh3 = plot_rtc_bf(df, kdam1, kdam2, :rh, "1", "ffd30cff", "RtcB")
rt1, rt2, rt3 = plot_rtc_bf(df, kdam1, kdam2, :rt, "2", "e96100ff", "rt")
rd1, rd2, rd3 = plot_rtc_bf(df, kdam1, kdam2, :rd, "3", "ac0606ff", "rd")

p1 = plot([rh1, rh2, rh3, rt1, rt2, rt3, rd1, rd2, rd3],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="Ribosomes (μM)", 
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),#yaxis_type="log",
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))

[p1_sig p1]

savefig(p, "/home/holliehindley/phd/rtc_model/paper_plots/tRNA/rtc_proteins.svg")
savefig(p1, "/home/holliehindley/phd/rtc_model/paper_plots/tRNA/ribosomes.svg")



kdam_range = range(0,1000,length=100)
kdam_range2 = range(1000,0,length=100)


res_trna1 = numerical_bistability_analysis(rtc_trna_model, params_trna, ssvals_trna, :trna, trna_species, kdam_range)
res_trna2 = numerical_bistability_analysis(rtc_trna_model, params_trna, ssvals_trna, :trna, trna_species, kdam_range2)
plot([scatter(x=kdam_range, y=res_trna1), scatter(x=kdam_range2, y=res_trna2)])

res_rtcb1 = numerical_bistability_analysis(rtc_trna_model, params_trna, ssvals_trna, :rtcb, trna_species, kdam_range)
res_rtcb2 = numerical_bistability_analysis(rtc_trna_model, params_trna, ssvals_trna, :rtcb, trna_species, kdam_range2)
plot([scatter(x=kdam_range, y=res_rtcb1), scatter(x=kdam_range2, y=res_rtcb2)])


res_rd1 = numerical_bistability_analysis(rtc_model_trna, params_trna, ssvals, :rd, trna_species, kdam_range)
res_rd2 = numerical_bistability_analysis(rtc_model_trna, params_trna, ssvals, :rd, trna_species, kdam_range2)
res_rt1 = numerical_bistability_analysis(rtc_model_trna, params_trna, ssvals, :rt, trna_species, kdam_range)
res_rt2 = numerical_bistability_analysis(rtc_model_trna, params_trna, ssvals, :rt, trna_species, kdam_range2)

ptrna1 = scatter(x=kdam_range, y=res_trna1, name="tRNA_h ON", line=attr(color="ffd30cff",width=6.5))
ptrna2 = scatter(x=kdam_range2, y=res_trna2, name="tRNA_h OFF", line=attr(color="ffd30cff",width=6.5))
prd1 = scatter(x=kdam_range, y=res_rd1, name="tRNA_d ON", line=attr(color="ac0606ff",width=6.5))
prd2 = scatter(x=kdam_range2, y=res_rd2, name="tRNA_d OFF", line=attr(color="ac0606ff",width=6.5))
prt1 = scatter(x=kdam_range, y=res_rt1, name="tRNA_t ON", line=attr(color="e96100ff",width=6.5))
prt2 = scatter(x=kdam_range2, y=res_rt2, name="tRNA_t OFF", line=attr(color="e96100ff",width=6.5))

p_trnas = plot([ptrna1, ptrna2, prd1, prd2, prt1, prt2], 
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", title = "Numerical bistability analysis",
yaxis_title="tRNAs (μM)", 
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))

res_rtcb1 = numerical_bistability_analysis(rtc_model_trna, params_trna, ssvals, :rtcb, trna_species, kdam_range)
res_rtcb2 = numerical_bistability_analysis(rtc_model_trna, params_trna, ssvals, :rtcb, trna_species, kdam_range2)
res_rtca1 = numerical_bistability_analysis(rtc_model_trna, params_trna, ssvals, :rtca, trna_species, kdam_range)
res_rtca2 = numerical_bistability_analysis(rtc_model_trna, params_trna, ssvals, :rtca, trna_species, kdam_range2)
res_rtcr1 = numerical_bistability_analysis(rtc_model_trna, params_trna, ssvals, :rtcr, trna_species, kdam_range)
res_rtcr2 = numerical_bistability_analysis(rtc_model_trna, params_trna, ssvals, :rtcr, trna_species, kdam_range2)

prtcb1 = scatter(x=kdam_range, y=res_rtcb1, name="RtcB ON", line=attr(color="b693ccff",width=6.5))
prtcb2 = scatter(x=kdam_range2, y=res_rtcb2, name="RtcB OFF", line=attr(color="b693ccff",width=6.5))
prtcr1 = scatter(x=kdam_range, y=res_rtcr1, name="RtcR ON", line=attr(color="4ca7a2ff",width=6.5))
prtcr2 = scatter(x=kdam_range2, y=res_rtcr2, name="RtcR OFF", line=attr(color="4ca7a2ff",width=6.5))
prtca1 = scatter(x=kdam_range, y=res_rtca1, name="RtcA ON", line=attr(color="e48080ff",width=6.5))
prtca2 = scatter(x=kdam_range2, y=res_rtca2, name="RtcA OFF", line=attr(color="e48080ff",width=6.5))

p_rtcs = plot([prtcb1, prtcb2, prtcr1, prtcr2, prtca1, prtca2], 
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", title = "Numerical bistability analysis",
yaxis_title="Rtc proteins (μM)", 
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))









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
kdeg_trna = 0.00001
trna_species = [:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :trna, :rd, :rt]
init_trna = [0,0,0,0,0,0,135.5,0,0] # tRNA initial conc = 135.5
params_trna = @LArray [10., c, kr*12, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 0.5, ktag, kdeg, kin_trna, 3578.9473684210525, na, nb, nr, 0.014, rh, thr_t] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam, :rh, :thr_t)

# params_trna = @LArray [10., c, kr*12, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 0.5, ktag, kdeg_trna, kin_trna, 3578.9473684210525, na, nb, nr, 0.014, rh, thr_t] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam, :rh, :thr_t)

kdam_range = range(0,100,length=50)
kdam_range2 = range(100,0,length=50)

params_trna2 = (L = 10., c = 0.001, kr = 0.125*12, Vmax_init = 39.51, Km_init = 250.,
θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923,
kdeg = 0.001, kin = kin_trna, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
kdam =  0.01, lam = 0.014, rh = rh, thr_t = thr_t)

# params_trna2 = (L = 10., c = 0.001, kr = 0.125*12, Vmax_init = 39.51, Km_init = 250.,
# θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
# krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923,
# kdeg = kdeg_trna, kin = kin_trna, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
# kdam =  0.01, lam = 0.014, rh = rh, thr_t = thr_t)

br = get_br(rtc_mod_trna, params_trna2, init_trna, 400.) # dsmax = 0.05, 0.08, 0.09 for to finish bifurcation at 0 
br1 = get_br(rtc_mod_trna, params_trna2, init_trna, 400.) # dsmax > 0.1 to finish bifurcation at kdam = 400 (see negative values)  
# Plots.plot(br, vars=(:param, :trna), linewidthstable=5)

# df = create_br_df(br)
# df1 = create_br_df(br1)
# bf = bf_point_df(br)
# bf1 = bf_point_df(br1)
# kdam1 = findall(x->x==bf.kdam[1],df.kdam)[1]
# kdam2 = findall(x->x==bf1.kdam[1],df1.kdam)[1]
# stable_trna = df.rh[1:kdam1]
# unstable_trna = df.rh[kdam1:end]
# stable_trna1 = df1.rh[1:kdam2]
# unstable_trna1 = df1.rh[kdam2:end]
# stable1 = scatter(x=df.kdam[1:kdam1],y=stable_trna, line=attr(color="#1f77b4"), name="BifurcationKit to 0", legendgroup=5)
# unstable1 = scatter(x=df.kdam[kdam1:end],y=unstable_trna, line=attr(color="#1f77b4", dash="dash"), name="BifurcationKit to 0", legendgroup=5)
# stable2 = scatter(x=df1.kdam[1:kdam2],y=stable_trna1, line=attr(color="#ff7f0e"), name="BifurcationKit to 400", legendgroup=6)
# unstable2 = scatter(x=df1.kdam[kdam2:end],y=unstable_trna1, line=attr(color="#ff7f0e", dash="dash"), name="BifurcationKit to 400", legendgroup=6)
# bf_p = (scatter(x=[df.kdam[kdam1]], y=[df.rh[kdam1]], legendgroup=5, line=attr(color="#1f77b4"), mode="markers", showlegend=false))
# bf_p1 = (scatter(x=[df1.kdam[kdam2]], y=[df1.rh[kdam2]], legendgroup=6, line=attr(color="#ff7f0e"), mode="markers", showlegend=false))
# trna_h_p1 = (scatter(x=df.kdam, y=df.rh, name="BifurcationKit", legendgroup=2))
# trna_h_p2 = (scatter(x=df1.kdam, y=df1.rh, name="BifurcationKit", legendgroup=7))
# plot([stable1, unstable1], Layout(title="θ tRNA = 30"))

# plot([stable2, unstable2])


function trna_plotting(specie, color, br, i)
    # res_trna1 = numerical_bistability_analysis(rtc_model_trna, params_trna, init_trna, :trna, trna_species, kdam_range)
    res_trna2 = numerical_bistability_analysis(rtc_model_trna, params_trna, init_trna, specie, all_species, kdam_range2)
    # ptrna1 = scatter(x=kdam_range, y=res_trna1, name="Healthy tRNA", legendgroup=3, yaxis="y2", line=attr(color=:gold,linewidth=3))
    ptrna2 = scatter(x=kdam_range2, y=res_trna2, name="", legendgroup=i, showlegend=false, line=attr(color=color,width=5))
# res_rtcb1 = numerical_bistability_analysis(rtc_model_trna, params_trna, init_trna, :rtcb, trna_species, kdam_range)
# res_rtcb2 = numerical_bistability_analysis(rtc_model_trna, params_trna, init_trna, :rtcb, trna_species, kdam_range2)
# prtcb1 = scatter(x=kdam_range, y=res_rtcb1, name="RtcB", legendgroup=3, line=attr(color=:plum,linewidth=3))
# prtcb2 = scatter(x=kdam_range2, y=res_rtcb2, name="", showlegend=false, legendgroup=3, line=attr(color=:plum,linewidth=3))
    df = create_br_df(br)
    # df1 = create_br_df(br1)
    bf = bf_point_df(br)
    # bf1 = bf_point_df(br1)
    kdam1 = findall(x->x==bf.kdam[1],df.kdam)[1]
    # kdam2 = findall(x->x==bf1.kdam[1],df1.kdam)[1]
    stable_trna = df[!,specie][1:kdam1]
    unstable_trna = df[!,specie][kdam1:end]
    stable1 = scatter(x=df.kdam[1:kdam1],y=stable_trna, line=attr(color=color, width=5), name="$specie", legendgroup=i)
    unstable1 = scatter(x=df.kdam[kdam1:end],y=unstable_trna, line=attr(color=color, width=5, dash="dash"), name="", showlegend=false, legendgroup=i)
    return ptrna2, stable1, unstable1
end

ptrna, stable_trna, unstable_trna = trna_plotting(:rh, "#ffd30cff", br, 1)
prd, stable_rd, unstable_rd = trna_plotting(:rd, "#ac0606ff", br, 2)
prt, stable_rt, unstable_rt = trna_plotting(:rt, "#e96100ff", br, 3)
p1 = plot([ptrna, stable_trna, unstable_trna, prt, stable_rt, unstable_rt, prd, stable_rd, unstable_rd],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="tRNAs (μM)",
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=22, color="black", family="sans-serif")))#, yaxis_type="log"))

plot([prd, stable_rd, unstable_rd])

prtcb, stable_rtcb, unstable_rtcb = trna_plotting(:rtcb, "#b693ccff", br, 1)
prtcr, stable_rtcr, unstable_rtcr = trna_plotting(:rtcr, "#4ca7a2ff", br, 2)
prtca, stable_rtca, unstable_rtca = trna_plotting(:rtca, "#e48080ff", br, 3)

p2 = plot([prtcb, stable_rtcb, unstable_rtcb, prtcr, stable_rtcr, unstable_rtcr, prtca, stable_rtca, unstable_rtca],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="Rtc protein (μM)",
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=22, color="black", family="sans-serif")))


# res_trna1 = numerical_bistability_analysis(rtc_model_trna, params_trna, init_trna, :rd, all_species, kdam_range)
# res_trna2 = numerical_bistability_analysis(rtc_model_trna, params_trna, init_trna, :rd, all_species, kdam_range2)
# ptrna1 = scatter(x=kdam_range, y=res_trna1, name="rd", legendgroup=3, line=attr(color=:red,width=5))
# ptrna2 = scatter(x=kdam_range2, y=res_trna2, name="", legendgroup=3, showlegend=false, line=attr(color=:red,width=5))
# plot([ptrna1, ptrna2])


savefig(p1, "/home/holliehindley/phd/may23_rtc/paper_plots/tRNA/trna_trnas.svg")
savefig(p2, "/home/holliehindley/phd/may23_rtc/paper_plots/tRNA/trna_proteins.svg")








res_trna1 = numerical_bistability_analysis(rtc_model_trna, params_trna, init_trna, :trna, trna_species, kdam_range)
res_trna2 = numerical_bistability_analysis(rtc_model_trna, params_trna, init_trna, :trna, trna_species, kdam_range2)
ptrna1 = scatter(x=kdam_range, y=res_trna1, name="Healthy tRNA", legendgroup=3, yaxis="y2", line=attr(color=:gold,linewidth=3))
ptrna2 = scatter(x=kdam_range2, y=res_trna2, name="", legendgroup=3, showlegend=false, yaxis="y2", line=attr(color=:gold,linewidth=3))
res_rtcb1 = numerical_bistability_analysis(rtc_model_trna, params_trna, init_trna, :rtcb, trna_species, kdam_range)
res_rtcb2 = numerical_bistability_analysis(rtc_model_trna, params_trna, init_trna, :rtcb, trna_species, kdam_range2)
prtcb1 = scatter(x=kdam_range, y=res_rtcb1, name="RtcB", legendgroup=3, line=attr(color=:plum,linewidth=3))
prtcb2 = scatter(x=kdam_range2, y=res_rtcb2, name="", showlegend=false, legendgroup=3, line=attr(color=:plum,linewidth=3))
res_rd1 = numerical_bistability_analysis(rtc_model_trna, params_trna, init_trna, :rd, trna_species, kdam_range)
res_rd2 = numerical_bistability_analysis(rtc_model_trna, params_trna, init_trna, :rd, trna_species, kdam_range2)
prd1 = scatter(x=kdam_range, y=res_rd1, name="Rd", legendgroup=3, line=attr(color=:plum,linewidth=3))
prd2 = scatter(x=kdam_range2, y=res_rd2, name="", showlegend=false, legendgroup=3, line=attr(color=:plum,linewidth=3))

df = create_br_df(br)
bf = bf_point_df(br)
kdam1 = findall(x->x==bf.kdam[1],df.kdam)[1]

stable_trna = df.rh[1:kdam1]
unstable_trna = df.rh[kdam1:end]
stable1 = scatter(x=df.kdam[1:kdam1],y=stable_trna, line=attr(color=:gold, linewidth=3), name="Healthy tRNA", legendgroup=5, yaxis="y2")
unstable1 = scatter(x=df.kdam[kdam1:end],y=unstable_trna, line=attr(color=:gold, linewidth=3, dash="dash"), name="", showlegend=false, legendgroup=5, yaxis="y2")

stable_rtcb = df.rtcb[1:kdam1]
unstable_rtcb = df.rtcb[kdam1:end]
stable2 = scatter(x=df.kdam[1:kdam1],y=stable_rtcb, line=attr(color=:plum, linewidth=3), name="RtcB", legendgroup=5)
unstable2 = scatter(x=df.kdam[kdam1:end],y=unstable_rtcb, line=attr(color=:plum, linewidth=3, dash="dash"), name="", showlegend=false, legendgroup=5)

bf_p = (scatter(x=[df.kdam[kdam1]], y=[df.rh[kdam1]], legendgroup=5, line=attr(color=:darkblue), mode="markers", showlegend=false, yaxis="y2"))
bf_p1 = (scatter(x=[df.kdam[kdam1]], y=[df.rtcb[kdam1]], legendgroup=5, line=attr(color=:darkblue), mode="markers", showlegend=false))

trna = plot([ptrna2, prtcb2, stable1, unstable1, stable2, unstable2, bf_p, bf_p1],
Layout(legend=attr(x=0.75,y=1),width=1000,height=750,yaxis2=attr(overlaying="y",side="right"), xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="RtcB (μM)", yaxis2_title="Healthy tRNA (μM)",
yaxis=attr(showline=true,linewidth=1,linecolor="black",mirror=true),xaxis=attr(showline=true,linewidth=1,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white"))

plot([ptrna1, ptrna2])
plot([prtcb1, prtcb2])
plot([prd1,prd2])

trna = plot([ptrna2, prtcb2],
Layout(legend=attr(x=0.75,y=1),width=1000,height=750,yaxis2=attr(overlaying="y",side="right"), xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="RtcB (μM)", yaxis2_title="Healthy tRNA (μM)",
yaxis=attr(showline=true,linewidth=1,linecolor="black",mirror=true),xaxis=attr(showline=true,linewidth=1,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white"))
