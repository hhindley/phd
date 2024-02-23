using PlotlyJS, DataFrames, DifferentialEquations, Parameters, LabelledArrays, StaticArrays, Colors

include("/home/holliehindley/phd/rtc_model/models/rtc_trna_model.jl")
include("/home/holliehindley/phd/rtc_model/parameters/trna_params.jl")
include("/home/holliehindley/phd/general_funcs/solving.jl")


colours =["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", :blue]

solu = sol(rtc_trna_model, init_trna, tspan, params_trna)
df = create_solu_df(solu, trna_species)
plot([scatter(x=df.time, y=col, name="$(names(df)[i])", legendgroup="$i", marker_color=colours[i]) for (col, i) in zip(eachcol(df[:,2:end]), range(2,length(names(df))))], Layout(xaxis_type="log", yaxis_tickformat=".2e", title="kdam = $(params_trna.kdam)"))


kdam_range = range(0,1, length=100)
res=[]
ps = deepcopy(params_trna)
for i in kdam_range
    ps[kdam] = i
    solu = sol(rtc_trna_model, ssvals_trna, tspan, ps)
    push!(res, get_all_ssvals(solu, trna_species))
end


df_ssvals = DataFrame(vcat(transpose(res)...), :auto)
rename!(df_ssvals, trna_species)

plot([scatter(x=kdam_range, y=col, name="$(names(df_ssvals)[i])") for (col,i) in zip(eachcol(df_ssvals), range(1,length(names(df_ssvals))))])









kdam_range = range(0,1,length=11)
kdam_range2 = range(1,0,length=11)


res_trna1 = numerical_bistability_analysis(rtc_trna_model, params_trna, ssvals_trna, :trna, trna_species, kdam_range)
res_trna2 = numerical_bistability_analysis(rtc_trna_model, params_trna, ssvals_trna, :trna, trna_species, kdam_range2)
res_rd1 = numerical_bistability_analysis(rtc_trna_model, params_trna, ssvals_trna, :rd, trna_species, kdam_range)
res_rd2 = numerical_bistability_analysis(rtc_trna_model, params_trna, ssvals_trna, :rd, trna_species, kdam_range2)
res_rt1 = numerical_bistability_analysis(rtc_trna_model, params_trna, ssvals_trna, :rt, trna_species, kdam_range)
res_rt2 = numerical_bistability_analysis(rtc_trna_model, params_trna, ssvals_trna, :rt, trna_species, kdam_range2)

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


res_rtcb1_m = numerical_bistability_analysis(rtc_model_trna, params_trna, ssvals, :rm_b, trna_species, kdam_range)
res_rtcb2_m = numerical_bistability_analysis(rtc_model_trna, params_trna, ssvals, :rm_b, trna_species, kdam_range2)
res_rtca1_m = numerical_bistability_analysis(rtc_model_trna, params_trna, ssvals, :rm_a, trna_species, kdam_range)
res_rtca2_m = numerical_bistability_analysis(rtc_model_trna, params_trna, ssvals, :rm_a, trna_species, kdam_range2)
res_rtcr1_m = numerical_bistability_analysis(rtc_model_trna, params_trna, ssvals, :rm_r, trna_species, kdam_range)
res_rtcr2_m = numerical_bistability_analysis(rtc_model_trna, params_trna, ssvals, :rm_r, trna_species, kdam_range2)

prtcb1_m = scatter(x=kdam_range, y=res_rtcb1_m, name="mRNA RtcB ON", line=attr(color="b693ccff",width=6.5))
prtcb2_m = scatter(x=kdam_range2, y=res_rtcb2_m, name="mRNA RtcB OFF", line=attr(color="b693ccff",width=6.5))
prtcr1_m = scatter(x=kdam_range, y=res_rtcr1_m, name="mRNA RtcR ON", line=attr(color="4ca7a2ff",width=6.5))
prtcr2_m = scatter(x=kdam_range2, y=res_rtcr2_m, name="mRNA RtcR OFF", line=attr(color="4ca7a2ff",width=6.5))
prtca1_m = scatter(x=kdam_range, y=res_rtca1_m, name="mRNA RtcA ON", line=attr(color="e48080ff",width=6.5))
prtca2_m = scatter(x=kdam_range2, y=res_rtca2_m, name="mRNA RtcA OFF", line=attr(color="e48080ff",width=6.5))

p_rtcs = plot([prtcb1_m, prtcb2_m, prtcr1_m, prtcr2_m, prtca1_m, prtca2_m], 
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", title = "Numerical bistability analysis",
yaxis_title="Rtc proteins (μM)", 
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))

plot([ptrna1, ptrna2, prd1, prd2, prt1, prt2, prtcb1, prtcb2, prtcr1, prtcr2, prtca1, prtca2, prtcb1_m, prtcb2_m, prtcr1_m, prtcr2_m, prtca1_m, prtca2_m])