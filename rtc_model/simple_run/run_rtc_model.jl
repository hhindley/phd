using ModelingToolkit, DifferentialEquations, PlotlyJS, LinearAlgebra, DataFrames, LabelledArrays, Printf, BifurcationKit, OrderedCollections, ProgressBars

include(joinpath(homedir(), "phd/general_funcs/all_model_funcs.jl"))
include(joinpath(homedir(), "phd/general_funcs/solving.jl"))
include(joinpath(homedir(), "phd/rtc_model/parameters/rtc_params.jl"))
include(joinpath(homedir(), "phd/rtc_model/parameters/rtc_params_molecs.jl"))

include(joinpath(homedir(), "phd/rtc_model/models/rtc_orig.jl"))
include(joinpath(homedir(), "phd/rtc_model/functions/bf_funcs/bf_funcs.jl"))

# colours =["#636EFA", "#EF553B", "#00CC96", "#AB63FA", "#FFA15A", "#19D3F3", "#FF6692", "#B6E880", "#FF97FF", "#FECB52", :blue]
colours =["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", :blue]

# concentration no damage 
solu_rtc = sol(rtc_model, init_rtc, tspan, params_rtc)
df = create_solu_df(solu_rtc, species_rtc)
p_rtc1 = plot([scatter(x=df.time, y=col, name="$(names(df)[i])", legendgroup="$i", marker_color=colours[i]) for (col, i) in zip(eachcol(df[:,2:end]), range(2,length(names(df))))], Layout(xaxis_type="log", title="kdam = $(params_rtc[kdam])", xaxis_title="Time (s)", yaxis_title="Concentration (μM)"))
# p_rtc1 = plot([scatter(x=df.time, y=col, name="$(names(df)[i])", legendgroup="$i", marker_color=colours[i]) for (col, i) in zip(eachcol(df[:,2:end]), range(2,length(names(df))))], Layout(yaxis_tickformat=".2e", title="kdam = $(params_rtc[kdam])"))

# open("/Users/s2257179/Desktop/solu_nodam_conc.html", "w") do io
#     PlotlyBase.to_html(io, p_rtc1.plot)
# end

# concentration with damage
params_kdam = deepcopy(params_rtc)
params_kdam[kdam] = 0.1
solu_rtc_kdam = sol(rtc_model, ssvals_rtc, tspan, params_kdam)
df_kdam = create_solu_df(solu_rtc_kdam, species_rtc)
p_kdam = plot([scatter(x=df_kdam.time, y=col, name="$(names(df_kdam)[i])", legendgroup="$i", marker_color=colours[i]) for (col, i) in zip(eachcol(df_kdam[:,2:end]), range(2,length(names(df_kdam))))], Layout(xaxis_type="log", title="kdam = $(params_kdam[kdam])", xaxis_title="Time (s)", yaxis_title="Concentration (μM)"))#, yaxis_tickformat=".2e"))

# p_conc = [p_rtc1 p_kdam]
# open("/Users/s2257179/Desktop/conc.html", "w") do io
#     PlotlyBase.to_html(io, p_conc.plot)
# end

# molecule number no damage
solu_rtc_molec = sol(rtc_model, init_rtc_molec, tspan, params_rtc_molec)
df_molec = create_solu_df(solu_rtc_molec, species_rtc)
p_rtc1_molec = plot([scatter(x=df_molec.time, y=col, name="$(names(df_molec)[i])", legendgroup="$i", marker_color=colours[i]) for (col, i) in zip(eachcol(df_molec[:,2:end]), range(2,length(names(df_molec))))], Layout(xaxis_type="log", yaxis_tickformat=".2e", title="kdam = $(params_rtc_molec[kdam])", xaxis_title="Time (s)", yaxis_title="Molecule number"))

# molecule number with damage 
params_kdam_molec = deepcopy(params_rtc_molec)
params_kdam_molec[kdam] = 0.1
solu_rtc_kdam_molec = sol(rtc_model, ssvals_rtc_molec, tspan, params_kdam_molec)
df_kdam_molec = create_solu_df(solu_rtc_kdam_molec, species_rtc)
p_kdam_molec = plot([scatter(x=df_kdam_molec.time, y=col, name="$(names(df_kdam_molec)[i])", legendgroup="$i", marker_color=colours[i]) for (col, i) in zip(eachcol(df_kdam_molec[:,2:end]), range(2,length(names(df_kdam_molec))))], Layout(xaxis_type="log", title="kdam = $(params_kdam[kdam])", xaxis_title="Time (s)", yaxis_title="Molecule number", yaxis_tickformat=".2e"))

# p_molec = [p_rtc1_molec p_kdam_molec]

p_all = [p_rtc1 p_kdam; p_rtc1_molec p_kdam_molec]
open("/Users/s2257179/Desktop/solutions.html", "w") do io
    PlotlyBase.to_html(io, p_all.plot)
end





[p_rtc p_kdam]


params_molec_kdam = deepcopy(params_rtc_molec)
params_molec_kdam[kdam] = 0.6

solu_rtc_kdam_molec = sol(rtc_model, ssvals_rtc_molec, (0,10000), params_molec_kdam)

df_kdam_molec = create_solu_df(solu_rtc_kdam_molec, species_rtc)

p_kdam_molec = plot([scatter(x=df_kdam_molec.time, y=col, name="$(names(df_kdam_molec)[i])", legendgroup="$i", marker_color=colours[i]) for (col, i) in zip(eachcol(df_kdam_molec[:,2:end]), range(2,length(names(df_kdam_molec))))], Layout(xaxis_type="log", yaxis_tickformat=".2e", title="kdam = $(params_molec_kdam[kdam])"))
p_kdam_molec = plot([scatter(x=df_kdam_molec.time, y=col, name="$(names(df_kdam_molec)[i])", legendgroup="$i", marker_color=colours[i]) for (col, i) in zip(eachcol(df_kdam_molec[:,2:end]), range(2,length(names(df_kdam_molec))))], Layout(xaxis_type="", yaxis_tickformat=".2e", title="kdam = $(params_molec_kdam[kdam])"))


# checking mRNA/protein ratios 
a_ratio = df.rtca[end]/df.rm_a[end]
b_ratio = df.rtcb[end]/df.rm_b[end]
r_ratio = df.rtcr[end]/df.rm_r[end]

df.rh[end]


br = get_br(rtc_model, ssvals_rtc, params_rtc, 1.5)
df = create_br_df(br)
plot(scatter(x=df.kdam, y=df.rtca), Layout(xaxis_title="kdam", yaxis_title="rtca"))

params_rtc1 = deepcopy(params_rtc)
params_rtc1[kdam] = 1.1
kdeg_range=range(0,0.1,length=100)
df = sweep_param(rtc_model, ssvals_rtc, params_rtc1, kdeg_range, kdeg, species_rtc)

plot(scatter(x=kdeg_range, y=df.rd))


kdam_range = range(0,1.5,length=50)
kdam_range1 = reverse(kdam_range)
res_ss = numerical_bistability_analysis(rtc_model, params_rtc, init_rtc, species_rtc, kdam_range, kdam)
res1_ss = numerical_bistability_analysis(rtc_model, params_rtc, ssvals_rtc, species_rtc, kdam_range1, kdam)

plot([scatter(x=kdam_range, y=res_ss.rh), scatter(x=kdam_range1, y=res1_ss.rh)])


include(joinpath(homedir(), "phd/rtc_model/paper/server_code/model_params_funcs_2024/solving.jl"))

df_ssvals = var_param(rtc_model, kdam, params_rtc, kdam_range, ssvals_rtc)
plot(scatter(x=kdam_range, y=df_ssvals.rtca), Layout(xaxis_title="kdam", yaxis_title="rtca"))
