using PlotlyJS, DataFrames, DifferentialEquations, Parameters, LabelledArrays, StaticArrays, Colors

include("/home/holliehindley/phd/rtc_model/models/rtc_orig.jl")
include("/home/holliehindley/phd/general_funcs/solving.jl")
include("/home/holliehindley/phd/rtc_model/parameters/init.jl")
include("/home/holliehindley/phd/rtc_model/parameters/params.jl")

colours =["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", :blue]

solu_rtc = sol(rtc_model, init_rtc, tspan, params_rtc)

df = create_solu_df(solu_rtc, species_rtc)

p_rtc1 = plot([scatter(x=df.time, y=col, name="$(names(df)[i])", legendgroup="$i", marker_color=colours[i]) for (col, i) in zip(eachcol(df[:,2:end]), range(2,length(names(df))))], Layout(xaxis_type="log", yaxis_tickformat=".2e", title="kdam = $(params_rtc.kdam)"))

open("./rtc.html", "w") do io
    PlotlyBase.to_html(io, p_rtc1.plot)
end

# with damage 
params_kdam = deepcopy(params_rtc)
params_kdam.kdam = 1.4

solu_rtc_kdam = sol(rtc_model, rtc_init, tspan, params_kdam)

df_kdam = create_solu_df(solu_rtc_kdam, species_rtc)

p_kdam = plot([scatter(x=df_kdam.time, y=col, name="$(names(df_kdam)[i])", legendgroup="$i", marker_color=colours[i]) for (col, i) in zip(eachcol(df_kdam[:,2:end]), range(2,length(names(df_kdam))))], Layout(xaxis_type="log", yaxis_tickformat=".2e", title="kdam = $(params_kdam.kdam)"))

[p_rtc p_kdam]

# checking mRNA/protein ratios 
a_ratio = df.rtca[end]/df.rm_a[end]
b_ratio = df.rtcb[end]/df.rm_b[end]
r_ratio = df.rtcr[end]/df.rm_r[end]

df.rh[end]

