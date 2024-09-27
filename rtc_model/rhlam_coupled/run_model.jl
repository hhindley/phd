model = "lamkin_coupled"
include(joinpath(homedir(), "phd/rtc_model/rhlam_coupled/models/$model.jl"))

colours =["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", :blue]

using GLMakie

# concentration no damage 
params_rtc1[kdam]
solu_rtc = sol(model, init_rtc, tspan, params_rtc1)
df = create_solu_df(solu_rtc, species_rtc)
lam = @. df.rh * tlr_el * lam_c_val
lines(df.time, lam)
lam[end]

p_rtc1 = plot([scatter(x=df.time, y=col, name="$(names(df)[i])", legendgroup="$i", marker_color=colours[i]) for (col, i) in zip(eachcol(df[:,2:end]), range(2,length(names(df))))], Layout(xaxis_type="log", title="kdam = $(params_rtc1[kdam])", xaxis_title="Time (s)", yaxis_title="Concentration (μM)"))

p = plot([scatter(x=df.time, y=df.rh*lam_c_val*tlr1, name="λ"),
scatter(x=df.time, y=df.rh*kin_c_val, yaxis="y2", name="kin")],
Layout(xaxis_title_text="Time (s)", yaxis_title_text="λ (min-1)", xaxis_type="log", 
yaxis2=attr(title="kin (μM aa-1)", overlaying="y", side="right", tickformat=".5f"),
legend=attr(x=0.1, y=0.9)))

open("/Users/s2257179/Desktop/dynamic_lam_kin.html", "w") do io
    PlotlyBase.to_html(io, p.plot)
end

# vary damage
fontsize_theme = Theme(fontsize = 25)
set_theme!(fontsize_theme)

kdam_range = range(0,500, length=200)
res = var_param(model, kdam, params_rtc1, kdam_range, ssvals_rtc)
f = Figure()
ax = Axis(f[1,1], xlabel="Damage rate (min⁻¹)", ylabel="RtcA (μM)")
lines!(ax, kdam_range, res.rtca, linewidth=5)

tlr_el = g_max_val*atp_val/(θtlr_val+atp_val)
lam_orig = @. res.rh * tlr_el * lam_c_val
# lam = @. 0.075 * (res.rh - 0.087)
lam_new = @. 0.001*(maximum(res.rh)[1] -res.rh)



f = Figure()
ax = Axis(f[1,1], xlabel="λ", ylabel="Rh")
lines!(ax, lam_orig, res.rh, linewidth=5)
lines!(ax, lam_new, res.rh, linewidth=5)


# hysteresis test
params1 = deepcopy(params_rtc1)
params1[kdam] = 0.01
solu_rtc_low = sol(test, ssvals_rtc, tspan, params1)
df_low = create_solu_df(solu_rtc_low, species_rtc)

params1[kdam] = 10
solu_rtc_high = sol(test, ssvals_rtc, tspan, params1)
df_high = create_solu_df(solu_rtc_high, species_rtc)

params1[kdam] = 1
solu_low = sol(test, [df_low[end,s] for s in species_rtc], tspan, params1)
solu_high = sol(test, [df_high[end,s] for s in species_rtc], tspan, params1)

df_low1 = create_solu_df(solu_low, species_rtc)
df_high1 = create_solu_df(solu_high, species_rtc)

plot([scatter(x=df_low.time, y=df_low.rtca, name="low: kdam = 0.01"), 
      scatter(x=df_high.time, y=df_high.rtca, name="high: kdam = 10"), 
      scatter(x=df_low1.time, y=df_low1.rtca, name="init low: kdam = 1"), 
      scatter(x=df_high1.time, y=df_high1.rtca, name="init high: kdam = 1")], 
      Layout(xaxis_type="log", yaxis_title="rtca", xaxis_title="Time (s)", title="hysteresis")
)



species = :rtca
p = plot([scatter(x=df.time, y=df[:,species], name="init"),
scatter(x=df1.time, y=df1[:,species], name="ss")],
Layout(xaxis_title_text="Time (min)", yaxis_title_text="$species", xaxis_type="log", 
legend=attr(x=0.1, y=0.9)))



init_rtc2 = [test.rm_a=>0.0,test.rtca=>1.0,test.rm_b=>0.0,test.rtcb=>0.0,test.rm_r=>0.0,test.rtcr=>0.0,test.rh=>11.29,test.rd=>0.0,test.rt=>0.0]
solu_rtc2 = sol(test, init_rtc2, tspan, params1)
df2 = create_solu_df(solu_rtc2, species_rtc)

init_rtc3 = [test.rm_a=>0.0,test.rtca=>1.0,test.rm_b=>0.0,test.rtcb=>0.0,test.rm_r=>0.0,test.rtcr=>0.0,test.rh=>19.29,test.rd=>0.0,test.rt=>0.0]
solu_rtc3 = sol(test, init_rtc3, tspan, params1)
df3 = create_solu_df(solu_rtc3, species_rtc)

species = :rm_a
p = plot([scatter(x=df.time, y=df[:,species], name="init"),
scatter(x=df1.time, y=df1[:,species], name="ss"),
scatter(x=df2.time, y=df2[:,species], name="test"),
scatter(x=df3.time, y=df3[:,species], name="test2")],
Layout(xaxis_title_text="Time (min)", yaxis_title_text="$species", xaxis_type="log", 
legend=attr(x=0.1, y=0.9)))


p = plot([scatter(x=df.time, y=df.rh*lam_c_val*tlr1, name="λ"),
scatter(x=df.time, y=df.rh*kin_c_val, yaxis="y2", name="kin"),
scatter(x=df1.time, y=df1.rh*lam_c_val*tlr1, name="λ ss"),
scatter(x=df1.time, y=df1.rh*kin_c_val, yaxis="y2", name="kin ss"),
scatter(x=df2.time, y=df2.rh*lam_c_val*tlr1, name="λ 15"),
scatter(x=df2.time, y=df2.rh*kin_c_val, yaxis="y2", name="kin 15")],
Layout(xaxis_title_text="Time (min)", yaxis_title_text="λ (min-1)", xaxis_type="log", 
yaxis2=attr(title="kin (μM aa-1)", overlaying="y", side="right", tickformat=".5f"),
legend=attr(x=0.1, y=0.9)))


p = plot([scatter(x=df.time, y=df.rtca, name="init"),
scatter(x=df1.time, y=df1.rtca, name="ss"),
scatter(x=df2.time, y=df2.rtca, name="test")],
Layout(xaxis_title_text="Time (min)", yaxis_title_text="rtca", xaxis_type="log", 
legend=attr(x=0.1, y=0.9)))

plot([scatter(x=df.time, y=df.rh, name="init"), scatter(x=df1.time, y=df1.rh, name="ss")], Layout(xaxis_type="log"))

