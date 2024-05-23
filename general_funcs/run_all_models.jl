using DifferentialEquations, PlotlyJS, DataFrames, Measures, LabelledArrays, BenchmarkTools, ModelingToolkit, TickTock

PATH = "/home/holliehindley/phd"
include("$PATH/general_funcs/solving.jl")

include("$PATH/growth_model/parameters/growth_model_params.jl")
include("$PATH/growth_model/parameters/gm_uM_parameters.jl")
include("$PATH/rtc_model/parameters/rtc_params.jl")
include("$PATH/combined_model/combined_params.jl")
include("$PATH/rtc_model/parameters/trna_params.jl")

include("$PATH/growth_model/model/growth_model.jl")
include("$PATH/combined_model/combined_model.jl")
include("$PATH/rtc_model/models/rtc_orig.jl")
include("$PATH/rtc_model/models/rtc_trna_model.jl")
include("$PATH/rtc_model/models/inhibition_models/rtc_inhibition_model.jl")
include("$PATH/rtc_model/models/inhibition_models/trna_inhib_models.jl")

solu_gm = sol(growth_model, init_gm, tspan, params_gm)
df_gm = create_solu_df(solu_gm, species_gm);
p_gm = plot_solu(df_gm)

params_dam = deepcopy(params_gm)
params_dam[abx] = 10
solu_gm_dam = sol(growth_model, ssvals_gm, tspan, params_dam)
df_gm_dam = create_solu_df(solu_gm_dam, species_gm);
p_gm_dam = plot_solu(df_gm_dam)

solu_gm_uM = sol(growth_model, init_gm_uM, tspan, params_gm_uM)
df_gm_uM = create_solu_df(solu_gm_uM, species_gm);
p_gm_uM = plot_solu(df_gm_uM)

params_dam_gm = deepcopy(params_gm_uM)
params_dam_gm[abx] = 10
solu_gm_uM_dam = sol(growth_model, ssvals_gm_uM, tspan, params_dam_gm)
df_gm_uM_dam = create_solu_df(solu_gm_uM_dam, species_gm);
p_gm_uM_dam = plot_solu(df_gm_uM_dam)

solu_comb = sol(combined_model, init_comb, tspan, params_comb)
df_comb = create_solu_df(solu_comb, species_comb);
p_comb = plot_solu(df_comb)

params_dam = deepcopy(params_comb)
params_dam[abx] = 4
solu_comb_dam = sol(combined_model, ssvals_comb, tspan, params_dam)
df_comb_dam = create_solu_df(solu_comb_dam, species_comb);
p_comb_dam = plot_solu(df_comb_dam)



solu_rtc = sol(rtc_model, init_rtc, tspan, params_rtc)
df_rtc = create_solu_df(solu_rtc, species_rtc)
p_rtc = plot_solu(df_rtc)

params_dam = deepcopy(params_rtc)
params_dam[kdam] = 0.1
solu_rtc_dam = sol(rtc_model, ssvals_rtc, tspan, params_dam)
df_rtc_dam = create_solu_df(solu_rtc_dam, species_rtc)
p_rtc_dam = plot_solu(df_rtc_dam)

# solu_rtc_inhib = sol(rtcb_inhib_model, init_inhib_rtca, tspan, params_inhib)
# df_rtc_inhib = create_solu_df(solu_rtc_inhib, species_inhib)
# p_rtc_inhib = plot_solu(df_rtc_inhib)


solu_trna = sol(rtc_trna_model, init_trna, tspan, params_trna)
df_trna = create_solu_df(solu_trna, species_trna)
p_trna = plot_solu(df_trna)

params_dam_trna = deepcopy(params_trna)
params_dam_trna[kdam] = 0.1
solu_trna_dam = sol(rtc_trna_model, ssvals_trna, tspan, params_dam_trna)
df_trna_dam = create_solu_df(solu_trna_dam, species_trna)
p_trna_dam = plot_solu(df_trna_dam)

# solu_trna_inhib = sol(rtcb_trna_inhib_model, init_trna_inhib_rtcb, tspan, params_trna_inhib)
# df_trna_inhib = create_solu_df(solu_trna_inhib, species_trna_inhib)
# p_trna = plot_solu(df_trna_inhib)


with_damage = [p_gm_dam, p_gm_uM_dam, p_comb_dam, p_rtc_dam, p_trna_dam]
without_damage = [p_gm, p_gm_uM, p_comb, p_rtc, p_trna]

model_names = ["growth_model", "growth_model_uM", "combined_model", "rtc_model", "trna_model"]

for (i,name) in zip(with_damage,model_names)
    open("$PATHgeneral_funcs/model_solutions/with_damage/$name.html", "w") do io
        PlotlyBase.to_html(io, i.plot)
    end
end

for (i,name) in zip(without_damage,model_names)
    open("$PATHgeneral_funcs/model_solutions/without_damage/$name.html", "w") do io
        PlotlyBase.to_html(io, i.plot)
    end
end

using ProgressBars

include("$PATH/rtc_model/models/rtc_orig.jl")

params1 = deepcopy(params_rtc)
params1[kin] = .0001


kdam_range = range(0,50,length=100)
res = sweep_param(rtc_model, ssvals_rtc, params1, kdam_range, kdam, species_rtc)

plot(scatter(x=kdam_range, y=res.rh), Layout(xaxis_title="kdam", yaxis_title="rh"))

# plot(scatter(x=kdam_range, y=repeat([params_rtc[lam]], length(kdam_range))))

lam1 = @. params_rtc[gr_c]*res.rh*params_rtc[g_max]*params_rtc[atp]/(params_rtc[θtlr]+params_rtc[atp])
plot(scatter(x=kdam_range, y=lam1))
