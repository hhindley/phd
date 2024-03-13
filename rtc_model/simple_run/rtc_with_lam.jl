using DifferentialEquations, PlotlyJS, DataFrames, Measures, LabelledArrays, BenchmarkTools, ModelingToolkit, TickTock

PATH = "/home/holliehindley/phd"
include("$PATH/general_funcs/solving.jl")

include("$PATH/rtc_model/parameters/rtc_params.jl")

include("$PATH/rtc_model/models/rtc_orig.jl")
include("$PATH/rtc_model/models/time_vary_params/vary_lam_2024.jl")

sol(model, )

solu_rtc = sol(rtc_model, init_rtc, tspan, params_rtc)
df_rtc = create_solu_df(solu_rtc, species_rtc)
p_rtc = plot_solu(df_rtc)


solu_rtc_lam = sol(rtc_model_lam, init_rtc, tspan, params_rtc)
df_rtc_lam = create_solu_df(solu_rtc_lam, species_rtc)
p_rtc = plot_solu(df_rtc_lam)

kdam_range=range(0,2,length=50)
kdam_rev = reverse(kdam_range)
res, res1 = full_numerical_bistab(rtc_model_lam, params_rtc, ssvals_rtc_lam, :rh, species_rtc, kdam_range, kdam_rev, kdam)

plot([scatter(x=kdam_range, y=res),scatter(x=kdam_rev, y=res1)])

plot([scatter(x=kdam_range, y=res*0.01),scatter(x=kdam_rev, y=res1*0.01)])


res, res1 = full_numerical_bistab(rtc_model, params_rtc, ssvals_rtc, :rh, species_rtc, kdam_range, kdam_rev, kdam)
plot([scatter(x=kdam_range, y=res),scatter(x=kdam_rev, y=res1)])
