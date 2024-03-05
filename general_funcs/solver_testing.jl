using DifferentialEquations, PlotlyJS, DataFrames, Measures, LabelledArrays, BenchmarkTools, ModelingToolkit, TickTock

include("/home/holliehindley/phd/general_funcs/solving.jl")

include("/home/holliehindley/phd/growth_model/parameters/growth_model_params.jl")
include("/home/holliehindley/phd/growth_model/parameters/gm_uM_parameters.jl")
include("/home/holliehindley/phd/rtc_model/parameters/rtc_params.jl")
include("/home/holliehindley/phd/combined_model/combined_params.jl")
include("/home/holliehindley/phd/rtc_model/parameters/trna_params.jl")

include("/home/holliehindley/phd/growth_model/model/growth_model.jl")
include("/home/holliehindley/phd/combined_model/combined_model.jl")
include("/home/holliehindley/phd/rtc_model/models/rtc_orig.jl")
include("/home/holliehindley/phd/rtc_model/models/rtc_trna_model.jl")
include("/home/holliehindley/phd/rtc_model/models/inhibition_models/rtc_inhibition_model.jl")
include("/home/holliehindley/phd/rtc_model/models/inhibition_models/trna_inhib_models.jl")

params_dam = deepcopy(params_gm)
params_dam[abx] = 4
prob_gm = ODEProblem(growth_model, ssvals_gm, tspan, params_dam; jac=true)

params_dam = deepcopy(params_comb)
params_dam[abx] = 4
prob_comb = ODEProblem(combined_model, ssvals_comb, tspan, params_dam; jac=true)

@time solu_gm = solve(prob_gm, Rosenbrock23());#, abstol=1e-12, reltol=1e-9);
@time solu_gm = solve(prob_gm, Rodas4());#, abstol=1e-12, reltol=1e-9);
@time solu_gm = solve(prob_gm, TRBDF2());#, abstol=1e-9, reltol=1e-6);
@time solu_gm = solve(prob_gm, QNDF());#, abstol=1e-12, reltol=1e-9);
@time solu_gm = solve(prob_gm, FBDF());#, abstol=1e-9, reltol=1e-6);

df_gm = create_solu_df(solu_gm, species_gm); plot_solu(df_gm)

@time solu_comb = solve(prob_comb, Rosenbrock23());#, abstol=1e-12, reltol=1e-9);
@time solu_comb = solve(prob_comb, Rodas4());#, abstol=1e-12, reltol=1e-9);
@time solu_comb = solve(prob_comb, TRBDF2());#, abstol=1e-9, reltol=1e-6);
@time solu_comb = solve(prob_comb, QNDF());#, abstol=1e-9, reltol=1e-6);
@time solu_comb = solve(prob_comb, FBDF());#, abstol=1e-9, reltol=1e-6);

df_comb = create_solu_df(solu_comb, species_comb); plot_solu(df_comb)


prob_gm = SteadyStateProblem(growth_model, init_gm, params_gm; jac=true)
println("QNDF")
@benchmark solu = solve(prob_gm, DynamicSS(QNDF()))
println("FBDF")
@benchmark solu = solve(prob_gm, DynamicSS(FBDF()))
println("TRBDF2")
@benchmark solu = solve(prob_gm, DynamicSS(TRBDF2()))
println("Rodas4")
@benchmark solu = solve(prob_gm, DynamicSS(Rodas4()))
println("Rosenbrock23")
@benchmark solu = solve(prob_gm, DynamicSS(Rosenbrock23()))

prob_comb = SteadyStateProblem(combined_model, init_comb, params_comb; jac=true)
println("QNDF")
@benchmark solu = solve(prob_comb, DynamicSS(QNDF()))
println("FBDF")
@benchmark solu = solve(prob_comb, DynamicSS(FBDF()))
println("TRBDF2")
@benchmark solu = solve(prob_comb, DynamicSS(TRBDF2()))
println("Rodas4")
@benchmark solu = solve(prob_comb, DynamicSS(Rodas4()))
println("Rosenbrock23")
@benchmark solu = solve(prob_comb, DynamicSS(Rosenbrock23()))


prob_rtc = SteadyStateProblem(rtc_model, init_rtc, params_rtc; jac=true)
println("QNDF")
@benchmark solu = solve(prob_rtc, DynamicSS(QNDF()))
println("FBDF")
@benchmark solu = solve(prob_rtc, DynamicSS(FBDF()))
println("TRBDF2")
@benchmark solu = solve(prob_rtc, DynamicSS(TRBDF2()))
println("Rodas4")
@benchmark solu = solve(prob_rtc, DynamicSS(Rodas4()))
println("Rosenbrock23")
@benchmark solu = solve(prob_rtc, DynamicSS(Rosenbrock23()))


prob_trna = SteadyStateProblem(rtc_trna_model, init_trna, params_trna; jac=true)
println("QNDF")
@benchmark solu = solve(prob_trna, DynamicSS(QNDF()))
println("FBDF")
@benchmark solu = solve(prob_trna, DynamicSS(FBDF()))
println("TRBDF2")
@benchmark solu = solve(prob_trna, DynamicSS(TRBDF2()))
println("Rodas4")
@benchmark solu = solve(prob_trna, DynamicSS(Rodas4()))
println("Rosenbrock23")
@benchmark solu = solve(prob_trna, DynamicSS(Rosenbrock23()))



solu_gm = sol(growth_model, init_gm, tspan, params_gm)
df_gm = create_solu_df(solu_gm, species_gm);
p_gm = plot_solu(df_gm)

params_dam = deepcopy(params_gm)
params_dam[abx] = 4
solu_gm_dam = sol(growth_model, init_gm, tspan, params_dam)
df_gm_dam = create_solu_df(solu_gm_dam, species_gm);
p_gm_dam = plot_solu(df_gm_dam)

solu_gm_uM = sol(growth_model, init_gm_uM, tspan, params_gm_uM)
df_gm_uM = create_solu_df(solu_gm_uM, species_gm);
p_gm_uM = plot_solu(df_gm_uM)

params_dam_gm = deepcopy(params_gm_uM)
params_dam_gm[abx] = 4
solu_gm_uM_dam = sol(growth_model, init_gm_uM, tspan, params_dam_gm)
df_gm_uM_dam = create_solu_df(solu_gm_uM_dam, species_gm);
p_gm_uM_dam = plot_solu(df_gm_uM_dam)

solu_comb = sol(combined_model, init_comb, tspan, params_comb)
df_comb = create_solu_df(solu_comb, species_comb);
p_comb = plot_solu(df_comb)

params_dam = deepcopy(params_comb)
params_dam[abx] = 4
solu_comb_dam = sol(combined_model, init_comb, tspan, params_dam)
df_comb_dam = create_solu_df(solu_comb_dam, species_comb);
p_comb_dam = plot_solu(df_comb_dam)



solu_rtc = sol(rtc_model, init_rtc, tspan, params_rtc)
df_rtc = create_solu_df(solu_rtc, species_rtc)
p_rtc = plot_solu(df_rtc)

params_dam = deepcopy(params_rtc)
params_dam[kdam] = 0.1
solu_rtc_dam = sol(rtc_model, init_rtc, tspan, params_dam)
df_rtc_dam = create_solu_df(solu_rtc_dam, species_rtc)
p_rtc_dam = plot_solu(df_rtc_dam)

solu_rtc_inhib = sol(rtcb_inhib_model, init_inhib_rtca, tspan, params_inhib)
df_rtc_inhib = create_solu_df(solu_rtc_inhib, species_inhib)
p_rtc_inhib = plot_solu(df_rtc_inhib)


solu_trna = sol(rtc_trna_model, init_trna, tspan, params_trna)
df_trna = create_solu_df(solu_trna, species_trna)
p_trna = plot_solu(df_trna)

params_dam_trna = deepcopy(params_trna)
params_dam_trna[kdam] = 0.1
solu_trna_dam = sol(rtc_trna_model, init_trna, tspan, params_dam_trna)
df_trna_dam = create_solu_df(solu_trna_dam, species_trna)
p_trna_dam = plot_solu(df_trna_dam)

solu_trna_inhib = sol(rtcb_trna_inhib_model, init_trna_inhib_rtcb, tspan, params_trna_inhib)
df_trna_inhib = create_solu_df(solu_trna_inhib, species_trna_inhib)
p_trna = plot_solu(df_trna_inhib)

open("/home/holliehindley/phd/general_funcs/model_solutions/with_damage/growth_model_uM.html", "w") do io
    PlotlyBase.to_html(io, p_gm_uM_dam.plot)
end