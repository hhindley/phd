using ModelingToolkit, DifferentialEquations, LinearAlgebra, DataFrames, LabelledArrays, Printf, BifurcationKit, OrderedCollections, ProgressBars

include(joinpath(homedir(), "phd/general_funcs/all_model_funcs.jl"))
include(joinpath(homedir(), "phd/general_funcs/solving.jl"))
include(joinpath(homedir(), "phd/rtc_model/parameters/rtc_params.jl"))
include(joinpath(homedir(), "phd/rtc_model/parameters/rtc_params_molecs.jl"))

include(joinpath(homedir(), "phd/rtc_model/models/rtc_orig.jl"))
include(joinpath(homedir(), "phd/rtc_model/functions/bf_funcs/bf_funcs.jl"))
include("/Users/s2257179/phd/rtc_model/paper/server_code/model_params_funcs_2024/solving.jl")


jac_sym = calculate_jacobian(rtc_model)


solu = sol(rtc_model, init_rtc, tspan, params_rtc)
df = create_solu_df(solu, species_rtc)

params1 = deepcopy(params_rtc)
kdam_range = [0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5]

df_ssvals = var_param(rtc_model, kdam, params1, kdam_range, ssvals_rtc)

cct = log(2)/lam_val

eigs_dict = Dict(kdam => [] for kdam in kdam_range)
for i in eachindex(kdam_range)
    eigs = calc_eigenvalues(i, params1, kdam_range[i], df_ssvals)
    inv_eigs = -1 ./ eigs
    eigs_dict[kdam_range[i]] = inv_eigs./cct 
end
eigs_dict
max_vals=[]
for kdam in kdam_range
    push!(max_vals, maximum(eigs_dict[kdam]))
end
maximum(max_vals)*cct # 350


