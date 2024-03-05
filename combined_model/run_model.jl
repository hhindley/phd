using Parameters, CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, Statistics, Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra, Printf, PlotlyJS, ProgressBars, ModelingToolkit

include("/home/holliehindley/phd/general_funcs/solving.jl")
include("/home/holliehindley/phd/rtc_model/functions/bf_funcs/bf_funcs.jl")
include("/home/holliehindley/phd/rtc_model/functions/sweep_params.jl")

include("/home/holliehindley/phd/growth_model/parameters/growth_model_params.jl")
include("/home/holliehindley/phd/growth_model/parameters/gm_uM_parameters.jl")
include("/home/holliehindley/phd/rtc_model/parameters/rtc_params.jl")
include("/home/holliehindley/phd/combined_model/combined_params.jl")

include("/home/holliehindley/phd/growth_model/model/growth_model.jl")
include("/home/holliehindley/phd/combined_model/combined_model.jl")
include("/home/holliehindley/phd/rtc_model/models/rtc_orig.jl")


solu = sol(combined_model, init_comb, tspan, params_comb)
df = create_solu_df(solu, species_comb); plot_solu(df)

abx_range = range(0,350,length=100)
abx_rev = reverse(abx_range)

new_params = deepcopy(params_comb)
new_params[ns] = 0.1
new_params[kdam_p] = 0.3

res, res2 = full_numerical_bistab(combined_model, new_params, ssvals_comb, :B, species_comb, abx_range, abx_rev, abx)
plot([scatter(x=abx_range, y=res), scatter(x=abx_rev, y=res2)])

kdam_p_range = range(0.3,0.6,length=3)
full_res=[]
for i in kdam_p_range
    # new_params = deepcopy(params_comb)
    new_params[kdam_p] = i
    res, res2 = full_numerical_bistab(combined_model, new_params, ssvals_comb, :B, species_comb, abx_range, abx_rev, abx)
    push!(full_res, [res,res2])
end

traces=[]
for i in range(1, length(full_res))
    t1 = scatter(x=abx_range, y=full_res[i][1], name="% dam $(kdam_p_range[i]) forw") #for i in range(1,length(full_res))
    t2 = scatter(x=abx_rev, y=full_res[i][2], name="% dam $(kdam_p_range[i]) rev") #for i in range(1,length(full_res))
    push!(traces, t1)
    push!(traces, t2)
end

plot([i for i in traces])


new_params = deepcopy(params_comb)
# new_params[kdam_p] = 0.3
solu = sol(combined_model, init_comb, tspan, new_params); df = create_solu_df(solu, species_comb); plot_solu(df)



abx_range = range(0,12,length=50)
sweep_abx = sweep_param(combined_model, ssvals_comb, new_params, abx_range, abx, species_comb, "with repair")

open("/home/holliehindley/phd/general_funcs/model_solutions/test.html", "w") do io
    PlotlyBase.to_html(io, sweep_abx.plot)
end



params_abx = deepcopy(params_comb)
params_abx[abx] = 2
kdamp_range = range(0,1,length=50)
sweep_kdamp = sweep_param(combined_model, ssvals_comb, params_abx, kdamp_range, kdam_p, species_comb, "")

ns_range = range(0.01,1,length=100)
sweep_ns = sweep_param(combined_model, ssvals_comb, params_abx, ns_range, ns, species_comb, "")

s0_range = 10 .^ range(log10(0.001),log10(1e4),length=1000)
sweep_s0 = sweep_param(combined_model, ssvals_comb, params_abx, s0_range, s0, species_comb, "")

wab_range = 10 .^ range(log10(1e-6),log10(1e-1),length=100)
sweep_wab = sweep_param(combined_model, ssvals_comb, params_abx, wab_range, w_BA, species_comb, "")

wr_range = 10 .^ range(log10(1e-6),log10(1e-1),length=100)
sweep_wr = sweep_param(combined_model, ssvals_comb, params_abx, wr_range, w_R, species_comb, "")




abx_range = range(0,12,length=50)
kdamp_range = range(0,1,length=50)



res = sweep_paramx2(combined_model, ssvals_comb, params_comb, species_comb, abx, kdam_p, abx_range, kdamp_range)

new_species = deepcopy(species_comb); push!(new_species, :lam, :rmf);
[display(plot_contour(res, i, abx_range, kdamp_range, "abx", "kdam_p")) for i in new_species]

