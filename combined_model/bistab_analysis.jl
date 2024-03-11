using Parameters, CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, Statistics, Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra, Printf, PlotlyJS, ProgressBars, ModelingToolkit

PATH = "/home/holliehindley/phd"

include("$PATH/general_funcs/solving.jl")
include("$PATH/rtc_model/functions/bf_funcs/bf_funcs.jl")
include("$PATH/rtc_model/functions/sweep_params.jl")
include("$PATH/growth_model/parameters/growth_model_params.jl")
include("$PATH/growth_model/parameters/gm_uM_parameters.jl")
include("$PATH/rtc_model/parameters/rtc_params.jl")
include("$PATH/combined_model/combined_params.jl")
include("$PATH/growth_model/model/growth_model.jl")
include("$PATH/combined_model/combined_model.jl")
include("$PATH/rtc_model/models/rtc_orig.jl")


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