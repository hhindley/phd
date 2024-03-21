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

# comparing lam under rtc and no rtc after introducing damage 
solu_dam = sol(combined_model, ssvals_comb, tspan, new_params(10, 0, 0)); df_dam = create_solu_df(solu_dam, species_comb); plot_solu(df_dam)
solu = sol(combined_model, init_comb, tspan, params_comb); df = create_solu_df(solu, species_comb); plot_solu(df)
solu_dam1 = sol(combined_model, ssvals_comb, tspan, new_params(10, ω_ab_val_comb, ω_r_val_comb)); df_dam1 = create_solu_df(solu_dam1, species_comb); plot_solu(df_dam1)

dict = Dict(pairs(eachcol(df))); lam = calc_lam(params_comb, dict)

dict_dam = Dict(pairs(eachcol(df_dam))); lam_dam = calc_lam(new_params(10,0,0), dict_dam)

dict_dam1 = Dict(pairs(eachcol(df_dam1))); lam_dam1 = calc_lam(new_params(10, ω_ab_val_comb, ω_r_val_comb), dict_dam1)

p = plot([scatter(x=df_dam1.time, y=lam_dam1, name="with dam/with rtc"), scatter(x=df_dam.time, y=lam_dam, name="with dam/no rtc")], Layout(xaxis_type="log"))
p1 = plot(scatter(x=df.time, y=lam, name="orig"), Layout(xaxis_type="log"))
[p1 p]


# comparing all species/variables with rtc and without rtc as abx concentration increases
abx_range = range(0,20,length=100)
w_rtc = sweep_param(combined_model, ssvals_comb, new_params(abx_val, ω_ab_val_comb, ω_r_val_comb), abx_range, abx, species_comb)
wo_rtc = sweep_param(combined_model, ssvals_comb, new_params(abx_val, 0, 0), abx_range, abx, species_comb)

function plot_comp(var)
    return plot([scatter(x=abx_range,y=w_rtc[!,var],name="with Rtc"), scatter(x=abx_range,y=wo_rtc[!,var],name="without Rtc")], Layout(xaxis_title="abx",yaxis_title="$var"))
end

new_species = deepcopy(species_comb); push!(new_species, :lam, :rmf);
[display(plot_comp(i)) for i in new_species]


# parameter sweeps
sweep_abx = sweep_param_with_plot(combined_model, ssvals_comb, params_comb, abx_range, abx, species_comb, "")

sweep_abx = sweep_param_with_plot(combined_model, ssvals_comb_new, params_comb_new, abx_range, abx, species_comb, "")

params_abx = deepcopy(params_comb)
params_abx[abx] = 10
kdamp_range = range(0,1,length=50)
sweep_kdamp = sweep_param_with_plot(combined_model, ssvals_comb, params_abx, kdamp_range, kdam_p, species_comb, "")

ns_range = range(0.01,1,length=50)
sweep_ns = sweep_param_with_plot(combined_model, ssvals_comb, params_abx, ns_range, ns, species_comb, "")

s0_range = 10 .^ range(log10(0.001),log10(1e4),length=50)
sweep_s0 = sweep_param_with_plot(combined_model, ssvals_comb, params_abx, s0_range, s0, species_comb, "")

wab_range = 10 .^ range(log10(1e-6),log10(1e-3),length=100)
sweep_wab = sweep_param_with_plot(combined_model, ssvals_comb, params_abx, wab_range, w_BA, species_comb, "")

wr_range = 10 .^ range(log10(1e-6),log10(1e-1),length=100)
sweep_wr = sweep_param_with_plot(combined_model, ssvals_comb, params_abx, wr_range, w_R, species_comb, "")



#double parameter sweep of abx and kdam_p with rtc and without rtc
abx_range = range(0,12,length=25)
kdamp_range = range(0,1,length=25)

res = sweep_paramx2(combined_model, ssvals_comb, params_comb, species_comb, abx, kdam_p, abx_range, kdamp_range)
res_nortc = sweep_paramx2(combined_model, ssvals_comb, new_params(0, 0, 0), species_comb, abx, kdam_p, abx_range, kdamp_range)

new_species = deepcopy(species_comb); push!(new_species, :lam, :rmf);

plot_contour(res, :rh, abx_range, kdamp_range, "abx", "kdam_p", "rtc active", 0.45, 1)
# [display(plot_contour(res, i, abx_range, kdamp_range, "abx", "kdam_p", "")) for i in new_species]

p2 = [(plot_contour(res, i, abx_range, kdamp_range, "abx", "kdam_p", "rtc active", 0.45)) for i in new_species];
p3 = [(plot_contour(res_nortc, i, abx_range, kdamp_range, "abx", "kdam_p", "no rtc", 1)) for i in new_species];

[p2[1] p3[1]]

# [[p2[i] for i in range(1,length(p2))] [p3[i] for i in range(1,length(p3))]]

for i in range(1,length(p2))
    display([p2[i] p3[i]])
end



# stability analysis 
abx_range = range(0,150,length=50)
abx_rev = reverse(abx_range)

n_params = deepcopy(params_comb_new)
# n_params[ns] = 0.5

res, res2 = full_numerical_bistab(combined_model, n_params, ssvals_comb_new, :B, species_comb, abx_range, abx_rev, abx)
plot([scatter(x=abx_range, y=res), scatter(x=abx_rev, y=res2)])

sweep_abx = sweep_param_with_plot(combined_model, ssvals_comb_new, params_comb_new, abx_range, abx, species_comb, "")

# n_params1 = deepcopy(params_comb)
# n_params1[ns] = 0.5
# n_params1[kdam_p] = 1

# res1, res21 = full_numerical_bistab(combined_model, n_params1, ssvals_comb, :B, species_comb, abx_range, abx_rev, abx)
# plot([scatter(x=abx_range, y=res1), scatter(x=abx_rev, y=res21)])






kdam_p_range = range(0.1,1,length=3)
full_res=[]
for i in kdam_p_range
    # new_params = deepcopy(params_comb)
    n_params[kdam_p] = i
    res, res2 = full_numerical_bistab(combined_model, n_params, ssvals_comb, :B, species_comb, abx_range, abx_rev, abx)
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




solu = sol(combined_model, init_comb, tspan, params_comb); df = create_solu_df(solu, species_comb); plot_solu(df)
solu_new = sol(combined_model, init_comb, tspan, params_comb_new); df_new = create_solu_df(solu_new, species_comb); plot_solu(df_new)

get_all_ssvals(solu_new, species_comb)