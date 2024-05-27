# using , Statistics, DataInterpolations, Printf
# using Distributions, BifurcationKit, ProgressBars, Suppressor
using ModelingToolkit, Distributions, Statistics, Suppressor, ProgressBars, DifferentialEquations, PlotlyJS, LinearAlgebra, DataFrames, LabelledArrays, Printf, BifurcationKit, OrderedCollections, LinearAlgebra, JLD2, Parameters, Revise

include("/home/hollie_hindley/Documents/paper/model_params_funcs_2024/solving.jl")
include("/home/hollie_hindley/Documents/paper/model_params_funcs_2024/rtc_orig.jl")
include("/home/hollie_hindley/Documents/paper/model_params_funcs_2024/params.jl")
include("/home/hollie_hindley/Documents/paper/model_params_funcs_2024/bf_funcs.jl")

kr_sd = kr_val*0.113 # estimated
Vmax_init_sd = 4.62
Km_init_sd = Km_init_val*0.1
d_sd = d_val*0.058
krep_sd = 31
ktag_sd = ktag_val*0.086
kma_sd = km_a_val*0.1
kmb_sd = km_b_val*0.1

kr_norm = Normal(kr_val, kr_sd)
Vmax_init_norm = Normal(Vmax_init_val, Vmax_init_sd)
Km_init_norm = Normal(Km_init_val, Km_init_sd)
d_norm = Normal(d_val, d_sd)
krep_norm = Normal(krep_val, krep_sd)
ktag_norm = Normal(ktag_val, ktag_sd)
kma_norm = Normal(km_a_val, kma_sd)
kmb_norm = Normal(km_b_val, kmb_sd)
krep_dist = Truncated(krep_norm, 0, 300)

kdeg_lnorm = LogNormal(kdeg_val)
# kin_lb = 0.01
# kin_ub = 2
# range_kin = kin_ub - kin_lb
# kin_std = 0.4
# kin_dist = Truncated(Normal(kin, kin_std), kin_lb, kin_ub)

# kdeg_lb = 0.0005
# kdeg_ub = 0.05
# range_kdeg = (kdeg_ub) - (kdeg_lb)
# kdeg_std = 0.015
# kdeg_dist = Truncated(Normal(kdeg, kdeg_std), kdeg_lb, kdeg_ub)

l=2000000
kr_new = rand(kr_norm,l)
Vmax_init_new = rand(Vmax_init_norm,l)
Km_init_new = rand(Km_init_norm,l)
d_new = rand(d_norm,l)
krep_new = rand(krep_dist,l)
ktag_new = rand(ktag_norm,l)
kdeg_new = rand(kdeg_lnorm,l)

kma_new = rand(kma_norm,l)
kmb_new = rand(kmb_norm,l)

function analysis()
    res_params=[]; res_df=[]; res_notbs_params=[]; res_notbs_df=[];
    for i in ProgressBar(range(1, length(kr_new)))
        params_new = deepcopy(params_rtc)
        params_new[kr] = kr_new[i]
        params_new[Vmax_init] = Vmax_init_new[i]
        params_new[Km_init] = Km_init_new[i]
        params_new[d] = d_new[i]
        params_new[krep] = krep_new[i]
        params_new[ktag] = ktag_new[i]
        params_new[km_a] = kma_new[i]
        params_new[km_b] = kmb_new[i]
        params_new[kdeg] = kdeg_new[i]

        br = get_br(rtc_model, ssvals_rtc, params_new, 3.)
        if length(br.specialpoint) == 4
            df = create_br_df(br)
            push!(res_df,df[!,[:kdam,:rh]])
            push!(res_params, @LArray [kr_new[i],Vmax_init_new[i],Km_init_new[i],d_new[i],krep_new[i],ktag_new[i],kma_new[i],kmb_new[i],kdeg_new[i]] (:kr,:Vmax_init,:Km_init,:d,:krep,:ktag,:km_a,:km_b,:kdeg))
        else
            df = create_br_df(br)
            push!(res_notbs_df,df[!,[:kdam,:rh]])
            push!(res_notbs_params, @LArray [kr_new[i],Vmax_init_new[i],Km_init_new[i],d_new[i],krep_new[i],ktag_new[i],kma_new[i],kmb_new[i],kdeg_new[i]] (:kr,:Vmax_init,:Km_init,:d,:krep,:ktag,:km_a,:km_b,:kdeg))
        end
    end
    return res_params, res_df, res_notbs_params, res_notbs_df
end

@time res_param, res_df, res_notbs_params, res_notbs_df = analysis();

# (length(res_param)/l)*100
# for i in range(1,length(res))
#     @show res[i].kin
# end

df_params = DataFrame(res_param)
df_solus = DataFrame(res_df)
df_params_notbs = DataFrame(res_notbs_params)
df_solus_notbs = DataFrame(res_notbs_df)


dict_solus=Dict(pairs(eachcol(df_solus)))
dict_params=Dict(pairs(eachcol(df_params)))
dict_solus_notbs=Dict(pairs(eachcol(df_solus_notbs)))
dict_params_notbs=Dict(pairs(eachcol(df_params_notbs)))

println((length(dict_solus[:kdam])/l)*100)
# plts=[]
# for (i,j) in zip(dict_solus[:kdam], dict_solus[:rh])
#     push!(plts, (scatter(x=i, y=j)))
# end

# plot([(i) for i in plts])

# plts=[]
# for (i,j) in zip(dict_solus_notbs[:kdam], dict_solus_notbs[:rh])
#     push!(plts, (scatter(x=i, y=j)))
# end

# plot([(i) for i in plts])


@suppress begin
    save("/home/hollie_hindley/Documents/paper/big_sensitivity_2024/solus_nokin.jld2",dict_solus)
    save("/home/hollie_hindley/Documents/paper/big_sensitivity_2024/params_nokin.jld2",dict_params)
    save("/home/hollie_hindley/Documents/paper/big_sensitivity_2024/solus_notbs_nokin.jld2",dict_solus_notbs)
    save("/home/hollie_hindley/Documents/paper/big_sensitivity_2024/params_notbs_nokin.jld2",dict_params_notbs)
    @warn("Ignoring the string warning")
end


# a1 = load("/home/hollie_hindley/Documents/data/solus_nokin.jld2")

# a2 = load("/home/hollie_hindley/Documents/data/solus_notbs_nokin.jld2")

# a1["kdam"]
# a2["kdam"]
# using Plots


# p1=plot();

# for (i,j) in zip(a1["kdam"], a1["rh"])
#     plot!(p1, i,j, legend=false)
# end

# p1



