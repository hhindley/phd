using Parameters, CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, Statistics
using Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra, Printf
# using Plots
using PlotlyJS, ProgressBars, Interpolations, QuadGK

include("/home/holliehindley/phd/may23_rtc/functions/bf_funcs/bf_funcs.jl");
include("/home/holliehindley/phd/may23_rtc/models/rtc_orig.jl");
include("/home/holliehindley/phd/may23_rtc/rtc_parameters/params.jl");
include("/home/holliehindley/phd/may23_rtc/rtc_parameters/init.jl");
include("/home/holliehindley/phd/may23_rtc/functions/solving.jl");

colours = ["#f4e5ffff","#e6c5ffff","#d7a1ffff","#c06affff","#a730ffff"]

atp_range = range(500,stop=5500,length=100)
lam_range = range(0.001,stop=0.04,length=100)

wab_range1 = 10 .^range(log10(1e-7),log10(1e-3),length=5)
results_wab, xvals = get_bs_region_results(atp_range, :atp, lam_range, :lam, wab_range1, :ω_ab, params_bf, 1.5)

function area_under_curve(xvals,yvals)
    x=xvals
    y=yvals
    int_orig = Interpolations.LinearInterpolation(x, y)
    f(x) = int_orig(x)
    a = minimum(x)
    b = maximum(x)

    result, error = quadgk(f, a, b)
    return @LArray [x,y,f,result] (:x,:y,:f1,:result)
end

auc1 = area_under_curve(xvals[1],results_wab[1][1])
auc2 = area_under_curve(xvals[2],results_wab[2][1])
auc3 = area_under_curve(xvals[3],results_wab[3][1])
auc4 = area_under_curve(xvals[4],results_wab[4][1])
auc5 = area_under_curve(xvals[5],results_wab[5][1])

wab_banana = plot([
scatter(x=xvals[1],y=results_wab[1][2], line=attr(color=colours[1])), scatter(x=auc1.x,y=auc1.f1.(auc1.x), fill="tonexty", mode="lines", line=attr(color=colours[1])),
scatter(x=xvals[2],y=results_wab[2][2], line=attr(color=colours[2])), scatter(x=auc2.x,y=auc2.f1.(auc2.x), fill="tonexty", mode="lines", line=attr(color=colours[2])),
scatter(x=xvals[3],y=results_wab[3][2], line=attr(color=colours[3])), scatter(x=auc3.x,y=auc3.f1.(auc3.x), fill="tonexty", mode="lines", line=attr(color=colours[3])),
scatter(x=xvals[4],y=results_wab[4][2], line=attr(color=colours[4])), scatter(x=auc4.x,y=auc4.f1.(auc4.x), fill="tonexty", mode="lines", line=attr(color=colours[4])),
scatter(x=xvals[5],y=results_wab[5][2], line=attr(color=colours[5])), scatter(x=auc5.x,y=auc5.f1.(auc5.x), fill="tonexty", mode="lines", line=attr(color=colours[5]))],
Layout(xaxis_title="ATP (μM)", yaxis_title="λ (min<sup>-1</sup>)",showlegend=false,
yaxis=attr(showline=true,linewidth=3,linecolor="black",range=(0.001,0.04)),xaxis=attr(showline=true,linewidth=3,linecolor="black", range=(500,5500)),# showlegend=false,
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))


savefig(wab_banana, "/home/holliehindley/phd/may23_rtc/paper_plots/wab_banana_final.svg")


wr_range1 = 10 .^ range(log10(1e-8),log10(1e-4),length=5)
results_wr, xvals_wr = get_bs_region_results(atp_range, :atp, lam_range, :lam, wr_range1, :ω_r, params_bf, 1.5)


auc1_wr = area_under_curve(xvals_wr[1],results_wr[1][1])
auc2_wr = area_under_curve(xvals_wr[2],results_wr[2][1])
auc3_wr = area_under_curve(xvals_wr[3],results_wr[3][1])
auc4_wr = area_under_curve(xvals_wr[4],results_wr[4][1])
auc5_wr = area_under_curve(xvals_wr[5],results_wr[5][1])

wr_banana = plot([
scatter(x=xvals_wr[1],y=results_wr[1][2], line=attr(color=colours[1])), scatter(x=auc1_wr.x,y=auc1_wr.f1.(auc1_wr.x), fill="tonexty", mode="lines", line=attr(color=colours[1])),
scatter(x=xvals_wr[2],y=results_wr[2][2], line=attr(color=colours[2])), scatter(x=auc2_wr.x,y=auc2_wr.f1.(auc2_wr.x), fill="tonexty", mode="lines", line=attr(color=colours[2])),
scatter(x=xvals_wr[3],y=results_wr[3][2], line=attr(color=colours[3])), scatter(x=auc3_wr.x,y=auc3_wr.f1.(auc3_wr.x), fill="tonexty", mode="lines", line=attr(color=colours[3])),
scatter(x=xvals_wr[4],y=results_wr[4][2], line=attr(color=colours[4])), scatter(x=auc4_wr.x,y=auc4_wr.f1.(auc4_wr.x), fill="tonexty", mode="lines", line=attr(color=colours[4])),
scatter(x=xvals_wr[5],y=results_wr[5][2], line=attr(color=colours[5])), scatter(x=auc5_wr.x,y=auc5_wr.f1.(auc5_wr.x), fill="tonexty", mode="lines", line=attr(color=colours[5]))],
Layout(xaxis_title="ATP (μM)", yaxis_title="λ (min<sup>-1</sup>)", showlegend=false,
yaxis=attr(showline=true,linewidth=3,linecolor="black",range=(0.001,0.04)),xaxis=attr(showline=true,linewidth=3,linecolor="black", range=(500,5500)),# showlegend=false,
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))

savefig(wr_banana, "/home/holliehindley/phd/may23_rtc/paper_plots/wr_banana_final.svg")








wr_range = 10 .^ range(log10(1e-7),log10(1e-4),length=50)
wab_range = 10 .^range(log10(1e-6),log10(1e-3),length=50)
atp_range = range(500,5000, length=2)
lam_range = range(0.01,0.04, length=2)

results_wab, xvals = get_bs_region_results(wab_range, :ω_ab, wr_range, :ω_r, lam_range, :lam, params_bf, 1.5)


function area_under_curve(xvals,yvals)
    x=xvals
    y=yvals
    int_orig = Interpolations.LinearInterpolation(x, y)
    f(x) = int_orig(x)
    a = minimum(x)
    b = maximum(x)

    result, error = quadgk(f, a, b)
    return @LArray [x,y,f,result] (:x,:y,:f1,:result)
end

auc1 = area_under_curve(xvals[1],results_wab[1][1])
auc2 = area_under_curve(xvals[2],results_wab[2][1])

wab_banana = plot([
scatter(x=xvals[1],y=results_wab[1][2]), scatter(x=auc1.x,y=auc1.f1.(auc1.x), fill="tonexty", mode="lines"),
scatter(x=xvals[2],y=results_wab[2][2]), scatter(x=auc2.x,y=auc2.f1.(auc2.x), fill="tonexty", mode="lines"),],
Layout(#xaxis_title="ATP (μM)", yaxis_title="λ (min<sup>-1</sup>)",showlegend=false,
yaxis=attr(showline=true,linewidth=3,linecolor="black",range=(1e-7,1e-4)),xaxis=attr(range=(1e-6,1e-3),showline=true,linewidth=3,linecolor="black"),# showlegend=false,
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=24, color="black", family="sans-serif")))


















df = double_param_vary(atp_range, :atp, lam_range, :lam, params_bf, 1.5)
bsp = df[df.bs .== :bp, :]
bsp[bsp[:,1] .== 500, :][:,1]



res_f, res_r, combos = numerical_param(rtc_model, params_rtc, rtc_init, species_rtc, atp_range, :atp, lam_range, :lam, kdam_range)
bistables = bistable_q(res_f, res_r, combos)

df_n = DataFrame(atp=[], lam=[])
for i in bistables
    push!(df_n.atp, i[1])
    push!(df_n.lam, i[2])
end

i=2
df_n[df_n.atp .== atp_range[i],:].lam == bsp[bsp.atp .== atp_range[i],:].wab
df_n[df_n.atp .== atp_range[i],:].lam
bsp[bsp.atp .== atp_range[i],:].wab

kdam_range=range(0,1.5,length=20)

params1 = deepcopy(params_rtc)
params1.atp = atp_range[i]
p = plot()
for i in lam_range
    params1.lam = i
    res1 = numerical_bistability_analysis(rtc_model, params1, rtc_init, :rtcb, species_rtc, kdam_range)
    res2 = numerical_bistability_analysis(rtc_model, params1, rtc_init, :rtcb, species_rtc, reverse(kdam_range))
    plot!(kdam_range,res1)
    plot!(reverse(kdam_range), res2)
end
p

diff = res1 - reverse(res2)
for i in diff 
    if i>0.001
        print("bistable")
    end
end

plot(kdam_range, res1)


















function numerical_param(rtc_model, params_, init_rtc, all_species, param_range1, param1, param_range2, param2, kdam_range)
    # params_trna1 = deepcopy(params_trna)
    res_f=[]
    res_r=[]
    param_combos=[]
    for i in ProgressBar(param_range1)
        # res1=[]
        params_[param1] = i
        for j in param_range2
            params_[param2] = j
            # @show params_[param1], params_[param2], params_[:ω_ab]
            res_rtcb_f = numerical_bistability_analysis(rtc_model, params_, init_rtc, :rtcb, all_species, kdam_range)
            res_rtcb_r = numerical_bistability_analysis(rtc_model, params_, init_rtc, :rtcb, all_species, reverse(kdam_range))
            push!(res_f, res_rtcb_f)
            push!(res_r, res_rtcb_r)
            push!(param_combos, (i,j))
        end   
        # push!(res, res1) 
    end
    return res_f, res_r, param_combos
end

# res_f, res_r, combos = numerical_param(rtc_model, params_rtc, rtc_init, species_rtc, atp_range, :atp, lam_range, :lam, kdam_range_f, kdam_range_r)
# bistables=[]; 
# for j in range(1, length(res_f))
#     # gap_before = res[j][i-1]-res[j][i]
#     # current_gap = res[j][i]-res[j][i+1]
#     diff = res_f[j]-reverse(res_r[j])
#     for i in diff
#         if i>0.001
#         # push!(bistables, ((res[j][i+1], res[j][i])))
#             push!(bistables, combos[j])
#         break 

#         end
#     end
# end
# bistables
# plot([scatter(x=kdam_range, y=res[i], name="$(param_combos[i])") for i in range(1,length(res))])

function bistable_q(res_f, res_r, param_combos)
    bistables=[]; 
    for j in range(1, length(res_f))
        # gap_before = res[j][i-1]-res[j][i]
        # current_gap = res[j][i]-res[j][i+1]
        diff = res_f[j]-reverse(res_r[j])
        for i in diff
            if i>0.001
            # push!(bistables, ((res[j][i+1], res[j][i])))
                push!(bistables, param_combos[j])
            break 

            end
        end
    end
    # bistables=[]; 
    # for j in range(1, length(res))
    #     for i in range(2,length(res[j])-1)
    #         # gap_before = res[j][i-1]-res[j][i]
    #         # current_gap = res[j][i]-res[j][i+1]
    #         if res[j][i] > 10 * res[j][i+1]
    #             # push!(bistables, ((res[j][i+1], res[j][i])))
    #             push!(bistables, param_combos[j])
    #             break 
    #         # else
    #         #     print("monostable")
    #         end
    #     end
    # end

    return bistables
end

bistables = bistable_q(res_n, combos)

# plot([scatter(x=[bistables[i][1]],y=[bistables[i][2]], mode="markers") for i in range(1, length(bistables))], Layout(xaxis_range=(400,5500), yaxis_range=(0.001,0.04)))


function get_max_min_yvals(param_range1, bistables)
    max_ = Float64[]
    min_ = Float64[]
    param1_vals1 = []
    for i in param_range1
        vals=[]
        param1_vals=[]
        for j in range(1, length(bistables))
            if i == bistables[j][1]
                push!(vals, bistables[j][2])
                push!(param1_vals, bistables[j][1])
            end
        end
        if length(vals) > 0 
            push!(max_, maximum(vals))
            push!(min_, minimum(vals))
            push!(param1_vals1, param1_vals[1])
        end
    end
    return max_,min_,param1_vals1
end

function numerical_triple_param_plotting(rtc_model, params_, init_rtc, all_species, param_range1, param1, param_range2, param2, kdam_range, w_range, w_change)
    p = Plots.plot()
    # max_s=[]
    # min_s=[]
    params_1 = deepcopy(params_)
    for i in w_range
        params_1[w_change] = i
        # @show params_1[w_change], params_1[param1], params_1[param2]
	    res_f, res_r, param_combos = numerical_param(rtc_model, params_1, init_rtc, all_species, param_range1, param1, param_range2, param2, kdam_range)
        bistables = bistable_q(res_f, res_r, param_combos)
        max_,min_,param1_vals1 = get_max_min_yvals(param_range1, bistables)
        # push!(max_s, max_)
        # push!(min_s, min_)
        p = Plots.plot!(param1_vals1, max_; fillrange=(min_), fillalpha = 0.45, 
        label="$(@sprintf "%g" i)", xlims=(minimum(param_range1), maximum(param_range1)), 
        ylims=(minimum(param_range2), maximum(param_range2)), xlabel="$param1", ylabel="$param2", tickfontsize=15,
        guidefontsize=18,guidefont="sans-serif",)
    end

    
    return p
end

 
kdam_range = range(0,1.5,length=10)

# atp_range = range(500,stop=5500,length=30)
# lam_range = range(0.01,stop=0.04,length=30)

# wab_range1 = 10 .^range(log10(2e-6),log10(2e-3),length=2)
p_atplam_wab = numerical_triple_param_plotting(rtc_model, params_rtc, rtc_init, species_rtc, atp_range, :atp, lam_range, :lam, kdam_range, wab_range1, :ω_ab)
Plots.savefig(p_atplam_wab, "/home/holliehindley/phd/may23_rtc/paper_plots/tRNA/p_atlam_wab.svg")

wr_range1 = 10 .^ range(log10(0.01/5e7),log10(0.01/1e4),length=5)

p_atplam_wr = numerical_triple_param_plotting(atp_range, :atp, lam_range, :lam, params_trna, kdam_range, wr_range1, :ω_r)
Plots.savefig(p_atplam_wr, "/home/holliehindley/phd/may23_rtc/paper_plots/tRNA/p_atlam_wr.svg")

println(wab_range1)


kdam_range1 = range(1.5,0,length=50)

params1 = deepcopy(params_rtc)
params1.ω_ab = 2e-3
res, combos = numerical_param(rtc_model, params1, rtc_init, species_rtc, atp_range, :atp, lam_range, :lam, kdam_range)
bistables = bistable_q(res, combos)
max_,min_ = get_max_min_yvals(atp_range, bistables)

max_ = Float64[]
min_ = Float64[]
param1_vals1 = []
for i in atp_range
    vals=[]
    param1_vals=[]
    for j in range(1, length(bistables))
        if i == bistables[j][1]
            push!(vals, bistables[j][2])
            push!(param1_vals, bistables[j][1])
        end
    end
    if length(vals) > 0 
        push!(max_, maximum(vals))
        push!(min_, minimum(vals))
        push!(param1_vals1, param1_vals[1])
    end
end

param1_vals1

plot(param1_vals1, max_; fillrange=(min_), fillalpha = 0.45,) 

bistables
max_ = Float64[]
min_ = Float64[]
for i in atp_range
    vals=[]
    for j in range(1, length(bistables))
        if i == bistables[j][1]
            push!(vals, bistables[j][2])
        end
    end
    if length(vals) > 0 
        push!(max_, maximum(vals))
        push!(min_, minimum(vals))
    end
end
max_
bistables

[i for i in atp_range]

params_bf = (L = L, c = c, kr = kr, Vmax_init = Vmax_init, Km_init = Km_init,
θtscr = θtscr, θtlr = θtlr, na = nA, nb = nB, nr = nR, d = d, 
krep = krep, ktag = ktag, atp = 500, km_a = km_a, km_b = km_b, g_max = g_max, 
kdeg = kdeg, kin = kin, ω_ab = ω_ab, ω_r = ω_r, 
kdam =  kdam, lam = 0.0147368, kc = kc, k_diss = k_diss) 

br = get_br(rtc_mod, params_bf, rtc_init, 1.5)
br.specialpoint[2].type

plot(br, vars=(:param, :rtcb))


kdam_range = range(0,1.5,length=10)
kdam_range1 = range(1.5,0,length=10)


res1_bs = numerical_bistability_analysis(rtc_model, params_rtc, rtc_init, :rtcb, species_rtc, kdam_range)
res2_bs = numerical_bistability_analysis(rtc_model, params_rtc, rtc_init, :rtcb, species_rtc, kdam_range1)

plot(kdam_range, res1_bs)
plot!(kdam_range1, res2_bs)

res1 = numerical_bistability_analysis(rtc_model, params1, rtc_init, :rtcb, species_rtc, kdam_range)
res2 = numerical_bistability_analysis(rtc_model, params1, rtc_init, :rtcb, species_rtc, kdam_range1)

plot(kdam_range, res1)
plot!(kdam_range1, res2)

diff = res1 - reverse(res2)


bs_diff = res1_bs - reverse(res2_bs)
for i in diff
    if i > 0.001
        print("bistable")
        break
    end
end



reverse(res2)
reverse(res2)*0.99 < res1 < 1.01*reverse(res2)

reverse(res2)*0.99
res1

params1 = deepcopy(params_rtc)
params1.atp = 500
params1.lam = 0.01631578947368421
res1 = numerical_bistability_analysis(rtc_model, params1, rtc_init, :rtcb, species_rtc, kdam_range)

res1[4:8]

i=5
gap_before = res1[i-1]-res1[i]
current_gap = res1[i]-res1[i+1]
current_gap > 10*gap_before

plot(kdam_range,res1)


wab_range1 = 10 .^range(log10(2e-6),log10(2e-3),length=2)

max_,min_ = numerical_triple_param_plotting(rtc_model, params_rtc, rtc_init, species_rtc, atp_range, :atp, lam_range, :lam, kdam_range, wab_range1, :ω_ab)

max_

plot(atp_range, max_[1])#; fillrange=(min_[1]), fillalpha = 0.45,) 
plot(atp_range, max_[2]; fillrange=(min_[2]), fillalpha = 0.45,) 
plot(atp_range, max_[3]; fillrange=(min_[3]), fillalpha = 0.45,) 

max_[1]
max_[2]
max_[3]


