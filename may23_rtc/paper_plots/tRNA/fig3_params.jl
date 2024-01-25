using Parameters, CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, Statistics
using Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra, Printf
# using Plots
using PlotlyJS, ProgressBars
include("/home/holliehindley/phd/may23_rtc/functions/solving.jl"); include("/home/holliehindley/phd/may23_rtc/functions/set_ups.jl"); include("/home/holliehindley/phd/may23_rtc/functions/plotting.jl"); 
include("/home/holliehindley/phd/may23_rtc/functions/sweep_params.jl"); include("/home/holliehindley/phd/may23_rtc/models/rtc_orig.jl"); include("/home/holliehindley/phd/may23_rtc/models/atp_lam_kin_t.jl"); 
include("/home/holliehindley/phd/may23_rtc/models/single_t.jl"); include("/home/holliehindley/phd/may23_rtc/models/combinations_t.jl"); 
include("/home/holliehindley/phd/may23_rtc/functions/bf_funcs/bf_funcs.jl");
include("/home/holliehindley/phd/colors_plotly.jl"); include("/home/holliehindley/phd/may23_rtc/models/inhibition_models/rtc_inhibition_model.jl");
include("/home/holliehindley/phd/may23_rtc/models/rtc_trna_model.jl")


@consts begin
    L = 10; #10 
    c = 0.001; 
    kr = 0.125; 
    Vmax_init = 39.51; 
    Km_init = 250; 
    θtscr = 160.01;  
    θtlr = 255.73; 
    # k_b = 17.7; 
    na = 338; 
    nb = 408; 
    nr = 532*6;
    d = 0.2; 
    krep = 137; 
    ktag = 9780;#0.1; 
    # atp = 4000;#2500; 
    km_a = 20; 
    km_b = 16;
    g_max = 2.0923; 
    gr_c = 0.0008856; # 0.000599; 
    kdeg = 0.001; 
    # kin = 0.054; #2.381 
    ω_ab = 4#4#0.093; #0.0828304057748932;#4; 
    ω_r = 0.0019*6#2e-7 #0.0019*6 #70.53; #0.0019*6#79.43865871861044; #0.0019*6;  
    ω_a = 4; 
    ω_b = 4;
    # kdam =  0.#0.000147;#0.05; 
    k = 2; # carrying capacity - changes depending on the data?
    # lam = 0.033;

    # rtca_0 = 0#0.00894; 
    # rtcb_0 = 0#0.0216; 
    # rh_0 = 11.29; #69.56; #69.4
    # rtcr_0 = 0# 0.0131 #0.04; # 8.67e-3; # change this based on keeping steady state level the whole time course (levels shouldn't really change)
    # rm_a_0 = 0; 
    # rm_b_0 = 0; 
    # rm_r_0 = 0#0.0131#0.04 # 0; 
    # rd_0 = 0; 
    # rt_0 = 0;
end

tspan = (0,1e9)
rh = 11.29 #75 # conc of ribosomes in exponential phase 
thr_t = 5#30 # was at 5 before to get saved plots # needs to be less than 30 
kin_trna = 1
kdeg_trna = 0.00001
trna_species = [:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :trna, :rd, :rt]
init_trna = [0,0,0,0,0,0,135.5,0,0] # tRNA initial conc = 135.5
params_trna = @LArray [10., c, kr*12, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 0.5, ktag, kdeg, kin_trna, 3578.9473684210525, na, nb, nr, 0.014, rh, thr_t] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam, :rh, :thr_t)

# params_trna = @LArray [10., c, kr*12, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 0.5, ktag, kdeg_trna, kin_trna, 3578.9473684210525, na, nb, nr, 0.014, rh, thr_t] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam, :rh, :thr_t)

kdam_range = range(0,400,length=1000)
kdam_range2 = range(400,0,length=1000)

params_trna2 = (L = 10., c = 0.001, kr = 0.125*12, Vmax_init = 39.51, Km_init = 250.,
θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923,
kdeg = 0.001, kin = kin_trna, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
kdam =  0.01, lam = 0.014, rh = rh, thr_t = thr_t)




# atp_range = range(1000,stop=5000,length=3)
# lam_range = range(0.01,stop=0.02,length=3)
# wr_range = range(0.001,0.01,length=3)


function trna_plotting(specie, color, br, i)
    # res_trna1 = numerical_bistability_analysis(rtc_model_trna, params_trna, init_trna, :trna, trna_species, kdam_range)
    res_trna2 = numerical_bistability_analysis(rtc_model_trna, params_trna, init_trna, specie, all_species, kdam_range2)
    # ptrna1 = scatter(x=kdam_range, y=res_trna1, name="Healthy tRNA", legendgroup=3, yaxis="y2", line=attr(color=:gold,linewidth=3))
    ptrna2 = scatter(x=kdam_range2, y=res_trna2, name="", legendgroup=i, showlegend=false, line=attr(color=color,width=5))
# res_rtcb1 = numerical_bistability_analysis(rtc_model_trna, params_trna, init_trna, :rtcb, trna_species, kdam_range)
# res_rtcb2 = numerical_bistability_analysis(rtc_model_trna, params_trna, init_trna, :rtcb, trna_species, kdam_range2)
# prtcb1 = scatter(x=kdam_range, y=res_rtcb1, name="RtcB", legendgroup=3, line=attr(color=:plum,linewidth=3))
# prtcb2 = scatter(x=kdam_range2, y=res_rtcb2, name="", showlegend=false, legendgroup=3, line=attr(color=:plum,linewidth=3))
    df = create_br_df(br)
    # df1 = create_br_df(br1)
    bf = bf_point_df(br)
    # bf1 = bf_point_df(br1)
    kdam1 = findall(x->x==bf.kdam[1],df.kdam)[1]
    # kdam2 = findall(x->x==bf1.kdam[1],df1.kdam)[1]
    stable_trna = df[!,specie][1:kdam1]
    unstable_trna = df[!,specie][kdam1:end]
    stable1 = scatter(x=df.kdam[1:kdam1],y=stable_trna, line=attr(color=color, width=5), name="$specie", legendgroup=i)
    unstable1 = scatter(x=df.kdam[kdam1:end],y=unstable_trna, line=attr(color=color, width=5, dash="dash"), name="", showlegend=false, legendgroup=i)
    return ptrna2, stable1, unstable1
end


lam_range = range(0.01,stop=0.02,length=3)

brs=[]; dfs=[]; res1=[]; res2=[];
for i in ProgressBar(lam_range)
    copyparams = deepcopy(params_trna2)
    params = merge(copyparams, (:lam=>i,))
    br = get_br(rtc_mod_trna, params, init_trna, 400.)
    # bf = bf_point_df(br)
    df = create_br_df(br)
    params_trna1 = deepcopy(params_trna)
    params_trna1.lam = i
    res_rtcb1 = numerical_bistability_analysis(rtc_model_trna, params_trna1, init_trna, :rtcb, trna_species, kdam_range)
    res_rtcb2 = numerical_bistability_analysis(rtc_model_trna, params_trna1, init_trna, :rtcb, trna_species, kdam_range2)
    push!(brs, br)
    push!(dfs, df)
    push!(res1, res_rtcb1)
    push!(res2, res_rtcb2)
end

prtcb1, stable_rtcb1, unstable_rtcb1 = trna_plotting(:rtcb, "#4ca7a2ff", brs[2], 2)
prtcb2, stable_rtcb2, unstable_rtcb2 = trna_plotting(:rtcb, "#e48080ff", brs[3], 3)
prtcb1a = plot(scatter(x=kdam_range, y=res1[1], name="RtcB", legendgroup=3, line=attr(color=:plum,width=5)))

lam = plot([prtcb1a, prtcb1, stable_rtcb1, unstable_rtcb1, prtcb2, stable_rtcb2, unstable_rtcb2],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="RtcB (μM)",
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=22, color="black", family="sans-serif")))

savefig(lam, "/home/holliehindley/phd/may23_rtc/paper_plots/tRNA/lam_rtcb.svg")


atp_range = range(1000,stop=5000,length=3)

brs_atp=[]; dfs_atp=[]; res1_atp=[]; res2_atp=[];
for i in ProgressBar(atp_range)
    copyparams = deepcopy(params_trna2)
    params = merge(copyparams, (:atp=>i,))
    br = get_br(rtc_mod_trna, params, init_trna, 400.)
    df = create_br_df(br)
    params_trna1 = deepcopy(params_trna)
    params_trna1.lam = i
    res_rtcb1 = numerical_bistability_analysis(rtc_model_trna, params_trna1, init_trna, :rtcb, trna_species, kdam_range)
    res_rtcb2 = numerical_bistability_analysis(rtc_model_trna, params_trna1, init_trna, :rtcb, trna_species, kdam_range2)
    push!(brs_atp, br)
    push!(dfs_atp, df)
    push!(res1_atp, res_rtcb1)
    push!(res2_atp, res_rtcb2)
end

prtcb_atp, stable_rtcb_atp, unstable_rtcb_atp = trna_plotting(:rtcb, "#b693ccff", brs_atp[1], 1)
prtcb1_atp, stable_rtcb1_atp, unstable_rtcb1_atp = trna_plotting(:rtcb, "#4ca7a2ff", brs_atp[2], 2)
prtcb2_atp, stable_rtcb2_atp, unstable_rtcb2_atp = trna_plotting(:rtcb, "#e48080ff", brs_atp[3], 3)

atp = plot([prtcb_atp, stable_rtcb_atp, unstable_rtcb_atp, prtcb1_atp, stable_rtcb1_atp, unstable_rtcb1_atp, prtcb2_atp, stable_rtcb2_atp, unstable_rtcb2_atp],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="RtcB (μM)",
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=22, color="black", family="sans-serif")))

savefig(atp, "/home/holliehindley/phd/may23_rtc/paper_plots/tRNA/atp_rtcb.svg")


wab_range = range(0.008,0.1,length=3)

brs_wab=[]; dfs_wab=[]; res1_wab=[]; res2_wab=[];
for i in ProgressBar(wab_range)
    copyparams = deepcopy(params_trna2)
    params = merge(copyparams, (:ω_ab=>i,))
    br = get_br(rtc_mod_trna, params, init_trna, 400.)
    df = create_br_df(br)
    params_trna1 = deepcopy(params_trna)
    params_trna1.lam = i
    res_rtcb1 = numerical_bistability_analysis(rtc_model_trna, params_trna1, init_trna, :rtcb, trna_species, kdam_range)
    res_rtcb2 = numerical_bistability_analysis(rtc_model_trna, params_trna1, init_trna, :rtcb, trna_species, kdam_range2)
    push!(brs_wab, br)
    push!(dfs_wab, df)
    push!(res1_wab, res_rtcb1)
    push!(res2_wab, res_rtcb2)
end

prtcb_wab, stable_rtcb_wab, unstable_rtcb_wab = trna_plotting(:rtcb, "#b693ccff", brs_wab[1], 1)
prtcb1_wab, stable_rtcb1_wab, unstable_rtcb1_wab = trna_plotting(:rtcb, "#4ca7a2ff", brs_wab[2], 2)

prtcb1a_wab = scatter(x=kdam_range, y=res1_wab[1], name="RtcB", legendgroup=3, line=attr(color=:plum,width=5))

wab = plot([prtcb_wab, stable_rtcb_wab, unstable_rtcb_wab, prtcb1_wab, stable_rtcb1_wab, unstable_rtcb1_wab, prtcb1a_wab],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="RtcB (μM)",
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=22, color="black", family="sans-serif")))

savefig(wab, "/home/holliehindley/phd/may23_rtc/paper_plots/tRNA/wab_rtcb.svg")


wr_range = range(0.001,0.01,length=3)

brs_wr=[];
for i in ProgressBar(wr_range)
    copyparams = deepcopy(params_trna2)
    params = merge(copyparams, (:ω_r=>i,))
    br = get_br(rtc_mod_trna, params, init_trna, 400.)
    # bf = bf_point_df(br)
    # df = create_br_df(br)
    push!(brs_wr, br)
end

prtcb_wr, stable_rtcb_wr, unstable_rtcb_wr = trna_plotting(:rtcb, "#b693ccff", brs_wr[1], 1)
prtcb1_wr, stable_rtcb1_wr, unstable_rtcb1_wr = trna_plotting(:rtcb, "#4ca7a2ff", brs_wr[2], 2)
prtcb2_wr, stable_rtcb2_wr, unstable_rtcb2_wr = trna_plotting(:rtcb, "#e48080ff", brs_wr[3], 3)

wr = plot([prtcb_wr, stable_rtcb_wr, unstable_rtcb_wr, prtcb1_wr, stable_rtcb1_wr, unstable_rtcb1_wr, prtcb2_wr, stable_rtcb2_wr, unstable_rtcb2_wr],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="RtcB (μM)",
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=22, color="black", family="sans-serif")))


savefig(wr, "/home/holliehindley/phd/may23_rtc/paper_plots/tRNA/wr_rtcb.svg")





# creating the banana plots 

function double_param_vary(param_range1, param1, param_range2, param2, params1)
    df = DataFrame(atp = Float64[], wab = Float64[], kdam1 = Float64[], kdam2 = Float64[], bs = Symbol[])
    params = deepcopy(params1)
    for i in param_range1
        params = merge(params, (param1=>i,))
        for j in param_range2
            params = merge(params, (param2=>j,))
            @show params[:ω_ab], params[:ω_r]
            br = get_br(rtc_mod_trna, params, init_trna, 400.)
            if length(br.specialpoint) == 2
                push!(df, (i, j, br.specialpoint[1].param, br.specialpoint[2].param, br.specialpoint[1].type))
            else
                push!(df, (i, j, br.specialpoint[2].param, br.specialpoint[3].param, br.specialpoint[2].type))
            end
        end    
    end
    return df
end

function plot_bistable_region(param_range1, param1, param_range2, param2, colour, title)
    df = double_param_vary(param_range1, param1, param_range2, param2, params1)
    bsp = df[df.bs .== :bp, :]
    # @show maximum(bsp[bsp[:,1] .== i, :][:,2])
    max_ = Float64[]
    min_ = Float64[]
    for i in param_range1
        if bsp[bsp[:,1] .== i, :][:,2] == Float64[]
            @show i
        else
            push!(max_, maximum(bsp[bsp[:,1] .== i, :][:,2]))
            push!(min_, minimum(bsp[bsp[:,1] .== i, :][:,2]))
        end
    end
    if length(max_) == length(param_range1)
        p = plot(param_range1, max_; fillrange=(min_), fillalpha = 0.35, fillcolor=colour, linecolor=:white, title=title, legend=false, xlims=(minimum(param_range1), maximum(param_range1)), ylims=(minimum(param_range2), maximum(param_range2)))
        # p = plot!(plotatpt, plotlamt, c=:black)
    else
        Nothing
    end

    return (p)
end







function bistable_region(param_range1, param1, param_range2, param2, params1)
    df = double_param_vary(param_range1, param1, param_range2, param2, params1)
    bsp = df[df.bs .== :bp, :]
    # @show maximum(bsp[bsp[:,1] .== i, :][:,2])
    max_ = Float64[]
    min_ = Float64[]
    for i in param_range1
        if bsp[bsp[:,1] .== i, :][:,2] == Float64[]
            @show i
        else
            push!(max_, maximum(bsp[bsp[:,1] .== i, :][:,2]))
            push!(min_, minimum(bsp[bsp[:,1] .== i, :][:,2]))
        end
    end
    return max_, min_
end
function get_bs_region_results_wr(param_range1, param1, param_range2, param2, wab)
    results_wr=[]
    for i in wr_range1
        
        params1 = (L = 10., c = 0.001, kr = 0.125*12, Vmax_init = 39.51, Km_init = 250.,
        θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
        krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923,
        kdeg = 0.001, kin = kin_trna, ω_ab = wab, ω_r = i, 
        kdam =  0.01, lam = 0.014, rh = rh, thr_t = thr_t)


        # @show (params_trna2[:ω_r])
        max_,min_ = bistable_region(param_range1, param1, param_range2, param2, params1)
        push!(results_wr, (max_,min_))
    end
    return results_wr
end

function get_bs_region_results_wab(param_range1, param1, param_range2, param2, wr)
    results_wab=[]
    for i in wab_range1
        
        params1 = (L = 10., c = 0.001, kr = 0.125*12, Vmax_init = 39.51, Km_init = 250.,
        θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
        krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923,
        kdeg = 0.001, kin = kin_trna, ω_ab = i, ω_r = wr, 
        kdam =  0.01, lam = 0.014, rh = rh, thr_t = thr_t)

	
        # @show (params1[:ω_ab], params1[:atp])
        max_,min_ = bistable_region(param_range1, param1, param_range2, param2, params1)
        push!(results_wab, (max_,min_))
    end
    return results_wab
end
using Plots
function plot_bs_region_same_plot(param_range1, param_range2, results, title, range1, param1, param2)
    colours =palette(:tab10)
    p = Plots.plot()
    for (i,j) in zip(range(1,length(results)), range(1,5))
        if length(results[i][1]) == length(param_range1)
            p = Plots.plot!(param_range1, results[i][1]; fillrange=(results[i][2]), fillalpha = 0.45, fillcolor=colours[j], 
            linecolor=colours[j], title=title, label="$(@sprintf "%g" (range1[j]))", xlims=(minimum(param_range1), maximum(param_range1)), 
            ylims=(minimum(param_range2), maximum(param_range2)), xlabel="$param1", ylabel="$param2", tickfontsize=15,
            guidefontsize=18,guidefont="sans-serif",)
            # display(p)
        else
            Nothing
        end
    end
    return p 
end

atp_range = range(400, stop=5500,length=50)
lam_range = range(0.001,stop=0.04,length=50)
# lam_range = range(0.01,stop=0.2,length=50)

wab_range1 = (10 .^ range(-4, stop=0, length = 5)) #0.05 normally

results_wab = get_bs_region_results_wab(atp_range, :atp, lam_range, :lam, 0.01)
p_atplam_wab = plot_bs_region_same_plot(atp_range, lam_range, results_wab, "changing ω_ab when ω_r = 0.01", wab_range1, "ATP", "λ")
Plots.savefig(p_atplam_wab, "/home/holliehindley/phd/may23_rtc/paper_plots/p_atlam_wab.svg")

wr_range1 = (10 .^ range(-5, stop=-1, length = 5))
results_wr = get_bs_region_results_wr(atp_range, :atp, lam_range, :lam, 0.05623413251903491)
p_atlam_wr = plot_bs_region_same_plot(atp_range, lam_range, results_wr, "changing ω_r when ω_ab = 0.056", wr_range1, "ATP", "λ")
Plots.savefig(p_atlam_wr, "/home/holliehindley/phd/may23_rtc/paper_plots/p_atlam_wr.svg")




# numerical parameter analysis because the bfkit is really sensitive and keeps giving errors and warnings, even though we get the same results from both analyses - used as a check but can use as normal too 
function numerical_param(param_range1, param1, param_range2, param2, params_trna1, kdam_range)
    # params_trna1 = deepcopy(params_trna)
    res1=[]
    param_combos=[]
    for i in ProgressBar(param_range1)
        # res1=[]
        params_trna1[param1] = i
        for j in param_range2
            params_trna1[param2] = j
            res_rtcb1 = numerical_bistability_analysis(rtc_model_trna, params_trna1, init_trna, :rtcb, trna_species, kdam_range)
            push!(res1, res_rtcb1)
            push!(param_combos, (i,j))
        end   
        # push!(res, res1) 
    end
    return res1, param_combos
end



# plot([scatter(x=kdam_range, y=res[i], name="$(param_combos[i])") for i in range(1,length(res))])

function bistable_q(res, param_combos)
    bistables=[]; 
    for j in range(1, length(res))
        for i in range(1,length(res[j])-1)
            if res[j][i+1] < 0.0001 * res[j][i]
                # push!(bistables, ((res[j][i+1], res[j][i])))
                push!(bistables, param_combos[j])
                break 
            # else
            #     print("monostable")
            end
        end
    end

    return bistables
end


# plot([scatter(x=[bistables[i][1]],y=[bistables[i][2]], mode="markers") for i in range(1, length(bistables))], Layout(xaxis_range=(400,5500), yaxis_range=(0.001,0.04)))


function get_max_min_yvals(param_range1, bistables)
    max_ = Float64[]
    min_ = Float64[]
    for i in param_range1
        vals=[]
        for j in range(1, length(bistables))
            if i == bistables[j][1]
                push!(vals, bistables[j][2])
            end
        end
        push!(max_, maximum(vals))
        push!(min_, minimum(vals))
    end
    return max_,min_
end


function numerical_triple_param_plotting(param_range1, param1, param_range2, param2, params_trna, kdam_range, w_range, w_change)
    p = Plots.plot()
    params_trna1 = deepcopy(params_trna)
    for i in w_range
        params_trna1[w_change] = i
        # @show params_trna1[w_change], params_trna1[param1], params_trna1[param2]
	    res, param_combos = numerical_param(param_range1, param1, param_range2, param2, params_trna1, kdam_range)
        bistables = bistable_q(res, param_combos)
        max_,min_ = get_max_min_yvals(param_range1, bistables)

        p = Plots.plot!(param_range1, max_; fillrange=(min_), fillalpha = 0.45, 
        label="$(@sprintf "%g" i)", xlims=(minimum(param_range1), maximum(param_range1)), 
        ylims=(minimum(param_range2), maximum(param_range2)), xlabel="$param1", ylabel="$param2", tickfontsize=15,
        guidefontsize=18,guidefont="sans-serif",)
    end

    
    return p 
end

using Plots 
kdam_range = range(0,400,length=50)
atp_range = range(400, stop=5500,length=100)
lam_range = range(0.001,stop=0.04,length=100)

wab_range1 = (10 .^ range(-4, stop=-1, length = 5)) #0.05 normally

p_atplam_wab = numerical_triple_param_plotting(atp_range, :atp, lam_range, :lam, params_trna, kdam_range, wab_range1, :ω_ab)
Plots.savefig(p_atplam_wab, "/home/holliehindley/phd/may23_rtc/paper_plots/tRNA/p_atlam_wab.svg")

wr_range1 = (10 .^ range(-5, stop=-2, length = 5))
p_atplam_wr = numerical_triple_param_plotting(atp_range, :atp, lam_range, :lam, params_trna, kdam_range, wr_range1, :ω_r)
Plots.savefig(p_atplam_wr, "/home/holliehindley/phd/may23_rtc/paper_plots/tRNA/p_atlam_wr.svg")

