using Parameters, CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, Statistics
using Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra, Printf
# using Plots
using PlotlyJS, ProgressBars
include("/home/holliehindley/phd/may23_rtc/functions/solving.jl"); include("/home/holliehindley/phd/may23_rtc/functions/set_ups.jl"); include("/home/holliehindley/phd/may23_rtc/functions/plotting.jl"); 
include("/home/holliehindley/phd/may23_rtc/functions/sweep_params.jl"); include("/home/holliehindley/phd/may23_rtc/models/rtc_orig.jl"); include("/home/holliehindley/phd/may23_rtc/models/atp_lam_kin_t.jl"); 
include("/home/holliehindley/phd/may23_rtc/models/single_t.jl"); include("/home/holliehindley/phd/may23_rtc/models/combinations_t.jl"); 
include("/home/holliehindley/phd/may23_rtc/functions/bf_funcs/bf_funcs.jl");
include("/home/holliehindley/phd/colors_plotly.jl"); include("/home/holliehindley/phd/may23_rtc/models/inhibition_models/rtc_inhibition_model.jl");

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

params2 = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
kdeg = 0.001, kin = 0.022222222, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
kdam =  0.01, lam = 0.014)

params_ = @LArray [L, c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 0.5, ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)

atp_range = range(1000,stop=5000,length=3)
lam_range = range(0.01,stop=0.02,length=3)
wr_range = range(0.001,0.01,length=3)

bfs=[]; dfs=[];
for i in ProgressBar(lam_range)
    copyparams = deepcopy(params2)
    params = merge(copyparams, (:lam=>i,))
    br = get_br(rtc_mod, params, initial, 3.)
    bf = bf_point_df(br)
    df = create_br_df(br)
    push!(bfs, bf)
    push!(dfs, df)
end

wr_rtcb1, wr_rtcb2, wr_rtcb3, wr_bf_rtcb = plot_rtc_bf(bfs[1], dfs[1], findall(x->x==bfs[1].kdam[1],dfs[1].kdam)[1], findall(x->x==bfs[1].kdam[2],dfs[1].kdam)[1], :rtcb)
wr_rtcb1a, wr_rtcb2a, wr_rtcb3a, wr_bf_rtcba = plot_rtc_bf(bfs[2], dfs[2], findall(x->x==bfs[2].kdam[1],dfs[2].kdam)[1], findall(x->x==bfs[2].kdam[2],dfs[2].kdam)[1], :rtcb)
wr_rtcb1b, wr_rtcb2b, wr_rtcb3b, wr_bf_rtcbb = plot_rtc_bf(bfs[3], dfs[3], findall(x->x==bfs[3].kdam[1],dfs[3].kdam)[1], findall(x->x==bfs[3].kdam[2],dfs[3].kdam)[1], :rtcb)

wr = plot([wr_rtcb1, wr_rtcb2, wr_rtcb3, wr_bf_rtcb,wr_rtcb1a, wr_rtcb2a, wr_rtcb3a, wr_bf_rtcba,wr_rtcb1b, wr_rtcb2b, wr_rtcb3b, wr_bf_rtcbb],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="RtcB (μM)",
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=22, color="black", family="sans-serif")))


PlotlyJS.savefig(wr, "/home/holliehindley/phd/may23_rtc/paper_plots/atp_rtcb.svg")



wab_range = range(0.008,0.1,length=3)

bfs=[]; dfs=[];
for i in ProgressBar(wab_range)
    copyparams = deepcopy(params2)
    params = merge(copyparams, (:ω_ab=>i,))
    br = get_br(rtc_mod, params, initial, 3.)
    bf = bf_point_df(br)
    df = create_br_df(br)
    push!(bfs, bf)
    push!(dfs, df)
end

wab_rtcb1, wab_rtcb2, wab_rtcb3, wab_bf_rtcb = plot_rtc_bf(bfs[1], dfs[1], findall(x->x==bfs[1].kdam[1],dfs[1].kdam)[1], findall(x->x==bfs[1].kdam[2],dfs[1].kdam)[1], :rtcb)
wab_rtcb1a, wab_rtcb2a, wab_rtcb3a, wab_bf_rtcba = plot_rtc_bf(bfs[2], dfs[2], findall(x->x==bfs[2].kdam[1],dfs[2].kdam)[1], findall(x->x==bfs[2].kdam[2],dfs[2].kdam)[1], :rtcb)
wab_rtcb1b, wab_rtcb2b, wab_rtcb3b, wab_bf_rtcbb = plot_rtc_bf(bfs[3], dfs[3], findall(x->x==bfs[3].kdam[1],dfs[3].kdam)[1], findall(x->x==bfs[3].kdam[2],dfs[3].kdam)[1], :rtcb)

wab = plot([wab_rtcb1, wab_rtcb2, wab_rtcb3, wab_bf_rtcb,wab_rtcb1a, wab_rtcb2a, wab_rtcb3a, wab_bf_rtcba,scatter(x=dfs[3].kdam,y=dfs[3].rtcb,showlegend=false,line=attr(width=5, color="#005356ff"))],
Layout(xaxis_title="Damage rate (min<sup>-1</sup>)", 
yaxis_title="RtcB (μM)",
yaxis=attr(showline=true,linewidth=3,linecolor="black"),xaxis=attr(showline=true,linewidth=3,linecolor="black"),
xaxis_showgrid=false,yaxis_showgrid=false,yaxis2_showgrid=false,plot_bgcolor="white",font=attr(size=22, color="black", family="sans-serif")))


PlotlyJS.savefig(wab, "/home/holliehindley/phd/may23_rtc/paper_plots/wab_inducibility.svg")







function double_param_vary(param_range1, param1, param_range2, param2, params1)
    df = DataFrame(atp = Float64[], wab = Float64[], kdam1 = Float64[], kdam2 = Float64[], bs = Symbol[])
    params = deepcopy(params1)
    for i in param_range1
        params = merge(params, (param1=>i,))
        for j in param_range2
            params = merge(params, (param2=>j,))
            @show params[:ω_ab], params[:ω_r]
            br = get_br(rtc_mod, params, initial, 1.)
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
        params1 = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
        θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
        krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
        kdeg = 0.001, kin = 0.022222222, ω_ab = wab, ω_r = i, 
        kdam =  0.01, lam = 0.014) 	
        @show (params1[:ω_r])
        max_,min_ = bistable_region(param_range1, param1, param_range2, param2, params1)
        push!(results_wr, (max_,min_))
    end
    return results_wr
end

function get_bs_region_results_wab(param_range1, param1, param_range2, param2, wr)
    results_wab=[]
    for i in wab_range1
        params1 = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
        θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
        krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
        kdeg = 0.001, kin = 0.022222222, ω_ab = i, ω_r = wr, 
        kdam =  0.01, lam = 0.014) 	
        @show (params1[:ω_ab], params1[:atp])
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

atp_range = range(400, stop=5500,length=100)
lam_range = range(0.001,stop=0.04,length=100)

wab_range1 = (10 .^ range(-4, stop=0, length = 5))
results_wab = get_bs_region_results_wab(atp_range, :atp, lam_range, :lam, 0.01)
p_atplam_wab = plot_bs_region_same_plot(atp_range, lam_range, results_wab, "changing ω_ab when ω_r = 0.01", wab_range1, "ATP", "λ")
Plots.savefig(p_atplam_wab, "/home/holliehindley/phd/may23_rtc/paper_plots/p_atlam_wab.svg")

wr_range1 = (10 .^ range(-5, stop=-1, length = 5))
results_wr = get_bs_region_results_wr(atp_range, :atp, lam_range, :lam, 0.05623413251903491)
p_atlam_wr = plot_bs_region_same_plot(atp_range, lam_range, results_wr, "changing ω_r when ω_ab = 0.056", wr_range1, "ATP", "λ")
Plots.savefig(p_atlam_wr, "/home/holliehindley/phd/may23_rtc/paper_plots/p_atlam_wr.svg")



ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002,