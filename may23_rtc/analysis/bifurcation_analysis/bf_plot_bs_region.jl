using Plots, Printf, Measures, CSV, DataInterpolations
using Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra, DataFrames

include("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/bf_funcs.jl"); include("/home/holliehindley/phd/may23_rtc/functions/set_ups.jl");

atp_range = range(500, stop=5500,length=50)
kin_range = range(0,stop=0.2,length=50)
lam_range = range(0.001,stop=0.04,length=50)
wr_range = (range(0.00001,stop=0.001,length=50))
wab_range = range(0.01, stop=4, length=50)

t, atp_t, lam_t, kin_t = set_time_vars("/home/holliehindley/phd/data/atp_for_rtcmodel.csv")
plotlamt = view(lam_t, 1:200066)
plotatpt = view(atp_t, 1:200066)
plotkint = view(kin_t, 1:200066)

plot(plotatpt, plotlamt)
plot(plotatpt, plotkint)
plot(plotlamt, plotkint)

function double_param_vary(param_range1, param1, param_range2, param2, params1)
    df = DataFrame(atp = Float64[], wab = Float64[], kdam1 = Float64[], kdam2 = Float64[], bs = Symbol[])
    params = deepcopy(params1)
    for i in param_range1
        params = merge(params, (param1=>i,))
        for j in param_range2
            params = merge(params, (param2=>j,))
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
        p = plot!(plotatpt, plotlamt, c=:black)
    else
        Nothing
    end

    return (p)
end

# atp vs wab
atp_wab = plot_bistable_region(atp_range, :atp, wab_range, :ω_ab, :cyan4, "")
# atp vs wr 
atp_wr = plot_bistable_region(atp_range, :atp, wr_range, :ω_r, :cyan4,"")
# atp vs kin
atp_kin = plot_bistable_region(atp_range, :atp, kin_range, :kin, :cyan4,"")
# atp vs lam
atp_lam = plot_bistable_region(atp_range, :atp, lam_range, :lam, :cyan4,"")
# wr vs wab
wr_wab = plot_bistable_region(wr_range, :ω_r, wab_range, :ω_ab, :cyan4,"")
# wr vs kin
wr_kin = plot_bistable_region(wr_range, :ω_r, kin_range, :kin, :cyan4,"")
# wr vs lam 
wr_lam = plot_bistable_region(wr_range, :ω_r, lam_range, :lam, :cyan4,"")
# wab vs kin
wab_kin = plot_bistable_region(wab_range, :ω_ab, kin_range, :kin, :cyan4,"")
# wab vs lam 
wab_lam = plot_bistable_region(wab_range, :ω_ab, lam_range, :lam, :cyan4,"")
# kin vs lam - get instability here 
kin_lam = plot_bistable_region(kin_range, :kin, lam_range, :lam, :cyan4,"")

l = @layout [a b c; d e f; g h i j];
all_ranges = plot(kin_lam, wab_lam, wab_kin, wr_lam, wr_wab, wr_kin, atp_kin, atp_lam, atp_wab, atp_wr, layout=l, size=(1200,800), left_margin=4Plots.mm, right_margin=2Plots.mm, plot_title="Parameter spaces where bistability is present (dark area = bs)")

# savefig(all_ranges, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/all_ranges.svg")
# savefig(atp_wab, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/wab_atp_bigger_range.svg")


plots = []
for i in (10 .^ range(-4, stop=-2, length = 5))
    for j in (10 .^ range(-2, stop=0, length = 5))
        params1 = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
        θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
        krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
        kdeg = 0.001, kin = 0.022222222, ω_ab = j, ω_r = i, 
        kdam =  0.01, lam = 0.014) 	
        @show (params1[:ω_ab], params1[:ω_r])
        push!(plots, plot_bistable_region(atp_range, :atp, lam_range, :lam, "ω_ab = $(@sprintf "%.2E" j), ω_r = $(@sprintf "%.2E" i)"))
    end
end
p_all = plot(plots[1], plots[2], plots[3], plots[4], plots[5], plots[6],
plots[7], plots[8], plots[9], plots[10], plots[11], plots[12],
plots[13], plots[14], plots[15], plots[16], plots[17], plots[18],
plots[19], plots[20], plots[21], plots[22], plots[23], plots[24], plots[25],
 margin=1mm, size=(2000,2000), layout=(5,5))


# savefig(p_all, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/p_all.svg")



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


function plot_bs_region_same_plot(param_range1, param_range2, plot_var1, plot_var2, results, title, range1, param1, param2)
    colours =palette(:tab10)
    p = plot()
    for (i,j) in zip(range(1,length(results)), range(1,5))
        if length(results[i][1]) == length(atp_range)
            p = plot!(param_range1, results[i][1]; fillrange=(results[i][2]), fillalpha = 0.45, fillcolor=colours[j], 
            linecolor=colours[j], title=title, label="$(@sprintf "%g" (range1[j]))", xlims=(minimum(param_range1), maximum(param_range1)), 
            ylims=(minimum(param_range2), maximum(param_range2)), xlabel="$param1", ylabel="$param2")
            # display(p)
        else
            Nothing
        end
    end
    return p = (plot!(plot_var1, plot_var2, c=:black, label=""))
end

function get_bs_region_results_wr(param_range1, param1, param_range2, param2)
    results_wr=[]
    for i in wr_range1
        params1 = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
        θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
        krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
        kdeg = 0.001, kin = 0.022222222, ω_ab = 1, ω_r = i, 
        kdam =  0.01, lam = 0.014) 	
        @show (params1[:ω_r])
        max_,min_ = bistable_region(param_range1, param1, param_range2, param2, params1)
        push!(results_wr, (max_,min_))
    end
    return results_wr
end

function get_bs_region_results_wab(param_range1, param1, param_range2, param2)
    results_wab=[]
    for i in wab_range1
        params1 = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
        θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
        krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
        kdeg = 0.001, kin = 0.022222222, ω_ab = i, ω_r = 0.0001, 
        kdam =  0.01, lam = 0.014) 	
        @show (params1[:ω_ab])
        max_,min_ = bistable_region(param_range1, param1, param_range2, param2, params1)
        push!(results_wab, (max_,min_))
    end
    return results_wab
end

atp_range = range(400, stop=5500,length=20)
lam_range = range(0.001,stop=0.04,length=20)
kin_range = range(0.002,stop=0.2,length=20)

wab_range1 = (10 .^ range(-3, stop=1, length = 5))

results_wab = get_bs_region_results_wab(atp_range, :atp, lam_range, :lam)
p_atplam_wab = plot_bs_region_same_plot(atp_range, lam_range, plotatpt, plotlamt, results_wab, "changing ω_ab when ω_r = 0.0001", wab_range1, "ATP", "λ")
# savefig(p_atplam_wab, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/p_atplam_wab.svg")

wab_range1 = (10 .^ range(-2.2, stop=0, length = 5))
results_kinatp_wab = get_bs_region_results_wab(atp_range, :atp, kin_range, :kin)
p_atpkin_wab = plot_bs_region_same_plot(atp_range, kin_range, plotatpt, plotkint, results_kinatp_wab, "changing ω_ab when ω_r = 0.001", wab_range1, "ATP", "kin")
# savefig(p_atpkin_wab, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/p_atpkin_wab.svg")

wab_range1 = (10 .^ range(-1, stop=1, length = 5))
results_kinlam_wab = get_bs_region_results_wab(lam_range, :lam, kin_range, :kin)
p_lamkin_wab = plot_bs_region_same_plot(lam_range, kin_range, plotlamt, plotkint, results_kinlam_wab, "changing ω_ab when ω_r = 0.001", wab_range1, "λ", "kin")
# savefig(p_atpkin_wab, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/p_atpkin_wab.svg")



wr_range1 = (10 .^ range(-7, stop=-2, length = 5))
results_wr = get_bs_region_results_wr(atp_range, :atp, lam_range, :lam)
p_atlam_wr = plot_bs_region_same_plot(atp_range, lam_range, plotatpt, plotlamt, results_wr, "changing ω_r when ω_ab = 0.01", wr_range1, "ATP", "λ")
# savefig(p_atlam_wr, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/p_atlam_wr.svg")


wr_range1 = (10 .^ range(-5, stop=-3, length = 5))
results_kinatp_wr = get_bs_region_results_wr(atp_range, :atp, kin_range, :kin)
p_kinatp_wr = plot_bs_region_same_plot(atp_range, kin_range, plotatpt, plotkint, results_kinatp_wr, "changing ω_r when ω_ab = 1", wr_range1, "ATP", "kin")
savefig(p_kinatp_wr, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/p_kinatp_wr.svg")




