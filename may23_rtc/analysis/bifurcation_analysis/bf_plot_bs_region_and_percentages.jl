using Plots, Printf, Measures, CSV, DataInterpolations
using Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra, DataFrames

include("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/bf_funcs.jl"); include("/home/holliehindley/phd/may23_rtc/functions/set_ups.jl");

atp_range = range(500, stop=5500,length=50)
kin_range = range(0,stop=0.2,length=50)
lam_range = range(0.001,stop=0.04,length=50)
wr_range = (range(0.00001,stop=0.001,length=50))
wab_range = range(0.01, stop=4, length=50)

# load time varying parameters and create object so they can be plotted 
t, atp_t, lam_t, kin_t = set_time_vars("/home/holliehindley/phd/data/atp_for_rtcmodel.csv")
plotlamt = view(lam_t, 1:200066)
plotatpt = view(atp_t, 1:200066)
plotkint = view(kin_t, 1:200066)
plotatpt = [(plotatpt...)...]
plotlamt = [(plotlamt...)...]
plot(plotatpt, plotlamt, plotkint)

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

#individual plots
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
        if length(results[i][1]) == length(param_range1)
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

function get_bs_region_results_both(param_range1, param1, param_range2, param2, wr_range, wab_range)
    results = []
    for (i,j) in zip(wr_range, wab_range)#(10 .^ range(-4, stop=-2, length = 5))
        # for j in wab_range#(10 .^ range(-2, stop=0, length = 5))
        params1 = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
        θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
        krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
        kdeg = 0.001, kin = 0.022222222, ω_ab = j, ω_r = i, 
        kdam =  0.01, lam = 0.014) 	
        @show (params1[:ω_ab], params1[:ω_r])
        max_,min_ = bistable_region(param_range1, param1, param_range2, param2, params1)
        push!(results, (max_,min_))
        # push!(results, plot_bistable_region(atp_range, :atp, lam_range, :lam, :red, "ω_ab = $(@sprintf "%.2E" j), ω_r = $(@sprintf "%.2E" i)"))
        # end
    end
    return results
end

# plots for each combination of time varying parameters 
atp_range = range(400, stop=5500,length=100)
lam_range = range(0.001,stop=0.04,length=100)
kin_range = range(0.003,stop=0.2,length=100)

wab_range1 = (10 .^ range(-3, stop=1, length = 5))
results_wab = get_bs_region_results_wab(atp_range, :atp, lam_range, :lam, 0.0001)
p_atplam_wab = plot_bs_region_same_plot(atp_range, lam_range, plotatpt, plotlamt, results_wab, "changing ω_ab when ω_r = 0.0001", wab_range1, "ATP", "λ")
# savefig(p_atplam_wab, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/p_atplam_wab.svg")


wab_range1 = (10 .^ range(-2.2, stop=0, length = 5))
results_kinatp_wab = get_bs_region_results_wab(atp_range, :atp, kin_range, :kin, 0.001)
p_atpkin_wab = plot_bs_region_same_plot(atp_range, kin_range, plotatpt, plotkint, results_kinatp_wab, "changing ω_ab when ω_r = 0.001", wab_range1, "ATP", "kin")
# savefig(p_atpkin_wab, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/p_atpkin_wab.svg")


lam_range = range(0.001,stop=0.04,length=20)
kin_range = range(0.003,stop=0.1,length=20)
wab_range1 = (10 .^ range(-2.05, stop=0, length = 5))
results_kinlam_wab = get_bs_region_results_wab(kin_range, :kin, lam_range, :lam, 0.00001)
p_lamkin_wab = plot_bs_region_same_plot(kin_range, lam_range, plotkint, plotlamt, results_kinlam_wab, "changing ω_ab when ω_r = 0.00001", wab_range1, "kin", "λ")
# savefig(p_lamkin_wab, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/p_lamkin_wab.svg")

plot(plotkint, plotlamt)

wr_range1 = (10 .^ range(-7, stop=-2, length = 5))
results_wr = get_bs_region_results_wr(atp_range, :atp, lam_range, :lam, 1)
p_atlam_wr = plot_bs_region_same_plot(atp_range, lam_range, plotatpt, plotlamt, results_wr, "changing ω_r when ω_ab = 1", wr_range1, "ATP", "λ")
# savefig(p_atlam_wr, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/p_atlam_wr.svg")


wr_range1 = (10 .^ range(-5, stop=-3, length = 5))
results_kinatp_wr = get_bs_region_results_wr(atp_range, :atp, kin_range, :kin, 1)
p_kinatp_wr = plot_bs_region_same_plot(atp_range, kin_range, plotatpt, plotkint, results_kinatp_wr, "changing ω_r when ω_ab = 1", wr_range1, "ATP", "kin")
# savefig(p_kinatp_wr, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/p_kinatp_wr.svg")



lam_range = range(0.001,stop=0.04,length=100)
kin_range = range(0.003,stop=0.1,length=100)
wr_range1 = (10 .^ range(-6.3, stop=-3, length = 5))
results_kinlam_wr = get_bs_region_results_wr(kin_range, :kin, lam_range, :lam, 0.1)
p_kinlam_wr = plot_bs_region_same_plot(kin_range, lam_range, plotkint, plotlamt, results_kinlam_wr, "changing ω_r when ω_ab = 0.1", wr_range1, "kin", "λ")
# savefig(p_kinlam_wr, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/p_kinlam_wr.svg")


# chaning both wab and wr at the same time to create a big dataframe for final plot
function plot_bs_region_same_plot(param_range1, param_range2, plot_var1, plot_var2, results, title, range1, range2, param1, param2)
    colours =palette(:tab20)
    p = plot()
    for (i,j) in zip(range(1,length(results)), range(1,20))
        if length(results[i][1]) == length(param_range1)
            p = plot!(param_range1, results[i][1]; fillrange=(results[i][2]), fillalpha = 0.45, fillcolor=colours[j], 
            linecolor=colours[j], title=title, label="wr = $(@sprintf "%g" (range1[j])), wab = $(@sprintf "%g" (range2[j]))", xlims=(minimum(param_range1), maximum(param_range1)), 
            ylims=(minimum(param_range2), maximum(param_range2)), xlabel="$param1", ylabel="$param2")
            # display(p)
        else
            Nothing
        end
    end
    return p = (plot!(plot_var1, plot_var2, c=:black, label=""))
end
function get_bs_region_results_both(param_range1, param1, param_range2, param2, wr_range, wab_range)
    results = []
    for i in wr_range#(10 .^ range(-4, stop=-2, length = 5))
        for j in wab_range#(10 .^ range(-2, stop=0, length = 5))
        params1 = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
        θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
        krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
        kdeg = 0.001, kin = 0.022222222, ω_ab = j, ω_r = i, 
        kdam =  0.01, lam = 0.014) 	
        @show (params1[:ω_ab], params1[:ω_r])
        max_,min_ = bistable_region(param_range1, param1, param_range2, param2, params1)
        push!(results, (max_,min_))
        # push!(results, plot_bistable_region(atp_range, :atp, lam_range, :lam, :red, "ω_ab = $(@sprintf "%.2E" j), ω_r = $(@sprintf "%.2E" i)"))
        end
    end
    return results
end

wr_range1 = (10 .^ range(-5, stop=-3, length = 20))
wab_range1 = (10 .^ range(-2, stop=1, length = 20))
results_both_atp_lam = get_bs_region_results_both(atp_range, :atp, lam_range, :lam, wr_range1, wab_range1)
plot_bs_region_same_plot(atp_range, lam_range, plotatpt, plotlamt, results_both_atp_lam, "changing both wr and wab", wr_range1, wab_range1, "ATP", "λ")


# working out the percentage area of the time varying parameter curve that is covered by the region of bistability for each parameter combination 
using PyCall
np = pyimport("numpy")
interp = pyimport("scipy.interpolate")

# splits the parameter curve into two and sorts into a left and right branch, interpolates to get the length of range1/2 number of points on each branch so a percentage can be calculated thereafter 
function param_curve(x, y, x_range) # need to change the ranges depeding on the number in atp_range when lenght=100 then [13:end-4] and [1:end-4] but when length=50 then [7:end-2] and [1:end-2]
    y1 = y[1:8]
    x1 = x[1:8]

    y2 = y[8:19]
    x2 = x[8:19]

    d1p = interp.interp1d(x1,y1)
    d2p = interp.interp1d(x2,y2)

    # test1 = d1p(x1)
    # test2 = d2p(x2)

    new_range = [(x_range...)...]
    range1 = new_range[13:end-4]
    range2 = new_range[1:end-4]

    rpoints=[]
    for i in range1
        push!(rpoints, d1p(i))
    end 

    lpoints=[]
    for i in range2
        push!(lpoints, d2p(i))
    end 

    rpoints = [(rpoints...)...]
    lpoints = [(lpoints...)...]
    return rpoints, lpoints
end
rpoints, lpoints = param_curve(plotatpt, plotlamt, atp_range)

# works out if the parameter curve is in the bistability area or not by considering the upper and lower bounds of the bistability region 
function inarea(points, upper, lower)
    inarea_y1 = []
    inarea_x1 = []
    for (i, up, low) in zip(points,upper,lower)
        if i<up 
            if i>low
                push!(inarea_y1, i)
                push!(inarea_x1, findfirst(isequal(i), points))
            end
        end
    end

    # xvals1=[]
    # for i in inarea_x1
    #     push!(xvals1, range1[i])
    # end
    return inarea_y1
end

# out of those points that are in the region works out the overall percentage of points that lie in the region 
function find_perc(rpoints, r_upper, r_lower, lpoints, l_upper, l_lower)
    inarea_y1 = inarea(rpoints, r_upper, r_lower)
    inarea_y2 = inarea(lpoints, l_upper, l_lower)

    tot = length(vcat(rpoints,lpoints))

    inside = length(vcat(inarea_y1,inarea_y2))

    perc = (inside/tot)*100
    return perc
end

upper = results_both_atp_lam[5][1]
lower = results_both_atp_lam[5][2]
r_upper = upper[13:end-4]
r_lower = lower[13:end-4]
l_upper = upper[1:end-4]
l_lower = lower[1:end-4]
find_perc(rpoints, r_upper, r_lower, lpoints, l_upper, l_lower)

# can use this to check that the points are correctly identified as inside the bistability area 
function checking_plot(x_range, rpoints, r_upper, r_lower, lpoints, l_upper, l_lower, upper, lower, x_plot, y_plot)

    new_range = [(x_range...)...]
    range1 = new_range[13:end-4]
    range2 = new_range[1:end-4]

    inarea_y1 = []
    inarea_x1 = []
    for (i, up, low) in zip(rpoints,r_upper,r_lower)
        if i<up 
            if i>low
                push!(inarea_y1, i)
                push!(inarea_x1, findfirst(isequal(i), rpoints))
            end
        end
    end

    xvals1=[]
    for i in inarea_x1
        push!(xvals1, range1[i])
    end

    inarea_y2 = []
    inarea_x2 = []
    for (i, up, low) in zip(lpoints,l_upper,l_lower)
        if i<up 
            if i>low
                push!(inarea_y2, i)
                push!(inarea_x2, findfirst(isequal(i), lpoints))
            end
        end
    end

    xvals2=[]
    for i in inarea_x2
        push!(xvals2, range2[i])
    end

    p = plot(atp_range, upper)
    p = plot!(atp_range, lower)
    p = plot!(x_plot, y_plot)
    p = scatter!(xvals1, inarea_y1)
    p = scatter!(xvals2, inarea_y2)
    return p
end

checking_plot(atp_range, rpoints, r_upper, r_lower, lpoints, l_upper, l_lower, upper, lower, plotatpt, plotlamt)





# creates array of each parameter that have been used in the bistability searches 
wabs=[]
wrs=[]
for i in wr_range1#(10 .^ range(-4, stop=-2, length = 5))
    for j in wab_range1
        push!(wrs, i)
        push!(wabs, j)
    end
end

# gets the percentage for each of the bistability regions 
percs=[]
for i in range(1, (length(wab_range1)*length(wr_range1)))
    upper = results_both_atp_lam[i][1]
    lower = results_both_atp_lam[i][2]
    r_upper = upper[13:end-4]
    r_lower = lower[13:end-4]
    l_upper = upper[1:end-4]
    l_lower = lower[1:end-4]
    push!(percs, find_perc(rpoints, r_upper, r_lower, lpoints, l_upper, l_lower))
end

# creates dataframe with all the data 
perc_df = DataFrame(wab=wabs, wr=wrs, percs=percs)
# CSV.write("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/PERCENTAGES.csv", perc_df)




# loads the dataframe 
perc_df = DataFrame(CSV.File("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/PERCENTAGES.csv"))
# reshapes the dataframe so the percentage column represents a 20x20 matrix 
plot_percs = reshape(perc_df.percs, (20,20))

# contour plots 
using Plots; pythonplot()

contour_logPlot = contour(wab_range1, wr_range1, plot_percs, scale=:log10, color=:plasma, fill=true, levels=20, xlabel="ω_ab range", ylabel="ω_r range", colorbar_title="%", title="% area of bistability region that crosses ATP/λ curve", titlefontsize=11)

contourPlot = contour(wab_range1, wr_range1, plot_percs, color=:plasma, fill=true, levels=20, xlabel="ω_ab range", ylabel="ω_r range", colorbar_title="%", title="% area of bistability region that crosses ATP/λ curve", titlefontsize=11)

savefig(contourPlot, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/contourPlot.svg")
savefig(contour_logPlot, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/contour_logPlot.svg")









































# sort!(perc_df, [order(:percs)])

# perc_df


# zeros = perc_df[perc_df.percs .== 0.,:]

# tens = perc_df[ (perc_df.percs .< 20),:]

# twentys = perc_df[(perc_df.percs .> 20) .& (perc_df.percs .< 30),:]

# thirtys = perc_df[(perc_df.percs .> 30) .& (perc_df.percs .< 40),:]

# fortys = perc_df[(perc_df.percs .> 40) .& (perc_df.percs .< 50),:]

# fiftys = perc_df[(perc_df.percs .> 50),:]

# # plot(xlims=(minimum(wab_range1), maximum(wab_range1)), ylims=(minimum(wr_range1), maximum(wr_range1)), xlabel="ω_ab", ylabel="ω_r", xaxis=:log, yaxis=:log)
# # scatter!(zeros.wab, zeros.wr, label="0%")#, fillrange=(zeros.wr))
# # scatter!(tens.wab, tens.wr, label="10-20%")#, fillrange=(tens.wr))
# # scatter!(twentys.wab, twentys.wr, label="20-30%")
# # scatter!(thirtys.wab, thirtys.wr, label="30-40%")
# # scatter!(fortys.wab, fortys.wr, label="40-50%")
# # scatter!(fiftys.wab, fiftys.wr, label="50+%")

# function get_boundaries(df)
#     min=[]
#     max=[]
#     x=[]
#     for i in wab_range1
#         if sizeof(df[df.wab .== i,:].wr) != 0
#             push!(min, minimum(df[df.wab .== i,:].wr))
#             push!(max, maximum(df[df.wab .== i,:].wr))  
#             push!(x, df[df.wab .== i,:].wab)     
#         end
#     end

#     x = ([i[1] for i in x])
#     min = convert(Array{Float64,1},min)
#     max = convert(Array{Float64,1},max)
#     return x, min, max
# end



# zero_x, zero_min, zero_max = get_boundaries(zeros)
# tens_x, tens_min, tens_max = get_boundaries(tens)
# twentys_x, twentys_min, twentys_max = get_boundaries(twentys)
# thirtys_x, thirtys_min, thirtys_max = get_boundaries(thirtys)
# fortys_x, fortys_min, fortys_max = get_boundaries(fortys)
# fiftys_x, fiftys_min, fiftys_max = get_boundaries(fiftys)


# #tens
# tens_x1 = tens_x[5:end]; tens_min1 = tens_min[5:end]; tens_max1 = tens_max[5:end];
# tens_x2 = tens_x[1:4]; tens_min2 = tens_min[1:4]; tens_max2 = tens_max[1:4];

# #twentys
# twentys_x_max2 = [twentys_x[1:9];twentys.wab[69];twentys_x[11]]; twentys_max2 = [twentys_max[1:9];twentys.wr[69];twentys_min[11]]
# twentys_x_min2 = twentys_x[1:11]; twentys_min2 = twentys_min[1:11];
# twentys_x_max1 = twentys_x[10:end]; twentys_max1 = twentys_max[10:end];
# twentys_x_min1 = [twentys_x[10];twentys.wab[93];twentys_x[12:end]]; twentys_min1 = [twentys_max[10];twentys.wr[93];twentys_min[12:end]];

# #thirtys
# thirtys_min2 = thirtys_min[1:end-3]
# thirtys_x_min2 = thirtys_x[1:end-3]
# thirtys_max1 = thirtys_max[9:end]
# thirtys_x_max1 = thirtys_x[9:end]
# thirtys_min1 = [thirtys.wr[140]; thirtys.wr[119]; thirtys.wr[138]; thirtys.wr[116]; thirtys.wr[136]; thirtys.wr[113]; thirtys.wr[134]; thirtys.wr[110]; thirtys.wr[132]; thirtys_min[end-2:end]]
# thirtys_x_min1 = [[thirtys.wab[140]]; thirtys.wab[119]; thirtys.wab[138]; thirtys.wab[116]; thirtys.wab[136]; thirtys.wab[113]; thirtys.wab[134]; thirtys.wab[110]; thirtys.wab[132]; thirtys_x[end-2:end]]
# thirtys_max2 = [thirtys_max[1:8]; thirtys.wr[135]; thirtys.wr[124]; thirtys.wr[133]; thirtys.wr[123]; thirtys.wr[131]; thirtys.wr[122]; thirtys.wr[129]; thirtys.wr[121]; thirtys_min[end-3]]
# thirtys_x_max2 = [[thirtys_x[1:8]; thirtys.wab[135]]; thirtys.wab[124]; thirtys.wab[133]; thirtys.wab[123]; thirtys.wab[131]; thirtys.wab[122]; thirtys.wab[129]; thirtys.wab[121]; thirtys_x[end-3]]

# #fortys
# fortys_x_min2 = fortys_x[1:end-3]
# fortys_min2 = fortys_min[1:end-3]
# fortys_x_max1 = [fortys_x[3:end];fortys.wab[27]]
# fortys_max1 = [fortys_max[3:end];fortys.wr[27]]
# fortys_x_min1 = [(reverse(fortys_x[end-2:end]));fortys.wab[35];fortys.wab[42];fortys.wab[36];fortys.wab[43];fortys.wab[37];fortys.wab[44];fortys.wab[38];fortys.wab[45];fortys.wab[39];fortys.wab[46];fortys.wab[40]]
# fortys_min1 = [(reverse(fortys_min[end-2:end]));fortys.wr[35];fortys.wr[42];fortys.wr[36];fortys.wr[43];fortys.wr[37];fortys.wr[44];fortys.wr[38];fortys.wr[45];fortys.wr[39];fortys.wr[46];fortys.wr[40]]


# scatter(fortys.wab, fortys.wr, color=:black, label="", axis=:log)
# scatter!(fortys.wab[1:40], fortys.wr[1:40])
# # plot!(fortys_x[1:end-3], fortys_min[1:end-3])
# plot!(fortys_x, fortys_min1, fillrange=(fortys_max1))
# plot(fortys_x[3:end], fortys_max[3:end])#, axis=:log)
# plot!(fortys_x_min1, fortys_min1)
# plot!(fortys_x[3:end], fortys_max[3:end], fillrange=(fortys_min1))

# fortys_min1
# fortys_x
# fortys_x[3:end]
# fortys_x_min1
# scatter!(fortys_x[3:end], fortys_max[3:end])
# scatter!(fortys_x_min1, fortys_min1)

# # plot(fortys_x_min2, fortys_min2, linewidth=lw, color=:maroon1)
# plot!(fortys_x_min1, fortys_min1, fillrange=(fortys_max1[3:end]))

# plot!(fortys_x_min2, fortys_min2)
# plot!(fortys_x_max1[2:end-1], fortys_max1[2:end-1])
# plot!(fortys_x_min1[1:end], fortys_min1[1:end])


# scatter(thirtys.wab, thirtys.wr, color=:dodgerblue1, label="", axis=:log)
# plot(thirtys_x_min2, thirtys_min2)
# plot!(thirtys_x_max2, thirtys_max2)
# plot!(thirtys_x_min1, thirtys_min1)
# plot!(thirtys_x_max1, thirtys_max1)
# plot!(thirtys_x_max1, thirtys_max1, fillrange=(thirtys_min1))
# plot!(thirtys_x_max2, thirtys_max2, fillrange=(thirtys_min2))


# lw=15

# plot(zero_x, zero_min, fillrange=(zero_max), axis=:log, linewidth=lw, color=:tomato1, label="0%", xlabel="ω_ab", ylabel="ω_r")

# plot!(tens_x1, tens_min1, linewidth=lw, color=:darkorange, label="10-20%")
# plot!(tens_x2, tens_max2, fillrange=(tens_min2), linewidth=lw, color=:darkorange, label="")

# plot!(twentys_x_max2, twentys_max2, fillrange=(twentys_min2), linewidth=lw, color=:goldenrod1, label="20-30%")
# plot!(twentys_x_max2, twentys_min2, linewidth=lw, color=:goldenrod1, label="")

# plot!(twentys_x_max1, twentys_max1, fillrange=(twentys_min1), linewidth=lw, color=:goldenrod1, label="")
# plot!(twentys_x_max1, twentys_min1, linewidth=lw, color=:goldenrod1, label="")

# plot!(thirtys_x_min1, thirtys_min1, fillrange=(thirtys_max1), linewidth=lw, color=:dodgerblue1, label="30-40%")
# plot!(thirtys_x_min1, thirtys_max1, linewidth=lw, color=:dodgerblue1, label="")

# plot!(thirtys_x_min2, thirtys_min2, fillrange=(thirtys_max2), linewidth=lw, color=:dodgerblue1, label="")
# plot!(thirtys_x_min2, thirtys_max2, linewidth=lw, color=:dodgerblue1, label="")

# plot!(fortys_x[1:end-3], fortys_min[1:end-3], color=:maroon1, linewidth=lw, label="")
# plot!(fortys_x[3:end], fortys_max[3:end], color=:maroon1)
# plot!([fortys.wab[27];fortys.wab[34];fortys.wab[41];fortys.wab[35];fortys.wab[42];fortys.wab[36];fortys.wab[43];fortys.wab[37];fortys.wab[44];fortys.wab[38];fortys.wab[45];fortys.wab[39];fortys.wab[46];fortys.wab[40]], [fortys.wr[27];fortys.wr[34];fortys.wr[41];fortys.wr[35];fortys.wr[42];fortys.wr[36];fortys.wr[43];fortys.wr[37];fortys.wr[44];fortys.wr[38];fortys.wr[45];fortys.wr[39];fortys.wr[46];fortys.wr[40]], color=:maroon1)

# plot!(fiftys_x, fiftys_min, fillrange=(fiftys_max), linewidth=lw, color=:purple1, label=">50%")
# plot!(fiftys_x, fiftys_max, linewidth=lw, color=:purple1, label="")

# plot!(fortys_x_max1, fortys_max1, fillrange=(fortys_min1), color=:maroon1)


# scatter(zeros.wab, zeros.wr, axis=:log, color=:tomato1, label="")
# scatter!(tens.wab, tens.wr, color=:darkorange, label="")
# scatter!(twentys.wab, twentys.wr, color=:goldenrod1, label="",axis=:log)
# scatter!(thirtys.wab, thirtys.wr, color=:dodgerblue1, label="")
# scatter!(fortys.wab, fortys.wr, color=:maroon1, label="")
# scatter!(fiftys.wab, fiftys.wr, color=:purple1, label="")

# perc_df