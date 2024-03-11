using PlotlyJS, Printf, Measures, CSV, DataInterpolations
using Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra, DataFrames

include("$PATHmay23_rtc/analysis/bifurcation_analysis/bf_funcs.jl"); include("$PATHmay23_rtc/functions/set_ups.jl");


# load time varying parameters and create object so they can be plotted 
t, atp_t, lam_t, kin_t = set_time_vars("$PATHdata/atp_for_rtcmodel.csv")
plotlamt = view(lam_t, 1:200066)
plotatpt = view(atp_t, 1:200066)
plotkint = view(kin_t, 1:200066)
plotatpt = [(plotatpt...)...]
plotlamt = [(plotlamt...)...]
plotkint = [(plotkint...)...]

plot(scatter(x=plotatpt, y=plotlamt))
plot(scatter(x=plotatpt[1:5], y=plotlamt[1:5]))

df_params = DataFrame(atp=plotatpt, lam=plotlamt, kin=plotkint)

# triple param vary and 3D plots

atp_range = range(500, stop=5500,length=20)
kin_range = range(0,stop=0.1,length=20)
lam_range = range(0.001,stop=0.04,length=20)

function triple_param_vary(param_range1, param1, param_range2, param2, param_range3, param3, params1)
    df = DataFrame(atp = Float64[], lam = Float64[], kin = Float64[], kdam1 = Float64[], kdam2 = Float64[], bs = Symbol[])
    params = deepcopy(params1)
    for i in param_range1
        params = merge(params, (param1=>i,))
        for j in param_range2
            params = merge(params, (param2=>j,))
            for k in param_range3
                params = merge(params, (param3=>k,))
                # @show params[:atp], params[:lam], params[:kin]
                br = get_br(rtc_mod, params, initial, 1.)
                if length(br.specialpoint) == 2
                    push!(df, (i, j, k, br.specialpoint[1].param, br.specialpoint[2].param, br.specialpoint[1].type))
                else
                    push!(df, (i, j, k, br.specialpoint[2].param, br.specialpoint[3].param, br.specialpoint[2].type))
                end
            end
        end    
    end
    return df
end

df = triple_param_vary(atp_range, :atp, lam_range, :lam, kin_range, :kin, params1)
# CSV.write("$PATHmay23_rtc/analysis/bifurcation_analysis/plots/df_triple_param_vary_100.csv", df)

bsp = df[df.bs .== :bp, :]

df = DataFrame(atp=atp_range, lam=lam_range, kin=kin_range)

bsp = DataFrame(CSV.File("$PATHmay23_rtc/analysis/bifurcation_analysis/plots/bsp_triple_param_vary_50.csv"))


trace1 = scatter(bsp, x=:atp, y=:lam, z=:kin, marker=attr(color=:kin, colorscale="Viridis"), type="scatter3d", mode="markers", name="")
trace2 = scatter(df_params, x=:atp, y=:lam, z=:kin, type="scatter3d", name="")
p = plot([trace1,trace2], Layout(scene=attr(xaxis_title="ATP", yaxis_title="λ", zaxis_title="kin"), showlegend=false))

PlotlyJS.savefig(p, "$PATHmay23_rtc/analysis/bifurcation_analysis/plots/3d_bs_region.png")

open("$PATHmay23_rtc/analysis/bifurcation_analysis/plots/3d_bs_region.html", "w") do io
    PlotlyBase.to_html(io, p.plot)
end




# bsp = DataFrame(CSV.File("$PATHmay23_rtc/analysis/bifurcation_analysis/plots/bsp_triple_param_vary_50.csv"))

wab_range1 = (10 .^ range(-3, stop=0, length = 3))

bsps = []
for i in wab_range1 
    params1 = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
    θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
    krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
    kdeg = 0.001, kin = 0.022222222, ω_ab = i, ω_r = 0.0001, 
    kdam =  0.01, lam = 0.014) 
    df = triple_param_vary(atp_range, :atp, lam_range, :lam, kin_range, :kin, params1)
    push!(bsps, df[df.bs .== :bp, :])
end

CSV.write("$PATHmay23_rtc/analysis/bifurcation_analysis/bf_3d_dfs/wab1e-3_20.csv", bsps[1])
CSV.write("$PATHmay23_rtc/analysis/bifurcation_analysis/bf_3d_dfs/wab0.0316_20.csv", bsps[2])
CSV.write("$PATHmay23_rtc/analysis/bifurcation_analysis/bf_3d_dfs/wab1_20.csv", bsps[3])

trace1 = scatter(bsps[1], x=:atp, y=:lam, z=:kin, marker=attr(color=:kin, colorscale=colors.viridis, opacity=1), type="scatter3d", mode="markers", name="ω_ab = 0.001")
trace2 = scatter(bsps[2], x=:atp, y=:lam, z=:kin, marker=attr(color=:kin, colorscale=colors.Wistia, opacity=1), type="scatter3d", mode="markers", name="ω_ab = 0.0316")
trace3 = scatter(bsps[3], x=:atp, y=:lam, z=:kin, marker=attr(color=:kin, colorscale=colors.cool, opacity=1), type="scatter3d", mode="markers", name="ω_ab = 1")

trace = scatter(df_params, x=:atp, y=:lam, z=:kin, type="scatter3d", name="")
p = plot([trace1,trace2,trace3], Layout(scene=attr(xaxis_title="ATP", yaxis_title="λ", zaxis_title="kin")))

p = plot([trace2,trace], Layout(scene=attr(xaxis_title="ATP", yaxis_title="λ", zaxis_title="kin")))#, showlegend=false))

open("$PATHmay23_rtc/analysis/bifurcation_analysis/plots/3dplots/3d_wab.html", "w") do io
    PlotlyBase.to_html(io, p.plot)
end

# wab001 = bsps[1]
# wab0056 = bsps[2]
# wab032 = bsps[3]
# wab178 = bsps[4]
# wab1 = bsps[5]

# CSV.write("$PATHmay23_rtc/analysis/bifurcation_analysis/bf_3d_dfs/wab001.csv", wab001)
# CSV.write("$PATHmay23_rtc/analysis/bifurcation_analysis/bf_3d_dfs/wab0056.csv", wab0056)
# CSV.write("$PATHmay23_rtc/analysis/bifurcation_analysis/bf_3d_dfs/wab032.csv", wab032)
# CSV.write("$PATHmay23_rtc/analysis/bifurcation_analysis/bf_3d_dfs/wab178.csv", wab178)
# CSV.write("$PATHmay23_rtc/analysis/bifurcation_analysis/bf_3d_dfs/wab1.csv", wab1)

wab001 = DataFrame(CSV.File("$PATHmay23_rtc/analysis/bifurcation_analysis/bf_3d_dfs/wab001.csv"))
wab0056 = DataFrame(CSV.File("$PATHmay23_rtc/analysis/bifurcation_analysis/bf_3d_dfs/wab0056.csv"))
wab032 = DataFrame(CSV.File("$PATHmay23_rtc/analysis/bifurcation_analysis/bf_3d_dfs/wab032.csv"))
wab178 = DataFrame(CSV.File("$PATHmay23_rtc/analysis/bifurcation_analysis/bf_3d_dfs/wab178.csv"))
wab1 = DataFrame(CSV.File("$PATHmay23_rtc/analysis/bifurcation_analysis/bf_3d_dfs/wab1.csv"))


# using Plots

# Plots.scatter(wab001.atp,wab001.lam,wab001.kin, camera=(-30,30))# xlims=(minimum(atp_range), maximum(atp_range)), ylims=(minimum(lam_range), maximum(lam_range)), zlims=(minimum(kin_range), maximum(kin_range)))
# # Plots.scatter!(wab0056.atp,wab0056.lam,wab0056.kin, xlims=(minimum(atp_range), maximum(atp_range)), ylims=(minimum(lam_range), maximum(lam_range)), zlims=(minimum(kin_range), maximum(kin_range)))
# Plots.scatter!(wab032.atp,wab032.lam,wab032.kin)#, xlims=(minimum(atp_range), maximum(atp_range)), ylims=(minimum(lam_range), maximum(lam_range)), zlims=(minimum(kin_range), maximum(kin_range)))
# # Plots.scatter!(wab178.atp,wab178.lam,wab178.kin, xlims=(minimum(atp_range), maximum(atp_range)), ylims=(minimum(lam_range), maximum(lam_range)), zlims=(minimum(kin_range), maximum(kin_range)))
# Plots.scatter!(wab1.atp,wab1.lam,wab1.kin)#, xlims=(minimum(atp_range), maximum(atp_range)), ylims=(minimum(lam_range), maximum(lam_range)), zlims=(minimum(kin_range), maximum(kin_range)))

# Plots.plot!(plotatpt,plotlamt,plotkint)



wr_range1 = (10 .^ range(-5, stop=-1, length = 3))

bsps = []
for i in wr_range1 
    params1 = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
    θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
    krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
    kdeg = 0.001, kin = 0.022222222, ω_ab = 0.1, ω_r = i, 
    kdam =  0.01, lam = 0.014) 
    df = triple_param_vary(atp_range, :atp, lam_range, :lam, kin_range, :kin, params1)
    push!(bsps, df[df.bs .== :bp, :])
end




trace1 = scatter(bsps[1], x=:atp, y=:lam, z=:kin, marker=attr(color=:kin, colorscale=colors.viridis, opacity=1), type="scatter3d", mode="markers", name="ω_r = 1e-5")
trace2 = scatter(bsps[2], x=:atp, y=:lam, z=:kin, marker=attr(color=:kin, colorscale=colors.Wistia, opacity=1), type="scatter3d", mode="markers", name="ω_r = 0.001")
trace3 = scatter(bsps[3], x=:atp, y=:lam, z=:kin, marker=attr(color=:kin, colorscale=colors.cool, opacity=1), type="scatter3d", mode="markers", name="ω_r = 0.1")

trace = scatter(df_params, x=:atp, y=:lam, z=:kin, type="scatter3d", name="")
p = plot([trace1,trace2,trace3, trace], Layout(scene=attr(xaxis_title="ATP", yaxis_title="λ", zaxis_title="kin")))

p = plot([trace1,trace], Layout(scene=attr(xaxis_title="ATP", yaxis_title="λ", zaxis_title="kin")))


open("$PATHmay23_rtc/analysis/bifurcation_analysis/plots/3dplots/3d_wr_range_params.html", "w") do io
    PlotlyBase.to_html(io, p.plot)
end