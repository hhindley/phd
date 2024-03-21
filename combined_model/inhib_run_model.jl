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
include("$PATH/combined_model/combined_inhib.jl")


abx_range = range(0,20,length=100)
w_rtc = sweep_param(combined_model, ssvals_comb, new_p, abx_range, abx, species_comb)
wo_rtca = sweep_param(combined_inhib_rtca, ssvals_inhib_comb_rtca, params_inhib_comb, abx_range, abx, species_inhib_comb)
wo_rtcb = sweep_param(combined_inhib_rtcb, ssvals_inhib_comb_rtcb, params_inhib_comb, abx_range, abx, species_inhib_comb)
wo_rtcr = sweep_param(combined_inhib_rtcr, ssvals_inhib_comb_rtcr, params_inhib_comb, abx_range, abx, species_inhib_comb)

function plot_comp(var)
    return plot([
        scatter(x=abx_range,y=w_rtc[!,var],name="with Rtc"), 
                scatter(x=abx_range,y=wo_rtca[!,var],name="without RtcA"),
                scatter(x=abx_range,y=wo_rtcb[!,var],name="without RtcB"),
                scatter(x=abx_range,y=wo_rtcr[!,var],name="without RtcR"),
                ], Layout(xaxis_title="abx",yaxis_title="$var"))
end

new_species = [:m_R, :m_A, :m_B, :m_rh, :R, :A, :B, :q, :rh, :rt, :rd, :a, :lam, :rmf]
[display(plot_comp(i)) for i in new_species]



#double parameter sweep of abx and kdam_p with rtc and without rtc
abx_range = range(0,12,length=25)
kdamp_range = range(0,1,length=25)

res = sweep_paramx2(combined_model, ssvals_comb, params_comb, species_comb, abx, kdam_p, abx_range, kdamp_range)
res_nortca = sweep_paramx2(combined_inhib_rtca, ssvals_inhib_comb_rtca, params_inhib_comb, species_inhib_comb, abx, kdam_p, abx_range, kdamp_range)
res_nortcb = sweep_paramx2(combined_inhib_rtcb, ssvals_inhib_comb_rtcb, params_inhib_comb, species_inhib_comb, abx, kdam_p, abx_range, kdamp_range)
res_nortcar = sweep_paramx2(combined_inhib_rtcr, ssvals_inhib_comb_rtcr, params_inhib_comb, species_inhib_comb, abx, kdam_p, abx_range, kdamp_range)

function plot_contour(res, specie, param_range1, param_range2, param1, param2, title, x, y)
    rh_res = reshape(res[!,specie], length(param_range1),length(param_range2))
    return plot(contour(x=param_range1, y=param_range2, z=rh_res, coloraxis="coloraxis"))#colorbar=attr(title="$specie", titleside="right", x=x, y=y, len=0.5)), Layout(xaxis_title=param1, yaxis_title=param2, title=title))
end

p2 = [(plot_contour(res, i, abx_range, kdamp_range, "abx", "kdam_p", "rtc active", 0.45, 0.2)) for i in new_species];
p3 = [(plot_contour(res_nortca, i, abx_range, kdamp_range, "abx", "kdam_p", "no rtcA", 1, 0.2)) for i in new_species];
p4 = [(plot_contour(res_nortcb, i, abx_range, kdamp_range, "abx", "kdam_p", "no rtcB", 1, 0.2)) for i in new_species];
p5 = [(plot_contour(res_nortcar, i, abx_range, kdamp_range, "abx", "kdam_p", "no rtcR", 1, 0.9)) for i in new_species];


# [[p2[i] for i in range(1,length(p2))] [p3[i] for i in range(1,length(p3))]]

for i in range(1,length(p2))
    display([p2[i] p3[i]; p4[i] p5[i]])
end





function plot_heatmap(specie)
    res1 = reshape(res[!,specie], length(abx_range),length(kdamp_range))
    res2 = reshape(res_nortcb[!,specie], length(abx_range),length(kdamp_range))
    res3 = reshape(res_nortca[!,specie], length(abx_range),length(kdamp_range))
    res4 = reshape(res_nortcar[!,specie], length(abx_range),length(kdamp_range))

    fig1 = make_subplots(rows=2,cols=2, shared_xaxes=true, shared_yaxes=true, vertical_spacing=0.05, horizontal_spacing=0.05,
    subplot_titles=["Rtc active" "RtcA inhib"; "RtcB inhib" "RtcR inhib"], y_title="abx", x_title="kdam_p");
    add_trace!(fig1, heatmap(x=abx_range, y=kdamp_range, z=res1, coloraxis="coloraxis"), row=1, col=1)
    add_trace!(fig1, heatmap(x=abx_range, y=kdamp_range, z=res2, coloraxis="coloraxis"), row=1, col=2)
    add_trace!(fig1, heatmap(x=abx_range, y=kdamp_range, z=res3, coloraxis="coloraxis"), row=2, col=1)
    add_trace!(fig1, heatmap(x=abx_range, y=kdamp_range, z=res4, coloraxis="coloraxis"), row=2, col=2)
    fig1.plot.layout.annotations[5].yshift=-20
    fig1.plot.layout.annotations[6].xshift=-20
    # relayout!(fig1, coloraxis_colorbar="inosenc")
    fig1.plot.layout.title="$specie"
    return (fig1)
end

heatmaps=[]
for i in new_species
    push!(heatmaps, plot_heatmap(i))
end

heatmaps[1]

for (i,name) in zip(heatmaps,new_species)
    open("$PATHcombined_model/plots/$name.html", "w") do io
        PlotlyBase.to_html(io, i.plot)
    end
end
# res1 = reshape(res[!,:lam], length(abx_range),length(kdamp_range))
# res2 = reshape(res_nortcb[!,:lam], length(abx_range),length(kdamp_range))
# res3 = reshape(res_nortca[!,:lam], length(abx_range),length(kdamp_range))
# res4 = reshape(res_nortcar[!,:lam], length(abx_range),length(kdamp_range))

# fig1 = make_subplots(rows=2,cols=2, shared_xaxes=true, shared_yaxes=true, vertical_spacing=0.05, horizontal_spacing=0.05,
# subplot_titles=["Rtc active" "RtcB inhib"; "RtcR inhib" "RtcR inhib"], y_title="abx", x_title="kdam_p");
# add_trace!(fig1, heatmap(x=abx_range, y=kdamp_range, z=res1, coloraxis="coloraxis"), row=1, col=1)
# add_trace!(fig1, heatmap(x=abx_range, y=kdamp_range, z=res2, coloraxis="coloraxis"), row=1, col=2)
# add_trace!(fig1, heatmap(x=abx_range, y=kdamp_range, z=res3, coloraxis="coloraxis"), row=2, col=1)
# add_trace!(fig1, heatmap(x=abx_range, y=kdamp_range, z=res4, coloraxis="coloraxis"), row=2, col=2)
# fig1.plot.layout.annotations[5].yshift=-20
# fig1.plot.layout.annotations[6].xshift=-20
# display(fig1)



# length(json(fig1.plot.layout.annotations, 2))