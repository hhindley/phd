using Plots
using Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra, DataFrames
include("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/bf_funcs.jl")


br = get_br(params1, initial)

plot(br, vars=(:param, :rh), putspecialptlegend=false, label="rfw")

p_rh = plot(br, vars = (:param, :rh), label="Healthy ribosomes", c=:palevioletred, putspecialptlegend = false)
p_rh1 = plot!(twinx(), br, vars = (:param, :rtca), c=:grey60, label="RtcBA")

# savefig(p_rm_a, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/p_rma.svg")
# savefig(p_rh, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/p_rh.svg")
# savefig(p_rh1, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/p_both.svg")



function plot_bf_against_kdam(param_range, param)
    df=DataFrame(param = Float64[], kdam1 = Float64[], kdam2 = Float64[])
    paramscopy = deepcopy(params1)
    for i in param_range
        params = merge(paramscopy, (param=>i,))
        # println("$param = $i")
        br = get_br(params, initial)
        for b in range(1,length(br.specialpoint))
            if br.specialpoint[b].type == :endpoint
                Nothing
                # push!(df_none, (i, br.specialpoint[1].param, br.specialpoint[2].param))
            else
                push!(df, (i, br.specialpoint[2].param, br.specialpoint[3].param))
            end
        end
    end
    df = df[1:2:end, :]
    p = plot(df.param, df.kdam1, xlabel="$param", ylabel="kdam", ylims=(0,4), label="bf1")
    p = plot!(df.param, df.kdam2, label="bf2")

end

atp_range = range(500,stop=5000,length=500)
kin_range = range(0,stop=0.05,length=500)
lam_range = range(0.001,stop=0.04,length=500)
wab_range = range(0.01, stop=4, length=500)
wr_range = (range(0.00001,stop=0.001,length=200))

p_atp = plot_bf_against_kdam(atp_range, :atp);
p_kin = plot_bf_against_kdam(kin_range, :kin);
p_lam = plot_bf_against_kdam(lam_range, :lam);
p_wab = plot_bf_against_kdam(wab_range, :ω_ab);
p_wr = plot_bf_against_kdam(wr_range, :ω_r);


l = @layout [a b; c d; e];
all_bf = plot(p_atp, p_kin, p_lam, p_wab, p_wr, layout=l, size=(1200,800), left_margin=4Plots.mm, right_margin=2Plots.mm, plot_title="Bifurcation points for kdam vs. parameter")
savefig(all_bf, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/all_bf.svg")









# using PlotlyJS
# p1 = PlotlyJS.scatter(
#     bsp,
#     x=:atp, y=:wab, z=:kdam1,
#     type="scatter3d", mode="markers", name="bf1"
# )
# p2 = PlotlyJS.scatter(
#     bsp,
#     x=:atp, y=:wab, z=:kdam2,
#     type="scatter3d", mode="markers"
# )

# layout = Layout(
#     scene=attr(
#         zaxis=attr(
#             range=[0,1],
#             title="kdam"),
#         xaxis=attr(
#             title="ATP"),    
#         yaxis=attr(
#             title="ω_ab"),)
# )
# plot([p1,p2], layout)




# function plot_3d_bf_against_kdam(param_range1, param1, param_range2, param2)
#     df = double_param_vary(atp_range, :atp, kin_range, :kin)
#     bsp = df[df.bs .== :bp, :]
#     p1 = PlotlyJS.scatter(
#         bsp,
#         x=:atp, y=:wab, z=:kdam1,
#         type="scatter3d", mode="markers", label="bf1"
#     )
#     p2 = PlotlyJS.scatter(
#         bsp,
#         x=:atp, y=:wab, z=:kdam2,
#         type="scatter3d", mode="markers"
#     )
    
#     layout = Layout(
#         scene=attr(
#             zaxis=attr(
#                 range=[0,1],
#                 title="kdam"),
#             xaxis=attr(
#                 title="ATP"),    
#             yaxis=attr(
#                 title="ω_ab"),)
#     )
#     plot([p1,p2], layout)

# end

