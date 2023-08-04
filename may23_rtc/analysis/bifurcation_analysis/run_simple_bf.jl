using Plots, Printf, Measures
using Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra, DataFrames

include("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/bf_funcs.jl")


params1 = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
kdeg = 0.001, kin = 0.022222222, ω_ab = 1, ω_r = 0.0001, 
kdam =  0.01, lam = 0.014) 	

params2 = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
kdeg = 0.001, kin = 0.022222222, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
kdam =  0.01, lam = 0.014) 	


br = get_br(rtc_mod, params1, initial, 1.)
br2 = get_br(rtc_mod, params2, initial, 3.)

plot(br, br2, vars = (:param, :rh), linewidthstable=2.5)
plot(br, br2, vars = (:param, :rd), linewidthstable=2.5)
plot(br, br2, vars = (:param, :rt), linewidthstable=2.5)
plot(br, br2, vars = (:param, :rm_a), linewidthstable=2.5)
plot(br, br2, vars = (:param, :rm_r), linewidthstable=2.5)
plot(br, br2, vars = (:param, :rtca), linewidthstable=2.5)
plot(br, br2, vars = (:param, :rtcr), linewidthstable=2.5)

p1 = Plots.plot(br2, vars=(:param, :rh), label="Rh")
p1 = Plots.plot!(br2, vars=(:param, :rt), label="Rt")
p1 = Plots.plot!(br2, vars=(:param, :rd), label="Rd")

savefig(p1, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/ribo_bf_plot.svg")

Plots.plot(br, vars=(:param, :rh))
Plots.plot!(br, vars=(:param, :rt))
Plots.plot!(br, vars=(:param, :rd))

Plots.plot(br, vars=(:param, :rm_a))
Plots.plot!(br, vars=(:param, :rm_r))
Plots.plot!(br, vars=(:param, :rtca))
Plots.plot!(br, vars=(:param, :rtcr))
 

p2 = Plots.plot(br2, vars=(:param, :rm_a), label="mRNA RtcBA")
p2 = Plots.plot!(br2, vars=(:param, :rtca), label="RtcA")
p2 = Plots.plot!(br2, vars=(:param, :rtcb), label="RtcB")
p2 = Plots.plot!(br2, vars=(:param, :rtcr), label="RtcR")

savefig(p2, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/mrnaprotein_bf_plot.svg")

p = Plots.plot(p2,p1,p3,p4,  size=(900,700), layout=@layout [a b; c d])
Plots.savefig(p, "/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/solve_vs_bfplot.svg")

function create_br_df(br)
    df = DataFrame(rm_a=[],rtca=[],rm_b=[],rtcb=[],rm_r=[],rtcr=[],rh=[],rd=[],rt=[],kdam=[])
    for i in eachcol(df)
        println(i)
    end
    for i in range(1,length(br.sol))
        for (s,d) in zip(range(1,9), eachcol(df))
            push!(d, br.sol[i][1][s])
        end
        push!(df.kdam, br.sol[i][2])
    end
    
    return df
end

df = create_br_df(br2)

using PlotlyJS
rma_p = scatter(x=df.kdam, y=df.rm_a, name="rm_a")
rtca_p = scatter(x=df.kdam, y=df.rtca, name="rtca")
rmb_p = scatter(x=df.kdam, y=df.rm_b, name="rm_b")
rtcb_p = scatter(x=df.kdam, y=df.rtcb, name="rtcb")
rmr_p = scatter(x=df.kdam, y=df.rm_r, name="rm_r")
rtcr_p = scatter(x=df.kdam, y=df.rtcr, name="rtcr")
rh_p = scatter(x=df.kdam, y=df.rh, name="rh")
rd_p = scatter(x=df.kdam, y=df.rd, name="rd")
rt_p = scatter(x=df.kdam, y=df.rt, name="rt")
plot([rma_p, rtca_p, rmb_p, rtcb_p, rmr_p, rtcr_p, rh_p, rd_p, rt_p])


atp = 3578.9473684210525 

function get_extra_vars(df, atp)
    kr=0.125
    L=10
    c=0.001
    Vmax_init = 39.51 
    Km_init = 250.
    fas=[]
    ras=[]
    alphas=[]
    vinits=[]
    for i in range(1,length(df.rt))
        alpha = df.rt[i]/kr
        push!(alphas, alpha)
        fa = (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6)
        push!(fas, fa)
        ra = fa*df.rtcr[i]
        push!(ras, ra)
        Vinit = ra*Vmax_init*atp/(Km_init+atp)
        push!(vinits,Vinit)
    end 
    return fas, ras, alphas
end



using PlotlyJS

fas, ras = get_extra_vars(df, atp)
fa = scatter(x=df.kdam, y=fas, mode="markers",  marker=attr(size=5, color=0:455), name="fa")#, Layout(yaxis_type=:log))
rtcr = scatter(x=df.kdam, y=df.rtcr, mode="markers",  marker=attr(size=5), name="rtcr")#, Layout(c="Viridis")))
ra = scatter(x=df.kdam, y=ras, mode="markers",  marker=attr(size=5), name="ra")#, Layout(yaxis_type=:log)
rm = scatter(x=df.kdam, y=df.rm_a, mode="markers",  marker=attr(size=5), name="mRNA RtcBA")#, Layout(c="Viridis")))
rt = scatter(x=df.kdam[2:end], y=df.rt[2:end], mode="markers",  marker=attr(size=5), name="rt")#, Layout(c="Viridis")))
rh = scatter(x=df.kdam, y=df.rh, mode="markers",  marker=attr(size=5), name="rh")#, Layout(c="Viridis")))

bf_line1 = scatter(x=[br2.specialpoint[2].param,br2.specialpoint[2].param],y=[1e-4,3])
bf_line2 = scatter(x=[br2.specialpoint[3].param,br2.specialpoint[3].param],y=[1e-4,3])
plot([fa,rtcr,ra, rm, rt, rh, bf_line1, bf_line2], Layout(yaxis_type=:log))


kr=0.125
L=10.
c=0.001

params1 = (L = L, c = c, kr = kr, Vmax_init = 39.51, Km_init = 250.,
θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
kdeg = 0.001, kin = 0.022222222, ω_ab = 1, ω_r = 0.0001, 
kdam =  0.01, lam = 0.014) 	

br = get_br(rtc_mod, params1, initial, 8.)

df1 = create_br_df(br)

fas1, ras1, alphas1 = get_extra_vars(df1, atp)

fa1 = scatter(x=df1.kdam, y=fas1, mode="markers",  marker=attr(size=5, color=0:length(df1.rtcr)), name="fa")
rtcr1 = scatter(x=df1.kdam, y=df1.rtcr, mode="markers",  marker=attr(size=5), name="rtcr")
ra1 = scatter(x=df1.kdam, y=ras1, mode="markers",  marker=attr(size=5), name="ra")
rm1 = scatter(x=df1.kdam, y=df1.rm_a, mode="markers",  marker=attr(size=5), name="mRNA RtcBA")
rt1 = scatter(x=df1.kdam[2:end], y=df1.rt[2:end], mode="markers",  marker=attr(size=5), name="rt")
rh1 = scatter(x=df1.kdam, y=df1.rh, mode="markers",  marker=attr(size=5), name="rh")
rtca1 = scatter(x=df1.kdam, y=df1.rtca, mode="markers",  marker=attr(size=5), name="RtcA")
alpha1 = scatter(x=df1.kdam[2:end], y=alphas1[2:end], mode="markers",  marker=attr(size=5), name="alpha")
rd1 = scatter(x=df1.kdam, y=df1.rd, mode="markers",  marker=attr(size=5), name="rd")

bf_line11 = scatter(x=[br.specialpoint[2].param,br.specialpoint[2].param],y=[1e-6,2.8])
bf_line21 = scatter(x=[br.specialpoint[3].param,br.specialpoint[3].param],y=[1e-6,2.8])
plot([fa1,rtcr1,ra1, rm1, rt1, rh1, rd1, rtca1, bf_line11, bf_line21], Layout(yaxis_type=:log))






plot(scatter(x=df.kdam, y=fas, mode="markers",  marker=attr(size=5, color=0:455)))#, Layout(yaxis_type=:log))#, Layout(c="Viridis")))
plot(scatter(x=df.kdam, y=ras, mode="markers",  marker=attr(size=5, color=0:455)), Layout(yaxis_type=:log))#, Layout(c="Viridis")))
plot(scatter(x=df.kdam, y=df.rtcr, mode="markers",  marker=attr(size=5, color=0:455)))#, Layout(c="Viridis")))
plot(scatter(x=df.kdam, y=vinits, mode="markers",  marker=attr(size=5, color=0:455)))#, Layout(c="Viridis")))


plot(df.kdam, fas, xlims=(0,3), ylims=(0.999999,1.000001), linewidth=0.5)
plot(df.kdam, fas, xlims=(0,3), ylims=(0.99999,1.00001), linewidth=0.5)
plot(df.kdam, fas, xlims=(0,3), ylims=(0.9999,1.0001), linewidth=0.5)
plot(df.kdam, fas, xlims=(0,3), ylims=(0.999,1.001), linewidth=0.5)
plot(df.kdam, fas, xlims=(0,3), ylims=(0.99,1.01), linewidth=0.5)
plot(df.kdam, fas, xlims=(0,3), ylims=(0.9,1.1), linewidth=0.5)


plots = []
for i in (10 .^ range(-4, stop=0, length = 5))
    for j in (10 .^ range(-2, stop=1, length = 5))
        params1 = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
        θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
        krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
        kdeg = 0.001, kin = 0.022222222, ω_ab = j, ω_r = i, 
        kdam =  0.01, lam = 0.014) 	
        @show (params1[:ω_ab], params1[:ω_r])
        push!(plots, plot(get_br(rtc_mod, params1, initial, 3.), title="ω_ab = $(@sprintf "%.2E" j), ω_r = $(@sprintf "%.2E" i)"))
    end
end

p_all = plot(plots[1], plots[2], plots[3], plots[4], plots[5], plots[6],
plots[7], plots[8], plots[9], plots[10], plots[11], plots[12],
plots[13], plots[14], plots[15], plots[16], plots[17], plots[18],
plots[19], plots[20], plots[21], plots[22], plots[23], plots[24], plots[25],
 margin=1mm, size=(2000,2000), layout=(5,5))

plots_br = []
for i in range(1,stop=5000,length=25)
    params1 = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
    θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
    krep = 137., ktag = 9780., atp = i, km_a = 20., km_b = 16., g_max = 2.0923, 
    kdeg = 0.001, kin = 0.022222222, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
    kdam =  0.01, lam = 0.014) 	
    
    push!(plots_br, plot(get_br(rtc_mod, params1, initial, 3.), title="atp=$i"))
end

plot(plots_br[1], plots_br[2], plots_br[3], plots_br[4], plots_br[5], plots_br[6],
plots_br[7], plots_br[8], plots_br[9], plots_br[10], plots_br[11], plots_br[12],
plots_br[13], plots_br[14], plots_br[15], plots_br[16], plots_br[17], plots_br[18],
plots_br[19], plots_br[20], plots_br[21], plots_br[22], plots_br[23], plots_br[24], plots_br[25],
 margin=1mm, size=(2000,2000), layout=(5,5))