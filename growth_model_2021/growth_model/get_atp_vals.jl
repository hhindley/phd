using DifferentialEquations, DataFrames, Measures, CSV, DiffEqCallbacks
using PlotlyJS

include("model.jl")
include("parameters.jl")
include("functions.jl")


tspan = (0,1e9)

solu = simple_solve!(odemodel!, init, tspan, params)

plotly_plot_sol(solu, "log", "")

species = [:cr, :em, :cp, :cq, :ct, :et, :cm, :mt, :mm, :q, :p, :si, :mq, :mp, :mr, :r, :a]
initsolDF = DataFrame([[j[i] for j in solu.u] for i=1:length(solu.u[1])], species)
sscr = initsolDF[end,:][:cr]
ssem = initsolDF[end,:][:em]
sscp = initsolDF[end,:][:cp]
sscq = initsolDF[end,:][:cq]
ssct = initsolDF[end,:][:ct]
sscm = initsolDF[end,:][:cm]
ssmt = initsolDF[end,:][:mt]
ssmm = initsolDF[end,:][:mm]
ssq = initsolDF[end,:][:q]
ssp = initsolDF[end,:][:p]
sssi = initsolDF[end,:][:si]
ssmq = initsolDF[end,:][:mq]
ssmp = initsolDF[end,:][:mp]
ssmr = initsolDF[end,:][:mr]
ssr = initsolDF[end,:][:r]
ssa = initsolDF[end,:][:a]
# sss = initsolDF[end,:][:s]
# ssN = initsolDF[end,:][:N]
sset = 1
ssinit = [sscr, ssem, sscp, sscq, ssct, sset, sscm, ssmt, ssmm, ssq, ssp, sssi, ssmq, ssmp, ssmr, ssr, ssa, s_0, N_0]

ns = 100;
pop_params = [dm, kb, ku, f, thetar, gmax, thetax, Kt, M, we, Km, vm, nx, Kq, vt, wr, wq, wp, nq, nr, ns, kin, d_s, d_n]

solu_ss = simple_solve!(pop_model!, ssinit, tspan, pop_params)

plotly_plot_sol_pop(solu_ss, "log", "")

M = 10^8
cr = get_curve(solu_ss, :cr); cp = get_curve(solu_ss, :cp); cq = get_curve(solu_ss, :cq); ct = get_curve(solu_ss, :ct); cm = get_curve(solu_ss, :cm);
Rt = @. cr+cp+cq+ct+cm
g_max = 1260
k_g = 7
a1 = get_curve(solu_ss, :a)
SF = 1e6/(6.022e23*1e-15)
a = a1*SF

lam = @. (Rt*(g_max*a/(k_g+a)))/(M)
p_atp1 = plot(scatter(x=solu_ss.t, y=a1), Layout(xaxis_type="log"))
p_atp = plot(scatter(x=solu_ss.t, y=a), Layout(xaxis_type="log"))
p_lam = plot(scatter(x=solu_ss.t, y=lam), Layout(xaxis_type="log"))
open("./atp_growth_model.html", "w") do io
    PlotlyBase.to_html(io, p_atp.plot)
end

plot(scatter(x=lam, y=a), Layout(yaxis_type="log", yaxis_title="ATP (μM)", xaxis_title="λ (/min)"))


plot(scatter(x=lam, y=lam_vs_a), Layout(yaxis_type="log"))

#load other growth rate 
csv_lam = DataFrame(CSV.File("/home/holliehindley/phd/data/results_colD_grfit.csv"))
csv_lam = select!(csv_lam, Not(["log(OD)", "log(OD) error", "gr error", "od"]))
csv_lam."gr"[csv_lam."gr".< 0] .= 0 #zero(eltype(lam_colD))


plot(scatter(x=csv_lam."t", y=csv_lam."gr"))
# for the actual data growth rate find the value of atp to use in rtc model 
atp_gr = []
for i in csv_lam."gr"
    push!(atp_gr, lam_vs_a(i))
end
atp_gr

csv_lam."gr"[1]
lam_vs_a(csv_lam."gr"[1])

plot(scatter(x=solu_ss.t, y=atp_gr))#, Layout(yaxis_type="log", xaxis_type="log"))

#make a dataframe with the atp values against time (for data) - now in rtc model script need to interpolate this so it can be used over time 
df_atp = DataFrame(t=csv_lam."t", atp=atp_gr)
CSV.write("/home/holliehindley/phd/data/atp_curve_from_growth_model.csv", df_atp)



using Dierckx
interp_func = Spline2D(lam, solu_ss.t, a, kx=3, ky=3, s=1e13)
#check Interpolation
atp_check = []
for i in 1:length(solu_ss.t)
    push!(atp_check, interp_func(lam[i], solu_ss.t[i]))
end
plot(scatter(x=solu_ss.t, y=atp_check), Layout(xaxis_type="log"))




atp_from_gr = []
for i in 1:length(csv_lam."gr")
    push!(atp_from_gr, interp_func(csv_lam."gr"[i], csv_lam."t"[i]))
end


plot(scatter(x=csv_lam."t", y=atp_from_gr))







plot(contour(x=lam, y=solu_ss.t, z=a))
scatterplot = scatter3d(
    x=a,
    y=lam,
    z=solu_ss.t,
    mode="markers",
    marker_color=y
)



layout=Layout(scene=attr(zaxis=attr( range=[z[1],650]), xaxis_title="ATP", yaxis_title="λ", zaxis_title="t"))
# Create the plot
plot(scatterplot, layout)


using Interpolations

itp_a = LinearInterpolation(solu_ss.t, a)

a_at_growth_rate = [itp_a(t) for t in solu_ss.t]

combined_data = hcat(lam, solu_ss.t, a_at_growth_rate)



