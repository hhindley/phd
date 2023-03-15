using DifferentialEquations, DataFrames, Measures, CSV
using PlotlyJS

include("model.jl")
include("parameters.jl")
include("functions.jl")


tspan = (0,1e9)

solu = simple_solve!(odemodel!, init, tspan, params)

plotly_plot_sol(solu, "", "")

species = [:cr, :em, :cp, :cq, :ct, :et, :cm, :mt, :mm, :q, :p, :si, :mq, :mp, :mr, :r, :a]
initsolDF = DataFrame([[j[i] for j in solu.u] for i=1:length(solu.u[1])], species)
sscr = initsolDF[end,:][:cr]
ssem = initsolDF[end,:][:em]
sscp = initsolDF[end,:][:cp]
sscq = initsolDF[end,:][:cq]
ssct = initsolDF[end,:][:ct]
sset = initsolDF[end,:][:et]
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

ssinit = [sscr, ssem, sscp, sscq, ssct, sset, sscm, ssmt, ssmm, ssq, ssp, sssi, ssmq, ssmp, ssmr, ssr, ssa]

solu_ss = simple_solve!(odemodel!, ssinit, tspan, params)

plotly_plot_sol(solu_ss, "", "")

M = 10^8
cr = get_curve(solu_ss, :cr); cp = get_curve(solu_ss, :cp); cq = get_curve(solu_ss, :cq); ct = get_curve(solu_ss, :ct); cm = get_curve(solu_ss, :cm);
Rt = @. cr+cp+cq+ct+cm
g_max = 1260
k_g = 7
a = get_curve(solu_ss, :a)



lam = @. (Rt*(g_max*a/(k_g+a)))/(M)

plot(scatter(x=solu_ss.t, y=lam), Layout(xaxis_type="log"))

using DataInterpolations#, Plots

# interpolate to get a function that will give the value of a at any value of lambda for the growth model 
lam_vs_a = QuadraticInterpolation(a, lam)

# Plots.scatter(lam_full)

plot(scatter(x=lam_vs_a, y=a), Layout(yaxis_type="log", xaxis_type="log"))

#load other growth rate 
csv_lam = DataFrame(CSV.File("/home/holliehindley/phd/data/results_colD_grfit.csv"))
csv_lam = select!(csv_lam, Not(["log(OD)", "log(OD) error", "gr error", "od"]))

# for the actual data growth rate find the value of atp to use in rtc model 
atp_gr = []
for i in csv_lam."gr"
    push!(atp_gr, lam_vs_a(i))
end
atp_gr

plot(scatter(x=solu_ss.t, y=atp_gr), Layout(yaxis_type="log", xaxis_type="log"))

#make a dataframe with the atp values against time (for data) - now in rtc model script need to interpolate this so it can be used over time 
df_atp = DataFrame(t=csv_lam."t", atp=atp_gr)
CSV.write("/home/holliehindley/phd/data/atp_curve_from_growth_model.csv", df_atp)

