using CSV, DataFrames, DataInterpolations, Plots# PlotlyJS

csv = DataFrame(CSV.File("/home/holliehindley/phd/data/results_colD_grfit.csv")) # read csv to a dataframe
gr = csv."gr"
t1 = csv."t"

# gr[gr.<0] .= 0




scatter(t1, gr, label="inut pdata",title="Initial fit", xlabel="time", ylabel="growth rate")
plot!(t1, gr, label="fit")

function plot_int(fit, title, t1, gr) # must be plotted with plots.jl
    p1 = Plots.scatter(t1, gr, label="inut pdata",title=title, xlabel="time", ylabel="growth rate")
    Plots.plot!(fit, label="fit")
    return display(p1)
end

plot_int(LinearInterpolation(gr,t1), "Linear Interpolation", t1, gr)
plot_int(QuadraticInterpolation(gr,t1), "Quadratic Interpolation", t1, gr) # quadratic is best
plot_int(LagrangeInterpolation(gr,t1), "Lagrange Interpolation", t1, gr)
plot_int(ConstantInterpolation(gr,t1), "Constant Interpolation", t1, gr)
plot_int(QuadraticSpline(gr,t1), "Quadratic Spline", t1, gr)
plot_int(CubicSpline(gr,t1), "Cubic Spline", t1, gr)
plot_int(BSplineInterpolation(gr,t1,2,:ArcLen,:Average), "BSpline Interpolation", t1, gr)
plot_int(BSplineApprox(gr,t1,1,2,:ArcLen,:Average), "BSpline Approx", t1, gr)

print(gr)
lam = QuadraticInterpolation(gr,ts)
print(lam)
lam(21.616666666666667)
lam[132]
lam[67]
print(gr[1:66])
plot()
typeof(lam)

lam_t = []
for i in t1
    push!(lam_t, lam(i))
end
lam_t-gr



csv_wt = DataFrame(CSV.File("/home/holliehindley/phd/data/results_rtcOFF_grfit.csv"))
gr_wt = csv_wt."gr"
t1_wt = csv_wt."t"


Plots.scatter(t1_wt, gr_wt, label="input data",title="Initial fit", xlabel="time", ylabel="growth rate")
Plots.plot!(t1_wt, gr_wt, label="fit")

plot_int(LinearInterpolation(gr_wt,t1_wt), "Linear Interpolation RtcOFF data", t1_wt, gr_wt)
plot_int(QuadraticInterpolation(gr_wt, t1_wt), "QuadraticInterpolation RtcOFF data", t1_wt, gr_wt)
plot_int(LagrangeInterpolation(gr_wt,t1_wt), "Lagrange Interpolation RtcOFF data", t1_wt, gr_wt)
plot_int(ConstantInterpolation(gr_wt,t1_wt), "Constant Interpolation RtcOFF data", t1_wt, gr_wt)
plot_int(QuadraticSpline(gr_wt,t1_wt), "Quadratic Spline RtcOFF data", t1_wt, gr_wt)
plot_int(CubicSpline(gr_wt,t1_wt), "Cubic Spline RtcOFF data", t1_wt, gr_wt)
plot_int(BSplineInterpolation(gr_wt,t1_wt,2,:ArcLen,:Average), "BSpline Interpolation RtcOFF data", t1_wt, gr_wt)
plot_int(BSplineApprox(gr_wt,t1_wt,3,5,:ArcLen,:Average), "BSpline Approx RtcOFF data", t1_wt, gr_wt)