using CSV, DataFrames, DataInterpolations, Plots# PlotlyJS
include("/home/holliehindley/phd/rtc_models/Oct2022_model/rtc_model.jl")
include("/home/holliehindley/phd/rtc_models/sol_species_funcs.jl")
include("/home/holliehindley/phd/rtc_models/params_init_tspan.jl")


csv = DataFrame(CSV.File("/home/holliehindley/phd/data/results_colD_grfit.csv")) # read csv to a dataframe
gr = csv."gr"
t1 = csv."t"*60

function plot_int(fit, title) # must be plotted with plots.jl
    p1 = scatter(t1, gr, label="input data",title=title, xlabel="time", ylabel="growth rate")
    plot!(fit, label="fit")
    return display(p1)
end

plot_int(lam, "Initial fit")
plot_int(LinearInterpolation(gr,ts), "Linear Interpolation")
plot_int(QuadraticInterpolation(gr,ts), "Quadratic Interpolation") # quadratic is best
plot_int(LagrangeInterpolation(gr,ts), "Lagrange Interpolation")
plot_int(ConstantInterpolation(gr,ts), "Constant Interpolation")
plot_int(QuadraticSpline(gr,ts), "Quadratic Spline")
plot_int(CubicSpline(gr,ts), "Cubic Spline")
plot_int(BSplineInterpolation(gr,ts,2,:ArcLen,:Average), "BSpline Interpolation")
plot_int(BSplineApprox(gr,ts,1,2,:ArcLen,:Average), "BSpline Approx")

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