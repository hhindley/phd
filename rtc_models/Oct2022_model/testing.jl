using CSV, DataFrames, DifferentialEquations, StaticArrays, BenchmarkTools, OrderedCollections, DataInterpolations, Plots# PlotlyJS
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

function rtc_model1!(initial, params, t) 
    L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam = params
    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt = initial
    # growth rate
    # lam = csv."gr"
    
    # dilution by growth and degradation 
    dil(species) = lam(t)*species
    deg(species) = d*species
    
    # MWC
    alpha = rt/kr 
    fa = (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6)
    ra = fa*rtcr
    
    # transcription
    Vinit = ra*Vmax_init*atp/(Km_init+atp)
    tscr_el_a = ω_ab*atp/(θtscr+atp)
    tscr_a = Vinit*tscr_el_a
    tscr_el_b = ω_ab*atp/(θtscr+atp)
    tscr_b = Vinit*tscr_el_b
    tscr_r = ω_r*atp/(θtscr+atp)

    # translation
    tlr_el = g_max*atp/(θtlr+atp)
    tlr(rm_x, nx) = (1/nx)*rh*rm_x*tlr_el # *1/nx nx = length of RtcA, RtcB and RtcR

    # ribosomes
    rtca1 = (atp*rtca)/(atp+(km_a*rd)) 
    rtcb1 = (atp*rtcb)/(atp+(km_b*rt)) 

    # rtcb1 = k_b*atp*rtcb/(k_b*atp+krep*rt)

    Vrep = krep*rtcb1*rt
    Vdam = kdam*rh
    Vinflux = kin*tlr_el
    Vtag = ktag*rtca1*rd


    # ODEs
    drm_a = tscr_a - dil(rm_a) - deg(rm_a)
    drtca = tlr(rm_a, na) - dil(rtca)    
    drm_b = tscr_b - dil(rm_b) - deg(rm_b)
    drtcb = tlr(rm_b, nb) - dil(rtcb)
    drm_r = tscr_r - dil(rm_r) - deg(rm_r)
    drtcr = tlr(rm_r, nr) - dil(rtcr)
    drh = Vrep - Vdam + Vinflux - dil(rh)
    drd = Vdam - Vtag - kdeg*rd - dil(rd)
    drt = Vtag - Vrep - dil(rt)

    # @show (params[23])

    @SVector [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt]
end

params = [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam]
initial = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, rh_0, rd_0, rt_0];

ts = csv."t"*60
tspan = (0, 1297);

prob = ODEProblem(rtc_model1!, initial, tspan, params)
solu = solve(prob, Rodas4())#, saveat=ts)
plot(solu)

plotly_plot_sol(solu)



