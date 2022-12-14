using CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, PlotlyJS
include("/home/holliehindley/phd/rtc_models/Oct2022_model/rtc_model.jl")
include("/home/holliehindley/phd/rtc_models/sol_species_funcs.jl")
include("/home/holliehindley/phd/rtc_models/params_init_tspan.jl")

csv = DataFrame(CSV.File("/home/holliehindley/phd/data/results_colD_grfit.csv")) # read csv to a dataframe
gr = csv."gr"
t1 = csv."t"

lam_colD = QuadraticInterpolation(gr,t1)
lam_colD(222)

sv_wt = DataFrame(CSV.File("/home/holliehindley/phd/data/results_rtcOFF_grfit.csv"))
gr_wt = csv_wt."gr"
t1_wt = csv_wt."t"

lam_wt = LinearInterpolation(gr_wt, t1_wt)
@show lam_wt

function rtc_model1!(initial, params, t) 
    L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam = params
    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt = initial
    # growth rate
    # lam = csv."gr"
    
    # dilution by growth and degradation 
    dil(species) = lam(t)*species
    deg(species) = d*species
    # @show lam(t)
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

    # @show (lam(t)), t

    @SVector [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt]
end

params = [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam_colD]
initial = [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, rh_0, rd_0, rt_0];

tspan = (0,1297)
prob = ODEProblem(rtc_model1!, initial, tspan, params)
solu = solve(prob, Rodas4())#, saveat=ts)


plotly_plot_sol(solu)

res = get_all_curves(solu, all_species)

dil_dict = OrderedDict(name => [] for name in all_species)
for (dict, species) in zip(values(dil_dict), all_species)
    for (t, i) in zip(solu.t, collect(1:length(solu.t)))
        push!(dict, lam_colD(t)*res[species][1][i])
    end
end

plot_from_dict(dil_dict, solu)


lambda = []
for t in solu.t
    push!(lambda, lam_colD(t))
end
alpha = res[:rt][1]/kr
fa = @. (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6)
ra = @. fa*res[:rtcr][1]
Vinit = @. ra*Vmax_init*atp/(Km_init+atp)
tscr_el_a = ω_ab*atp/(θtscr+atp)
tscr_a = Vinit*tscr_el_a
tscr_el_b = ω_ab*atp/(θtscr+atp)
tscr_b = Vinit*tscr_el_b
tscr_r = ω_r*atp/(θtscr+atp)
tlr_el = g_max*atp/(θtlr+atp)
tlr(rm_x, nx) = @. (1/nx)*res[:rh][1]*rm_x*tlr_el
tlr_r = tlr(res[:rm_r][1], nr); tlr_a = tlr(res[:rm_a][1], na); tlr_b = tlr(res[:rm_b][1], nb);
rtca1 = @. (atp*res[:rtca][1])/(atp+(km_a*res[:rd][1])) 
rtcb1 = @. (atp*res[:rtcb][1])/(atp+(km_b*res[:rt][1])) 
Vrep = @. krep*rtcb1*res[:rt][1]
Vdam = @. kdam*res[:rh][1]
Vinflux = kin*tlr_el
Vtag = @. ktag*rtca1*res[:rd][1]


p1 = plot(scatter(x=solu.t, y=alpha), Layout(title="alpha"));
p2 = plot(scatter(x=solu.t, y=fa), Layout(title="fa"));
p3 = plot(scatter(x=solu.t, y=ra), Layout(title="ra"));
p4 = plot(scatter(x=solu.t, y=Vinit), Layout(title="Vinit"));
p5 = plot(scatter(x=solu.t, y=tscr_a), Layout(title="tscr_a"));
p6 = plot(scatter(x=solu.t, y=tscr_b), Layout(title="tscr_b"));
p7 = plot(scatter(x=solu.t, y=tlr_r), Layout(title="tlr_r"));
p8 = plot(scatter(x=solu.t, y=tlr_a), Layout(title="tlr_a"));
p9 = plot(scatter(x=solu.t, y=tlr_b), Layout(title="tlr_b"));
p10 = plot(scatter(x=solu.t, y=rtca1), Layout(title="rtca1"));
p11 = plot(scatter(x=solu.t, y=rtcb1), Layout(title="rtcb1"));
p12 = plot(scatter(x=solu.t, y=Vrep), Layout(title="Vrep"));
p13 = plot(scatter(x=solu.t, y=Vdam), Layout(title="Vdam"));
p14 = plot(scatter(x=solu.t, y=Vtag), Layout(title="Vtag"));
p15 = plot(scatter(x=solu.t, y=lambda), Layout(title="λ"))

[p1 p2 p3; p4 p5 p6; p7 p8 p9; p10 p11 p12; p13 p14 p15]