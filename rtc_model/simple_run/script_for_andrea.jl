using CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, PlotlyJS, Statistics

function rtc_model1!(initial, params, t) 
    L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam = params
    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt = initial
    
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
    tlr(rm_x, nx) = (1/nx)*rh*rm_x*tlr_el 

    # ribosomes
    rtca1 = (atp*rtca)/(atp+(km_a*rd)) 
    rtcb1 = (atp*rtcb)/(atp+(km_b*rt)) 

    Vrep = krep*rtcb1*rt
    Vdam = kdam*rh
    Vinflux = kin*tlr_el #lam(t)*rh 
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

function sol(model, init, tspan, params)
    # params, init = choose_param_vector(model)
    prob = ODEProblem(model, init, tspan, params)
    solu = solve(prob, Rodas4())
    return solu
end

function extend_gr_curve(csv)
    mean_gr = mean((csv."gr"[Int64((length(csv."t")*2/3)+1):end]))
    df = DataFrame(t=Float64[], gr=Float64[])
    for t in collect(csv."t"[end]+10:5000:1e9)
        push!(df, [t, mean_gr])
    end    
    new_df = vcat(csv, df)
    lam = QuadraticInterpolation(new_df."gr",new_df."t")
    return lam, new_df
end 

function get_curve(sol, species)
    df = DataFrame(sol)
    if length(sol[1]) == 9
        rename!(df, [:time, :rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt])
    else
        rename!(df, [:time, :rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt, :atp])
    end
    species = df[:, species]
    return species
end

function plotly_plot_sol(sol, log, log1)
    rm_a = get_curve(sol, :rm_a); rm_b = get_curve(sol, :rm_b); rm_r = get_curve(sol, :rm_r); rtca = get_curve(sol, :rtca); rtcb = get_curve(sol, :rtcb); rtcr = get_curve(sol, :rtcr); rh = get_curve(sol, :rh); rt = get_curve(sol, :rt); rd = get_curve(sol, :rd);

    rma_curve = scatter(x=sol.t, y=rm_a, name="mRNA RtcA")
    rmb_curve = scatter(x=sol.t, y=rm_b, name="mRNA RtcB")
    rmr_curve = scatter(x=sol.t, y=rm_r, name="mRNA RtcR")
    rtca_curve = scatter(x=sol.t, y=rtca, name="RtcA")
    rtcb_curve = scatter(x=sol.t, y=rtcb, name="RtcB")
    rtcr_curve = scatter(x=sol.t, y=rtcr, name="RtcR")
    rh_curve = scatter(x=sol.t, y=rh, name="Rh")
    rt_curve = scatter(x=sol.t, y=rt, name="Rt")
    rd_curve = scatter(x=sol.t, y=rd, name="Rd")
    return (plot([rma_curve, rmb_curve, rmr_curve, rtca_curve, rtcb_curve, rtcr_curve, rh_curve, rt_curve, rd_curve] ,Layout(xaxis_type=log, yaxis_type=log1)))
end

L = 10; 
c = 0.001; 
kr = 0.125; 
Vmax_init = 39.51; 
Km_init = 250; 
θtscr = 160.01;  
θtlr = 255.73; 
k_b = 17.7; 
na = 338; 
nb = 408; 
nr = 532*6;
d = 0.2; 
krep = 137; 
ktag = 0.1; 
atp = 2500;
km_a = 20; 
km_b = 16;
g_max = 2.0923; 
kdeg = 0.001; 
kin = 0.054; 
ω_ab = 4;
ω_r = 0.0019*6;
ω_a = 4; 
ω_b = 4;
kdam = 0.000147;

rtca_0 = 0.00894; 
rtcb_0 = 0.0216; 
rh_0 = 11.29; 
rtcr_0 = 0.04; 
rm_a_0 = 0; 
rm_b_0 = 0; 
rm_r_0 = 0.04;
rd_0 = 0; 
rt_0 = 0;
initial = @SVector [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, rh_0, rd_0, rt_0];


csv = DataFrame(CSV.File("/home/holliehindley/phd/data/results_colD_grfit.csv")) 
csv = select!(csv, Not(["log(OD)", "log(OD) error", "gr error", "od"]))

lam_colD, new_df = extend_gr_curve(csv)
lam_colD[lam_colD.< 0] .= 0 

tspan = (0, lam_colD[end])
params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam_colD] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
solu = sol(rtc_model1!, initial, tspan, params)
p = plotly_plot_sol(solu, "log", "")



