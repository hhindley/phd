using DifferentialEquations, StaticArrays, BenchmarkTools, DataFrames #, PlotlyJS

function rtc_model(initial, params, t)
    L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θ, max, thr, K1_tag, K2_tag, K1_rep, K2_rep, gr_c, d, krep, kdam, ktag, kdeg, kin, atp = params
    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt = initial

    # growth rate
    lam = gr_c*rh
    
    # dilution by growth and degradation 
    dil(species) = lam*species
    deg(species) = d*species
    
    # MWC
    alpha = rt/kr 
    fa = (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6)
    ra = fa*rtcr
    
    # transcription 
    Vinit = ra*Vmax_init*atp/(Km_init+atp)
    tscr_el_ab = ω_ab*atp/(θ+atp)
    tscr_ab = Vinit*tscr_el_ab
    tscr_r = ω_r*atp/(θ+atp)

    # translation
    tlr_el = max*atp/(thr+atp)
    tlr(rm_x) = rh*rm_x*tlr_el

    # ribosomes
    rdrtca = rtca*rd/(K1_tag+rd+K2_tag*atp)
    rtrtcb = rtcb*rt/(K1_rep+rt+K2_rep*atp)
    Vrep = krep*rtrtcb*atp
    Vdam = kdam*rh
    Vinflux = kin*tlr_el
    Vtag = ktag*rdrtca*atp

    # ODEs
    drm_a = tscr_ab - dil(rm_a) - deg(rm_a)
    drtca = tlr(rm_a) - dil(rtca)    
    drm_b = tscr_ab - dil(rm_b) - deg(rm_b)
    drtcb = tlr(rm_b) - dil(rtcb)
    drm_r = tscr_r - dil(rm_r) - deg(rm_r)
    drtcr = tlr(rm_r) - dil(rtcr)
    drh = Vrep - Vdam + Vinflux - dil(rh)
    drd = Vdam - Vtag - kdeg*rd - dil(rd)
    drt = Vtag - Vrep - dil(rt)

    @SVector [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt]
end

species = [:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt]

L = 100; c = 0.01; kr = 10; Vmax_init = 5; Km_init = 55.829; ω_ab = 4.14; ω_r = 4.14; 
θ = 20; max = 4; thr = 20; K1_tag = 1; K2_tag = 10; K1_rep = 1; K2_rep = 10; gr_c = 0.01; 
d = 0.01; krep = 1;  kdam = 0.05; ktag = 1; kdeg = 0.001; kin = 0.4; atp = 10;

rm_a_0 = 0; rtca_0 = 1; rm_b_0 = 0; rtcb_0 = 1; rm_r_0 = 0; rtcr_0 = 0;
rh_0 = 10; rd_0 = 0; rt_0 = 0;



params = @SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θ, max, thr, K1_tag, K2_tag, K1_rep, K2_rep, gr_c, d, krep, kdam, ktag, kdeg, kin, atp]
init = @SVector [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, rh_0, rd_0, rt_0]

tspan = (0, 1e9)



# prob = ODEProblem(rtc_model, init, tspan, params)
# solu = solve(prob, Rodas4())
# Plots.plot(solu[2:end], xaxis=(:log10, (1,Inf)), yaxis=(:log10, (1,Inf)), ylabel="[species]", labels=["rm_a" "rtca" "rm_b" "rtcb" "rm_r" "rtcr" "rh" "rd" "rt"])

# function run()
#     for i in collect(1:10000)
#         sol = solve(prob, Rodas4())
#     end
# end

# @time run()

# df = DataFrame(sol)
# rename!(df, [:time, :rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt])

# print(df)

# plot(df, x=:time, y=:rtcr)

# data = AbstractTrace[]
# for col in names(df)[2:end]
#     push!(data, scatter(x=:time, y=df[!,col], mode="lines", name=col))
# end

# layout = Layout(yaxis_type="log10", xaxis_type="log10")

# plot(data)


