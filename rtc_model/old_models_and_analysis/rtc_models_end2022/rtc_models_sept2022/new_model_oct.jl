
function rtc_model_new(initial, params, t) # change KM in params back to k_a
    L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θ, max, thr, km, k_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp = params
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
    rtca1 = (atp*rtca)/(atp+(km*rd)) 
    # rtca1 = k_a*atp*rtca/(k_a*atp+ktag*rd)
    # rtca1 = (atp*rtca)/(atp+(ktag/k_a)*rd) 
    rtcb1 = k_b*atp*rtcb/(k_b*atp+krep*rt)
    # rtcb1 = (atp*rtcb)/(atp+(krep/k_b)*rt) 

    Vrep = krep*rtcb1*rt
    Vdam = kdam*rh
    Vinflux = kin*tlr_el
    Vtag = ktag*rtca1*rd


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
θ = 20; max = 4; thr = 20; gr_c = 0.01; k_a = 10e-5; k_b = 10e-5;
d = 0.01; krep = 10;  kdam = 0.05; ktag = 10; kdeg = 0.001; kin = 0.4; atp = 10;
km = 5;
# L = 1000; c = 0.001; kr = 125; Vmax_init = 39.51; Km_init = 250; ω_ab = 4.14; ω_r = 4.14; 
# θ = 4.38; max = 1260; thr = 7; gr_c = 0.01; k_a = 0.18; k_b = 17.7;
# d = 5; krep = 137;  kdam = 0.05; ktag = 9680; kdeg = 0.001; kin = 0.4; atp = 2500;
# KM = 5;


rm_a_0 = 0; rtca_0 = 1; rm_b_0 = 0; rtcb_0 = 1; rm_r_0 = 0; rtcr_0 = 0;
rh_0 = 10; rd_0 = 0; rt_0 = 0;



params = @SVector [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θ, max, thr, km, k_b, gr_c, d, krep, kdam, ktag, kdeg, kin, atp]
init = @SVector [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, rh_0, rd_0, rt_0]



tspan = (0, 1e9)

