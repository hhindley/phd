include("$PATH/general_funcs/all_model_funcs.jl")
include("$PATH/rtc_model/parameters/rtc_params.jl")
include("$PATH/rtc_model/parameters/trna_params.jl")

indexof(sym,syms) = findfirst(isequal(sym),syms)
@variables t 
@parameters L c kr Vmax_init Km_init ω_ab ω_r θtscr g_max θtlr km_a km_b d krep kdam ktag kdeg kin atp na nb nr lam kc k_diss rh thr_t
species_trna1 = @syms rm_a(t) rtca(t) rm_b(t) rtcb(t) rm_r(t) rtcr(t) trna(t) rd(t) rt(t) 
species_trna = [Symbol(i) for i in species_trna1]

D = Differential(t)

@mtkmodel RTC_TRNA begin
    @parameters begin
        L 
        c 
        kr
        Vmax_init 
        Km_init 
        ω_ab  
        ω_r 
        θtscr
        g_max
        θtlr 
        km_a 
        km_b 
        d 
        krep 
        kdam 
        ktag 
        kdeg 
        kin 
        atp 
        na 
        nb 
        nr 
        lam 
        kc 
        k_diss 
        rh
        thr_t
    end
    @variables begin
        rm_a(t) 
        rtca(t) 
        rm_b(t) 
        rtcb(t) 
        rm_r(t) 
        rtcr(t) 
        trna(t) 
        rd(t) 
        rt(t)

        rhs_rm_a(t) 
        rhs_rtca(t) 
        rhs_rm_b(t) 
        rhs_rtcb(t) 
        rhs_rm_r(t) 
        rhs_rtcr(t) 
        rhs_trna(t) 
        rhs_rd(t) 
        rhs_rt(t) 

        alpha(t)
        fa(t)
        ra(t)
        Voc(t)
        sig_o(t)
        tscr_ab(t)
        tscr_r(t)
        tlr_el(t)
        Vrep(t)
        Vdam(t)
        Vinflux(t)
        Vtag(t)
        

    end

    @equations begin
        alpha ~ rt/kr # unitless
        fa ~ (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6) # unitless 
        ra ~ fa*rtcr # uM 
        
        # transcription
        Voc ~ Vmax_init*atp/(Km_init+atp) # uM min-1 
        sig_o ~ ra*Voc/k_diss # uM
    
        tscr_ab ~ sig_o*ω_ab*atp/(θtscr+atp) # uM min-1
        tscr_r ~ ω_r*atp/(θtscr+atp) # uM min-1

        tlr_el ~ (g_max*atp/(θtlr+atp)) * trna/(thr_t+trna) 

        # # ribosomes
        Vrep ~ rtcb*rt*krep/(rt+km_b) # uM min-1 
        Vdam ~ kdam*trna # uM min-1
        Vinflux ~ kin* g_max*atp/(θtlr+atp) #* trna/(thr_t+trna)  # uM min-1 
        Vtag ~ rtca*rd*ktag/(rd+km_a) # uM min-1 

        rhs_rm_a ~ tscr_ab - dil(rm_a,lam) - deg(rm_a)
        rhs_rm_b ~ tscr_ab - lam*(rm_b) - deg(rm_b)
        rhs_rm_r ~ tscr_r - lam*(rm_r) - deg(rm_r)

        rhs_rtca ~ tlr(rm_a, na, rh, tlr_el) - dil(rtca,lam)     
        rhs_rtcb ~ tlr(rm_b, nb, rh, tlr_el) - dil(rtcb,lam)
        rhs_rtcr ~ tlr(rm_r, nr, rh, tlr_el) - dil(rtcr,lam)

        rhs_trna ~ Vrep - Vdam + Vinflux - dil(trna,lam)
        rhs_rd ~ Vdam - Vtag - kdeg*rd - dil(rd,lam)
        rhs_rt ~ Vtag - Vrep - dil(rt,lam)

        D(rm_a) ~ rhs_rm_a
        D(rtca) ~ rhs_rtca
        D(rm_b) ~ rhs_rm_b
        D(rtcb) ~ rhs_rtcb 
        D(rm_r) ~ rhs_rm_r
        D(rtcr) ~ rhs_rtcr
        D(trna) ~ rhs_trna
        D(rd) ~ rhs_rd
        D(rt) ~ rhs_rt
    end
end


@mtkbuild rtc_trna_model = RTC_TRNA()

init_trna = [rtc_trna_model.rm_a=>0.0,rtc_trna_model.rtca=>0.0,rtc_trna_model.rm_b=>0.0,rtc_trna_model.rtcb=>0.0,rtc_trna_model.rm_r=>0.0,rtc_trna_model.rtcr=>0.0,rtc_trna_model.trna=>135.5,rtc_trna_model.rd=>0.0,rtc_trna_model.rt=>0.0] # tRNA initial conc = 135.5

tspan = (0, 1e9);

params_trna = Dict(L=>L_val, c=>c_val, kr=>kr_val*12, Vmax_init=>Vmax_init_val, Km_init=>Km_init_val, θtscr=>θtscr_val, θtlr=>θtlr_val, na=>nA_val, nb=>nB_val, nr=>nR_val, d=>d_val, krep=>krep_val, ktag=>ktag_val,
atp=>atp_val, km_a=>km_a_val, km_b=>km_b_val, g_max=>g_max_val, kdeg=>kdeg_trna_val, kin=>kin_trna_val, ω_ab=>ω_ab_val, ω_r=>ω_r_val, kdam=>kdam_val, lam=>lam_val, kc=>kc_val, k_diss=>k_diss_val, rh=>rh_val, thr_t=>thr_t_val)

ssvals_trna = steady_states(rtc_trna_model, init_trna, params_trna)

# function rtc_model_trna(initial, params, t) 
#     L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, kc, k_diss, rh, thr_t = params
#     rm_a, rtca, rm_b, rtcb, rm_r, rtcr, trna, rd, rt = initial


#     # dilution by growth and degradation 
#     dil(species) = lam*species
#     deg(species) = d*species
    
#     # MWC
#     alpha = rt/kr # unitless
#     fa = (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6) # unitless 
#     ra = fa*rtcr # uM 
    
#     # transcription
#     Voc = Vmax_init*atp/(Km_init+atp) # uM min-1 
#     sig_o = ra*Voc/k_diss # uM

#     tscr_el_a = ω_ab*atp/(θtscr+atp) # min-1
#     tscr_a = sig_o*tscr_el_a # uM min-1
#     tscr_el_b = ω_ab*atp/(θtscr+atp) # min-1
#     tscr_b = sig_o*tscr_el_b # uM min-1
#     tscr_r = ω_r*atp/(θtscr+atp) # uM min-1

#     # translation
#     tlr_el = (g_max*atp/(θtlr+atp)) * trna/(thr_t+trna) 
#     tlr(rm_x, nx) = (1/nx)*kc*rh*rm_x*tlr_el # uM min-1

#     # # ribosomes
#     Vrep = rtcb*rt*krep/(rt+km_b) # uM min-1 
#     Vdam = kdam*trna # uM min-1
#     Vinflux = kin*(g_max*atp/(θtlr+atp)) # uM min-1 
#     Vtag = rtca*rd*ktag/(rd+km_a) # uM min-1 


#     # ODEs
#     drm_a = tscr_a - dil(rm_a) - deg(rm_a)
#     drtca = tlr(rm_a, na) - dil(rtca)    
#     drm_b = tscr_b - dil(rm_b) - deg(rm_b)
#     drtcb = tlr(rm_b, nb) - dil(rtcb) 
#     drm_r = tscr_r - dil(rm_r) - deg(rm_r)
#     drtcr = tlr(rm_r, nr) - dil(rtcr)
#     dtrna = Vrep - Vdam + Vinflux - dil(trna)
#     drd = Vdam - Vtag - kdeg*rd - dil(rd)
#     drt = Vtag - Vrep - dil(rt)

#     # @SVector [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt]
#     [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, dtrna, drd, drt]

# end





# function rtc_mod_trna!(dz, z, p, t)
#     @unpack L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, kc, k_diss, rh, thr_t = p
#     rm_a, rtca, rm_b, rtcb, rm_r, rtcr, trna, rd, rt = z


#     # dilution by growth and degradation 
#     dil(species) = lam*species
#     deg(species) = d*species
    
#     # MWC
#     alpha = rt/kr # unitless
#     fa = (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6) # unitless 
#     ra = fa*rtcr # uM 
    
#     # transcription
#     Voc = Vmax_init*atp/(Km_init+atp) # uM min-1 
#     sig_o = ra*Voc/k_diss # uM

#     tscr_el_a = ω_ab*atp/(θtscr+atp) # min-1
#     tscr_a = sig_o*tscr_el_a # uM min-1
#     tscr_el_b = ω_ab*atp/(θtscr+atp) # min-1
#     tscr_b = sig_o*tscr_el_b # uM min-1
#     tscr_r = ω_r*atp/(θtscr+atp) # uM min-1

#     # translation
#     tlr_el = (g_max*atp/(θtlr+atp)) * trna/(thr_t+trna) 
#     tlr(rm_x, nx) = (1/nx)*kc*rh*rm_x*tlr_el # uM min-1

#     # # ribosomes
#     Vrep = rtcb*rt*krep/(rt+km_b) # uM min-1 
#     Vdam = kdam*trna # uM min-1
#     Vinflux = kin*(g_max*atp/(θtlr+atp)) # uM min-1 
#     Vtag = rtca*rd*ktag/(rd+km_a) # uM min-1 


#     # ODEs
#     dz[1] = tscr_a - dil(rm_a) - deg(rm_a)
#     dz[2] = tlr(rm_a, na) - dil(rtca)    
#     dz[3] = tscr_b - dil(rm_b) - deg(rm_b)
#     dz[4] = tlr(rm_b, nb) - dil(rtcb) 
#     dz[5] = tscr_r - dil(rm_r) - deg(rm_r)
#     dz[6] = tlr(rm_r, nr) - dil(rtcr)
#     dz[7] = Vrep - Vdam + Vinflux - dil(trna)
#     dz[8] = Vdam - Vtag - kdeg*rd - dil(rd)
#     dz[9] = Vtag - Vrep - dil(rt)
    
#     dz
    
# end

# rtc_mod_trna(z, p) = rtc_mod_trna!(similar(z), z, p, 0)
