include("/home/holliehindley/phd/general_funcs/all_model_funcs.jl")
include("/home/holliehindley/phd/rtc_model/parameters/rtc_params.jl")

indexof(sym, syms) = findfirst(isequal(sym),syms)

@variables t 
@parameters L c kr Vmax_init Km_init ω_ab ω_r θtscr g_max θtlr km_a km_b d krep kdam ktag kdeg kin atp na nb nr lam kc k_diss k_inhib1 k_inhib2 inhib
species_inhib1 = @syms rm_a(t) rtca(t) rm_b(t) rtcb(t) rm_r(t) rtcr(t) rh(t) rd(t) rt(t) rtc_i(t)
species_inhib = [Symbol(i) for i in species_inhib1]
    
D = Differential(t)

function build_inhib_model(inhib_protein1)
    inhib_protein = inhib_protein1
    @mtkmodel RTC_INHIB begin
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
            k_inhib1
            k_inhib2
            inhib
        end
        @variables begin
            rm_a(t) 
            rtca(t) 
            rm_b(t) 
            rtcb(t) 
            rm_r(t) 
            rtcr(t) 
            rh(t) 
            rd(t) 
            rt(t) 
            rtc_i(t)

            rhs_rm_a(t) 
            rhs_rtca(t) 
            rhs_rm_b(t) 
            rhs_rtcb(t) 
            rhs_rm_r(t) 
            rhs_rtcr(t) 
            rhs_rh(t) 
            rhs_rd(t) 
            rhs_rt(t)
            rhs_rtc_i(t) 

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

            tlr_el ~ g_max*atp/(θtlr+atp)

            # # ribosomes
            Vrep ~ rtcb*rt*krep/(rt+km_b) # uM min-1 
            Vdam ~ kdam*rh # uM min-1
            Vinflux ~ kin* g_max*atp/(θtlr+atp) # uM min-1 
            Vtag ~ rtca*rd*ktag/(rd+km_a) # uM min-1 

            rhs_rm_a ~ tscr_ab - dil(rm_a,lam) - deg(rm_a)
            rhs_rm_b ~ tscr_ab - dil(rm_b,lam) - deg(rm_b)
            rhs_rm_r ~ tscr_r - dil(rm_r,lam) - deg(rm_r)

            rhs_rh ~ Vrep - Vdam + Vinflux - dil(rh,lam)
            rhs_rd ~ Vdam - Vtag - kdeg*rd - dil(rd,lam)
            rhs_rt ~ Vtag - Vrep - dil(rt,lam)

            if inhib_protein == :rtca_inhib
                rhs_rtc_i ~ k_inhib1*rtca*inhib - k_inhib2*rtc_i - dil(rtc_i,lam)            
                rhs_rtca ~ tlr(rm_a, na, rh, tlr_el) - dil(rtca,lam) - k_inhib1*rtca*inhib + k_inhib2*rtc_i   
                rhs_rtcb ~ tlr(rm_b, nb, rh, tlr_el) - dil(rtcb,lam)
                rhs_rtcr ~ tlr(rm_r, nr, rh, tlr_el) - dil(rtcr,lam)
            elseif inhib_protein == :rtcb_inhib
                rhs_rtc_i ~ k_inhib1*rtcb*inhib - k_inhib2*rtc_i - dil(rtc_i,lam)            
                rhs_rtca ~ tlr(rm_a, na, rh, tlr_el) - dil(rtca,lam)
                rhs_rtcb ~ tlr(rm_b, nb, rh, tlr_el) - dil(rtcb,lam) - k_inhib1*rtcb*inhib + k_inhib2*rtc_i   
                rhs_rtcr ~ tlr(rm_r, nr, rh, tlr_el) - dil(rtcr,lam)
            else 
                rhs_rtc_i ~ k_inhib1*rtcr*inhib - k_inhib2*rtc_i - dil(rtc_i,lam)            
                rhs_rtca ~ tlr(rm_a, na, rh, tlr_el) - dil(rtca,lam)
                rhs_rtcb ~ tlr(rm_b, nb, rh, tlr_el) - dil(rtcb,lam) 
                rhs_rtcr ~ tlr(rm_r, nr, rh, tlr_el) - dil(rtcr,lam) - k_inhib1*rtcr*inhib + k_inhib2*rtc_i   
            end

            D(rm_a) ~ rhs_rm_a
            D(rtca) ~ rhs_rtca
            D(rm_b) ~ rhs_rm_b
            D(rtcb) ~ rhs_rtcb 
            D(rm_r) ~ rhs_rm_r
            D(rtcr) ~ rhs_rtcr
            D(rh) ~ rhs_rh
            D(rd) ~ rhs_rd
            D(rt) ~ rhs_rt
            D(rtc_i) ~ rhs_rtc_i
        end
    end

    return @mtkbuild rtc_inhib_model = RTC_INHIB()
end 

rtca_inhib_model = build_inhib_model(:rtca_inhib)
rtcb_inhib_model = build_inhib_model(:rtcb_inhib)
rtcr_inhib_model = build_inhib_model(:rtcr_inhib)


params_inhib = Dict(L=>L_val, c=>c_val, kr=>kr_val, Vmax_init=>Vmax_init_val, Km_init=>Km_init_val, θtscr=>θtscr_val, θtlr=>θtlr_val, na=>nA_val, nb=>nB_val, nr=>nR_val, d=>d_val, krep=>krep_val, ktag=>ktag_val,
atp=>atp_val, km_a=>km_a_val, km_b=>km_b_val, g_max=>g_max_val, kdeg=>kdeg_val, kin=>kin_val, ω_ab=>ω_ab_val, ω_r=>ω_r_val, kdam=>kdam_val, lam=>lam_val, kc=>kc_val, k_diss=>k_diss_val, k_inhib1=>k_inhib1_val, k_inhib2=>k_inhib2_val, inhib=>inhib_val)

init_inhib_rtca = [rtca_inhib_model.rm_a=>0.0,rtca_inhib_model.rtca=>0.0,rtca_inhib_model.rm_b=>0.0,rtca_inhib_model.rtcb=>0.0,rtca_inhib_model.rm_r=>0.0,rtca_inhib_model.rtcr=>0.0,rtca_inhib_model.rh=>11.29,rtca_inhib_model.rd=>0.0,rtca_inhib_model.rt=>0.0,rtca_inhib_model.rtc_i=>0.0]
init_inhib_rtcb = [rtcb_inhib_model.rm_a=>0.0,rtcb_inhib_model.rtca=>0.0,rtcb_inhib_model.rm_b=>0.0,rtcb_inhib_model.rtcb=>0.0,rtcb_inhib_model.rm_r=>0.0,rtcb_inhib_model.rtcr=>0.0,rtcb_inhib_model.rh=>11.29,rtcb_inhib_model.rd=>0.0,rtcb_inhib_model.rt=>0.0,rtcb_inhib_model.rtc_i=>0.0]
init_inhib_rtcr = [rtcr_inhib_model.rm_a=>0.0,rtcr_inhib_model.rtca=>0.0,rtcr_inhib_model.rm_b=>0.0,rtcr_inhib_model.rtcb=>0.0,rtcr_inhib_model.rm_r=>0.0,rtcr_inhib_model.rtcr=>0.0,rtcr_inhib_model.rh=>11.29,rtcr_inhib_model.rd=>0.0,rtcr_inhib_model.rt=>0.0,rtcr_inhib_model.rtc_i=>0.0]

ssvals_rtca_inhib = steady_states(rtca_inhib_model, init_inhib_rtca, params_inhib)
ssvals_rtcb_inhib = steady_states(rtcb_inhib_model, init_inhib_rtcb, params_inhib)
ssvals_rtcr_inhib = steady_states(rtcr_inhib_model, init_inhib_rtcr, params_inhib)











# function rtc_inhib_model_rtcb(initial, params, t) 
#     L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, kc, k_inhib1, k_inhib2, inhib = params
#     rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt, rtc_i = initial


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
#     tlr_el = g_max*atp/(θtlr+atp) # aa min-1 molec-1
#     tlr(rm_x, nx) = (1/nx)*kc*rh*rm_x*tlr_el # uM min-1

#     # # ribosomes
#     Vrep = rtcb*rt*krep/(rt+km_b) # uM min-1 
#     Vdam = kdam*rh # uM min-1
#     Vinflux = kin*tlr_el # uM min-1 
#     Vtag = rtca*rd*ktag/(rd+km_a) # uM min-1 

#     # ODEs
#     drm_a = tscr_a - dil(rm_a) - deg(rm_a)
#     drtca = tlr(rm_a, na) - dil(rtca)    
#     drm_b = tscr_b - dil(rm_b) - deg(rm_b)
#     drtcb = tlr(rm_b, nb) - dil(rtcb) - k_inhib1*rtcb*inhib + k_inhib2*rtc_i
#     drm_r = tscr_r - dil(rm_r) - deg(rm_r)
#     drtcr = tlr(rm_r, nr) - dil(rtcr)
#     drh = Vrep - Vdam + Vinflux - dil(rh)
#     drd = Vdam - Vtag - kdeg*rd - dil(rd)
#     drt = Vtag - Vrep - dil(rt)
#     drtc_i = k_inhib1*rtcb*inhib - k_inhib2*rtc_i - dil(rtc_i)
#     # @SVector [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt]
#     [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt, drtc_i]

# end


# function rtc_inhib_mod_rtcb!(dz, z, p, t)
#     @unpack L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, kc, k_inhib1, k_inhib2, inhib = p
#     rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt, rtc_i = z


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
#     tlr_el = g_max*atp/(θtlr+atp) # aa min-1 molec-1
#     tlr(rm_x, nx) = (1/nx)*kc*rh*rm_x*tlr_el # uM min-1

#     # # ribosomes
#     Vrep = rtcb*rt*krep/(rt+km_b) # uM min-1 
#     Vdam = kdam*rh # uM min-1
#     Vinflux = kin*tlr_el # uM min-1 
#     Vtag = rtca*rd*ktag/(rd+km_a) # uM min-1 

#     # ODEs
#     dz[1] = tscr_a - dil(rm_a) - deg(rm_a)
#     dz[2] = tlr(rm_a, na) - dil(rtca)    
#     dz[3] = tscr_b - dil(rm_b) - deg(rm_b)
#     dz[4] = tlr(rm_b, nb) - dil(rtcb) - k_inhib1*rtcb*inhib + k_inhib2*rtc_i
#     dz[5] = tscr_r - dil(rm_r) - deg(rm_r)
#     dz[6] = tlr(rm_r, nr) - dil(rtcr)
#     dz[7] = Vrep - Vdam + Vinflux - dil(rh)
#     dz[8] = Vdam - Vtag - kdeg*rd - dil(rd)
#     dz[9] = Vtag - Vrep - dil(rt)
#     dz[10] = k_inhib1*rtcb*inhib - k_inhib2*rtc_i - dil(rtc_i)

#     dz
# end


# rtc_inhib_mod_rtcb(z, p) = rtc_inhib_mod_rtcb!(similar(z), z, p, 0)








# function rtc_inhib_model_rtca(initial, params, t) 
#     L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, kc, k_inhib1, k_inhib2, inhib = params
#     rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt, rtc_i = initial


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
#     tlr_el = g_max*atp/(θtlr+atp) # aa min-1 molec-1
#     tlr(rm_x, nx) = (1/nx)*kc*rh*rm_x*tlr_el # uM min-1

#     # # ribosomes
#     Vrep = rtcb*rt*krep/(rt+km_b) # uM min-1 
#     Vdam = kdam*rh # uM min-1
#     Vinflux = kin*tlr_el # uM min-1 
#     Vtag = rtca*rd*ktag/(rd+km_a) # uM min-1 

#     # ODEs
#     drm_a = tscr_a - dil(rm_a) - deg(rm_a)
#     drtca = tlr(rm_a, na) - dil(rtca) - k_inhib1*rtca*inhib + k_inhib2*rtc_i
#     drm_b = tscr_b - dil(rm_b) - deg(rm_b)
#     drtcb = tlr(rm_b, nb) - dil(rtcb) 
#     drm_r = tscr_r - dil(rm_r) - deg(rm_r)
#     drtcr = tlr(rm_r, nr) - dil(rtcr)
#     drh = Vrep - Vdam + Vinflux - dil(rh)
#     drd = Vdam - Vtag - kdeg*rd - dil(rd)
#     drt = Vtag - Vrep - dil(rt)
#     drtc_i = k_inhib1*rtca*inhib - k_inhib2*rtc_i - dil(rtc_i)
#     # @SVector [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt]
#     [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt, drtc_i]

# end


# function rtc_inhib_mod_rtca!(dz, z, p, t)
#     @unpack L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, kc, k_inhib1, k_inhib2, inhib = p
#     rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt, rtc_i = z


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
#     tlr_el = g_max*atp/(θtlr+atp) # aa min-1 molec-1
#     tlr(rm_x, nx) = (1/nx)*kc*rh*rm_x*tlr_el # uM min-1

#     # # ribosomes
#     Vrep = rtcb*rt*krep/(rt+km_b) # uM min-1 
#     Vdam = kdam*rh # uM min-1
#     Vinflux = kin*tlr_el # uM min-1 
#     Vtag = rtca*rd*ktag/(rd+km_a) # uM min-1 

#     # ODEs
#     dz[1] = tscr_a - dil(rm_a) - deg(rm_a)
#     dz[2] = tlr(rm_a, na) - dil(rtca) - k_inhib1*rtca*inhib + k_inhib2*rtc_i
#     dz[3] = tscr_b - dil(rm_b) - deg(rm_b)
#     dz[4] = tlr(rm_b, nb) - dil(rtcb) 
#     dz[5] = tscr_r - dil(rm_r) - deg(rm_r)
#     dz[6] = tlr(rm_r, nr) - dil(rtcr)
#     dz[7] = Vrep - Vdam + Vinflux - dil(rh)
#     dz[8] = Vdam - Vtag - kdeg*rd - dil(rd)
#     dz[9] = Vtag - Vrep - dil(rt)
#     dz[10] = k_inhib1*rtca*inhib - k_inhib2*rtc_i - dil(rtc_i)
#     dz
# end


# rtc_inhib_mod_rtca(z, p) = rtc_inhib_mod_rtca!(similar(z), z, p, 0)





# function rtc_inhib_model_rtcr(initial, params, t) 
#     L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, kc, k_inhib1, k_inhib2, inhib = params
#     rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt, rtc_i = initial


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
#     tlr_el = g_max*atp/(θtlr+atp) # aa min-1 molec-1
#     tlr(rm_x, nx) = (1/nx)*kc*rh*rm_x*tlr_el # uM min-1

#     # # ribosomes
#     Vrep = rtcb*rt*krep/(rt+km_b) # uM min-1 
#     Vdam = kdam*rh # uM min-1
#     Vinflux = kin*tlr_el # uM min-1 
#     Vtag = rtca*rd*ktag/(rd+km_a) # uM min-1 

#     # ODEs
#     drm_a = tscr_a - dil(rm_a) - deg(rm_a)
#     drtca = tlr(rm_a, na) - dil(rtca)    
#     drm_b = tscr_b - dil(rm_b) - deg(rm_b)
#     drtcb = tlr(rm_b, nb) - dil(rtcb) 
#     drm_r = tscr_r - dil(rm_r) - deg(rm_r)
#     drtcr = tlr(rm_r, nr) - dil(rtcr) - k_inhib1*rtcr*inhib + k_inhib2*rtc_i
#     drh = Vrep - Vdam + Vinflux - dil(rh)
#     drd = Vdam - Vtag - kdeg*rd - dil(rd)
#     drt = Vtag - Vrep - dil(rt)
#     drtc_i = k_inhib1*rtcr*inhib - k_inhib2*rtc_i - dil(rtc_i)
#     # @SVector [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt]
#     [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt, drtc_i]

# end


# function rtc_inhib_mod_rtcr!(dz, z, p, t)
#     @unpack L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, kc, k_inhib1, k_inhib2, inhib = p
#     rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt, rtc_i = z


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
#     tlr_el = g_max*atp/(θtlr+atp) # aa min-1 molec-1
#     tlr(rm_x, nx) = (1/nx)*kc*rh*rm_x*tlr_el # uM min-1

#     # # ribosomes
#     Vrep = rtcb*rt*krep/(rt+km_b) # uM min-1 
#     Vdam = kdam*rh # uM min-1
#     Vinflux = kin*tlr_el # uM min-1 
#     Vtag = rtca*rd*ktag/(rd+km_a) # uM min-1 

#     # ODEs
#     dz[1] = tscr_a - dil(rm_a) - deg(rm_a)
#     dz[2] = tlr(rm_a, na) - dil(rtca)    
#     dz[3] = tscr_b - dil(rm_b) - deg(rm_b)
#     dz[4] = tlr(rm_b, nb) - dil(rtcb)
#     dz[5] = tscr_r - dil(rm_r) - deg(rm_r)
#     dz[6] = tlr(rm_r, nr) - dil(rtcr) - k_inhib1*rtcr*inhib + k_inhib2*rtc_i
#     dz[7] = Vrep - Vdam + Vinflux - dil(rh)
#     dz[8] = Vdam - Vtag - kdeg*rd - dil(rd)
#     dz[9] = Vtag - Vrep - dil(rt)
#     dz[10] = k_inhib1*rtcr*inhib - k_inhib2*rtc_i - dil(rtc_i)

#     dz
# end


# rtc_inhib_mod_rtcr(z, p) = rtc_inhib_mod_rtcr!(similar(z), z, p, 0)











# function rtc_inhib_model_rt(initial, params, t) 
#     L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, kc, k_inhib1, k_inhib2, inhib = params
#     rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt, rtc_i = initial


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
#     tlr_el = g_max*atp/(θtlr+atp) # aa min-1 molec-1
#     tlr(rm_x, nx) = (1/nx)*kc*rh*rm_x*tlr_el # uM min-1

#     # # ribosomes
#     Vrep = rtcb*rt*krep/(rt+km_b) # uM min-1 
#     Vdam = kdam*rh # uM min-1
#     Vinflux = kin*tlr_el # uM min-1 
#     Vtag = rtca*rd*ktag/(rd+km_a) # uM min-1 

#     # ODEs
#     drm_a = tscr_a - dil(rm_a) - deg(rm_a)
#     drtca = tlr(rm_a, na) - dil(rtca)    
#     drm_b = tscr_b - dil(rm_b) - deg(rm_b)
#     drtcb = tlr(rm_b, nb) - dil(rtcb) 
#     drm_r = tscr_r - dil(rm_r) - deg(rm_r)
#     drtcr = tlr(rm_r, nr) - dil(rtcr) 
#     drh = Vrep - Vdam + Vinflux - dil(rh)
#     drd = Vdam - Vtag - kdeg*rd - dil(rd)
#     drt = Vtag - Vrep - dil(rt) - k_inhib1*rt*inhib + k_inhib2*rtc_i
#     drtc_i = k_inhib1*rt*inhib - k_inhib2*rtc_i - dil(rtc_i)
#     # @SVector [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt]
#     [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt, drtc_i]

# end


# function rtc_inhib_mod_rt!(dz, z, p, t)
#     @unpack L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, kc, k_inhib1, k_inhib2, inhib = p
#     rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt, rtc_i = z


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
#     tlr_el = g_max*atp/(θtlr+atp) # aa min-1 molec-1
#     tlr(rm_x, nx) = (1/nx)*kc*rh*rm_x*tlr_el # uM min-1

#     # # ribosomes
#     Vrep = rtcb*rt*krep/(rt+km_b) # uM min-1 
#     Vdam = kdam*rh # uM min-1
#     Vinflux = kin*tlr_el # uM min-1 
#     Vtag = rtca*rd*ktag/(rd+km_a) # uM min-1 

#     # ODEs
#     dz[1] = tscr_a - dil(rm_a) - deg(rm_a)
#     dz[2] = tlr(rm_a, na) - dil(rtca)    
#     dz[3] = tscr_b - dil(rm_b) - deg(rm_b)
#     dz[4] = tlr(rm_b, nb) - dil(rtcb)
#     dz[5] = tscr_r - dil(rm_r) - deg(rm_r)
#     dz[6] = tlr(rm_r, nr) - dil(rtcr)
#     dz[7] = Vrep - Vdam + Vinflux - dil(rh)
#     dz[8] = Vdam - Vtag - kdeg*rd - dil(rd)
#     dz[9] = Vtag - Vrep - dil(rt) - k_inhib1*rt*inhib + k_inhib2*rtc_i
#     dz[10] =k_inhib1*rt*inhib - k_inhib2*rtc_i - dil(rtc_i)

#     dz
# end


# rtc_inhib_mod_rt(z, p) = rtc_inhib_mod_rt!(similar(z), z, p, 0)










# function rtc_inhib_model_rtcab(initial, params, t) 
#     L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, kc, k_inhib1, k_inhib2, inhib = params
#     rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt, rtca_i, rtcb_i = initial


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
#     tlr_el = g_max*atp/(θtlr+atp) # aa min-1 molec-1
#     tlr(rm_x, nx) = (1/nx)*kc*rh*rm_x*tlr_el # uM min-1

#     # # ribosomes
#     Vrep = rtcb*rt*krep/(rt+km_b) # uM min-1 
#     Vdam = kdam*rh # uM min-1
#     Vinflux = kin*tlr_el # uM min-1 
#     Vtag = rtca*rd*ktag/(rd+km_a) # uM min-1 

#     # ODEs
#     drm_a = tscr_a - dil(rm_a) - deg(rm_a)
#     drtca = tlr(rm_a, na) - dil(rtca) - k_inhib1*rtca*inhib + k_inhib2*rtca_i
#     drm_b = tscr_b - dil(rm_b) - deg(rm_b)
#     drtcb = tlr(rm_b, nb) - dil(rtcb) - k_inhib1*rtcb*inhib + k_inhib2*rtcb_i
#     drm_r = tscr_r - dil(rm_r) - deg(rm_r)
#     drtcr = tlr(rm_r, nr) - dil(rtcr)
#     drh = Vrep - Vdam + Vinflux - dil(rh)
#     drd = Vdam - Vtag - kdeg*rd - dil(rd)
#     drt = Vtag - Vrep - dil(rt)
#     drtca_i = k_inhib1*rtca*inhib - k_inhib2*rtca_i - dil(rtca_i)
#     drtcb_i = k_inhib1*rtcb*inhib - k_inhib2*rtcb_i - dil(rtcb_i)
#     # @SVector [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt]
#     [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt, drtca_i, drtcb_i]

# end


# function rtc_inhib_mod_rtcab!(dz, z, p, t)
#     @unpack L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, kc, k_inhib1, k_inhib2, inhib = p
#     rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt, rtca_i, rtcb_i = z


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
#     tlr_el = g_max*atp/(θtlr+atp) # aa min-1 molec-1
#     tlr(rm_x, nx) = (1/nx)*kc*rh*rm_x*tlr_el # uM min-1

#     # # ribosomes
#     Vrep = rtcb*rt*krep/(rt+km_b) # uM min-1 
#     Vdam = kdam*rh # uM min-1
#     Vinflux = kin*tlr_el # uM min-1 
#     Vtag = rtca*rd*ktag/(rd+km_a) # uM min-1 

#     # ODEs
#     dz[1] = tscr_a - dil(rm_a) - deg(rm_a)
#     dz[2] = tlr(rm_a, na) - dil(rtca) - k_inhib1*rtca*inhib + k_inhib2*rtca_i
#     dz[3] = tscr_b - dil(rm_b) - deg(rm_b)
#     dz[4] = tlr(rm_b, nb) - dil(rtcb) - k_inhib1*rtcb*inhib + k_inhib2*rtcb_i
#     dz[5] = tscr_r - dil(rm_r) - deg(rm_r)
#     dz[6] = tlr(rm_r, nr) - dil(rtcr)
#     dz[7] = Vrep - Vdam + Vinflux - dil(rh)
#     dz[8] = Vdam - Vtag - kdeg*rd - dil(rd)
#     dz[9] = Vtag - Vrep - dil(rt)
#     dz[10] = k_inhib1*rtca*inhib - k_inhib2*rtca_i - dil(rtca_i)
#     dz[11] = k_inhib1*rtcb*inhib - k_inhib2*rtcb_i - dil(rtcb_i)

#     dz
# end


# rtc_inhib_mod_rtcab(z, p) = rtc_inhib_mod_rtcab!(similar(z), z, p, 0)












# function rtc_inhib_model_rtcbr(initial, params, t) 
#     L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, kc, k_inhib1, k_inhib2, inhib = params
#     rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt, rtcr_i, rtcb_i = initial


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
#     tlr_el = g_max*atp/(θtlr+atp) # aa min-1 molec-1
#     tlr(rm_x, nx) = (1/nx)*kc*rh*rm_x*tlr_el # uM min-1

#     # # ribosomes
#     Vrep = rtcb*rt*krep/(rt+km_b) # uM min-1 
#     Vdam = kdam*rh # uM min-1
#     Vinflux = kin*tlr_el # uM min-1 
#     Vtag = rtca*rd*ktag/(rd+km_a) # uM min-1 

#     # ODEs
#     drm_a = tscr_a - dil(rm_a) - deg(rm_a)
#     drtca = tlr(rm_a, na) - dil(rtca) 
#     drm_b = tscr_b - dil(rm_b) - deg(rm_b)
#     drtcb = tlr(rm_b, nb) - dil(rtcb) - k_inhib1*rtcb*inhib + k_inhib2*rtcb_i
#     drm_r = tscr_r - dil(rm_r) - deg(rm_r)
#     drtcr = tlr(rm_r, nr) - dil(rtcr)- k_inhib1*rtcr*inhib + k_inhib2*rtcr_i
#     drh = Vrep - Vdam + Vinflux - dil(rh)
#     drd = Vdam - Vtag - kdeg*rd - dil(rd)
#     drt = Vtag - Vrep - dil(rt)
#     drtcr_i = k_inhib1*rtcr*inhib - k_inhib2*rtcr_i - dil(rtcr_i)
#     drtcb_i = k_inhib1*rtcb*inhib - k_inhib2*rtcb_i - dil(rtcb_i)
#     # @SVector [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt]
#     [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt, drtcr_i, drtcb_i]

# end


# function rtc_inhib_mod_rtcbr!(dz, z, p, t)
#     @unpack L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, kc, k_inhib1, k_inhib2, inhib = p
#     rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt, rtcr_i, rtcb_i = z


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
#     tlr_el = g_max*atp/(θtlr+atp) # aa min-1 molec-1
#     tlr(rm_x, nx) = (1/nx)*kc*rh*rm_x*tlr_el # uM min-1

#     # # ribosomes
#     Vrep = rtcb*rt*krep/(rt+km_b) # uM min-1 
#     Vdam = kdam*rh # uM min-1
#     Vinflux = kin*tlr_el # uM min-1 
#     Vtag = rtca*rd*ktag/(rd+km_a) # uM min-1 

#     # ODEs
#     dz[1] = tscr_a - dil(rm_a) - deg(rm_a)
#     dz[2] = tlr(rm_a, na) - dil(rtca)
#     dz[3] = tscr_b - dil(rm_b) - deg(rm_b)
#     dz[4] = tlr(rm_b, nb) - dil(rtcb) - k_inhib1*rtcb*inhib + k_inhib2*rtcb_i
#     dz[5] = tscr_r - dil(rm_r) - deg(rm_r)
#     dz[6] = tlr(rm_r, nr) - dil(rtcr) - k_inhib1*rtcr*inhib + k_inhib2*rtcr_i
#     dz[7] = Vrep - Vdam + Vinflux - dil(rh)
#     dz[8] = Vdam - Vtag - kdeg*rd - dil(rd)
#     dz[9] = Vtag - Vrep - dil(rt)
#     dz[10] = k_inhib1*rtcr*inhib - k_inhib2*rtcr_i - dil(rtcr_i)
#     dz[11] = k_inhib1*rtcb*inhib - k_inhib2*rtcb_i - dil(rtcb_i)

#     dz
# end


# rtc_inhib_mod_rtcbr(z, p) = rtc_inhib_mod_rtcbr!(similar(z), z, p, 0)














# function rtc_inhib_model_rtcar(initial, params, t) 
#     L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, kc, k_inhib1, k_inhib2, inhib = params
#     rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt, rtcr_i, rtca_i = initial


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
#     tlr_el = g_max*atp/(θtlr+atp) # aa min-1 molec-1
#     tlr(rm_x, nx) = (1/nx)*kc*rh*rm_x*tlr_el # uM min-1

#     # # ribosomes
#     Vrep = rtcb*rt*krep/(rt+km_b) # uM min-1 
#     Vdam = kdam*rh # uM min-1
#     Vinflux = kin*tlr_el # uM min-1 
#     Vtag = rtca*rd*ktag/(rd+km_a) # uM min-1 

#     # ODEs
#     drm_a = tscr_a - dil(rm_a) - deg(rm_a)
#     drtca = tlr(rm_a, na) - dil(rtca) - k_inhib1*rtca*inhib + k_inhib2*rtca_i
#     drm_b = tscr_b - dil(rm_b) - deg(rm_b)
#     drtcb = tlr(rm_b, nb) - dil(rtcb) 
#     drm_r = tscr_r - dil(rm_r) - deg(rm_r)
#     drtcr = tlr(rm_r, nr) - dil(rtcr)- k_inhib1*rtcr*inhib + k_inhib2*rtcr_i
#     drh = Vrep - Vdam + Vinflux - dil(rh)
#     drd = Vdam - Vtag - kdeg*rd - dil(rd)
#     drt = Vtag - Vrep - dil(rt)
#     drtcr_i = k_inhib1*rtcr*inhib - k_inhib2*rtcr_i - dil(rtcr_i)
#     drtca_i = k_inhib1*rtca*inhib - k_inhib2*rtca_i - dil(rtca_i)
#     # @SVector [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt]
#     [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt, drtcr_i, drtca_i]

# end


# function rtc_inhib_mod_rtcar!(dz, z, p, t)
#     @unpack L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, kc, k_inhib1, k_inhib2, inhib = p
#     rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt, rtcr_i, rtca_i = z


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
#     tlr_el = g_max*atp/(θtlr+atp) # aa min-1 molec-1
#     tlr(rm_x, nx) = (1/nx)*kc*rh*rm_x*tlr_el # uM min-1

#     # # ribosomes
#     Vrep = rtcb*rt*krep/(rt+km_b) # uM min-1 
#     Vdam = kdam*rh # uM min-1
#     Vinflux = kin*tlr_el # uM min-1 
#     Vtag = rtca*rd*ktag/(rd+km_a) # uM min-1 

#     # ODEs
#     dz[1] = tscr_a - dil(rm_a) - deg(rm_a)
#     dz[2] = tlr(rm_a, na) - dil(rtca)- k_inhib1*rtca*inhib + k_inhib2*rtca_i
#     dz[3] = tscr_b - dil(rm_b) - deg(rm_b)
#     dz[4] = tlr(rm_b, nb) - dil(rtcb) 
#     dz[5] = tscr_r - dil(rm_r) - deg(rm_r)
#     dz[6] = tlr(rm_r, nr) - dil(rtcr) - k_inhib1*rtcr*inhib + k_inhib2*rtcr_i
#     dz[7] = Vrep - Vdam + Vinflux - dil(rh)
#     dz[8] = Vdam - Vtag - kdeg*rd - dil(rd)
#     dz[9] = Vtag - Vrep - dil(rt)
#     dz[10] = k_inhib1*rtcr*inhib - k_inhib2*rtcr_i - dil(rtcr_i)
#     dz[11] = k_inhib1*rtca*inhib - k_inhib2*rtca_i - dil(rtca_i)

#     dz
# end


# rtc_inhib_mod_rtcar(z, p) = rtc_inhib_mod_rtcar!(similar(z), z, p, 0)










# function rtc_inhib_model_rtcat(initial, params, t) 
#     L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, kc, k_inhib1, k_inhib2, inhib = params
#     rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt, rt_i, rtca_i = initial


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
#     tlr_el = g_max*atp/(θtlr+atp) # aa min-1 molec-1
#     tlr(rm_x, nx) = (1/nx)*kc*rh*rm_x*tlr_el # uM min-1

#     # # ribosomes
#     Vrep = rtcb*rt*krep/(rt+km_b) # uM min-1 
#     Vdam = kdam*rh # uM min-1
#     Vinflux = kin*tlr_el # uM min-1 
#     Vtag = rtca*rd*ktag/(rd+km_a) # uM min-1 

#     # ODEs
#     drm_a = tscr_a - dil(rm_a) - deg(rm_a)
#     drtca = tlr(rm_a, na) - dil(rtca) - k_inhib1*rtca*inhib + k_inhib2*rtca_i
#     drm_b = tscr_b - dil(rm_b) - deg(rm_b)
#     drtcb = tlr(rm_b, nb) - dil(rtcb) 
#     drm_r = tscr_r - dil(rm_r) - deg(rm_r)
#     drtcr = tlr(rm_r, nr) - dil(rtcr)
#     drh = Vrep - Vdam + Vinflux - dil(rh)
#     drd = Vdam - Vtag - kdeg*rd - dil(rd)
#     drt = Vtag - Vrep - dil(rt) - k_inhib1*rt*inhib + k_inhib2*rt_i
#     drt_i = k_inhib1*rt*inhib - k_inhib2*rt_i - dil(rt_i)
#     drtca_i = k_inhib1*rtca*inhib - k_inhib2*rtca_i - dil(rtca_i)
#     # @SVector [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt]
#     [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt, drt_i, drtca_i]

# end


# function rtc_inhib_mod_rtcat!(dz, z, p, t)
#     @unpack L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, kc, k_inhib1, k_inhib2, inhib = p
#     rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt, rt_i, rtca_i = z


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
#     tlr_el = g_max*atp/(θtlr+atp) # aa min-1 molec-1
#     tlr(rm_x, nx) = (1/nx)*kc*rh*rm_x*tlr_el # uM min-1

#     # # ribosomes
#     Vrep = rtcb*rt*krep/(rt+km_b) # uM min-1 
#     Vdam = kdam*rh # uM min-1
#     Vinflux = kin*tlr_el # uM min-1 
#     Vtag = rtca*rd*ktag/(rd+km_a) # uM min-1 

#     # ODEs
#     dz[1] = tscr_a - dil(rm_a) - deg(rm_a)
#     dz[2] = tlr(rm_a, na) - dil(rtca)- k_inhib1*rtca*inhib + k_inhib2*rtca_i
#     dz[3] = tscr_b - dil(rm_b) - deg(rm_b)
#     dz[4] = tlr(rm_b, nb) - dil(rtcb) 
#     dz[5] = tscr_r - dil(rm_r) - deg(rm_r)
#     dz[6] = tlr(rm_r, nr) - dil(rtcr) 
#     dz[7] = Vrep - Vdam + Vinflux - dil(rh)
#     dz[8] = Vdam - Vtag - kdeg*rd - dil(rd)
#     dz[9] = Vtag - Vrep - dil(rt)- k_inhib1*rt*inhib + k_inhib2*rt_i
#     dz[10] = k_inhib1*rt*inhib - k_inhib2*rt_i - dil(rt_i)
#     dz[11] = k_inhib1*rtca*inhib - k_inhib2*rtca_i - dil(rtca_i)

#     dz
# end


# rtc_inhib_mod_rtcat(z, p) = rtc_inhib_mod_rtcat!(similar(z), z, p, 0)






# function rtc_inhib_model_rtcbt(initial, params, t) 
#     L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, kc, k_inhib1, k_inhib2, inhib = params
#     rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt, rt_i, rtcb_i = initial


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
#     tlr_el = g_max*atp/(θtlr+atp) # aa min-1 molec-1
#     tlr(rm_x, nx) = (1/nx)*kc*rh*rm_x*tlr_el # uM min-1

#     # # ribosomes
#     Vrep = rtcb*rt*krep/(rt+km_b) # uM min-1 
#     Vdam = kdam*rh # uM min-1
#     Vinflux = kin*tlr_el # uM min-1 
#     Vtag = rtca*rd*ktag/(rd+km_a) # uM min-1 

#     # ODEs
#     drm_a = tscr_a - dil(rm_a) - deg(rm_a)
#     drtca = tlr(rm_a, na) - dil(rtca) 
#     drm_b = tscr_b - dil(rm_b) - deg(rm_b)
#     drtcb = tlr(rm_b, nb) - dil(rtcb) - k_inhib1*rtcb*inhib + k_inhib2*rtcb_i
#     drm_r = tscr_r - dil(rm_r) - deg(rm_r)
#     drtcr = tlr(rm_r, nr) - dil(rtcr)
#     drh = Vrep - Vdam + Vinflux - dil(rh)
#     drd = Vdam - Vtag - kdeg*rd - dil(rd)
#     drt = Vtag - Vrep - dil(rt) - k_inhib1*rt*inhib + k_inhib2*rt_i
#     drt_i = k_inhib1*rt*inhib - k_inhib2*rt_i - dil(rt_i)
#     drtcb_i = k_inhib1*rtcb*inhib - k_inhib2*rtcb_i - dil(rtcb_i)
#     # @SVector [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt]
#     [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt, drt_i, drtcb_i]

# end


# function rtc_inhib_mod_rtcbt!(dz, z, p, t)
#     @unpack L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, kc, k_inhib1, k_inhib2, inhib = p
#     rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt, rt_i, rtcb_i = z


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
#     tlr_el = g_max*atp/(θtlr+atp) # aa min-1 molec-1
#     tlr(rm_x, nx) = (1/nx)*kc*rh*rm_x*tlr_el # uM min-1

#     # # ribosomes
#     Vrep = rtcb*rt*krep/(rt+km_b) # uM min-1 
#     Vdam = kdam*rh # uM min-1
#     Vinflux = kin*tlr_el # uM min-1 
#     Vtag = rtca*rd*ktag/(rd+km_a) # uM min-1 

#     # ODEs
#     dz[1] = tscr_a - dil(rm_a) - deg(rm_a)
#     dz[2] = tlr(rm_a, na) - dil(rtca)
#     dz[3] = tscr_b - dil(rm_b) - deg(rm_b)
#     dz[4] = tlr(rm_b, nb) - dil(rtcb) - k_inhib1*rtcb*inhib + k_inhib2*rtcb_i
#     dz[5] = tscr_r - dil(rm_r) - deg(rm_r)
#     dz[6] = tlr(rm_r, nr) - dil(rtcr) 
#     dz[7] = Vrep - Vdam + Vinflux - dil(rh)
#     dz[8] = Vdam - Vtag - kdeg*rd - dil(rd)
#     dz[9] = Vtag - Vrep - dil(rt)- k_inhib1*rt*inhib + k_inhib2*rt_i
#     dz[10] = k_inhib1*rt*inhib - k_inhib2*rt_i - dil(rt_i)
#     dz[11] = k_inhib1*rtcb*inhib - k_inhib2*rtcb_i - dil(rtcb_i)

#     dz
# end


# rtc_inhib_mod_rtcbt(z, p) = rtc_inhib_mod_rtcbt!(similar(z), z, p, 0)










# function rtc_inhib_model_rtcrt(initial, params, t) 
#     L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, kc, k_inhib1, k_inhib2, inhib = params
#     rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt, rt_i, rtc_i = initial


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
#     tlr_el = g_max*atp/(θtlr+atp) # aa min-1 molec-1
#     tlr(rm_x, nx) = (1/nx)*kc*rh*rm_x*tlr_el # uM min-1

#     # # ribosomes
#     Vrep = rtcb*rt*krep/(rt+km_b) # uM min-1 
#     Vdam = kdam*rh # uM min-1
#     Vinflux = kin*tlr_el # uM min-1 
#     Vtag = rtca*rd*ktag/(rd+km_a) # uM min-1 

#     # ODEs
#     drm_a = tscr_a - dil(rm_a) - deg(rm_a)
#     drtca = tlr(rm_a, na) - dil(rtca) 
#     drm_b = tscr_b - dil(rm_b) - deg(rm_b)
#     drtcb = tlr(rm_b, nb) - dil(rtcb) 
#     drm_r = tscr_r - dil(rm_r) - deg(rm_r)
#     drtcr = tlr(rm_r, nr) - dil(rtcr) - k_inhib1*rtcr*inhib + k_inhib2*rtc_i
#     drh = Vrep - Vdam + Vinflux - dil(rh)
#     drd = Vdam - Vtag - kdeg*rd - dil(rd)
#     drt = Vtag - Vrep - dil(rt) - k_inhib1*rt*inhib + k_inhib2*rt_i
#     drt_i = k_inhib1*rt*inhib - k_inhib2*rt_i - dil(rt_i)
#     drtcr_i = k_inhib1*rtcr*inhib - k_inhib2*rtc_i - dil(rtc_i)
#     # @SVector [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt]
#     [drm_a, drtca, drm_b, drtcb, drm_r, drtcr, drh, drd, drt, drt_i, drtc_i]

# end


# function rtc_inhib_mod_rtcrt!(dz, z, p, t)
#     @unpack L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam, kc, k_inhib1, k_inhib2, inhib = p
#     rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rh, rd, rt, rt_i, rtc_i = z


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
#     tlr_el = g_max*atp/(θtlr+atp) # aa min-1 molec-1
#     tlr(rm_x, nx) = (1/nx)*kc*rh*rm_x*tlr_el # uM min-1

#     # # ribosomes
#     Vrep = rtcb*rt*krep/(rt+km_b) # uM min-1 
#     Vdam = kdam*rh # uM min-1
#     Vinflux = kin*tlr_el # uM min-1 
#     Vtag = rtca*rd*ktag/(rd+km_a) # uM min-1 

#     # ODEs
#     dz[1] = tscr_a - dil(rm_a) - deg(rm_a)
#     dz[2] = tlr(rm_a, na) - dil(rtca)
#     dz[3] = tscr_b - dil(rm_b) - deg(rm_b)
#     dz[4] = tlr(rm_b, nb) - dil(rtcb) 
#     dz[5] = tscr_r - dil(rm_r) - deg(rm_r)
#     dz[6] = tlr(rm_r, nr) - dil(rtcr) - k_inhib1*rtcr*inhib + k_inhib2*rtc_i
#     dz[7] = Vrep - Vdam + Vinflux - dil(rh)
#     dz[8] = Vdam - Vtag - kdeg*rd - dil(rd)
#     dz[9] = Vtag - Vrep - dil(rt)- k_inhib1*rt*inhib + k_inhib2*rt_i
#     dz[10] = k_inhib1*rt*inhib - k_inhib2*rt_i - dil(rt_i)
#     dz[11] = k_inhib1*rtcr*inhib - k_inhib2*rtc_i - dil(rtc_i)

#     dz
# end


# rtc_inhib_mod_rtcrt(z, p) = rtc_inhib_mod_rtcrt!(similar(z), z, p, 0)











