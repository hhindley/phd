include("$PATH/general_funcs/all_model_funcs.jl")
include("$PATH/rtc_model/parameters/rtc_params.jl")

indexof(sym, syms) = findfirst(isequal(sym),syms)

@variables t 
@parameters L c kr Vmax_init Km_init ω_ab ω_r θtscr g_max θtlr km_a km_b d krep kdam ktag kdeg kin atp na nb nr lam kc k_diss 

D = Differential(t)

@mtkmodel RTC_LAM begin
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

        rhs_rm_a(t) 
        rhs_rtca(t) 
        rhs_rm_b(t) 
        rhs_rtcb(t) 
        rhs_rm_r(t) 
        rhs_rtcr(t) 
        rhs_rh(t) 
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
        lam_t(t)
        

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

        lam_t ~ 0.01*rh

        rhs_rm_a ~ tscr_ab - dil(rm_a,lam_t) - deg(rm_a)
        rhs_rtca ~ tlr(rm_a, na, rh, tlr_el) - dil(rtca,lam_t)     
        rhs_rm_b ~ tscr_ab - dil(rm_b,lam_t) - deg(rm_b)
        rhs_rtcb ~ tlr(rm_b, nb, rh, tlr_el) - dil(rtcb,lam_t)
        rhs_rm_r ~ tscr_r - dil(rm_r,lam_t) - deg(rm_r)
        rhs_rtcr ~ tlr(rm_r, nr, rh, tlr_el) - dil(rtcr,lam_t)
        rhs_rh ~ Vrep - Vdam + Vinflux - dil(rh,lam_t)
        rhs_rd ~ Vdam - Vtag - kdeg*rd - dil(rd,lam_t)
        rhs_rt ~ Vtag - Vrep - dil(rt,lam_t)

        D(rm_a) ~ rhs_rm_a
        D(rtca) ~ rhs_rtca
        D(rm_b) ~ rhs_rm_b
        D(rtcb) ~ rhs_rtcb 
        D(rm_r) ~ rhs_rm_r
        D(rtcr) ~ rhs_rtcr
        D(rh) ~ rhs_rh
        D(rd) ~ rhs_rd
        D(rt) ~ rhs_rt
    end
end

@mtkbuild rtc_model_lam = RTC_LAM()

ssvals_rtc_lam = steady_states(rtc_model_lam, init_rtc, params_rtc)



