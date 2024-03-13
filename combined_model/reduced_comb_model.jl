include("$PATH/general_funcs/all_model_funcs.jl")
indexof(sym, syms) = findfirst(isequal(sym),syms)

@variables t
@parameters d kr L c w_rh w_t w_m w_q w_BA w_R θ_rh θ_nr Kq nq Vmax_init Km_init nrh nx nR nA nB gmax Kgamma M kb ku vt vm s0 Kt Km km_a km_b krep kdam ktag kdeg k_diss ns
species_comb1 = @syms m_rh(t) m_t(t) m_m(t) m_q(t) m_R(t) m_A(t) m_B(t) c_rh(t) c_t(t) c_m(t) c_q(t) c_R(t) c_A(t) c_B(t) et(t) em(t) q(t) R(t) A(t) B(t) rh(t) rt(t) rd(t) si(t) a(t) 
species_red_comb = [Symbol(i) for i in species_comb1]

D = Differential(t)

@mtkmodel REDUCED_COMBINED_MODEL begin
    @parameters begin
        d 
        kr 
        L 
        c 
        w_rh 
        w_t 
        w_m 
        w_q 
        w_BA 
        w_R 
        θ_rh 
        θ_nr 
        Kq 
        nq 
        Vmax_init 
        Km_init 
        nrh 
        nx 
        nR 
        nA 
        nB 
        gmax 
        Kgamma 
        M 
        kb 
        ku 
        vt 
        vm 
        s0 
        Kt 
        Km 
        km_a 
        km_b 
        krep 
        kdam
        ktag 
        kdeg 
        k_diss 
        ns
    end
    @variables begin
        m_rh(t) 
        m_t(t) 
        m_m(t) 
        m_q(t) 
        m_R(t) 
        m_A(t) 
        m_B(t) 
        c_rh(t) 
        c_t(t) 
        c_m(t) 
        c_q(t) 
        c_R(t) 
        c_A(t) 
        c_B(t) 
        et(t) 
        em(t) 
        q(t) 
        R(t) 
        A(t) 
        B(t) 
        rh(t) 
        rt(t) 
        rd(t) 
        si(t) 
        a(t)

        rhs_m_rh(t) 
        rhs_m_t(t) 
        rhs_m_m(t) 
        rhs_m_q(t) 
        rhs_m_R(t) 
        rhs_m_A(t) 
        rhs_m_B(t) 
        rhs_c_rh(t) 
        rhs_c_t(t) 
        rhs_c_m(t) 
        rhs_c_q(t) 
        rhs_c_R(t) 
        rhs_c_A(t) 
        rhs_c_B(t) 
        rhs_et(t) 
        rhs_em(t) 
        rhs_q(t) 
        rhs_R(t) 
        rhs_A(t) 
        rhs_B(t) 
        rhs_rh(t) 
        rhs_rt(t) 
        rhs_rd(t)  
        rhs_si(t) 
        rhs_a(t)
        
        alpha(t)
        fa(t)
        ra(t)
        Voc(t)
        sig_o(t)

        gamma(t)
        ttrate(t)
        lam(t)

        vimp(t)
        vcat(t)

        Vrep(t)
        Vdam(t)
        Vtag(t)
    end
    @equations begin
        # MWC
        alpha ~ rt/kr # unitless 
        fa ~ (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6) # unitless
        ra ~ fa*R # uM 

        # transcription 
        Voc ~ Vmax_init*a/(Km_init+a) # uM min-1 
        sig_o ~ ra*Voc/k_diss

        # translation 
        gamma ~ gmax*a/(Kgamma+a) # aa uM-1 min-1 
        ttrate ~ (c_q+c_rh+c_t+c_m+c_R+c_A+c_B)*gamma # aa min-1

        # growth rate 
        lam ~ ttrate/M # min-1

        # metabolism and import 
        vimp ~ (et*vt*s0/(Kt+s0)) # uM min-1 
        vcat ~ (em*vm*si/(Km+si)) # uM min-1 

        # rtc
        Vrep ~ B*rt*krep/(rt+km_b) # uM min-1
        Vdam ~ rh*kdam # uM min-1
        Vtag ~ A*rd*ktag/(rd+km_a) # uM min-1

        # ODEs
        # mRNA
        rhs_m_rh ~ ω_p(w_rh, θ_rh, a) - dil(m_rh, lam) - deg(m_rh) + v_x(c_rh, nrh, gamma) - rh_bind(m_rh, rh) + rh_unbind(c_rh)
        rhs_m_t ~ ω_p(w_t, θ_nr, a) - dil(m_t, lam) - deg(m_t) + v_x(c_t, nx, gamma) - rh_bind(m_t, rh) + rh_unbind(c_t)
        rhs_m_m ~ ω_p(w_m, θ_nr, a) - dil(m_m, lam) - deg(m_m) + v_x(c_m, nx, gamma) - rh_bind(m_m, rh) + rh_unbind(c_m)
        rhs_m_q ~ ω_q(w_q, θ_nr, a, q) - dil(m_q, lam) - deg(m_q) + v_x(c_q, nx, gamma) - rh_bind(m_q, rh) + rh_unbind(c_q) 
        rhs_m_R ~ ω_p(w_R, θ_nr, a) - dil(m_R, lam) - deg(m_R) + v_x(c_R, nR, gamma) - rh_bind(m_R, rh) + rh_unbind(c_R) 
        rhs_m_A ~ ω_rtcBA(θ_nr, a, sig_o) - dil(m_A, lam) - deg(m_A) + v_x(c_A, nA, gamma) - rh_bind(m_A, rh) + rh_unbind(c_A) 
        rhs_m_B ~ ω_rtcBA(θ_nr, a, sig_o) - dil(m_B, lam) - deg(m_B) + v_x(c_B, nB, gamma) - rh_bind(m_B, rh) + rh_unbind(c_B) 

        # mRNA:ribosome complexes 
        rhs_c_rh ~ rh_bind(m_rh, rh) - dil(c_rh, lam) - rh_unbind(c_rh) - v_x(c_rh, nrh, gamma)
        rhs_c_t ~ rh_bind(m_t, rh) - dil(c_t, lam) - rh_unbind(c_t) - v_x(c_t, nx, gamma)
        rhs_c_m ~ rh_bind(m_m, rh) - dil(c_m, lam) - rh_unbind(c_m) - v_x(c_m, nx, gamma)
        rhs_c_q ~ rh_bind(m_q, rh) - dil(c_q, lam) - rh_unbind(c_q) - v_x(c_q, nx, gamma)
        rhs_c_R ~ rh_bind(m_R, rh) - dil(c_R, lam) - rh_unbind(c_R) - v_x(c_R, nR, gamma)
        rhs_c_A ~ rh_bind(m_A, rh) - dil(c_A, lam) - rh_unbind(c_A) - v_x(c_A, nA, gamma)
        rhs_c_B ~ rh_bind(m_B, rh) - dil(c_B, lam) - rh_unbind(c_B) - v_x(c_B, nB, gamma)

        # proteins 
        rhs_et ~ v_x(c_t, nx, gamma) - dil(et, lam)
        rhs_em ~ v_x(c_m, nx, gamma) - dil(em, lam)
        rhs_q ~ v_x(c_q, nx, gamma) - dil(q, lam)
        rhs_R ~ v_x(c_R, nR, gamma) - dil(R, lam)
        rhs_A ~ v_x(c_A, nA, gamma) - dil(A, lam)
        rhs_B ~ v_x(c_B, nB, gamma) - dil(B, lam)

        # ribosomes
        rhs_rh ~ v_x(c_rh, nrh, gamma) - dil(rh, lam) - (rh_bind(m_rh, rh) + rh_bind(m_t, rh) + rh_bind(m_m, rh) + rh_bind(m_q, rh) + rh_bind(m_R, rh) + rh_bind(m_A, rh) + rh_bind(m_B, rh)) + (rh_unbind(c_rh) + rh_unbind(c_t) + rh_unbind(c_m) + rh_unbind(c_q) + rh_unbind(c_R) + rh_unbind(c_A) + rh_unbind(c_B)) + (v_x(c_rh, nrh, gamma) + v_x(c_t, nx, gamma) + v_x(c_m, nx, gamma) + v_x(c_q, nx, gamma) + v_x(c_R, nR, gamma) + v_x(c_A, nA, gamma) + v_x(c_B, nB, gamma)) + Vrep - Vdam
        rhs_rt ~ Vtag - Vrep - dil(rt, lam)
        rhs_rd ~ Vdam - Vtag - rd*kdeg - dil(rd, lam)

        # nutrient
        rhs_si ~ vimp - vcat - dil(si, lam)

        # energy
        rhs_a ~ ns*vcat - ttrate - dil(a, lam)

        D(m_rh) ~ rhs_m_rh
        D(m_t) ~ rhs_m_t
        D(m_m) ~ rhs_m_m
        D(m_q) ~ rhs_m_q
        D(m_R) ~ rhs_m_R
        D(m_A) ~ rhs_m_A
        D(m_B) ~ rhs_m_B
        D(c_rh) ~ rhs_c_rh
        D(c_t) ~ rhs_c_t
        D(c_m) ~ rhs_c_m
        D(c_q) ~ rhs_c_q
        D(c_R) ~ rhs_c_R
        D(c_A) ~ rhs_c_A
        D(c_B) ~ rhs_c_B
        D(et) ~ rhs_et
        D(em) ~ rhs_em
        D(q) ~ rhs_q
        D(R) ~ rhs_R
        D(A) ~ rhs_A
        D(B) ~ rhs_B
        D(rh) ~ rhs_rh
        D(rt) ~ rhs_rt
        D(rd) ~ rhs_rd
        D(si) ~ rhs_si
        D(a) ~ rhs_a
    end
end

@mtkbuild reduced_combined_model = REDUCED_COMBINED_MODEL()

init_red_comb = [reduced_combined_model.m_rh=>0.0,reduced_combined_model.m_t=>0.0,reduced_combined_model.m_m=>0.0,reduced_combined_model.m_q=>0.0,reduced_combined_model.m_R=>0.0,reduced_combined_model.m_A=>0.0,reduced_combined_model.m_B=>0.0,reduced_combined_model.c_rh=>0.0,reduced_combined_model.c_t=>0.0,reduced_combined_model.c_m=>0.0,reduced_combined_model.c_q=>0.0,reduced_combined_model.c_R=>0.0,reduced_combined_model.c_A=>0.0,reduced_combined_model.c_B=>0.0,reduced_combined_model.et=>0.0,reduced_combined_model.em=>0.0,reduced_combined_model.q=>0.0,reduced_combined_model.R=>0.0,reduced_combined_model.A=>0.0,reduced_combined_model.B=>0.0,reduced_combined_model.rh=>0.0166,reduced_combined_model.rt=>0.0,reduced_combined_model.rd=>0.0,reduced_combined_model.si=>0.0,reduced_combined_model.a=>1.66]

params_red_comb = Dict(
    # rtc params
    d=>d_val, kr=>kr_val, L=>L_val, c=>c_val, Vmax_init=>Vmax_init_val, Km_init=>Km_init_val,
    w_BA=>ω_ab_val, w_R=>ω_r_val, nR=>nR_val, nA=>nA_val, nB=>nB_val, km_a=>km_a_val, km_b=>km_b_val, 
    krep=>krep_val, kdam=>0, ktag=>ktag_val, kdeg=>kdeg_val, k_diss=>k_diss_val,
    # growth model params
    w_rh=>wr_uM_val, w_t=>we_uM_val, w_m=>we_uM_val, w_q=>wq_uM_val, gmax=>gmax_val, 
    θ_rh=>thetar_uM_val, θ_nr=>thetax_uM_val, Kq=>Kq_uM_val, Kt=>Kt_uM_val, Km=>Km_uM_val, Kgamma=>Kgamma_uM_val, 
    nq=>nq_val, nrh=>nr_val, nx=>nx_val, 
    kb=>kb_uM_val, ku=>ku_val, 
    s0=>s0_uM_val, M=>M_uM_val, ns=>ns_val, vt=>vt_val, vm=>vm_val)


params_red_comb_new = Dict(
    # rtc params
    d=>d_val, kr=>kr_val, L=>L_val, c=>c_val, Vmax_init=>Vmax_init_val, Km_init=>Km_init_val,
    w_BA=>ω_ab_val_comb_new, w_R=>ω_r_val_comb_new, nR=>nR_val, nA=>nA_val, nB=>nB_val, km_a=>km_a_val, km_b=>km_b_val, 
    krep=>krep_val, kdam=>0, ktag=>ktag_val, kdeg=>kdeg_val, k_diss=>k_diss_val,
    # growth model params
    w_rh=>wr_val_comb_new, w_t=>we_val_comb_new, w_m=>we_val_comb_new, w_q=>wq_val_comb_new, gmax=>gmax_val_comb_new, 
    θ_rh=>thetar_uM_val, θ_nr=>thetax_uM_val, Kq=>Kq_uM_val, Kt=>Kt_uM_val, Km=>Km_uM_val, Kgamma=>Kgamma_uM_val, 
    nq=>nq_val, nrh=>nr_val, nx=>nx_val, 
    kb=>kb_uM_val, ku=>ku_val, 
    s0=>s0_val_comb_new, M=>M_uM_val, ns=>ns_val, vt=>vt_val_comb_new, vm=>vm_val)

ssvals_red_comb = steady_states(reduced_combined_model, init_red_comb, params_red_comb)
ssvals_red_comb_new = steady_states(reduced_combined_model, init_red_comb, params_red_comb_new)

calc_lam(params_red_comb_new, ssvals_red_comb_new)