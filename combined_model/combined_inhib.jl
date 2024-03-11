include("$PATH/general_funcs/all_model_funcs.jl")
indexof(sym, syms) = findfirst(isequal(sym),syms)

@variables t
@parameters d kr L c w_rh w_t w_m w_q w_BA w_R θ_rh θ_nr Kq nq Vmax_init Km_init nrh nx nR nA nB gmax Kgamma M kb ku abx kon koff vt vm s0 Kt Km km_a km_b krep kdam_p ktag kdeg k_diss ns k_inhib1 k_inhib2 inhib
species_inhib_comb1 = @syms m_rh(t) m_t(t) m_m(t) m_q(t) m_R(t) m_A(t) m_B(t) c_rh(t) c_t(t) c_m(t) c_q(t) c_R(t) c_A(t) c_B(t) et(t) em(t) q(t) R(t) A(t) B(t) rh(t) rt(t) rd(t) z_rh(t) z_t(t) z_m(t) z_q(t) z_R(t) z_A(t) z_B(t) si(t) a(t) rtc_i(t)
species_inhib_comb = [Symbol(i) for i in species_inhib_comb1]

D = Differential(t)

function build_combined_inhib(inhib_protein1)
    inhib_protein = inhib_protein1

    @mtkmodel COMBINED_INHIB_MODEL begin
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
            abx 
            kon 
            koff 
            vt 
            vm 
            s0 
            Kt 
            Km 
            km_a 
            km_b 
            krep 
            kdam_p 
            ktag 
            kdeg 
            k_diss 
            ns
            k_inhib1
            k_inhib2
            inhib
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
            z_rh(t) 
            z_t(t) 
            z_m(t) 
            z_q(t) 
            z_R(t) 
            z_A(t) 
            z_B(t) 
            si(t) 
            a(t)
            rtc_i(t)

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
            rhs_z_rh(t) 
            rhs_z_t(t) 
            rhs_z_m(t) 
            rhs_z_q(t) 
            rhs_z_R(t) 
            rhs_z_A(t) 
            rhs_z_B(t) 
            rhs_si(t) 
            rhs_a(t)
            rhs_rtc_i(t) 
            
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
            sig_o ~ ra*Voc/k_diss # uM

            # translation 
            gamma ~ gmax*a/(Kgamma+a) # aa min-1 
            ttrate ~ (c_q+c_rh+c_t+c_m+c_R+c_A+c_B)*gamma # aa uM min-1

            # growth rate 
            lam ~ ttrate/M # min-1

            # metabolism and import 
            vimp ~ (et*vt*s0/(Kt+s0)) # uM min-1 
            vcat ~ (em*vm*si/(Km+si)) # uM min-1 

            # rtc
            Vrep ~ B*rt*krep/(rt+km_b) # uM min-1
            Vdam ~ (z_rh + z_t + z_m + z_q + z_R + z_A + z_B)*koff*kdam_p*ku #(z_rh + z_t + z_m + z_q + z_R + z_A + z_B)*koff*kdam # uM min-1
            Vtag ~ A*rd*ktag/(rd+km_a) # uM min-1

            # ODEs
            # mRNA
            rhs_m_rh ~ ω_p(w_rh, θ_rh, a) - dil(m_rh, lam) - deg(m_rh) + v_x(c_rh, nrh, gamma) - rh_bind(m_rh, rh) + rh_unbind(c_rh) + zm_diss(z_rh)*ku*kdam_p
            rhs_m_t ~ ω_p(w_t, θ_nr, a) - dil(m_t, lam) - deg(m_t) + v_x(c_t, nx, gamma) - rh_bind(m_t, rh) + rh_unbind(c_t) + zm_diss(z_t)*ku*kdam_p
            rhs_m_m ~ ω_p(w_m, θ_nr, a) - dil(m_m, lam) - deg(m_m) + v_x(c_m, nx, gamma) - rh_bind(m_m, rh) + rh_unbind(c_m) + zm_diss(z_m)*ku*kdam_p
            rhs_m_q ~ ω_q(w_q, θ_nr, a, q) - dil(m_q, lam) - deg(m_q) + v_x(c_q, nx, gamma) - rh_bind(m_q, rh) + rh_unbind(c_q) + zm_diss(z_q)*ku*kdam_p
            rhs_m_R ~ ω_p(w_R, θ_nr, a) - dil(m_R, lam) - deg(m_R) + v_x(c_R, nR, gamma) - rh_bind(m_R, rh) + rh_unbind(c_R) + zm_diss(z_R)*ku*kdam_p
            rhs_m_A ~ ω_rtcBA(θ_nr, a, sig_o) - dil(m_A, lam) - deg(m_A) + v_x(c_A, nA, gamma) - rh_bind(m_A, rh) + rh_unbind(c_A) + zm_diss(z_A)*ku*kdam_p
            rhs_m_B ~ ω_rtcBA(θ_nr, a, sig_o) - dil(m_B, lam) - deg(m_B) + v_x(c_B, nB, gamma) - rh_bind(m_B, rh) + rh_unbind(c_B) + zm_diss(z_B)*ku*kdam_p

            # mRNA:ribosome complexes 
            rhs_c_rh ~ rh_bind(m_rh, rh) - dil(c_rh, lam) - rh_unbind(c_rh) - v_x(c_rh, nrh, gamma) - zm(c_rh) + zm_diss(z_rh)*(1-kdam_p)
            rhs_c_t ~ rh_bind(m_t, rh) - dil(c_t, lam) - rh_unbind(c_t) - v_x(c_t, nx, gamma) - zm(c_t) + zm_diss(z_t)*(1-kdam_p)
            rhs_c_m ~ rh_bind(m_m, rh) - dil(c_m, lam) - rh_unbind(c_m) - v_x(c_m, nx, gamma) - zm(c_m) + zm_diss(z_m)*(1-kdam_p)
            rhs_c_q ~ rh_bind(m_q, rh) - dil(c_q, lam) - rh_unbind(c_q) - v_x(c_q, nx, gamma) - zm(c_q) + zm_diss(z_q)*(1-kdam_p)
            rhs_c_R ~ rh_bind(m_R, rh) - dil(c_R, lam) - rh_unbind(c_R) - v_x(c_R, nR, gamma) - zm(c_R) + zm_diss(z_R)*(1-kdam_p)
            rhs_c_A ~ rh_bind(m_A, rh) - dil(c_A, lam) - rh_unbind(c_A) - v_x(c_A, nA, gamma) - zm(c_A) + zm_diss(z_A)*(1-kdam_p)
            rhs_c_B ~ rh_bind(m_B, rh) - dil(c_B, lam) - rh_unbind(c_B) - v_x(c_B, nB, gamma) - zm(c_B) + zm_diss(z_B)*(1-kdam_p)

            # proteins 
            rhs_et ~ v_x(c_t, nx, gamma) - dil(et, lam)
            rhs_em ~ v_x(c_m, nx, gamma) - dil(em, lam)
            rhs_q ~ v_x(c_q, nx, gamma) - dil(q, lam)
            # rhs_R ~ v_x(c_R, nR, gamma) - dil(R, lam)
            # rhs_A ~ v_x(c_A, nA, gamma) - dil(A, lam)
            # rhs_B ~ v_x(c_B, nB, gamma) - dil(B, lam)

            if inhib_protein == :rtca_inhib
                rhs_rtc_i ~ k_inhib1*A*inhib - k_inhib2*rtc_i - dil(rtc_i,lam)            
                rhs_A ~ v_x(c_A, nA, gamma) - dil(A, lam) - k_inhib1*A*inhib + k_inhib2*rtc_i   
                rhs_B ~ v_x(c_B, nB, gamma) - dil(B, lam)
                rhs_R ~ v_x(c_R, nR, gamma) - dil(R, lam)
            elseif inhib_protein == :rtcb_inhib
                rhs_rtc_i ~ k_inhib1*B*inhib - k_inhib2*rtc_i - dil(rtc_i,lam)            
                rhs_A ~ v_x(c_A, nA, gamma) - dil(A, lam)
                rhs_B ~ v_x(c_B, nB, gamma) - dil(B, lam) - k_inhib1*B*inhib + k_inhib2*rtc_i   
                rhs_R ~ v_x(c_R, nR, gamma) - dil(R, lam)
            else 
                rhs_rtc_i ~ k_inhib1*R*inhib - k_inhib2*rtc_i - dil(rtc_i,lam)            
                rhs_A ~ v_x(c_A, nA, gamma) - dil(A, lam)
                rhs_B ~ v_x(c_B, nB, gamma) - dil(B, lam)
                rhs_R ~ v_x(c_R, nR, gamma) - dil(R, lam) - k_inhib1*R*inhib + k_inhib2*rtc_i   
            end

            # ribosomes
            rhs_rh ~ v_x(c_rh, nrh, gamma) - dil(rh, lam) - (rh_bind(m_rh, rh) + rh_bind(m_t, rh) + rh_bind(m_m, rh) + rh_bind(m_q, rh) + rh_bind(m_R, rh) + rh_bind(m_A, rh) + rh_bind(m_B, rh)) + (rh_unbind(c_rh) + rh_unbind(c_t) + rh_unbind(c_m) + rh_unbind(c_q) + rh_unbind(c_R) + rh_unbind(c_A) + rh_unbind(c_B)) + (v_x(c_rh, nrh, gamma) + v_x(c_t, nx, gamma) + v_x(c_m, nx, gamma) + v_x(c_q, nx, gamma) + v_x(c_R, nR, gamma) + v_x(c_A, nA, gamma) + v_x(c_B, nB, gamma)) + Vrep
            rhs_rt ~ Vtag - Vrep - dil(rt, lam)
            rhs_rd ~ Vdam - Vtag - rd*kdeg - dil(rd, lam)

            # zombie complexes 
            rhs_z_rh ~ zm(c_rh) - zm_diss(z_rh)*(1-kdam_p) - zm_diss(z_rh)*kdam_p - dil(z_rh, lam)
            rhs_z_t ~ zm(c_t) - zm_diss(z_t)*(1-kdam_p) - zm_diss(z_t)*kdam_p - dil(z_t, lam)
            rhs_z_m ~ zm(c_m) - zm_diss(z_m)*(1-kdam_p) - zm_diss(z_m)*kdam_p - dil(z_m, lam)
            rhs_z_q ~ zm(c_q) - zm_diss(z_q)*(1-kdam_p) - zm_diss(z_q)*kdam_p - dil(z_q, lam)
            rhs_z_R ~ zm(c_R) - zm_diss(z_R)*(1-kdam_p) - zm_diss(z_R)*kdam_p - dil(z_R, lam)
            rhs_z_A ~ zm(c_A) - zm_diss(z_A)*(1-kdam_p) - zm_diss(z_A)*kdam_p - dil(z_A, lam)
            rhs_z_B ~ zm(c_B) - zm_diss(z_B)*(1-kdam_p) - zm_diss(z_B)*kdam_p - dil(z_B, lam)

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
            D(z_rh) ~ rhs_z_rh
            D(z_t) ~ rhs_z_t
            D(z_m) ~ rhs_z_m
            D(z_q) ~ rhs_z_q
            D(z_R) ~ rhs_z_R
            D(z_A) ~ rhs_z_A
            D(z_B) ~ rhs_z_B
            D(si) ~ rhs_si
            D(a) ~ rhs_a
            D(rtc_i) ~ rhs_rtc_i

        end
    end

    @mtkbuild combined_inhib_model = COMBINED_INHIB_MODEL()
end

combined_inhib_rtca = build_combined_inhib(:rtca_inhib)
combined_inhib_rtcb = build_combined_inhib(:rtcb_inhib)
combined_inhib_rtcr = build_combined_inhib(:rtcr_inhib)

init_inhib_comb_rtca = [combined_inhib_rtca.m_rh=>0.0,combined_inhib_rtca.m_t=>0.0,combined_inhib_rtca.m_m=>0.0,combined_inhib_rtca.m_q=>0.0,combined_inhib_rtca.m_R=>0.0,combined_inhib_rtca.m_A=>0.0,combined_inhib_rtca.m_B=>0.0,combined_inhib_rtca.c_rh=>0.0,combined_inhib_rtca.c_t=>0.0,combined_inhib_rtca.c_m=>0.0,combined_inhib_rtca.c_q=>0.0,combined_inhib_rtca.c_R=>0.0,combined_inhib_rtca.c_A=>0.0,combined_inhib_rtca.c_B=>0.0,combined_inhib_rtca.et=>0.0,combined_inhib_rtca.em=>0.0,combined_inhib_rtca.q=>0.0,combined_inhib_rtca.R=>0.0,combined_inhib_rtca.A=>0.0,combined_inhib_rtca.B=>0.0,combined_inhib_rtca.rh=>0.0166,combined_inhib_rtca.rt=>0.0,combined_inhib_rtca.rd=>0.0,combined_inhib_rtca.z_rh=>0.0,combined_inhib_rtca.z_t=>0.0,combined_inhib_rtca.z_m=>0.0,combined_inhib_rtca.z_q=>0.0,combined_inhib_rtca.z_R=>0.0,combined_inhib_rtca.z_A=>0.0,combined_inhib_rtca.z_B=>0.0,combined_inhib_rtca.si=>0.0,combined_inhib_rtca.a=>1.66,combined_inhib_rtca.rtc_i=>0.0]
init_inhib_comb_rtcb = [combined_inhib_rtcb.m_rh=>0.0,combined_inhib_rtcb.m_t=>0.0,combined_inhib_rtcb.m_m=>0.0,combined_inhib_rtcb.m_q=>0.0,combined_inhib_rtcb.m_R=>0.0,combined_inhib_rtcb.m_A=>0.0,combined_inhib_rtcb.m_B=>0.0,combined_inhib_rtcb.c_rh=>0.0,combined_inhib_rtcb.c_t=>0.0,combined_inhib_rtcb.c_m=>0.0,combined_inhib_rtcb.c_q=>0.0,combined_inhib_rtcb.c_R=>0.0,combined_inhib_rtcb.c_A=>0.0,combined_inhib_rtcb.c_B=>0.0,combined_inhib_rtcb.et=>0.0,combined_inhib_rtcb.em=>0.0,combined_inhib_rtcb.q=>0.0,combined_inhib_rtcb.R=>0.0,combined_inhib_rtcb.A=>0.0,combined_inhib_rtcb.B=>0.0,combined_inhib_rtcb.rh=>0.0166,combined_inhib_rtcb.rt=>0.0,combined_inhib_rtcb.rd=>0.0,combined_inhib_rtcb.z_rh=>0.0,combined_inhib_rtcb.z_t=>0.0,combined_inhib_rtcb.z_m=>0.0,combined_inhib_rtcb.z_q=>0.0,combined_inhib_rtcb.z_R=>0.0,combined_inhib_rtcb.z_A=>0.0,combined_inhib_rtcb.z_B=>0.0,combined_inhib_rtcb.si=>0.0,combined_inhib_rtcb.a=>1.66,combined_inhib_rtcb.rtc_i=>0.0]
init_inhib_comb_rtcr = [combined_inhib_rtcr.m_rh=>0.0,combined_inhib_rtcr.m_t=>0.0,combined_inhib_rtcr.m_m=>0.0,combined_inhib_rtcr.m_q=>0.0,combined_inhib_rtcr.m_R=>0.0,combined_inhib_rtcr.m_A=>0.0,combined_inhib_rtcr.m_B=>0.0,combined_inhib_rtcr.c_rh=>0.0,combined_inhib_rtcr.c_t=>0.0,combined_inhib_rtcr.c_m=>0.0,combined_inhib_rtcr.c_q=>0.0,combined_inhib_rtcr.c_R=>0.0,combined_inhib_rtcr.c_A=>0.0,combined_inhib_rtcr.c_B=>0.0,combined_inhib_rtcr.et=>0.0,combined_inhib_rtcr.em=>0.0,combined_inhib_rtcr.q=>0.0,combined_inhib_rtcr.R=>0.0,combined_inhib_rtcr.A=>0.0,combined_inhib_rtcr.B=>0.0,combined_inhib_rtcr.rh=>0.0166,combined_inhib_rtcr.rt=>0.0,combined_inhib_rtcr.rd=>0.0,combined_inhib_rtcr.z_rh=>0.0,combined_inhib_rtcr.z_t=>0.0,combined_inhib_rtcr.z_m=>0.0,combined_inhib_rtcr.z_q=>0.0,combined_inhib_rtcr.z_R=>0.0,combined_inhib_rtcr.z_A=>0.0,combined_inhib_rtcr.z_B=>0.0,combined_inhib_rtcr.si=>0.0,combined_inhib_rtcr.a=>1.66,combined_inhib_rtcr.rtc_i=>0.0]

params_inhib_comb = Dict(d=>d_val, kr=>kr_val, L=>L_val, c=>c_val, w_rh=>wr_uM_val, w_t=>we_uM_val, w_m=>we_uM_val, w_q=>wq_uM_val, w_BA=>ω_ab_val_comb, w_R=>ω_r_val_comb, θ_rh=>thetar_uM_val, θ_nr=>thetax_uM_val, Kq=>Kq_uM_val, nq=>nq_val, Vmax_init=>Vmax_init_val, Km_init=>Km_init_val, nrh=>nr_val, nx=>nx_val, nR=>nR_val, nA=>nA_val, nB=>nB_val, gmax=>gmax_val, Kgamma=>Kgamma_uM_val, M=>M_uM_val, kb=>kb_uM_val, ku=>ku_val, abx=>abx_val, kon=>kon_val, koff=>koff_val, vt=>vt_val, vm=>vm_val, s0=>s0_uM_val, Kt=>Kt_uM_val, Km=>Km_uM_val, km_a=>km_a_val, km_b=>km_b_val, krep=>krep_val, kdam_p=>kdam_p_val, ktag=>ktag_val, kdeg=>kdeg_val, k_diss=>k_diss_val, ns=>ns_val, k_inhib1=>k_inhib1_val_comb, k_inhib2=>k_inhib2_val_comb, inhib=>inhib_val_comb)

ssvals_inhib_comb_rtca = steady_states(combined_inhib_rtca, init_inhib_comb_rtca, params_inhib_comb)
ssvals_inhib_comb_rtcb = steady_states(combined_inhib_rtcb, init_inhib_comb_rtcb, params_inhib_comb)
ssvals_inhib_comb_rtcr = steady_states(combined_inhib_rtcr, init_inhib_comb_rtcr, params_inhib_comb)

# prob2 = ODEProblem(growth_model, init_gm, tspan, params_gm; jac=true);
# solu2 = solve(prob2, Rodas4(), abstol=1e-12, reltol=1e-9);



# function combined_model(initial, params, t)
#     d, kr, L, c, w_rh, w_t, w_m, w_q, w_BA, w_R, θ_rh, θ_nr, Kq, hq, Vmax_init, Km_init, nrh, nx, nR, nA, nB, gmax, Kgamma, M, kb, ku, abx, kon, koff, vt, vm, s0, Kt, Km, km_a, km_b, krep, kdam_p, ktag, kdeg, k_diss, ns = params
#     m_rh, m_t, m_m, m_q, m_R, m_A, m_B, c_rh, c_t, c_m, c_q, c_R, c_A, c_B, et, em, q, R, A, B, rh, rt, rd, z_rh, z_t, z_m, z_q, z_R, z_A, z_B, si, a = initial
    
#     # MWC
#     alpha = rt/kr # unitless 
#     fa = (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6) # unitless
#     ra = fa*R # uM 

#     # transcription 
#     Voc = Vmax_init*a/(Km_init+a) # uM min-1 
#     sig_o = ra*Voc/k_diss
#     ω_rtcBA(θ_x) = sig_o * w_BA*a/(θ_x + a) # transcription rate for RtcBA # uM min-1

#     ω_p(w_x, θ_x) = w_x*a/(θ_x + a) # transcription rate for rh, et, em, rtcR # uM min-1
#     ω_q(θ_x) = ((w_q*a/(θ_x + a))/(1+(q/Kq)^hq)) # transcription rate for q # uM min-1

#     # translation 
#     gamma = gmax*a/(Kgamma+a) # aa min-1 uM-1
#     v_x(c_x, nx) = c_x*gamma/nx # min-1
#     ttrate = (c_q+c_rh+c_t+c_m+c_R+c_A+c_B)*gamma # aa min-1

#     # growth rate 
#     λ = ttrate/M # min-1

#     # dilution by growth and degradation 
#     dil(x) = λ*x # uM min-1 
#     deg(x) = d*x # uM min-1 

#     # ribosome binding/unbinding
#     rh_bind(m_x) = kb*rh*m_x # uM min-1 
#     rh_unbind(c_x) = ku*c_x # uM min-1 

#     # zombie complex formation 
#     zm(c_x) = c_x*abx*kon # uM min-1
#     zm_diss(z_x) = koff*z_x # uM min-1 

#     # metabolism and import 
#     vimp = (et*vt*s0/(Kt+s0)) # uM min-1 
#     vcat = (em*vm*si/(Km+si)) # uM min-1 

#     # rtc
#     Vrep = B*rt*krep/(rt+km_b) # uM min-1
#     Vdam = (zm_diss(z_rh)+zm_diss(z_t)+zm_diss(z_m)+zm_diss(z_q)+zm_diss(z_R)+zm_diss(z_A)+zm_diss(z_B))*kdam_p*ku #(z_rh + z_t + z_m + z_q + z_R + z_A + z_B)*koff*kdam # uM min-1
#     Vtag = A*rd*ktag/(rd+km_a) # uM min-1

#     # ODEs
#     # mRNA
#     dm_rh = ω_p(w_rh, θ_rh) - dil(m_rh) - deg(m_rh) + v_x(c_rh, nrh) - rh_bind(m_rh) + rh_unbind(c_rh) + zm_diss(z_rh)*ku*kdam_p
#     dm_t = ω_p(w_t, θ_nr) - dil(m_t) - deg(m_t) + v_x(c_t, nx) - rh_bind(m_t) + rh_unbind(c_t) + zm_diss(z_t)*ku*kdam_p
#     dm_m = ω_p(w_m, θ_nr) - dil(m_m) - deg(m_m) + v_x(c_m, nx) - rh_bind(m_m) + rh_unbind(c_m) + zm_diss(z_m)*ku*kdam_p
#     dm_q = ω_q(θ_nr) - dil(m_q) - deg(m_q) + v_x(c_q, nx) - rh_bind(m_q) + rh_unbind(c_q) + zm_diss(z_q)*ku*kdam_p
#     dm_R = ω_p(w_R, θ_nr) - dil(m_R) - deg(m_R) + v_x(c_R, nR) - rh_bind(m_R) + rh_unbind(c_R) + zm_diss(z_R)*ku*kdam_p
#     dm_A = ω_rtcBA(θ_nr) - dil(m_A) - deg(m_A) + v_x(c_A, nA) - rh_bind(m_A) + rh_unbind(c_A) + zm_diss(z_A)*ku*kdam_p
#     dm_B = ω_rtcBA(θ_nr) - dil(m_B) - deg(m_B) + v_x(c_B, nB) - rh_bind(m_B) + rh_unbind(c_B) + zm_diss(z_B)*ku*kdam_p

#     # mRNA:ribosome complexes 
#     dc_rh = rh_bind(m_rh) - dil(c_rh) - rh_unbind(c_rh) - v_x(c_rh, nrh) - zm(c_rh) + zm_diss(z_rh)*(1-kdam_p)
#     dc_t = rh_bind(m_t) - dil(c_t) - rh_unbind(c_t) - v_x(c_t, nx) - zm(c_t) + zm_diss(z_t)*(1-kdam_p)
#     dc_m = rh_bind(m_m) - dil(c_m) - rh_unbind(c_m) - v_x(c_m, nx) - zm(c_m) + zm_diss(z_m)*(1-kdam_p)
#     dc_q = rh_bind(m_q) - dil(c_q) - rh_unbind(c_q) - v_x(c_q, nx) - zm(c_q) + zm_diss(z_q)*(1-kdam_p)
#     dc_R = rh_bind(m_R) - dil(c_R) - rh_unbind(c_R) - v_x(c_R, nR) - zm(c_R) + zm_diss(z_R)*(1-kdam_p)
#     dc_A = rh_bind(m_A) - dil(c_A) - rh_unbind(c_A) - v_x(c_A, nA) - zm(c_A) + zm_diss(z_A)*(1-kdam_p)
#     dc_B = rh_bind(m_B) - dil(c_B) - rh_unbind(c_B) - v_x(c_B, nB) - zm(c_B) + zm_diss(z_B)*(1-kdam_p)

#     # proteins 
#     det = v_x(c_t, nx) - dil(et)
#     dem = v_x(c_m, nx) - dil(em)
#     dq = v_x(c_q, nx) - dil(q)
#     dR = v_x(c_R, nR) - dil(R)
#     dA = v_x(c_A, nA) - dil(A)
#     dB = v_x(c_B, nB) - dil(B)

#     # ribosomes
#     drh = v_x(c_rh, nrh) - dil(rh) - (rh_bind(m_rh) + rh_bind(m_t) + rh_bind(m_m) + rh_bind(m_q) + rh_bind(m_R) + rh_bind(m_A) + rh_bind(m_B)) + (rh_unbind(c_rh) + rh_unbind(c_t) + rh_unbind(c_m) + rh_unbind(c_q) + rh_unbind(c_R) + rh_unbind(c_A) + rh_unbind(c_B)) + (v_x(c_rh, nrh) + v_x(c_t, nx) + v_x(c_m, nx) + v_x(c_q, nx) + v_x(c_R, nR) + v_x(c_A, nA) + v_x(c_B, nB)) + Vrep
#     drt = Vtag - Vrep - dil(rt)
#     drd = Vdam - Vtag - rd*kdeg - dil(rd)

#     # zombie complexes 
#     dz_rh = zm(c_rh) - zm_diss(z_rh)*(1-kdam_p) - zm_diss(z_rh)*kdam_p - dil(z_rh)
#     dz_t = zm(c_t) - zm_diss(z_t)*(1-kdam_p) - zm_diss(z_t)*kdam_p - dil(z_t)
#     dz_m = zm(c_m) - zm_diss(z_m)*(1-kdam_p) - zm_diss(z_m)*kdam_p - dil(z_m)
#     dz_q = zm(c_q) - zm_diss(z_q)*(1-kdam_p) - zm_diss(z_q)*kdam_p - dil(z_q)
#     dz_R = zm(c_R) - zm_diss(z_R)*(1-kdam_p) - zm_diss(z_R)*kdam_p - dil(z_R)
#     dz_A = zm(c_A) - zm_diss(z_A)*(1-kdam_p) - zm_diss(z_A)*kdam_p - dil(z_A)
#     dz_B = zm(c_B) - zm_diss(z_B)*(1-kdam_p) - zm_diss(z_B)*kdam_p - dil(z_B)

#     # nutrient
#     dsi = vimp - vcat - dil(si)

#     # energy
#     da = ns*vcat - ttrate - dil(a)

#     [dm_rh, dm_t, dm_m, dm_q, dm_R, dm_A, dm_B, dc_rh, dc_t, dc_m, dc_q, dc_R, dc_A, dc_B, det, dem, dq, dR, dA, dB, drh, drt, drd, dz_rh, dz_t, dz_m, dz_q, dz_R, dz_A, dz_B, dsi, da]

# end



# function combined_model_v2(initial, params, t)
#     d, kr, L, c, w_rh, w_t, w_m, w_q, w_BA, w_R, θ_rh, θ_nr, Kq, hq, Vmax_init, Km_init, nrh, nx, nR, nA, nB, gmax, Kgamma, M, kb, ku, abx, kon, koff, vt, vm, s0, Kt, Km, km_a, km_b, krep, kdam, ktag, kdeg, kdiss, ns = params
#     m_rh, m_t, m_m, m_q, m_R, m_A, m_B, c_rh, c_t, c_m, c_q, c_R, c_A, c_B, et, em, q, R, A, B, rh, rt, rd, z_rh, z_t, z_m, z_q, z_R, z_A, z_B, c_rhD, c_tD, c_mD, c_qD, c_RD, c_AD, c_BD, si, a = initial
#     # m_rh, m_t, m_m, m_q, m_R, m_A, m_B, c_rh, c_t, c_m, c_q, c_R, c_A, c_B, et, em, q, R, A, B, rh, rt, rd, z_rh, z_t, z_m, z_q, z_R, z_A, z_B, si, a = initial
    
#     # MWC
#     alpha = rt/kr # unitless 
#     fa = (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6) # unitless
#     ra = fa*R # uM 

#     # transcription 
#     Voc = Vmax_init*a/(Km_init+a) # uM min-1 
#     sig_o = ra*Voc/kdiss
#     ω_rtcBA(θ_x) = (w_BA*a/(θ_x + a)) * sig_o # transcription rate for RtcBA # uM min-1

#     ω_p(w_x, θ_x) = w_x*a/(θ_x + a) # transcription rate for rh, et, em, rtcR # uM min-1
#     ω_q(θ_x) = ((w_q*a/(θ_x + a))/(1+(q/Kq)^hq)) # transcription rate for q # uM min-1

#     # translation 
#     gamma = gmax*a/(Kgamma+a) # aa min-1 uM-1
#     v_x(c_x, nx) = c_x*gamma/nx # min-1
#     ttrate = (c_q+c_rh+c_t+c_m+c_R+c_A+c_B)*gamma # aa min-1

#     # growth rate 
#     λ = ttrate/M # min-1

#     # dilution by growth and degradation 
#     dil(x) = λ*x # uM min-1 
#     deg(x) = d*x # uM min-1 

#     # ribosome binding/unbinding
#     rh_bind(m_x) = kb*rh*m_x # uM min-1 
#     rh_unbind(c_x) = ku*c_x # uM min-1 

#     # zombie complex formation 
#     zm(c_x) = c_x*abx*kon # uM min-1
#     zm_diss(z_x) = koff*z_x # uM min-1 

#     # metabolism and import 
#     vimp = (et*vt*s0/(Kt+s0)) # uM min-1 
#     vcat = (em*vm*si/(Km+si)) # uM min-1 

#     # rtc
#     Vrep = B*rt*krep/(rt+km_b) # uM min-1
#     Vdam = (c_rhD + c_tD + c_mD + c_qD + c_RD + c_AD + c_BD) * ku 
#     # Vdam = (zm_diss(z_rh)+zm_diss(z_t)+zm_diss(z_m)+zm_diss(z_q)+zm_diss(z_R)+zm_diss(z_A)+zm_diss(z_B))*kdam*ku #(z_rh + z_t + z_m + z_q + z_R + z_A + z_B)*koff*kdam # uM min-1
#     Vtag = A*rd*ktag/(rd+km_a) # uM min-1

#     # ODEs
#     # mRNA
#     dm_rh = ω_p(w_rh, θ_rh) - dil(m_rh) - deg(m_rh) + v_x(c_rh, nrh) - rh_bind(m_rh) + rh_unbind(c_rh) + rh_unbind(c_rhD)
#     dm_t = ω_p(w_t, θ_nr) - dil(m_t) - deg(m_t) + v_x(c_t, nx) - rh_bind(m_t) + rh_unbind(c_t) + rh_unbind(c_tD)
#     dm_m = ω_p(w_m, θ_nr) - dil(m_m) - deg(m_m) + v_x(c_m, nx) - rh_bind(m_m) + rh_unbind(c_m) + rh_unbind(c_mD)
#     dm_q = ω_q(θ_nr) - dil(m_q) - deg(m_q) + v_x(c_q, nx) - rh_bind(m_q) + rh_unbind(c_q) + rh_unbind(c_qD)
#     dm_R = ω_p(w_R, θ_nr) - dil(m_R) - deg(m_R) + v_x(c_R, nR) - rh_bind(m_R) + rh_unbind(c_R) + rh_unbind(c_RD)
#     dm_A = ω_rtcBA(θ_nr) - dil(m_A) - deg(m_A) + v_x(c_A, nA) - rh_bind(m_A) + rh_unbind(c_A) + rh_unbind(c_AD)
#     dm_B = ω_rtcBA(θ_nr) - dil(m_B) - deg(m_B) + v_x(c_B, nB) - rh_bind(m_B) + rh_unbind(c_B) + rh_unbind(c_BD)

#     # mRNA:ribosome complexes 
#     dc_rh = rh_bind(m_rh) - dil(c_rh) - rh_unbind(c_rh) - v_x(c_rh, nrh) - zm(c_rh) + zm_diss(z_rh)*(1-kdam)
#     dc_t = rh_bind(m_t) - dil(c_t) - rh_unbind(c_t) - v_x(c_t, nx) - zm(c_t) + zm_diss(z_t)*(1-kdam)
#     dc_m = rh_bind(m_m) - dil(c_m) - rh_unbind(c_m) - v_x(c_m, nx) - zm(c_m) + zm_diss(z_m)*(1-kdam)
#     dc_q = rh_bind(m_q) - dil(c_q) - rh_unbind(c_q) - v_x(c_q, nx) - zm(c_q) + zm_diss(z_q)*(1-kdam)
#     dc_R = rh_bind(m_R) - dil(c_R) - rh_unbind(c_R) - v_x(c_R, nR) - zm(c_R) + zm_diss(z_R)*(1-kdam)
#     dc_A = rh_bind(m_A) - dil(c_A) - rh_unbind(c_A) - v_x(c_A, nA) - zm(c_A) + zm_diss(z_A)*(1-kdam)
#     dc_B = rh_bind(m_B) - dil(c_B) - rh_unbind(c_B) - v_x(c_B, nB) - zm(c_B) + zm_diss(z_B)*(1-kdam)

#     # proteins 
#     det = v_x(c_t, nx) - dil(et)
#     dem = v_x(c_m, nx) - dil(em)
#     dq = v_x(c_q, nx) - dil(q)
#     dR = v_x(c_R, nR) - dil(R)
#     dA = v_x(c_A, nA) - dil(A)
#     dB = v_x(c_B, nB) - dil(B)

#     # ribosomes
#     drh = v_x(c_rh, nrh) - dil(rh) - (rh_bind(m_rh) + rh_bind(m_t) + rh_bind(m_m) + rh_bind(m_q)) + (rh_unbind(c_rh) + rh_unbind(c_t) + rh_unbind(c_m) + rh_unbind(c_q)) + (v_x(c_rh, nrh) + v_x(c_t, nx) + v_x(c_m, nx) + v_x(c_q, nx)) + Vrep
#     drt = Vtag - Vrep - dil(rt)
#     drd = Vdam - Vtag - rd*kdeg - dil(rd)

#     # zombie complexes 
#     dz_rh = zm(c_rh) - zm_diss(z_rh)*(1-kdam) - zm_diss(z_rh)*kdam - dil(z_rh)
#     dz_t = zm(c_t) - zm_diss(z_t)*(1-kdam) - zm_diss(z_t)*kdam - dil(z_t)
#     dz_m = zm(c_m) - zm_diss(z_m)*(1-kdam) - zm_diss(z_m)*kdam - dil(z_m)
#     dz_q = zm(c_q) - zm_diss(z_q)*(1-kdam) - zm_diss(z_q)*kdam - dil(z_q)
#     dz_R = zm(c_R) - zm_diss(z_R)*(1-kdam) - zm_diss(z_R)*kdam - dil(z_R)
#     dz_A = zm(c_A) - zm_diss(z_A)*(1-kdam) - zm_diss(z_A)*kdam - dil(z_A)
#     dz_B = zm(c_B) - zm_diss(z_B)*(1-kdam) - zm_diss(z_B)*kdam - dil(z_B)

#     dc_rhD = zm_diss(z_rh)*kdam - rh_unbind(c_rhD) - dil(c_rhD)
#     dc_tD = zm_diss(z_t)*kdam - rh_unbind(c_tD) - dil(c_tD)
#     dc_mD = zm_diss(z_m)*kdam - rh_unbind(c_mD) - dil(c_mD)
#     dc_qD = zm_diss(z_q)*kdam - rh_unbind(c_qD) - dil(c_qD)
#     dc_RD = zm_diss(z_R)*kdam - rh_unbind(c_RD) - dil(c_RD)
#     dc_AD = zm_diss(z_A)*kdam - rh_unbind(c_AD) - dil(c_AD)
#     dc_BD = zm_diss(z_B)*kdam - rh_unbind(c_BD) - dil(c_BD)

#     # nutrient
#     dsi = vimp - vcat - dil(si)

#     # energy
#     da = ns*vcat - ttrate - dil(a)

#     [dm_rh, dm_t, dm_m, dm_q, dm_R, dm_A, dm_B, dc_rh, dc_t, dc_m, dc_q, dc_R, dc_A, dc_B, det, dem, dq, dR, dA, dB, drh, drt, drd, dz_rh, dz_t, dz_m, dz_q, dz_R, dz_A, dz_B, dc_rhD, dc_tD, dc_mD, dc_qD, dc_RD, dc_AD, dc_BD, dsi, da]
#     # [dm_rh, dm_t, dm_m, dm_q, dm_R, dm_A, dm_B, dc_rh, dc_t, dc_m, dc_q, dc_R, dc_A, dc_B, det, dem, dq, dR, dA, dB, drh, drt, drd, dz_rh, dz_t, dz_m, dz_q, dz_R, dz_A, dz_B, dsi, da]

# end