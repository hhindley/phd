using ModelingToolkit, DifferentialEquations, PlotlyJS, LinearAlgebra, DataFrames, LabelledArrays, Printf, BifurcationKit, OrderedCollections, ProgressBars, Combinatorics

include(joinpath(homedir(), "phd/general_funcs/all_model_funcs.jl"))
include(joinpath(homedir(), "phd/general_funcs/solving.jl"))
include(joinpath(homedir(), "phd/rtc_model/parameters/rtc_params.jl"))
include(joinpath(homedir(), "phd/rtc_model/parameters/rtc_params_molecs.jl"))

include(joinpath(homedir(), "phd/rtc_model/models/rtc_orig.jl"))
include(joinpath(homedir(), "phd/rtc_model/functions/bf_funcs/bf_funcs.jl"))
include(joinpath(homedir(), "phd/rtc_model/paper/server_code/model_params_funcs_2024/solving.jl"))

# colours =["#636EFA", "#EF553B", "#00CC96", "#AB63FA", "#FFA15A", "#19D3F3", "#FF6692", "#B6E880", "#FF97FF", "#FECB52", :blue]
colours =["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", :blue]



indexof(sym, syms) = findfirst(isequal(sym),syms)

@variables t 
@parameters L c kr Vmax_init Km_init ω_ab ω_r θtscr g_max θtlr km_a km_b d krep kdam ktag kdeg kin_c atp na nb nr lam_c kc k_diss 
species_rtc1 = @syms rm_a(t) rtca(t) rm_b(t) rtcb(t) rm_r(t) rtcr(t) rh(t) rd(t) rt(t) 
species_rtc = [Symbol(i) for i in species_rtc1]
  
D = Differential(t)

@mtkmodel TEST begin
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
        kin_c
        atp 
        na 
        nb 
        nr 
        lam_c
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
        lam(t)
        kin(t)
        

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

        lam ~ rh*tlr_el*lam_c
        kin ~ rh*kin_c
        # kin ~ lam*rh

        

        # # ribosomes
        Vrep ~ rtcb*rt*krep/(rt+km_b) # uM min-1 
        Vdam ~ kdam*rh # uM min-1
        Vinflux ~ kin* g_max*atp/(θtlr+atp) # uM min-1 
        Vtag ~ rtca*rd*ktag/(rd+km_a) # uM min-1 

        rhs_rm_a ~ tscr_ab - dil(rm_a,lam) - deg(rm_a)
        rhs_rtca ~ tlr(rm_a, na, rh, tlr_el) - dil(rtca,lam)     
        rhs_rm_b ~ tscr_ab - dil(rm_b,lam) - deg(rm_b)
        rhs_rtcb ~ tlr(rm_b, nb, rh, tlr_el) - dil(rtcb,lam)
        rhs_rm_r ~ tscr_r - dil(rm_r,lam) - deg(rm_r)
        rhs_rtcr ~ tlr(rm_r, nr, rh, tlr_el) - dil(rtcr,lam)
        rhs_rh ~ Vrep - Vdam + Vinflux - dil(rh,lam)
        rhs_rd ~ Vdam - Vtag - kdeg*rd - dil(rd,lam)
        rhs_rt ~ Vtag - Vrep - dil(rt,lam)

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

@mtkbuild test = TEST()

tlr1 = g_max_val*atp_val/(θtlr_val+atp_val)


lam_c_val = 8e-7 #8e-7
kin_c_val = 1.5e-5 #1.5e-5 # this and the above value give pretty much same model concs as other model but kin is a bit high so need to adjust 
ω_ab_val = 1.3e-5
# L_val = 100000
# c_val = 0.1

params_rtc1 = OrderedDict(L=>L_val, c=>c_val, kr=>kr_val, Vmax_init=>Vmax_init_val, Km_init=>Km_init_val, θtscr=>θtscr_val, θtlr=>θtlr_val, na=>nA_val, nb=>nB_val, nr=>nR_val, d=>d_val, krep=>krep_val, ktag=>ktag_val,
atp=>atp_val, km_a=>km_a_val, km_b=>km_b_val, g_max=>g_max_val, kdeg=>kdeg_val, kin_c=>kin_c_val, ω_ab=>ω_ab_val, ω_r=>ω_r_val, kdam=>kdam_val, lam_c=>lam_c_val, kc=>kc_val, k_diss=>k_diss_val)

init_rtc = [test.rm_a=>0.0,test.rtca=>0.0,test.rm_b=>0.0,test.rtcb=>0.0,test.rm_r=>0.0,test.rtcr=>0.0,test.rh=>11.29,test.rd=>0.0,test.rt=>0.0]

ssvals_rtc = steady_states(test, init_rtc, params_rtc1)



