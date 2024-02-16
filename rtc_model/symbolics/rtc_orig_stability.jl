using ModelingToolkit, DifferentialEquations, Plots, Latexify, LinearAlgebra

@variables t 
@parameters L c kr Vmax_init Km_init ω_ab ω_r θtscr g_max θtlr km_a km_b d krep kdam ktag kdeg kin atp na nb nr lam kc k_diss 
@syms rm_a(t) rtca(t) rm_b(t) rtcb(t) rm_r(t) rtcr(t) rh(t) rd(t) rt(t) 
    
D = Differential(t)

@mtkmodel RTC begin
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

        # # ribosomes
        Vrep ~ rtcb*rt*krep/(rt+km_b) # uM min-1 
        Vdam ~ kdam*rh # uM min-1
        Vinflux ~ kin* g_max*atp/(θtlr+atp) # uM min-1 
        Vtag ~ rtca*rd*ktag/(rd+km_a) # uM min-1 

        rhs_rm_a ~ tscr_ab - lam*(rm_a) - d*(rm_a)
        rhs_rtca ~ (1/na)*kc*rh*rm_a*(g_max*atp/(θtlr+atp)) - lam*(rtca)     
        rhs_rm_b ~ tscr_ab - lam*(rm_b) - d*(rm_b)
        rhs_rtcb ~ (1/nb)*kc*rh*rm_b*(g_max*atp/(θtlr+atp)) - lam*(rtcb)
        rhs_rm_r ~ tscr_r - lam*(rm_r) - d*(rm_r)
        rhs_rtcr ~ (1/nr)*kc*rh*rm_r*(g_max*atp/(θtlr+atp)) - lam*(rtcr)
        rhs_rh ~ Vrep - Vdam + Vinflux - lam*(rh)
        rhs_rd ~ Vdam - Vtag - kdeg*rd - lam*(rd)
        rhs_rt ~ Vtag - Vrep - lam*(rt)

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


@mtkbuild rtc = RTC()

rtc = ode_order_lowering(rtc)
init_ = [rtc.rm_a=>0.0,rtc.rtca=>0.0,rtc.rm_b=>0.0,rtc.rtcb=>0.0,rtc.rm_r=>0.0,rtc.rtcr=>0.0,rtc.rh=>11.29,rtc.rd=>0.0,rtc.rt=>0.0]

param_dict = Dict(L=>10, c=>0.001, kr=>0.125, Vmax_init=>39.51, Km_init=>250, θtscr=>160.01, θtlr=>255.73, na=>338, nb=>408, nr=>532*6, d=>0.2, krep=>137, ktag=>9780,
atp=>3000, km_a=>20, km_b=>16, g_max=>1260, kdeg=>0.001, kin=>0.022/100, ω_ab=>1e-5, ω_r=>1e-6, kdam=>0, lam=>0.014, kc=>0.6, k_diss=>0.006)

prob = ODEProblem(rtc, init_, (0.0,1e9), param_dict; jac=true)

# latexify(equations(rtc)) |> render

equations(rtc)
observed(rtc)
# states(rtc)
# parameters(rtc)
solu = solve(prob)

jac_sym=calculate_jacobian(rtc)


var_dict = Dict(rm_a(t)=>solu[rtc.rm_a, end], rtca(t)=>solu[rtc.rtca, end], rm_b(t)=>solu[rtc.rm_b, end], rtcb(t)=>solu[rtc.rtcb, end], rm_r(t)=>solu[rtc.rm_r, end], rtcr(t)=>solu[rtc.rtcr, end],
rh(t)=>solu[rtc.rh, end], rd(t)=>solu[rtc.rd, end], rt(t)=>solu[rtc.rt, end])

param_replace = Symbolics.value.(substitute(jac_sym, param_dict))

vals=[]
for i in 1:size(param_replace, 1)
    for j in 1:size(param_replace, 2)
        new::String = repr(param_replace[i,j])
        param_replace[i,j] = eval(Meta.parse(new))
        push!(vals, substitute(param_replace[i,j], var_dict))
    end
end

jac_mat = Float64.(transpose(reshape(vals, (9,9))))
eigvals(jac_mat)

latexify(param_replace) |> render
