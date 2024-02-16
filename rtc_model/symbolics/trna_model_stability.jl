using ModelingToolkit, DifferentialEquations, PlotlyJS, Latexify, LinearAlgebra

@variables t 
@parameters L c kr Vmax_init Km_init ω_ab ω_r θtscr g_max θtlr km_a km_b d krep kdam ktag kdeg kin atp na nb nr lam kc k_diss rh thr_t
@syms rm_a(t) rtca(t) rm_b(t) rtcb(t) rm_r(t) rtcr(t) trna(t) rd(t) rt(t) 
    
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
        Vinflux ~ kin* g_max*atp/(θtlr+atp) # uM min-1 
        Vtag ~ rtca*rd*ktag/(rd+km_a) # uM min-1 

        rhs_rm_a ~ tscr_ab - lam*(rm_a) - d*(rm_a)
        rhs_rtca ~ (1/na)*kc*rh*rm_a*tlr_el - lam*(rtca)     
        rhs_rm_b ~ tscr_ab - lam*(rm_b) - d*(rm_b)
        rhs_rtcb ~ (1/nb)*kc*rh*rm_b*tlr_el - lam*(rtcb)
        rhs_rm_r ~ tscr_r - lam*(rm_r) - d*(rm_r)
        rhs_rtcr ~ (1/nr)*kc*rh*rm_r*tlr_el - lam*(rtcr)
        rhs_trna ~ Vrep - Vdam + Vinflux - lam*(trna)
        rhs_rd ~ Vdam - Vtag - kdeg*rd - lam*(rd)
        rhs_rt ~ Vtag - Vrep - lam*(rt)

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


@mtkbuild rtc_trna = RTC_TRNA()

jac_sym=calculate_jacobian(rtc_trna)

init_ = [rtc_trna.rm_a=>0.0,rtc_trna.rtca=>0.0,rtc_trna.rm_b=>0.0,rtc_trna.rtcb=>0.0,rtc_trna.rm_r=>4e-6,rtc_trna.rtcr=>0.0,rtc_trna.trna=>0,rtc_trna.rd=>0.0,rtc_trna.rt=>0.0]

param_dict = Dict(L=>10, c=>0.001, kr=>0.125*12, Vmax_init=>39.51, Km_init=>250, θtscr=>160.01, θtlr=>255.73, na=>338, nb=>408, nr=>532*6, d=>0.2, krep=>137, ktag=>9780,
atp=>3000, km_a=>20, km_b=>16, g_max=>1260, kdeg=>0.00001, kin=>0.00175, ω_ab=>1e-5, ω_r=>1e-6, kdam=>0.1, lam=>0.014, kc=>0.6, k_diss=>0.006, rh=>11.29, thr_t=>10)

prob = ODEProblem(rtc_trna, init_, (0.0,1e9), param_dict; jac=true)

# latexify(equations(rtc)) |> render

equations(rtc_trna)
observed(rtc_trna)
# states(rtc)
# parameters(rtc)
solu = solve(prob)
df = create_solu_df(solu, trna_species)
plot([scatter(x=df.time, y=col, name="$(names(df)[i])", legendgroup="$i") for (col, i) in zip(eachcol(df[:,2:end]), range(2,length(names(df))))], Layout(xaxis_type="log", yaxis_tickformat=".2e", title="kdam = $(params_trna.kdam)"))

var_dict = Dict(rm_a(t)=>solu[rtc_trna.rm_a, end], rtca(t)=>solu[rtc_trna.rtca, end], rm_b(t)=>solu[rtc_trna.rm_b, end], rtcb(t)=>solu[rtc_trna.rtcb, end], rm_r(t)=>solu[rtc_trna.rm_r, end], rtcr(t)=>solu[rtc_trna.rtcr, end],
trna(t)=>solu[rtc_trna.trna, end], rd(t)=>solu[rtc_trna.rd, end], rt(t)=>solu[rtc_trna.rt, end])

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
eigs = eigvals(jac_mat)

for i in range(1,length(eigs))
    if eigs[i] < 0 
        println("negative")
    else
        push!(i)
    end
end




params_trna.kdam = 0.0
solu = sol(rtc_model_trna, init_trna, tspan, params_trna)
df = create_solu_df(solu, trna_species)
ssvals = ss_init_vals(df, trna_species)

br = get_br(rtc_mod_trna, params_trna_bf, ssvals, 1.)
bf = bf_point_df(br)
df = create_br_df(br)
kdam1 = findall(x->x==bf.kdam[1],df.kdam)[1]
kdam2 = findall(x->x==bf.kdam[2],df.kdam)[1]

unstable_df = df[kdam1:kdam2,:]
rename!(unstable_df, :rh=>:trna)
unstable_df_n = unstable_df[1:10:end,:]
length(unstable_df.kdam)


kdam_range = range(0,1,length=455)
kdam_range2 = range(1,0,length=455)

res_trna1 = numerical_bistability_analysis(rtc_model_trna, params_trna, ssvals, :trna, trna_species, kdam_range)
res_trna2 = numerical_bistability_analysis(rtc_model_trna, params_trna, ssvals, :trna, trna_species, kdam_range2)
res_rd1 = numerical_bistability_analysis(rtc_model_trna, params_trna, ssvals, :rd, trna_species, kdam_range)
res_rd2 = numerical_bistability_analysis(rtc_model_trna, params_trna, ssvals, :rd, trna_species, kdam_range2)
res_rt1 = numerical_bistability_analysis(rtc_model_trna, params_trna, ssvals, :rt, trna_species, kdam_range)
res_rt2 = numerical_bistability_analysis(rtc_model_trna, params_trna, ssvals, :rt, trna_species, kdam_range2)

res_rtcb1 = numerical_bistability_analysis(rtc_model_trna, params_trna, ssvals, :rtcb, trna_species, kdam_range)
res_rtcb2 = numerical_bistability_analysis(rtc_model_trna, params_trna, ssvals, :rtcb, trna_species, kdam_range2)
res_rtca1 = numerical_bistability_analysis(rtc_model_trna, params_trna, ssvals, :rtca, trna_species, kdam_range)
res_rtca2 = numerical_bistability_analysis(rtc_model_trna, params_trna, ssvals, :rtca, trna_species, kdam_range2)
res_rtcr1 = numerical_bistability_analysis(rtc_model_trna, params_trna, ssvals, :rtcr, trna_species, kdam_range)
res_rtcr2 = numerical_bistability_analysis(rtc_model_trna, params_trna, ssvals, :rtcr, trna_species, kdam_range2)

res_rtcb1_m = numerical_bistability_analysis(rtc_model_trna, params_trna, ssvals, :rm_b, trna_species, kdam_range)
res_rtcb2_m = numerical_bistability_analysis(rtc_model_trna, params_trna, ssvals, :rm_b, trna_species, kdam_range2)
res_rtca1_m = numerical_bistability_analysis(rtc_model_trna, params_trna, ssvals, :rm_a, trna_species, kdam_range)
res_rtca2_m = numerical_bistability_analysis(rtc_model_trna, params_trna, ssvals, :rm_a, trna_species, kdam_range2)
res_rtcr1_m = numerical_bistability_analysis(rtc_model_trna, params_trna, ssvals, :rm_r, trna_species, kdam_range)
res_rtcr2_m = numerical_bistability_analysis(rtc_model_trna, params_trna, ssvals, :rm_r, trna_species, kdam_range2)

ss_vals_kdam_ON = DataFrame("kdam"=>kdam_range,"trna"=>res_trna1,"rd"=>res_rd1, "rt"=>res_rt1, "rtca"=>res_rtca1, "rtcb"=>res_rtcb1, "rtcr"=>res_rtcr1, "rm_a"=>res_rtca1_m, "rm_b"=>res_rtcb1_m, "rm_r"=>res_rtcr1_m)

ss_vals_kdam_OFF = DataFrame("kdam"=>kdam_range2,"trna"=>res_trna2,"rd"=>res_rd2, "rt"=>res_rt2, "rtca"=>res_rtca2, "rtcb"=>res_rtcb2, "rtcr"=>res_rtcr2, "rm_a"=>res_rtca2_m, "rm_b"=>res_rtcb2_m, "rm_r"=>res_rtcr2_m)

function stability(ss_vals_kdam, kdam_length)
    unstable=[]
    for i in range(1, length(kdam_length))
        param_dict = Dict(L=>10, c=>0.001, kr=>0.125*12, Vmax_init=>39.51, Km_init=>250, θtscr=>160.01, θtlr=>255.73, na=>338, nb=>408, nr=>532*6, d=>0.2, krep=>137, ktag=>9780,
        atp=>3000, km_a=>20, km_b=>16, g_max=>1260, kdeg=>0.00001, kin=>0.00175, ω_ab=>1e-5, ω_r=>1e-6, kdam=>ss_vals_kdam[i,:kdam], lam=>0.014, kc=>0.6, k_diss=>0.006, rh=>11.29, thr_t=>10)

        var_dict = Dict(rm_a(t)=>ss_vals_kdam[i,:rm_a], rtca(t)=>ss_vals_kdam[i,:rtca], rm_b(t)=>ss_vals_kdam[i,:rm_b], rtcb(t)=>ss_vals_kdam[i,:rtcb], rm_r(t)=>ss_vals_kdam[i,:rm_r], rtcr(t)=>ss_vals_kdam[i,:rtcr],
        trna(t)=>ss_vals_kdam[i,:trna], rd(t)=>ss_vals_kdam[i,:rd], rt(t)=>ss_vals_kdam[i,:rt])

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
        eigs = eigvals(jac_mat)
        if all(<(0), real.(eigs)) == true
            println("stable")
        else
            println("unstable")
            push!(unstable, i)
            @show eigs
        end
    end
    return unstable
end

unstable_ON = stability(ss_vals_kdam_ON, kdam_range)
unstable_OFF = stability(ss_vals_kdam_OFF, kdam_range2)
unstable = stability(unstable_df, unstable_df.kdam)


unstab = scatter(x=unstable_df.kdam[unstable[1]:unstable[end]], y=unstable_df.trna[unstable[1]:unstable[end]], name="unstable")
stab3 = scatter(x=unstable_df.kdam[1:unstable[1]], y=unstable_df.trna[1:unstable[1]])
stab4 = scatter(x=unstable_df.kdam[unstable[end]:end], y=unstable_df.trna[unstable[end]:end])
stab1 = scatter(x=kdam_range, y=res_trna1, name="tRNA_h ON", line=attr(color="ffd30cff"))
stab2 = scatter(x=kdam_range2, y=res_trna2, name="tRNA_h OFF", line=attr(color="ffd30cff"))
plot([unstab, stab1, stab2, stab3, stab4])