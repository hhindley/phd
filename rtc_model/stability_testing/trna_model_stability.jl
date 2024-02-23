using ModelingToolkit, DifferentialEquations, PlotlyJS, Latexify, LinearAlgebra, DataFrames, LabelledArrays, Printf

include("/home/holliehindley/phd/general_funcs/solving.jl")
include("/home/holliehindley/phd/rtc_model/parameters/trna_params.jl")
include("/home/holliehindley/phd/rtc_model/models/rtc_trna_model.jl")
include("/home/holliehindley/phd/rtc_model/functions/bf_funcs/bf_funcs.jl")

indexof(sym,syms) = findfirst(isequal(sym),syms)
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
        Vinflux ~ kin* g_max*atp/(θtlr+atp) #* trna/(thr_t+trna)  # uM min-1 
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

init_ = [rtc_trna.rm_a=>0.0,rtc_trna.rtca=>0.0,rtc_trna.rm_b=>0.0,rtc_trna.rtcb=>0.0,rtc_trna.rm_r=>0.0,rtc_trna.rtcr=>0.0,rtc_trna.trna=>135.5,rtc_trna.rd=>0.0,rtc_trna.rt=>0.0]

param_dict = Dict(L=>10, c=>0.001, kr=>0.125*12, Vmax_init=>39.51, Km_init=>250, θtscr=>160.01, θtlr=>255.73, na=>338, nb=>408, nr=>532*6, d=>0.2, krep=>137, ktag=>9780,
atp=>3000, km_a=>20, km_b=>16, g_max=>1260, kdeg=>0.00001, kin=>0.00175, ω_ab=>1e-5, ω_r=>1e-6, kdam=>0., lam=>0.014, kc=>0.6, k_diss=>0.006, rh=>11.29, thr_t=>10)

solu = sol(rtc_trna, init_, (0,1e9), param_dict)

prob_trna = ODEProblem(rtc_trna, init_, (0.0,1e9), param_dict; jac=true)
odefun = prob_trna.f
F = (u,p) -> odefun(u,p,0)
J = (u,p) -> odefun.jac(u,p,0)
id_kdam = indexof(kdam, parameters(rtc_trna))
par_tm = prob_trna.p
prob_bf = BifurcationProblem(F, prob_trna.u0, setproperties(par_tm), (@lens _[id_kdam]); record_from_solution = (x, p) -> (rm_a = x[1], rtca = x[2], rm_b = x[3], rtcb = x[4], rm_r = x[5], rtcr = x[6], trna = x[7], rd = x[8], rt = x[9]),)
opts_br = ContinuationPar(p_min = 0., p_max = 4000., ds = 0.000161026,# a=0.1,
        dsmax = 0.00573615, dsmin = 0.000161026,# 0.15
        # options to detect bifurcations
        detect_bifurcation = 3, n_inversion = 2, max_bisection_steps = 10, #3,2,10
        # number of eigenvalues
        # nev =100, #tolParamBisectionEvent=1e-30, 
        # maximum number of continuation steps
        max_steps = 1000,)# dsmin_bisection=1e-30)#, tol_bisection_eigenvalue=1e-10)# a=0.9, )

    
br = continuation(prob_bf, PALC(θ=0.000529832), opts_br; plot = false, bothside=true, normC = norminf)

br = get_br(rtc_trna, param_dict, init_, 4000.)

# latexify(equations(rtc)) |> render

equations(rtc_trna)
observed(rtc_trna)
# states(rtc)
# parameters(rtc)
solu = solve(prob_trna)
df = create_solu_df(solu, trna_species)
plot([scatter(x=df.time, y=col, name="$(names(df)[i])", legendgroup="$i") for (col, i) in zip(eachcol(df[:,2:end]), range(2,length(names(df))))], Layout(xaxis_type="log", yaxis_tickformat=".2e", title="kdam = $(params_trna.kdam)"))

var_dict = Dict(rm_a(t)=>solu[rtc_trna.rm_a, end], rtca(t)=>solu[rtc_trna.rtca, end], rm_b(t)=>solu[rtc_trna.rm_b, end], rtcb(t)=>solu[rtc_trna.rtcb, end], rm_r(t)=>solu[rtc_trna.rm_r, end], rtcr(t)=>solu[rtc_trna.rtcr, end],
trna(t)=>solu[rtc_trna.trna, end], rd(t)=>solu[rtc_trna.rd, end], rt(t)=>solu[rtc_trna.rt, end])


ssvals = [var_dict[rm_a(t)], var_dict[rtca(t)], var_dict[rm_b(t)], var_dict[rtcb(t)], var_dict[rm_r(t)], var_dict[rtcr(t)], var_dict[trna(t)], var_dict[rd(t)], var_dict[rt(t)]]















br = get_br(rtc_mod_trna, params_trna_bf, ssvals, 4000.)

# using Plots
# Plots.plot(brs[17], linewidthstable=5)
# Plots.plot(brs[18], linewidthstable=5)
Plots.plot(br, vars = (:param, :trna), linewidthstable=5)
bf = bf_point_df(br)
df_bfk = create_br_df(br)
kdam1 = findall(x->x==bf.kdam[1],df_bfk.kdam)[1]
kdam2 = findall(x->x==bf.kdam[2],df_bfk.kdam)[1]


rename!(df_bfk, :rh=>:trna)

plot(scatter(x=df_bfk.kdam, y=df_bfk.trna))

# unstable_df = df[kdam1:kdam2,:]
# stable1_df = df[1:kdam1,:]
# stable2_df = df[kdam2:end,:]

# unstable_df_n = unstable_df[1:10:end,:]
# length(unstable_df.kdam)


kdam_range = range(0,4000,length=100)
kdam_range2 = range(4000,0,length=100)

res_trna1 = numerical_bistability_analysis(rtc_model_trna, params_trna, ssvals, :trna, trna_species, kdam_range)
res_trna2 = numerical_bistability_analysis(rtc_model_trna, params_trna, ssvals, :trna, trna_species, kdam_range2)
plot([scatter(x=kdam_range, y=res_trna1), scatter(x=kdam_range2, y=res_trna2)])

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
        param_dict_dam = deepcopy(param_dict)
        param_dict_dam[kdam] = ss_vals_kdam[i,:kdam]

        var_dict = Dict(rm_a(t)=>ss_vals_kdam[i,:rm_a], rtca(t)=>ss_vals_kdam[i,:rtca], rm_b(t)=>ss_vals_kdam[i,:rm_b], rtcb(t)=>ss_vals_kdam[i,:rtcb], rm_r(t)=>ss_vals_kdam[i,:rm_r], rtcr(t)=>ss_vals_kdam[i,:rtcr],
        trna(t)=>ss_vals_kdam[i,:trna], rd(t)=>ss_vals_kdam[i,:rd], rt(t)=>ss_vals_kdam[i,:rt])

        all_vals = merge(var_dict, param_dict_dam)

        vals=[]
        for i in 1:size(jac_sym, 1)
            for j in 1:size(jac_sym, 2)
                new::String = repr(jac_sym[i,j])
                jac_sym[i,j] = eval(Meta.parse(new))
                push!(vals, substitute(jac_sym[i,j], all_vals))
            end
        end
        
        jac_mat = parse.(Float64,string.(transpose(reshape(vals, (9,9)))))
        
        eigs = eigvals(jac_mat)
        
        if all(<(0), real.(eigs)) == true
            # println("stable")
            @show eigs
        else
            # println("unstable")
            push!(unstable, i)
            @show eigs
        end
    end
    if length(unstable) > 0 
        Nothing
    else
        push!(unstable, "stable")
    end
    return unstable
end

unstable_ON = stability(ss_vals_kdam_ON, kdam_range)
unstable_OFF = stability(ss_vals_kdam_OFF, kdam_range2)
bfk = stability(df_bfk, df_bfk.kdam)

function plot_stability(trna, res_trna1, res_trna2)
    unstab = scatter(x=df_bfk.kdam[bfk[1]:bfk[end]], y=df_bfk[!,trna][bfk[1]:bfk[end]], name="bf-kit dash unstable", legendgroup=1)
    stab3 = scatter(x=df_bfk.kdam[1:bfk[1]], y=df_bfk[!,trna][1:bfk[1]], name="bf-kit dash stable", legendgroup=1, showlegend=true)
    stab4 = scatter(x=df_bfk.kdam[bfk[end]:end], y=df_bfk[!,trna][bfk[end]:end], name="bf-kit dash stable", legendgroup=1, showlegend=true)
    stab1 = scatter(x=kdam_range, y=res_trna1, name="numerical $(unstable_ON[1])", line=attr(color="ffd30cff"), legendgroup=3)
    stab2 = scatter(x=kdam_range2, y=res_trna2, name="numerical $(unstable_OFF[1])", line=attr(color="ffd30cff"), legendgroup=3)

    return plot([unstab, stab1, stab2, stab3, stab4])
end

trna_p = plot_stability(:trna, res_trna1, res_trna2)
trna_p = plot_stability(:rtcb, res_rtcb1, res_rtcb2)
trna_p = plot_stability(:rm_b, res_rtcb1_m, res_rtcb2_m)




x=40
function check_stab(x,n, model, df_bfk, param_dict)
    ss = [df_bfk[x,:rm_a],df_bfk[x,:rtca],df_bfk[x,:rm_b],df_bfk[x,:rtcb],df_bfk[x,:rm_r],df_bfk[x,:rtcr],df_bfk[x,:trna],df_bfk[x,:rd],df_bfk[x,:rt]]
    ss = @. ss*n
    ps = deepcopy(param_dict)
    ps[kdam] = df_bfk[x,:kdam]

    prob = ODEProblem(rtc_trna, ss, (0.0,1e12), ps; jac=true)

    solu = solve(prob)

    var_dict = Dict(rm_a(t)=>solu[rtc_trna.rm_a, end], rtca(t)=>solu[rtc_trna.rtca, end], rm_b(t)=>solu[rtc_trna.rm_b, end], rtcb(t)=>solu[rtc_trna.rtcb, end], rm_r(t)=>solu[rtc_trna.rm_r, end], rtcr(t)=>solu[rtc_trna.rtcr, end],
    trna(t)=>solu[rtc_trna.trna, end], rd(t)=>solu[rtc_trna.rd, end], rt(t)=>solu[rtc_trna.rt, end])


    return [scatter(x=[df_bfk.kdam[x]], y=[df_bfk.rtcb[x]], name="old", legendgroup=x, mode="markers", marker=attr(color=:purple)), scatter(x=[df_bfk.kdam[x]], y=[ss[4]], name="new", legendgroup=x, mode="markers", marker=attr(color=:green)), scatter(x=[df_bfk.kdam[x]], y=[var_dict[rtcb(t)]], name="solu", mode="markers", legendgroup=x, marker=attr(color=:red))]
    
end

scatters=[]
for i in range(1,length(br.sol))
    push!(scatters, check_stab(i,1))
end

scatters = reduce(vcat, scatters)
pushfirst!(scatters, scatter(x=df_bfk.kdam, y=df_bfk.rtcb))
# pushfirst!(scatters, scatter(x=ss_vals_kdam_ON.kdam, y=ss_vals_kdam_ON.rtcb))
plot([i for i in scatters])

plot([scatter(x=df_bfk.kdam, y=df_bfk.rtcb), scatter(x=[df_bfk.kdam[x]], y=[df_bfk.rtcb[x]], name="old"), scatter(x=[df_bfk.kdam[x]], y=[ss[4]], name="new", mode="markers", marker=attr(color=:blue)), scatter(x=[df_bfk.kdam[x]], y=[var_dict[rtcb(t)]], name="solu")])