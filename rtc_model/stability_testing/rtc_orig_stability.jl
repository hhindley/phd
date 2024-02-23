using ModelingToolkit, DifferentialEquations, PlotlyJS, Latexify, LinearAlgebra, DataFrames, LabelledArrays, Printf

include("/home/holliehindley/phd/general_funcs/solving.jl")
include("/home/holliehindley/phd/rtc_model/parameters/params.jl")
include("/home/holliehindley/phd/rtc_model/parameters/init.jl")
include("/home/holliehindley/phd/rtc_model/models/rtc_orig.jl")
include("/home/holliehindley/phd/rtc_model/functions/bf_funcs/bf_funcs.jl")


include("/home/holliehindley/phd/rtc_model/functions/bf_funcs/bf_funcs.jl");
include("/home/holliehindley/phd/rtc_model/models/rtc_orig.jl");
include("/home/holliehindley/phd/rtc_model/parameters/params.jl")
include("/home/holliehindley/phd/rtc_model/parameters/init.jl")


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

jac_sym=calculate_jacobian(rtc)

init_ = [rtc.rm_a=>0.0,rtc.rtca=>0.0,rtc.rm_b=>0.0,rtc.rtcb=>0.0,rtc.rm_r=>0.0,rtc.rtcr=>0.0,rtc.rh=>11.29,rtc.rd=>0.0,rtc.rt=>0.0]

param_dict = Dict(L=>10, c=>0.001, kr=>0.125, Vmax_init=>39.51, Km_init=>250, θtscr=>160.01, θtlr=>255.73, na=>338, nb=>408, nr=>532*6, d=>0.2, krep=>137, ktag=>9780,
atp=>3000, km_a=>20, km_b=>16, g_max=>1260, kdeg=>0.001, kin=>0.022/100, ω_ab=>1e-5, ω_r=>1e-6, kdam=>0, lam=>0.014, kc=>0.6, k_diss=>0.006)

prob = ODEProblem(rtc, init_, (0.0,1e9), param_dict; jac=true)

odefun = prob.f
F = (u,p) -> odefun(u,p,0)
J = (u,p) -> odefun.jac(u,p,0)
id_kdam = indexof(kdam, parameters(rtc))
par_tm = prob.p
prob_bf = BifurcationProblem(F, prob.u0, setproperties(par_tm), (@lens _[id_kdam]); record_from_solution = (x, p) -> (rm_a = x[1], rtca = x[2], rm_b = x[3], rtcb = x[4], rm_r = x[5], rtcr = x[6], trna = x[7], rd = x[8], rt = x[9]),)
opts_br = ContinuationPar(p_min = 0., p_max = 1.5, ds = 0.001, a=0.1,
        dsmax = 0.05, # 0.15
        # options to detect bifurcations
        detect_bifurcation = 3, n_inversion = 4, max_bisection_steps = 20, #3,2,10
        # number of eigenvalues
        nev = 2, 
        # maximum number of continuation steps
        max_steps = 50000,)
    
br = continuation(prob_bf, PALC(θ=0.5), opts_br; plot = false, bothside=true, normC = norminf)
plot(scatter(x=df_bfk.kdam, y=df_bfk.rh))
# latexify(equations(rtc)) |> render

# equations(rtc)
# observed(rtc)
# states(rtc)
# parameters(rtc)
solu = solve(prob)



var_dict = Dict(rm_a(t)=>solu[rtc.rm_a, end], rtca(t)=>solu[rtc.rtca, end], rm_b(t)=>solu[rtc.rm_b, end], rtcb(t)=>solu[rtc.rtcb, end], rm_r(t)=>solu[rtc.rm_r, end], rtcr(t)=>solu[rtc.rtcr, end],
rh(t)=>solu[rtc.rh, end], rd(t)=>solu[rtc.rd, end], rt(t)=>solu[rtc.rt, end])

ssvals = [var_dict[rm_a(t)], var_dict[rtca(t)], var_dict[rm_b(t)], var_dict[rtcb(t)], var_dict[rm_r(t)], var_dict[rtcr(t)], var_dict[rh(t)], var_dict[rd(t)], var_dict[rt(t)]]

br = get_br(rtc_mod, params_bf, init_rtc, 1.5)

# using Plots
Plots.plot(br, linewidthstable=5)

bf = bf_point_df(br)
df_bfk = create_br_df(br)
kdam1 = findall(x->x==bf.kdam[1],df_bfk.kdam)[1]
kdam2 = findall(x->x==bf.kdam[2],df_bfk.kdam)[1]

# unstable_df = df[kdam1:kdam2,:]
# stable1_df = df[1:kdam1,:]
# stable2_df = df[kdam2:end,:]

# unstable_df_n = unstable_df[1:10:end,:]
# length(unstable_df.kdam)


kdam_range = range(0,1.5,length=50)
kdam_range2 = range(1.5,0,length=50)

res_trna1 = numerical_bistability_analysis(rtc_model, params_rtc, ssvals, :rh, species_rtc, kdam_range)
res_trna2 = numerical_bistability_analysis(rtc_model, params_rtc, ssvals, :rh, species_rtc, kdam_range2)
res_rd1 = numerical_bistability_analysis(rtc_model, params_rtc, ssvals, :rd, species_rtc, kdam_range)
res_rd2 = numerical_bistability_analysis(rtc_model, params_rtc, ssvals, :rd, species_rtc, kdam_range2)
res_rt1 = numerical_bistability_analysis(rtc_model, params_rtc, ssvals, :rt, species_rtc, kdam_range)
res_rt2 = numerical_bistability_analysis(rtc_model, params_rtc, ssvals, :rt, species_rtc, kdam_range2)

res_rtcb1 = numerical_bistability_analysis(rtc_model, params_rtc, ssvals, :rtcb, species_rtc, kdam_range)
res_rtcb2 = numerical_bistability_analysis(rtc_model, params_rtc, ssvals, :rtcb, species_rtc, kdam_range2)
res_rtca1 = numerical_bistability_analysis(rtc_model, params_rtc, ssvals, :rtca, species_rtc, kdam_range)
res_rtca2 = numerical_bistability_analysis(rtc_model, params_rtc, ssvals, :rtca, species_rtc, kdam_range2)
res_rtcr1 = numerical_bistability_analysis(rtc_model, params_rtc, ssvals, :rtcr, species_rtc, kdam_range)
res_rtcr2 = numerical_bistability_analysis(rtc_model, params_rtc, ssvals, :rtcr, species_rtc, kdam_range2)

res_rtcb1_m = numerical_bistability_analysis(rtc_model, params_rtc, ssvals, :rm_b, species_rtc, kdam_range)
res_rtcb2_m = numerical_bistability_analysis(rtc_model, params_rtc, ssvals, :rm_b, species_rtc, kdam_range2)
res_rtca1_m = numerical_bistability_analysis(rtc_model, params_rtc, ssvals, :rm_a, species_rtc, kdam_range)
res_rtca2_m = numerical_bistability_analysis(rtc_model, params_rtc, ssvals, :rm_a, species_rtc, kdam_range2)
res_rtcr1_m = numerical_bistability_analysis(rtc_model, params_rtc, ssvals, :rm_r, species_rtc, kdam_range)
res_rtcr2_m = numerical_bistability_analysis(rtc_model, params_rtc, ssvals, :rm_r, species_rtc, kdam_range2)

ss_vals_kdam_ON = DataFrame("kdam"=>kdam_range,"rh"=>res_trna1,"rd"=>res_rd1, "rt"=>res_rt1, "rtca"=>res_rtca1, "rtcb"=>res_rtcb1, "rtcr"=>res_rtcr1, "rm_a"=>res_rtca1_m, "rm_b"=>res_rtcb1_m, "rm_r"=>res_rtcr1_m)

ss_vals_kdam_OFF = DataFrame("kdam"=>kdam_range2,"rh"=>res_trna2,"rd"=>res_rd2, "rt"=>res_rt2, "rtca"=>res_rtca2, "rtcb"=>res_rtcb2, "rtcr"=>res_rtcr2, "rm_a"=>res_rtca2_m, "rm_b"=>res_rtcb2_m, "rm_r"=>res_rtcr2_m)

function stability(ss_vals_kdam, kdam_length)
    unstable=[]
    for i in range(1, length(kdam_length))
        param_dict_dam = deepcopy(param_dict)
        param_dict_dam[kdam] = ss_vals_kdam[i,:kdam]

        var_dict = Dict(rm_a(t)=>ss_vals_kdam[i,:rm_a], rtca(t)=>ss_vals_kdam[i,:rtca], rm_b(t)=>ss_vals_kdam[i,:rm_b], rtcb(t)=>ss_vals_kdam[i,:rtcb], rm_r(t)=>ss_vals_kdam[i,:rm_r], rtcr(t)=>ss_vals_kdam[i,:rtcr],
        rh(t)=>ss_vals_kdam[i,:rh], rd(t)=>ss_vals_kdam[i,:rd], rt(t)=>ss_vals_kdam[i,:rt])

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

trna_p = plot_stability(:rh, res_trna1, res_trna2)
trna_p = plot_stability(:rtcb, res_rtcb1, res_rtcb2)
trna_p = plot_stability(:rm_b, res_rtcb1_m, res_rtcb2_m)




x=40
function check_stab(x,n)
    ss = [df_bfk[x,:rm_a],df_bfk[x,:rtca],df_bfk[x,:rm_b],df_bfk[x,:rtcb],df_bfk[x,:rm_r],df_bfk[x,:rtcr],df_bfk[x,:rh],df_bfk[x,:rd],df_bfk[x,:rt]]
    ss = @. ss*n
    ps = deepcopy(param_dict)
    ps[kdam] = df_bfk[x,:kdam]

    prob = ODEProblem(rtc, ss, (0.0,1e9), ps; jac=true)

    solu = solve(prob)

    var_dict = Dict(rm_a(t)=>solu[rtc.rm_a, end], rtca(t)=>solu[rtc.rtca, end], rm_b(t)=>solu[rtc.rm_b, end], rtcb(t)=>solu[rtc.rtcb, end], rm_r(t)=>solu[rtc.rm_r, end], rtcr(t)=>solu[rtc.rtcr, end],
    rh(t)=>solu[rtc.rh, end], rd(t)=>solu[rtc.rd, end], rt(t)=>solu[rtc.rt, end])


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