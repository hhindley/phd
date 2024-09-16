using ModelingToolkit, DifferentialEquations, PlotlyJS, LinearAlgebra, DataFrames, LabelledArrays, Printf, BifurcationKit, OrderedCollections, ProgressBars

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

kdam_range = range(0,5, length=100) 
df_ssvals = var_param(test, kdam, params_rtc1, kdam_range, ssvals_rtc)
plot(scatter(x=kdam_range, y=df_ssvals.rtca), Layout(xaxis_title="kdam", yaxis_title="rtca"))

df_init = var_param(test, kdam, params_rtc1, kdam_range, init_rtc)
plot(scatter(x=kdam_range, y=df_init.rtca), Layout(xaxis_title="kdam", yaxis_title="rtca"))

ssvals_diff = deepcopy(ssvals_rtc)
ssvals_diff[2] = 0
df_diff = var_param(test, kdam, params_rtc1, kdam_range, ssvals_diff)
plot(scatter(x=kdam_range, y=df_diff.rtca), Layout(xaxis_title="kdam", yaxis_title="rtca"))

res=[]
minus_num = 20
for i in eachindex(ssvals_rtc)
    ssvals_diff = deepcopy(ssvals_rtc)
    if ssvals_diff[i]-minus_num < 0
        ssvals_diff[i] = 0
    else 
        ssvals_diff[i] = ssvals_diff[i]-minus_num
    end
    println(ssvals_diff)
    df_diff = var_param(test, kdam, params_rtc1, kdam_range, ssvals_diff)
    push!(res, df_diff)
end
plot([scatter(x=kdam_range, y=res[i].rtca, name="$(species_rtc[i]) - $minus_num") for i in eachindex(res)])

res1=[]
minus_num = 20
for i in eachindex(ssvals_rtc)
    ssvals_diff = deepcopy(ssvals_rtc)
    if ssvals_diff[i]-minus_num < 0 
        ssvals_diff[i] = 0
    else 
        ssvals_diff[i] = ssvals_diff[i]-minus_num
    end
    if i+1 < length(ssvals_diff)
        if ssvals_diff[i+1]-minus_num < 0 
            ssvals_diff[i+1] = 0
        else 
            ssvals_diff[i+1] = ssvals_diff[i+1]-minus_num
        end
        if ssvals_diff[i+2]-minus_num < 0 
            ssvals_diff[i+2] = 0
        else 
            ssvals_diff[i+2] = ssvals_diff[i+2]-minus_num
        end
        println(ssvals_diff)
        df_diff = var_param(test, kdam, params_rtc1, kdam_range, ssvals_diff)
        push!(res1, df_diff)
    end
end
res1
plot([scatter(x=kdam_range, y=res1[i], name="$(species_rtc[i]) - $minus_num") for i in eachindex(res1)])

using Combinatorics

res1 = []
combs = []
minus_num = 2
num_conditions = length(ssvals_rtc)
for k in eachindex(ssvals_rtc)
    println(k)
    for comb in combinations(1:num_conditions, k)
        # println(comb)
        ssvals_diff = deepcopy(ssvals_rtc)
        for i in comb
            # println(i)
            if ssvals_diff[i] - minus_num < 0
                ssvals_diff[i] = 0
            else
                ssvals_diff[i] = ssvals_diff[i] - minus_num
            end
        end

        df_diff = var_param(test, kdam, params_rtc1, kdam_range, ssvals_diff)

        if round(df_diff.rtca[end], digits=3) != round(df_ssvals.rtca[end], digits=3)
            # println(ssvals_diff)
            push!(res1, df_diff)
            push!(combs, comb)
        end
        
    end
end

p = plot([scatter(x=kdam_range, y=res1[i].rtca, name="change $([species_rtc[s] for s in combs[i]])") for i in eachindex(res1)], Layout(xaxis_title="kdam",yaxis_title="rtca",title="$([round(ssvals_rtc[i], digits=6) for i in eachindex(ssvals_rtc)]) - $minus_num"))
open("/Users/s2257179/Desktop/minus2.html", "w") do io
    PlotlyBase.to_html(io, p.plot)
end


L_range = [10, 100, 1000, 10000, 100000]
c_range = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2]

L_ssvals = []
for i in L_range
    L_params = deepcopy(params_rtc1)
    L_params[L] = i
    push!(L_ssvals, var_param(test, kdam, L_params, kdam_range, ssvals_rtc))
end
plot([scatter(x=kdam_range, y=i.rtca, name="L = $(L_range[j])", line=attr(color=colours[j])) for (i,j) in zip(L_ssvals, 1:length(L_range))], Layout(xaxis_title="kdam", yaxis_title="rtca"))

c_ssvals = []
for i in c_range
    c_params = deepcopy(params_rtc1)
    c_params[c] = i
    push!(c_ssvals, var_param(test, kdam, c_params, kdam_range, ssvals_rtc))
end
plot([scatter(x=kdam_range, y=i.rtca, name="c = $(c_range[j])", line=attr(color=colours[j])) for (i,j) in zip(c_ssvals, 1:length(c_range))], Layout(xaxis_title="kdam", yaxis_title="rtca"))

Lc_ssvals = []
Lc_vals = []
for i in L_range
    for j in c_range
        Lc_params = deepcopy(params_rtc1)
        Lc_params[L] = i
        Lc_params[c] = j
        push!(Lc_vals, (i,j))
        push!(Lc_ssvals, var_param(test, kdam, Lc_params, kdam_range, ssvals_rtc))
    end
end

plot([scatter(x=kdam_range, y=Lc_ssvals[i].rtca, name="L,c = $(Lc_vals[i])") for i in eachindex(Lc_vals)], Layout(xaxis_title="kdam", yaxis_title="rtca"))

# concentration no damage 
solu_rtc = sol(test, init_rtc, tspan, params_rtc1)
df = create_solu_df(solu_rtc, species_rtc)
p_rtc1 = plot([scatter(x=df.time, y=col, name="$(names(df)[i])", legendgroup="$i", marker_color=colours[i]) for (col, i) in zip(eachcol(df[:,2:end]), range(2,length(names(df))))], Layout(xaxis_type="log", title="kdam = $(params_rtc1[kdam])", xaxis_title="Time (s)", yaxis_title="Concentration (μM)"))

p = plot([scatter(x=df.time, y=df.rh*lam_c_val*tlr1, name="λ"),
scatter(x=df.time, y=df.rh*kin_c_val, yaxis="y2", name="kin")],
Layout(xaxis_title_text="Time (s)", yaxis_title_text="λ (min-1)", xaxis_type="log", 
yaxis2=attr(title="kin (μM aa-1)", overlaying="y", side="right", tickformat=".5f"),
legend=attr(x=0.1, y=0.9)))

open("/Users/s2257179/Desktop/dynamic_lam_kin.html", "w") do io
    PlotlyBase.to_html(io, p.plot)
end

params1 = deepcopy(params_rtc1)
params1[kdam] = 1
solu_rtc = sol(test, init_rtc, tspan, params1)
df = create_solu_df(solu_rtc, species_rtc)

solu_rtc1 = sol(test, ssvals_rtc, tspan, params1)
df1 = create_solu_df(solu_rtc1, species_rtc)

ssvals_rtc


species = :rtca
p = plot([scatter(x=df.time, y=df[:,species], name="init"),
scatter(x=df1.time, y=df1[:,species], name="ss")],
Layout(xaxis_title_text="Time (min)", yaxis_title_text="$species", xaxis_type="log", 
legend=attr(x=0.1, y=0.9)))

init_rtc2 = [test.rm_a=>0.0,test.rtca=>1.0,test.rm_b=>0.0,test.rtcb=>0.0,test.rm_r=>0.0,test.rtcr=>0.0,test.rh=>11.29,test.rd=>0.0,test.rt=>0.0]
solu_rtc2 = sol(test, init_rtc2, tspan, params1)
df2 = create_solu_df(solu_rtc2, species_rtc)

init_rtc3 = [test.rm_a=>0.0,test.rtca=>1.0,test.rm_b=>0.0,test.rtcb=>0.0,test.rm_r=>0.0,test.rtcr=>0.0,test.rh=>19.29,test.rd=>0.0,test.rt=>0.0]
solu_rtc3 = sol(test, init_rtc3, tspan, params1)
df3 = create_solu_df(solu_rtc3, species_rtc)

species = :rm_a
p = plot([scatter(x=df.time, y=df[:,species], name="init"),
scatter(x=df1.time, y=df1[:,species], name="ss"),
scatter(x=df2.time, y=df2[:,species], name="test"),
scatter(x=df3.time, y=df3[:,species], name="test2")],
Layout(xaxis_title_text="Time (min)", yaxis_title_text="$species", xaxis_type="log", 
legend=attr(x=0.1, y=0.9)))


p = plot([scatter(x=df.time, y=df.rh*lam_c_val*tlr1, name="λ"),
scatter(x=df.time, y=df.rh*kin_c_val, yaxis="y2", name="kin"),
scatter(x=df1.time, y=df1.rh*lam_c_val*tlr1, name="λ ss"),
scatter(x=df1.time, y=df1.rh*kin_c_val, yaxis="y2", name="kin ss"),
scatter(x=df2.time, y=df2.rh*lam_c_val*tlr1, name="λ 15"),
scatter(x=df2.time, y=df2.rh*kin_c_val, yaxis="y2", name="kin 15")],
Layout(xaxis_title_text="Time (min)", yaxis_title_text="λ (min-1)", xaxis_type="log", 
yaxis2=attr(title="kin (μM aa-1)", overlaying="y", side="right", tickformat=".5f"),
legend=attr(x=0.1, y=0.9)))


p = plot([scatter(x=df.time, y=df.rtca, name="init"),
scatter(x=df1.time, y=df1.rtca, name="ss"),
scatter(x=df2.time, y=df2.rtca, name="test")],
Layout(xaxis_title_text="Time (min)", yaxis_title_text="rtca", xaxis_type="log", 
legend=attr(x=0.1, y=0.9)))

plot([scatter(x=df.time, y=df.rh, name="init"), scatter(x=df1.time, y=df1.rh, name="ss")], Layout(xaxis_type="log"))


kin_val
2.2e-4

init_rtc2 = [test.rm_a=>0.0,test.rtca=>0.000001,test.rm_b=>0.0,test.rtcb=>0.0,test.rm_r=>0.0,test.rtcr=>0.0,test.rh=>20.29,test.rd=>0.0,test.rt=>0.0]

new_params = deepcopy(params_rtc1)
kdam_range = range(0,5, length=100) 
ssvals=[]
for i in kdam_range
    new_params[kdam] = i
    solu_rtc = sol(test, init_rtc2, tspan, new_params)
    df = create_solu_df(solu_rtc, species_rtc)
    ss = [i[end] for i in eachcol(df[:,2:end])]
    # ss = steady_states(test, init_rtc, new_params)
    push!(ssvals, ss)
end
df_ssvals = DataFrame(vcat(transpose(ssvals)...), :auto)
rename!(df_ssvals, species_rtc)
plot(scatter(x=kdam_range, y=df_ssvals.rh), Layout(xaxis_title="kdam", yaxis_title="rh"))


ssvals1=[]
for i in kdam_range
    new_params[kdam] = i
    solu_rtc = sol(test, ssvals_rtc, tspan, new_params)
    df = create_solu_df(solu_rtc, species_rtc)
    ss = [i[end] for i in eachcol(df[:,2:end])]
    # ss = steady_states(test, init_rtc, new_params)
    push!(ssvals1, ss)
end
df_ssvals1 = DataFrame(vcat(transpose(ssvals1)...), :auto)
rename!(df_ssvals1, species_rtc)
plot(scatter(x=kdam_range, y=df_ssvals1.rh), Layout(xaxis_title="kdam", yaxis_title="rh"))

plot([scatter(x=kdam_range, y=df_ssvals.rh, name="init"), scatter(x=kdam_range, y=df_ssvals1.rh, name="ss")], Layout(xaxis_title="kdam", yaxis_title="rh"))


kdam_range = range(0,10, length=100) 
kdam_range1 = reverse(kdam_range)
res = numerical_bistability_analysis(test, params_rtc1, init_rtc, :rtca, species_rtc, kdam_range, kdam)
res1 = numerical_bistability_analysis(test, params_rtc1, init_rtc, :rtca, species_rtc, kdam_range1, kdam)

plot([scatter(x=kdam_range, y=res), scatter(x=kdam_range1, y=res1)])


br = get_br(test, ssvals_rtc, params_rtc1, 1.5)
bf = bf_point_df(br)
df = create_br_df(br)
kdam1 = findall(x->x==bf.kdam[1],df.kdam)[1]
kdam2 = findall(x->x==bf.kdam[2],df.kdam)[1]
p_conc = plot([scatter(x=df.kdam, y=df.rtcb, name="RtcB", line=attr(color=:blue)), scatter(x=df.kdam, y=df.rtca, name="RtcA", line=attr(color=:red))], Layout(xaxis_title="Damage rate (min<sup>-1</sup>)",yaxis_title="Concentration μM"))


prob = ODEProblem(test, ssvals_rtc, tspan, params_rtc1; jac=true)
odefun = prob.f
F = (u,p) -> odefun(u,p,0)
J = (u,p) -> odefun.jac(u,p,0)
par_tm = prob.p[1]
# id_kdam = indexof(kdam, parameters(model))
id_kdam = indexof(0.0, par_tm)
# Bifurcation Problem

prob = BifurcationProblem(F, prob.u0, (par_tm), (@lens _[id_kdam]); J=J,
record_from_solution = (x, p) -> (rm_a = x[1], rtca = x[2], rm_b = x[3], rtcb = x[4], rm_r = x[5], rtcr = x[6], rh = x[7], rd = x[8], rt = x[9]),)
opts_br = ContinuationPar(p_min = 0., p_max = 1., ds = 0.001, #a=0.1,
dsmax = 0.15, dsmin = 0.0001, # 0.15
# options to detect bifurcations
detect_bifurcation = 3, n_inversion = 2, max_bisection_steps = 20, #3,2,10
# number of eigenvalues
nev = 2, 
# maximum number of continuation steps
max_steps = 50000,)# dsminBisection=1e-30, tolBisectionEigenvalue=1e-30)# a=0.9, )
# tolStability=1e-10, tolBisectionEigenvalue=1e-10)#,tolParamBisectionEvent=1e-1)
# only using parameters that make a difference to solution
# continuation of equilibria
br = continuation(prob, PALC(θ=0.1), opts_br; plot = false, bothside=true, normC = norminf)

df1 = create_br_df(br)

plot([scatter(x=df1.kdam, y=df1.rtca, name="RtcA", line=attr(color=:red)), scatter(x=df1.kdam, y=df1.rtcb, name="RtcB", line=attr(color=:blue))], Layout(xaxis_title="Damage rate (min<sup>-1</sup>)",yaxis_title="Concentration μM"))