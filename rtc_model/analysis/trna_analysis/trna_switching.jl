using Parameters, CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, Statistics
using Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra, Printf
# using Plots
using PlotlyJS, ProgressBars
include("/home/holliehindley/phd/may23_rtc/functions/solving.jl"); include("/home/holliehindley/phd/may23_rtc/functions/set_ups.jl"); include("/home/holliehindley/phd/may23_rtc/functions/plotting.jl"); 
include("/home/holliehindley/phd/may23_rtc/functions/sweep_params.jl"); include("/home/holliehindley/phd/may23_rtc/models/rtc_orig.jl"); include("/home/holliehindley/phd/may23_rtc/models/atp_lam_kin_t.jl"); 
include("/home/holliehindley/phd/may23_rtc/models/single_t.jl"); include("/home/holliehindley/phd/may23_rtc/models/combinations_t.jl"); 
include("/home/holliehindley/phd/may23_rtc/functions/bf_funcs/bf_funcs.jl");
include("/home/holliehindley/phd/colors_plotly.jl"); include("/home/holliehindley/phd/may23_rtc/models/rtc_trna_model.jl")
include("/home/holliehindley/phd/may23_rtc/functions/bf_funcs/init_switch_funcs.jl"); include("/home/holliehindley/phd/may23_rtc/models/inhibition_models/rtc_inhibition_model.jl");

@consts begin
    L = 10; #10 
    c = 0.001; 
    kr = 0.125; 
    Vmax_init = 39.51; 
    Km_init = 250; 
    θtscr = 160.01;  
    θtlr = 255.73; 
    # k_b = 17.7; 
    na = 338; 
    nb = 408; 
    nr = 532*6;
    d = 0.2; 
    krep = 137; 
    ktag = 9780;#0.1; 
    # atp = 4000;#2500; 
    km_a = 20; 
    km_b = 16;
    g_max = 2.0923; 
    gr_c = 0.0008856; # 0.000599; 
    kdeg = 0.001; 
    kin = 0.022222222 #0.054; #2.381 
    ω_ab = 4#4#0.093; #0.0828304057748932;#4; 
    ω_r = 0.0019*6#2e-7 #0.0019*6 #70.53; #0.0019*6#79.43865871861044; #0.0019*6;  
    ω_a = 4; 
    ω_b = 4;
    # kdam =  0.#0.000147;#0.05; 
    k = 2; # carrying capacity - changes depending on the data?
    # lam = 0.033;

    # rtca_0 = 0#0.00894; 
    # rtcb_0 = 0#0.0216; 
    # rh_0 = 11.29; #69.56; #69.4
    # rtcr_0 = 0# 0.0131 #0.04; # 8.67e-3; # change this based on keeping steady state level the whole time course (levels shouldn't really change)
    # rm_a_0 = 0; 
    # rm_b_0 = 0; 
    # rm_r_0 = 0#0.0131#0.04 # 0; 
    # rd_0 = 0; 
    # rt_0 = 0;
end

tspan = (0,1e9)
rh = 11.29 #75 # conc of ribosomes in exponential phase 
thr_t = 5#30 # was at 5 before to get saved plots # needs to be less than 30 
kin_trna = 1

init_trna = [0,0,0,0,0,0,135.5,0,0] # tRNA initial conc = 135.5

function checking_bistability(model, params, init, specie, all_species, kdam_range)
    param_init = deepcopy(params)
    new_params = deepcopy(params)
    first_params = deepcopy(params)
    first_params[:kdam]=kdam_range[1]
    solu = sol(model, init, tspan, first_params)
    ss = get_ssval(solu, specie, all_species)
    init_first = ss_init_vals(solu, all_species)
    res =[]
    for i in range(2, length(kdam_range))
        param_init[:kdam]=kdam_range[i-1]
        solu_init = sol(model, init_first, tspan, param_init)
        init_ss = ss_init_vals(solu_init, all_species)
        new_params[:kdam] = kdam_range[i]
        solu_new = sol(model, init_ss, tspan, new_params)
        push!(res, get_ssval(solu_new, specie, all_species))
    end
    pushfirst!(res, ss)
    return res
end

trna_species = [:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :trna, :rd, :rt]

params_trna2 = (L = 10., c = 0.001, kr = 0.125, Vmax_init = 39.51, Km_init = 250.,
θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923,
kdeg = 0.001, kin = kin_trna, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
kdam =  0.01, lam = 0.014, rh = rh, thr_t = thr_t)
params_trna = @LArray [10., c, kr*12, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 0.5, ktag, kdeg, kin_trna, 3578.9473684210525, na, nb, nr, 0.014, rh, thr_t] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam, :rh, :thr_t)

br = get_br(rtc_mod_trna, params_trna2, init_trna, 400.)

df = create_br_df(br)
bf = bf_point_df(br)
kdam1 = findall(x->x==bf.kdam[1],df.kdam)[1]
unstable_firsthalf = df[!,:rh][kdam1:end][1:2094]
unstable1 = plot(scatter(x=df.kdam[kdam1:end],y=unstable_firsthalf, line=attr(color="#1f77b4", dash="dash"),showlegend=false))


kdam_range = range(0,400,length=1000)
kdam_range2 = range(400,0,length=1000)
kdam_range_short = range(20,350,length=10)
new_ps = deepcopy(params_trna)
first_init = DataFrame(rm_a=[],rtca=[],rm_b=[],rtcb=[],rm_r=[],rtcr=[],trna=[],rd=[],rt=[]); 
new_init = DataFrame(rm_a=[],rtca=[],rm_b=[],rtcb=[],rm_r=[],rtcr=[],trna=[],rd=[],rt=[]); 
final_ss = DataFrame(rm_a=[],rtca=[],rm_b=[],rtcb=[],rm_r=[],rtcr=[],trna=[],rd=[],rt=[]);
for i in kdam_range_short
    for (s, f, n, l) in zip(range(1,9),eachcol(first_init), eachcol(new_init), eachcol(final_ss))
        new_ps.kdam = i    
        solu = sol(rtc_model_trna, init_trna, tspan, new_ps)
        ssvals = ss_init_vals(solu, trna_species)
        add_vals=[]
        for spec in all_species
            push!(add_vals, 0.1*maximum(df[!,spec]))
        end
        ssvalsA = ssvals+add_vals
        solu2 = sol(rtc_model_trna, ssvalsA, tspan, new_ps)
        new_ssvals = ss_init_vals(solu2, trna_species)
        push!(f, ssvals[s])
        push!(n, ssvalsA[s])
        push!(l, new_ssvals[s])
    end
end

function plot_bf_species_switch(br, specie, specie2, kdam_range2, kdam_range_short, first_init, new_init, final_ss)
    df = create_br_df(br)
    bf = bf_point_df(br)
    kdam1 = findall(x->x==bf.kdam[1],df.kdam)[1]
    stable_trna = df[!,specie][1:kdam1]
    unstable_trna = df[!,specie][kdam1:end]
    stable1 = (scatter(x=df.kdam[1:kdam1],y=stable_trna, line=attr(color="#1f77b4"),showlegend=false))
    unstable1 = (scatter(x=df.kdam[kdam1:end],y=unstable_trna, line=attr(color="#1f77b4", dash="dash"),showlegend=false))
    res_trna2 = checking_bistability(rtc_model_trna, params_trna, init_trna, specie2, trna_species, kdam_range2)
    zeross = scatter(x=kdam_range2, y=res_trna2, name="↓ kdam", legendgroup=4, line=attr(color="#9467bd"),showlegend=false)
    diffs_ss = []
    for (i,j, kdam_val) in ProgressBar(zip(first_init[!,specie2], new_init[!,specie2], kdam_range_short))
        push!(diffs_ss,scatter(x=[kdam_val,kdam_val], y=[i,j],line=attr(color="#e377c2"),showlegend=false))
    end
    final_ss_p = scatter(x=kdam_range_short, y=final_ss[!,specie2],showlegend=false)
    plot([stable1, unstable1, zeross, diffs_ss[1], diffs_ss[2],
     diffs_ss[3], diffs_ss[4], diffs_ss[5], diffs_ss[6], diffs_ss[7], diffs_ss[8], 
    diffs_ss[9], diffs_ss[10], final_ss_p], Layout(title="$specie"))
end

prma = plot_bf_species_switch(br, :rm_a, :rm_a, kdam_range2, kdam_range_short, first_init, new_init, final_ss);
prtca = plot_bf_species_switch(br, :rtca, :rtca, kdam_range2, kdam_range_short, first_init, new_init, final_ss);    
prmb = plot_bf_species_switch(br, :rm_b, :rm_b, kdam_range2, kdam_range_short, first_init, new_init, final_ss);
prtcb = plot_bf_species_switch(br, :rtcb, :rtcb, kdam_range2, kdam_range_short, first_init, new_init, final_ss);     
prtcr = plot_bf_species_switch(br, :rtcr, :rtcr, kdam_range2, kdam_range_short, first_init, new_init, final_ss) ;
ptrna = plot_bf_species_switch(br, :rh, :trna, kdam_range2, kdam_range_short, first_init, new_init, final_ss)  ;
prd = plot_bf_species_switch(br, :rd, :rd, kdam_range2, kdam_range_short, first_init, new_init, final_ss)    ;
prt = plot_bf_species_switch(br, :rt, :rt, kdam_range2, kdam_range_short, first_init, new_init, final_ss) ;
    
p = [prma prtca prmb prtcb; prtcr ptrna prd prt]



# switching from on to off 
res_trna1 = numerical_bistability_analysis(rtc_model_trna, params_trna, init_trna, :rm_a, trna_species, kdam_range)
ptrna1 = plot(scatter(x=kdam_range, y=res_trna1, name="↑ kdam", legendgroup=3, line=attr(color="#e377c2")))
a=QuadraticInterpolation(Float64.(res_trna1),Float64.(collect(kdam_range)))
a(20)

maximum(df[!,:rh])

kdam_range = range(0,400,length=1000)
kdam_range2 = range(400,0,length=1000)
kdam_range_short = range(20,290,length=10)
new_ps = deepcopy(params_trna)
first_init = DataFrame(rm_a=[],rtca=[],rm_b=[],rtcb=[],rm_r=[],rtcr=[],trna=[],rd=[],rt=[]); 
new_init = DataFrame(rm_a=[],rtca=[],rm_b=[],rtcb=[],rm_r=[],rtcr=[],trna=[],rd=[],rt=[]); 
final_ss = DataFrame(rm_a=[],rtca=[],rm_b=[],rtcb=[],rm_r=[],rtcr=[],trna=[],rd=[],rt=[]);
for i in ProgressBar(kdam_range_short)
    ssvals=[]
    new_ps.kdam = i   
    for spec in trna_species
        res1 = numerical_bistability_analysis(rtc_model_trna, params_trna, init_trna, spec, trna_species, kdam_range)
        rQ = QuadraticInterpolation(Float64.(res1), Float64.(collect(kdam_range)))
        push!(ssvals, rQ(i))
    end
    add_vals=[]
    for (spec,i) in zip(all_species, range(1,9))
        push!(add_vals, 0.5*ssvals[i])
    end

    ssvalsA = ssvals-add_vals
    @show ssvals, ssvalsA
    solu2 = sol(rtc_model_trna, ssvalsA, tspan, new_ps)
    new_ssvals = ss_init_vals(solu2, trna_species)
    for (s, f, n, l) in zip(range(1,9),eachcol(first_init), eachcol(new_init), eachcol(final_ss))
        push!(f, ssvals[s])
        push!(n, ssvalsA[s])
        push!(l, new_ssvals[s])
    end
end

function plot_bf_species_switch(br, specie, specie2, kdam_range2, kdam_range_short, first_init, new_init, final_ss)
    df = create_br_df(br)
    bf = bf_point_df(br)
    kdam1 = findall(x->x==bf.kdam[1],df.kdam)[1]
    stable_trna = df[!,specie][1:kdam1]
    unstable_trna = df[!,specie][kdam1:end]
    stable1 = (scatter(x=df.kdam[1:kdam1],y=stable_trna, line=attr(color="#1f77b4"),showlegend=false))
    unstable1 = (scatter(x=df.kdam[kdam1:end],y=unstable_trna, line=attr(color="#1f77b4", dash="dash"),showlegend=false))
    res_trna2 = checking_bistability(rtc_model_trna, params_trna, init_trna, specie2, trna_species, kdam_range2)
    zeross = scatter(x=kdam_range2, y=res_trna2, name="↓ kdam", legendgroup=4, line=attr(color="#9467bd"),showlegend=false)
    diffs_ss = []
    for (i,j, kdam_val) in ProgressBar(zip(first_init[!,specie2], new_init[!,specie2], kdam_range_short))
        push!(diffs_ss,scatter(x=[kdam_val,kdam_val], y=[i,j],line=attr(color="#e377c2"),showlegend=false))
    end
    final_ss_p = scatter(x=kdam_range_short, y=final_ss[!,specie2],showlegend=false)
    plot([stable1, unstable1, zeross, diffs_ss[1], diffs_ss[2],
     diffs_ss[3], diffs_ss[4], diffs_ss[5], diffs_ss[6], diffs_ss[7], diffs_ss[8], 
    diffs_ss[9], diffs_ss[10], final_ss_p], Layout(title="$specie"))
end

prma = plot_bf_species_switch(br, :rm_a, :rm_a, kdam_range2, kdam_range_short, first_init, new_init, final_ss);
prtca = plot_bf_species_switch(br, :rtca, :rtca, kdam_range2, kdam_range_short, first_init, new_init, final_ss);    
prmb = plot_bf_species_switch(br, :rm_b, :rm_b, kdam_range2, kdam_range_short, first_init, new_init, final_ss);
prtcb = plot_bf_species_switch(br, :rtcb, :rtcb, kdam_range2, kdam_range_short, first_init, new_init, final_ss);     
prtcr = plot_bf_species_switch(br, :rtcr, :rtcr, kdam_range2, kdam_range_short, first_init, new_init, final_ss) ;
ptrna = plot_bf_species_switch(br, :rh, :trna, kdam_range2, kdam_range_short, first_init, new_init, final_ss)  ;
prd = plot_bf_species_switch(br, :rd, :rd, kdam_range2, kdam_range_short, first_init, new_init, final_ss)    ;
prt = plot_bf_species_switch(br, :rt, :rt, kdam_range2, kdam_range_short, first_init, new_init, final_ss) ;
    
p = [prma prtca prmb prtcb; prtcr ptrna prd prt]


