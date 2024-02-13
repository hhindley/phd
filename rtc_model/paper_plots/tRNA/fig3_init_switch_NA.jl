using Parameters, CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, Statistics
using Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra, Printf
# using Plots
using PlotlyJS, ProgressBars
include("/home/holliehindley/phd/may23_rtc/functions/solving.jl"); include("/home/holliehindley/phd/may23_rtc/functions/set_ups.jl"); include("/home/holliehindley/phd/may23_rtc/functions/plotting.jl"); 
include("/home/holliehindley/phd/may23_rtc/functions/sweep_params.jl"); include("/home/holliehindley/phd/may23_rtc/models/rtc_orig.jl"); include("/home/holliehindley/phd/may23_rtc/models/atp_lam_kin_t.jl"); 
include("/home/holliehindley/phd/may23_rtc/models/single_t.jl"); include("/home/holliehindley/phd/may23_rtc/models/combinations_t.jl"); 
include("/home/holliehindley/phd/may23_rtc/functions/bf_funcs/bf_funcs.jl");
include("/home/holliehindley/phd/colors_plotly.jl"); include("/home/holliehindley/phd/may23_rtc/models/inhibition_models/rtc_inhibition_model.jl");
include("/home/holliehindley/phd/may23_rtc/models/rtc_trna_model.jl")

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
    # kin = 0.054; #2.381 
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
kdeg_trna = 0.00001
trna_species = [:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :trna, :rd, :rt]
init_trna = [0,0,0,0,0,0,135.5,0,0] # tRNA initial conc = 135.5
params_trna = @LArray [10., c, kr*12, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 0.5, ktag, kdeg, kin_trna, 3578.9473684210525, na, nb, nr, 0.014, rh, thr_t] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam, :rh, :thr_t)

# params_trna = @LArray [10., c, kr*12, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 0.5, ktag, kdeg_trna, kin_trna, 3578.9473684210525, na, nb, nr, 0.014, rh, thr_t] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam, :rh, :thr_t)

kdam_range = range(0,400,length=1000)
kdam_range2 = range(400,0,length=1000)

params_trna2 = (L = 10., c = 0.001, kr = 0.125*12, Vmax_init = 39.51, Km_init = 250.,
θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923,
kdeg = 0.001, kin = kin_trna, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
kdam =  0.01, lam = 0.014, rh = rh, thr_t = thr_t)

# params_trna2 = (L = 10., c = 0.001, kr = 0.125*12, Vmax_init = 39.51, Km_init = 250.,
# θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = 0.2, 
# krep = 137., ktag = 9780., atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923,
# kdeg = kdeg_trna, kin = kin_trna, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
# kdam =  0.01, lam = 0.014, rh = rh, thr_t = thr_t)

br = get_br(rtc_mod_trna, params_trna2, init_trna, 400.)


bf0 = bf_point_df(br)
df0 = create_br_df(br)
kdam01 = findall(x->x==bf0.kdam[1],df0.kdam)[1]
kdam02 = findall(x->x==bf0.kdam[2],df0.kdam)[1]




params1 = @LArray [L, c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 1., ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)

kdam_range_onoff = range(0, df0.kdam[kdam01]-0.00008*df0.kdam[kdam01], length=100)
tspan=(0,1e9)

svals_onoff = DataFrame(rm_a=[],rtca=[],rm_b=[],rtcb=[],rm_r=[],rtcr=[],rh=[],rd=[],rt=[])

numerical_df = DataFrame(rm_a=[i for i in numerical_bistability_analysis(rtc_model_trna, params_trna, init_trna, :rm_a, all_species, kdam_range2)],
rtca=[i for i in numerical_bistability_analysis(rtc_model_trna, params_trna, init_trna, :rtca, all_species, kdam_range2)],
rm_b=[i for i in numerical_bistability_analysis(rtc_model_trna, params_trna, init_trna, :rm_b, all_species, kdam_range2)],
rtcb=[i for i in numerical_bistability_analysis(rtc_model_trna, params_trna, init_trna, :rtcb, all_species, kdam_range2)],
rm_r=[i for i in numerical_bistability_analysis(rtc_model_trna, params_trna, init_trna, :rm_r, all_species, kdam_range2)],
rtcr=[i for i in numerical_bistability_analysis(rtc_model_trna, params_trna, init_trna, :rtcr, all_species, kdam_range2)],
rh=[i for i in numerical_bistability_analysis(rtc_model_trna, params_trna, init_trna, :rh, all_species, kdam_range2)],
rd=[i for i in numerical_bistability_analysis(rtc_model_trna, params_trna, init_trna, :rd, all_species, kdam_range2)],
rt=[i for i in numerical_bistability_analysis(rtc_model_trna, params_trna, init_trna, :rt, all_species, kdam_range2)], 
kdam=[i for i in kdam_range2])

function setup_ssvals_trna(rtc_mod_trna, kdam_val, params_trna2, init_trna, numerical_df)
    br2 = get_br(rtc_mod_trna, params_trna2, init_trna, 400.)
    df = create_br_df(br2)
    df_bf = bf_point_df(br2)

    kdam1 = findall(x->x==df_bf.kdam[1],df.kdam)[1]
    # kdam2 = findall(x->x==df_bf.kdam[2],df.kdam)[1]
    first=DataFrame(kdam=[],rm_a=[],rtca=[],rm_b=[],rtcb=[],rm_r=[],rtcr=[],rh=[],rd=[],rt=[])
    last=DataFrame(kdam=[],rm_a=[],rtca=[],rm_b=[],rtcb=[],rm_r=[],rtcr=[],rh=[],rd=[],rt=[])

    for i in df.kdam[1:kdam1]
        push!(first.kdam, i)
    end
    
    for i in reverse(numerical_df.kdam)
        push!(last.kdam, i)
    end

    for (col,col1) in zip(eachcol(df)[1:9],eachcol(first)[2:end])
        for i in col[1:kdam1]
            push!(col1, i)
        end
    end

    for (col,col1) in zip(eachcol(numerical_df)[1:9],eachcol(last)[2:end])
        for i in reverse(col)
            push!(col1, i)
        end
    end

    first=Float64.(first)
    last=Float64.(last)
    
    int_rma1 = QuadraticInterpolation(first.rm_a, first.kdam)
    int_rma2 = QuadraticInterpolation(last.rm_a, last.kdam)
    int_rtca1 = QuadraticInterpolation(first.rtca, first.kdam)
    int_rtca2 = QuadraticInterpolation(last.rtca, last.kdam)
    int_rmb1 = QuadraticInterpolation(first.rm_b, first.kdam)
    int_rmb2 = QuadraticInterpolation(last.rm_b, last.kdam)
    int_rtcb1 = QuadraticInterpolation(first.rtcb, first.kdam)
    int_rtcb2 = QuadraticInterpolation(last.rtcb, last.kdam)
    int_rmr1 = QuadraticInterpolation(first.rm_r, first.kdam)
    int_rmr2 = QuadraticInterpolation(last.rm_r, last.kdam)
    int_rtcr1 = QuadraticInterpolation(first.rtcr, first.kdam)
    int_rtcr2 = QuadraticInterpolation(last.rtcr, last.kdam)
    int_rh1 = QuadraticInterpolation(first.rh, first.kdam)
    int_rh2 = QuadraticInterpolation(last.rh, last.kdam)
    int_rd1 = QuadraticInterpolation(first.rd, first.kdam)
    int_rd2 = QuadraticInterpolation(last.rd, last.kdam)
    int_rt1 = QuadraticInterpolation(first.rt, first.kdam)
    int_rt2 = QuadraticInterpolation(last.rt, last.kdam)

    return DataFrame(species=all_species,ss_val_on=[int_rma1(kdam_val),int_rtca1(kdam_val),int_rmb1(kdam_val),int_rtcb1(kdam_val),int_rmr1(kdam_val),int_rtcr1(kdam_val),int_rh1(kdam_val),int_rd1(kdam_val),int_rt1(kdam_val)],ss_val_off=[int_rma2(kdam_val),int_rtca2(kdam_val),int_rmb2(kdam_val),int_rtcb2(kdam_val),int_rmr2(kdam_val),int_rtcr2(kdam_val),int_rh2(kdam_val),int_rd2(kdam_val),int_rt2(kdam_val)]);
end

branches1 = setup_ssvals_trna(rtc_mod_trna, 200, params_trna2, init_trna, numerical_df)




res=[]
for kdam_val in ProgressBar(kdam_range_onoff)
    psm = deepcopy(params_trna)
    psm.kdam = kdam_val
    branches1 = setup_ssvals_trna(rtc_mod_trna, kdam_val, params_trna2, init_trna, numerical_df)
    # @show psm
    # @show branches1
    n = 600; l = 10;
    upper_ranges = get_all_ranges(set_ss_range_zerotossval, branches1, "ss_val_on", n, l)
    # @show upper_ranges[4]
    all, init_vals = get_rh_init_switch_all_ranges(rtc_model_trna, upper_ranges, branches1.ss_val_on,:rh,l,psm,9, all_species)
    binary = upper_or_lower(all, branches1.ss_val_off[7], l, 9)
    
    push!(res, binary)
    
    # inds = get_switch_ind(binary, l)
    # vals = get_switch_vals(inds, init_vals)
    # push!(svals_onoff.rm_a, vals[1])
    # push!(svals_onoff.rtca, vals[2])
    # push!(svals_onoff.rm_b, vals[3])
    # push!(svals_onoff.rtcb, vals[4])
    # push!(svals_onoff.rm_r, vals[5])
    # push!(svals_onoff.rtcr, vals[6])
    # push!(svals_onoff.rh, vals[7])
    # push!(svals_onoff.rd, vals[8])
    # push!(svals_onoff.rt, vals[9])
end

res

switch_ind=[]
for col in eachcol(res)
    @show col
    # print(length(findall(x->x==1,col)))

    if length(findall(x->x==1,col)) == 10
        push!(switch_ind,NaN)
    else
        # push!(switch_ind,findall(x->x==round(branches1.ss_val_on[7];digits=3),col)[1])
        push!(switch_ind,findall(x->x==1,col)[1])

    end
end
switch_ind


#THERE IS NO SWITCHING FROM ON TO OFF IN THIS MODEL WITH CHANGING ONE SPECIES AT A TIME 