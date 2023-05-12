using Parameters, LabelledArrays, StaticArrays, CSV, DataFrames, DifferentialEquations, BenchmarkTools, OrderedCollections, DataInterpolations, PlotlyJS, Statistics
include("/home/holliehindley/phd/may23_rtc/functions/solving.jl"); include("/home/holliehindley/phd/may23_rtc/functions/plotting.jl"); include("/home/holliehindley/phd/may23_rtc/functions/sweep_params.jl"); include("/home/holliehindley/phd/may23_rtc/models/rtc_orig.jl"); include("/home/holliehindley/phd/may23_rtc/models/atp_lam_kin_t.jl"); include("/home/holliehindley/phd/may23_rtc/analysis/t_param_setup.jl");

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
    atp = 4000;#2500; 
    km_a = 20; 
    km_b = 16;
    g_max = 2.0923; 
    gr_c = 0.0008856; # 0.000599; 
    kdeg = 0.001; 
    kin = 0.054; #2.381 
    ω_ab = 4#4#0.093; #0.0828304057748932;#4; 
    ω_r = 2e-7 #0.0019*6 #70.53; #0.0019*6#79.43865871861044; #0.0019*6;  
    ω_a = 4; 
    ω_b = 4;
    kdam = 0;#0.000147;#0.05; 
    k = 2; # carrying capacity - changes depending on the data?
    lam = 0.033;

    rtca_0 = 0#0.00894; 
    rtcb_0 = 0#0.0216; 
    rh_0 = 11.29; #69.56; #69.4
    rtcr_0 = 0# 0.0131 #0.04; # 8.67e-3; # change this based on keeping steady state level the whole time course (levels shouldn't really change)
    rm_a_0 = 0; 
    rm_b_0 = 0; 
    rm_r_0 = 0#0.0131#0.04 # 0; 
    rd_0 = 0; 
    rt_0 = 0;

    tspan = (0, 1e9);
end
params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
initial = @SVector [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, rh_0, rd_0, rt_0]

kin_range = collect(0:0.1:10)
kdam_range = collect(0:0.01:1)


params[:kin] = 3
kdam_list = collect(0:0.1:1)#[0.6,0.68,0.75,0.8,1]
param_comparison(rtc_model, initial, tspan, params, :rh, :kdam, kdam_list, "kin=$(params[:kin])", "log")
param_comparison(rtc_model, initial, tspan, params, :rt, :kdam, kdam_list, "kin=$(params[:kin])", "log")
param_comparison(rtc_model, initial, tspan, params, :rd, :kdam, kdam_list, "kin=$(params[:kin])", "log")
param_comparison(rtc_model, initial, tspan, params, :rm_a, :kdam, kdam_list, "kin=$(params[:kin])", "log")
param_comparison(rtc_model, initial, tspan, params, :rtca, :kdam, kdam_list, "kin=$(params[:kin])", "log")
param_comparison(rtc_model, initial, tspan, params, :rtcr, :kdam, kdam_list, "kin=$(params[:kin])", "log")


params[:kdam] = 0.71
kin_list = collect(2.5:0.2:4.5)
param_comparison(rtc_model, initial, tspan, params, :rh, :kin, kin_list, "kdam=$(params[:kdam])", "log")
param_comparison(rtc_model, initial, tspan, params, :rt, :kin, kin_list, "kdam=$(params[:kdam])", "log")
param_comparison(rtc_model, initial, tspan, params, :rd, :kin, kin_list, "kdam=$(params[:kdam])", "log")
param_comparison(rtc_model, initial, tspan, params, :rm_a, :kin, kin_list, "kdam=$(params[:kdam])", "log")
param_comparison(rtc_model, initial, tspan, params, :rtca, :kin, kin_list, "kdam=$(params[:kdam])", "log")
param_comparison(rtc_model, initial, tspan, params, :rtcr, :kin, kin_list, "kdam=$(params[:kdam])", "log")

# changing both parameters to match the border of rh plot 
param_list = [[0.95,0.73,0.56,0.1],[5,3.6,2.7,1]]
param = [:kdam, :kin]
param_comparison(rtc_model, initial, tspan, params, :rh, param, param_list, "", "log")
param_comparison(rtc_model, initial, tspan, params, :rt, param, param_list, "", "log")
param_comparison(rtc_model, initial, tspan, params, :rd, param, param_list, "", "log")
param_comparison(rtc_model, initial, tspan, params, :rm_a, param, param_list, "", "log")
param_comparison(rtc_model, initial, tspan, params, :rtca, param, param_list, "", "log")
param_comparison(rtc_model, initial, tspan, params, :rtcr, param, param_list, "", "log")


kdam_range = collect(0:0.001:0.04)
params[:kin] = 0.054#4
params
results = change_param(kdam_range, :kdam, rtc_model, initial, all_species, lam, atp, params[:kin])
p = (plot_change_param_sols(kdam_range, results, "kdam", ""))








function sweep_paramx2_new(model, lam, atp, kin, species, func, param1, param2, param_range1, param_range2)
    all_res = []
    params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
    for i in param_range1
        params[param1] = i
        res1 = []
        for val in param_range2 
            params[param2] = val
            solu = sol(model, initial, tspan, params)
            push!(res1, func(solu, species))
        end
        push!(all_res, res1)

    end
    # @show size(all_res)
    vec = []
    for i in (1:length(param_range1))
        append!(vec, values(all_res[i]))
    end
    # @show (vec)
    vec = reshape(vec, (length(param_range1),length(param_range1)))
    @show (vec)
    return (contour(x=param_range1, y=param_range2, z=vec, colorbar=attr(title="$species", titleside="right"), opacity=0.5))
end

p_rh = (sweep_paramx2_new(rtc_model, lam, atp, kin, :rh, get_ssval, :kin, :kdam, kin_range, kdam_range))
p_rd = (sweep_paramx2_new(rtc_model, lam, atp, kin, :rd, get_ssval, :kin, :kdam, kin_range, kdam_range))
p_rt = (sweep_paramx2_new(rtc_model, lam, atp, kin, :rt, get_ssval, :kin, :kdam, kin_range, kdam_range))

fig = plot([p_rh, p_rd, p_rt])#, Layout(grid=false))#xaxis_visible=false, yaxis_visible=false, showticklabels=true))
relayout!(fig, template=:plotly_white)