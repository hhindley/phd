using Parameters, LabelledArrays, StaticArrays, CSV, DataFrames, DifferentialEquations, BenchmarkTools, OrderedCollections, DataInterpolations, PlotlyJS, Statistics
include("$PATHmay23_rtc/functions/solving.jl"); include("$PATHmay23_rtc/functions/plotting.jl"); include("$PATHmay23_rtc/functions/sweep_params.jl"); include("$PATHmay23_rtc/models/rtc_orig.jl"); include("$PATHmay23_rtc/models/atp_lam_kin_t.jl"); include("$PATHmay23_rtc/analysis/t_param_setup.jl");

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
    ω_ab = 2#0.093; #0.0828304057748932;#4; 
    ω_r = 2e-7 #0.0019*6 #70.53; #0.0019*6#79.43865871861044; #0.0019*6;  
    ω_a = 4; 
    ω_b = 4;
    kdam = 0.01;#0.000147;#0.05; 
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

L_range = 10 .^(range(1,stop=3,length=101))
c_range = 10 .^(range(-4,stop=0,length=101))
wab_range = range(0,stop=2,length=101)
wr_range = 10 .^(range(-7,stop=-2,length=101))
atp_range = range(500,stop=5000,length=101)
kin_range = range(0,stop=0.2,length=101)
lam_range = range(0.001,stop=0.04,length=101)
kdam_range = 10 .^ range(-4,stop=0,length=101)


results = change_param(kdam_range, :kdam, rtc_model, initial, all_species, params)
plot(scatter(x=kdam_range, y=results[:rm_a]), Layout(xaxis_title="kdam", yaxis_title="rm_a"))


# peaks
sweep_paramx2_new(rtc_model, :rm_a, get_ssval, :atp, :kdam, atp_range, kdam_range, "", "log")
sweep_paramx2_new(rtc_model, :rm_a, get_ssval, :ω_ab, :kdam, wab_range, kdam_range, "", "log")
sweep_paramx2_new(rtc_model, :rm_a, get_ssval, :kin, :kdam, kin_range, kdam_range, "", "log")

range1 = range(500, stop=5000, length=5)
param2x_plot_same_species(range1, :atp, kdam_range, :kdam, :rm_a)
res_atp = param2x_same_species(range1, :atp, kdam_range, :kdam, :rm_a)
max_atp = [(argmax(i)) for i in res_atp]


wab_range1 = range(0.001, stop=2, length=10)
param2x_plot_same_species(wab_range1, :ω_ab, kdam_range, :kdam, :rm_a)
res_wab = param2x_same_species(wab_range1, :ω_ab, kdam_range, :kdam, :rm_a)
max_wab = [(argmax(i)) for i in res_wab]


kin_range1 = range(0.0001, stop=0.2, length=10)
param2x_plot_same_species(kin_range1, :kin, kdam_range, :kdam, :rm_a)
res_kin = param2x_same_species(kin_range1, :kin, kdam_range, :kdam, :rm_a)
max_kin = [(argmax(i)) for i in res_kin]


kdam_range1 = 10 .^range(-3, stop=0, length=10)
param2x_plot_same_species(kdam_range1, :kdam, lam_range, :lam, :rm_a)



solu = sol(rtc_model, initial, tspan, params)
plotly_plot_sol(solu, "log", "", "")

params1 = deepcopy(params)

ss_init = ss_init_vals(solu)
ss_init
tspan1 = (0,1000)
solu1 = sol(rtc_model, ss_init, tspan1, params)
plotly_plot_sol(solu1, "", "", "")

# params1[:c] = 0.01
params1[:kdam] = 0.0
solu2 = sol(rtc_model, ss_init, tspan, params1)
plotly_plot_sol(solu2, "log", "", "")





L_range = 10 .^(range(1,stop=3,length=10))
c_range = 10 .^(range(-4,stop=0,length=10))
wab_range = range(0,stop=2,length=10)
wr_range = 10 .^(range(-7,stop=-2,length=10))
atp_range = range(500,stop=5000,length=10)
kin_range = range(0,stop=0.2,length=10)
lam_range = range(0.001,stop=0.04,length=10)
kdam_range = 10 .^ range(-4,stop=0,length=10)

param_change = [:L, :c, :ω_ab, :ω_r, :atp, :kin, :lam, :kdam]
param_ranges = [L_range, c_range, wab_range, wr_range, atp_range, kin_range, lam_range, kdam_range]

for (i, j) in zip(param_change, param_ranges)
    cuv = []
    params1 = deepcopy(params)
    for param in j
        params1[i] = param
        solu = sol(rtc_model, ss_init, (180,1440), params1)
        push!(cuv, scatter(x=solu.t/60, y=get_curve(solu, :rm_a), name="$i = $param"))
    end
    display(plot([i for i in cuv], Layout(title="rm_a", yaxis_tickformat=".2e")))
end

