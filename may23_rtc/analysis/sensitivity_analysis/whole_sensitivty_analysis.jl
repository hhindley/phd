using Parameters, CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, Statistics, DataInterpolations, Printf 
using Distributions
using StatsPlots
include("/home/holliehindley/phd/may23_rtc/functions/solving.jl"); include("/home/holliehindley/phd/may23_rtc/functions/set_ups.jl"); include("/home/holliehindley/phd/may23_rtc/functions/plotting.jl"); 
include("/home/holliehindley/phd/may23_rtc/functions/sweep_params.jl"); include("/home/holliehindley/phd/may23_rtc/models/rtc_orig.jl"); include("/home/holliehindley/phd/may23_rtc/models/atp_lam_kin_t.jl"); 
include("/home/holliehindley/phd/may23_rtc/models/single_t.jl"); include("/home/holliehindley/phd/may23_rtc/models/combinations_t.jl"); 
include("/home/holliehindley/phd/may23_rtc/functions/bf_funcs/bf_funcs.jl");
include("/home/holliehindley/phd/colors_plotly.jl");

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

# kr_range = range(kr-(kr*0.113), kr*1.113, length=100) # 11.3%
# Vmax_init_range = range(Vmax_init-4.62, Vmax_init+4.62, length=100) # +-4.62
# Km_init_range = range(Km_init-(Km_init*0.1), Km_init*1.1, length=100) # 10%
# d_range = range(d-(d*0.058), d*1.058, length=100) # 5.8%
# krep_range = range(krep-31, krep+31, length=100) # +-31
# ktag_range = range(ktag-(ktag*0.086), ktag*1.086, length=100) # 8.6%

# kdeg_range = range(kdeg-(kdeg*0.1), kdeg*1.1, length=100)
# kin_range = range(kin-(kin*0.1), kin*1.1, length=100)

kr_sd = kr*0.113 # estimated
Vmax_init_sd = 4.62
Km_init_sd = Km_init*0.1
d_sd = d*0.058
krep_sd = 31
ktag_sd = ktag*0.086
# kdeg_sd = kdeg*0.1 # estimated
# kin_sd = kin*0.5 # estimated 
kma_sd = km_a*0.1
kmb_sd = km_b*0.1

kr_norm = Normal(kr, kr_sd)
Vmax_init_norm = Normal(Vmax_init, Vmax_init_sd)
Km_init_norm = Normal(Km_init, Km_init_sd)
d_norm = Normal(d, d_sd)
krep_norm = Normal(krep, krep_sd)
ktag_norm = Normal(ktag, ktag_sd)
kma_norm = Normal(km_a, kma_sd)
kmb_norm = Normal(km_b, kmb_sd)

kin_lb = 0.01
kin_ub = 2
range_kin = kin_ub - kin_lb
kin_std = 0.4
kin_dist = Truncated(Normal(kin, kin_std), kin_lb, kin_ub)

kin_new = rand(kin_dist,l)
histogram(kin_new, bins=50, label="kin")

kdeg_lb = 0.0005
kdeg_ub = 0.05
range_kdeg = (kdeg_ub) - (kdeg_lb)
kdeg_std = 0.015
kdeg_dist = Truncated(Normal(kdeg, kdeg_std), kdeg_lb, kdeg_ub)

krep_dist = Truncated(krep_norm, 0, 300)

kdeg_lnorm = LogNormal(kdeg)
histogram(rand(kdeg_lnorm, l), bins=100)

l=1000000
kr_new = rand(kr_norm,l)
Vmax_init_new = rand(Vmax_init_norm,l)
Km_init_new = rand(Km_init_norm,l)
d_new = rand(d_norm,l)
krep_new = rand(krep_dist,l)
ktag_new = rand(ktag_norm,l)
kdeg_new = rand(kdeg_dist,l)
kin_new = rand(kin_dist,l)
kma_new = rand(kma_norm,l)
kmb_new = rand(kmb_norm,l)



p_kr = histogram(kr_new, bins=100, label="kr", ylabel="f(x)", xlabel="x")#, formatter=y->y/1e4)#, x->x*1e6)
p_Vmax = histogram(Vmax_init_new, bins=100, label="Vmax_init", ylabel="f(x)", xlabel="x")#, formatter=y->y/1e4)
p_Km = histogram(Km_init_new, bins=100, label="Km_init", ylabel="f(x)", xlabel="x")#, formatter=y->y/1e4)
p_d = histogram(d_new, bins=100, label="d", ylabel="f(x)", xlabel="x")#, formatter=y->y/1e4)#, x->x*1e5)
p_krep = histogram(krep_new, bins=100, label="krep", ylabel="f(x)", xlabel="x")#, formatter=y->y/1e4)
p_ktag = histogram(ktag_new, bins=100, label="ktag", ylabel="f(x)", xlabel="x")#, formatter=y->y/1e4)
p_kdeg = histogram(kdeg_new, bins=100, label="kdeg", ylabel="f(x)", xlabel="x")#, formatter=y->y/1e4)#, x->x*1e6)
p_kin = histogram(kin_new, bins=100, label="kin", ylabel="f(x)", xlabel="x")#, formatter=y->y/1e4)
p_kma = histogram(kma_new, bins=100, label="km_a", ylabel="f(x)", xlabel="x")#, formatter=y->y/1e4)
p_kmb = histogram(kmb_new, bins=100, label="km_b", ylabel="f(x)", xlabel="x")#, formatter=y->y/1e4)#, x->x*1e3)

p=plot(p_kr, p_Vmax, p_Km, p_d, p_krep, p_ktag, p_kdeg, p_kma, p_kmb, size=(1000,1000), yaxis=(formatter=y->y/1e4))

savefig(p, "/home/holliehindley/phd/may23_rtc/analysis/sensitivity_analysis/dists.svg")

plot(x->pdf(kr_new, x))

pdf(kr_new, l)
# plot(plot(krep_norm, label="krep"), plot(kr_norm, label="kr"), plot(Vmax_init_norm, label="Vmax_init"), plot(Km_init_norm, label="Km_init"),
# plot(d_norm, label="d"), plot(ktag_norm, label="ktag"), plot(kdeg_norm, label="kdeg"), plot(kin_dist, label="kin"), 
# plot(kma_norm, label="km_a"), plot(kmb_norm, label="km_b"), size=(1000,400))

# kin_sd = kin*0.5 # estimated 
# kin_norm = Normal(kin, kin_sd)
# kin_new = rand(kin_norm,l)


res=[]; new_ps = []; p = Plots.plot();
for i in range(1, length(kr_new))
    params_new = (L = 10., c = 0.001, kr = kr_new[i], Vmax_init = Vmax_init_new[i], Km_init = Km_init_new[i],
    θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = d_new[i], 
    krep = krep_new[i], ktag = ktag_new[i], atp = 3578.9473684210525, km_a = kma_new[i], km_b = kmb_new[i], g_max = 2.0923, 
    kdeg = kdeg_new[i], kin = kin_new[i], ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
    kdam =  0., lam = 0.014)
    # params_new = (L = 10., c = 0.001, kr = kr_new[i], Vmax_init = Vmax_init_new[i], Km_init = Km_init_new[i],
    # θtscr = 160.01, θtlr = 255.73, na = 338., nb = 408., nr = 532. *6, d = d_new[i], 
    # krep = krep_new[i], ktag = ktag_new[i], atp = 3578.9473684210525, km_a = 20., km_b = 16., g_max = 2.0923, 
    # kdeg = kdeg_new[i], kin = kin, ω_ab = 0.05623413251903491, ω_r = 0.010000000000000002, 
    # kdam =  0., lam = 0.014)
    br = get_br(rtc_mod, params_new, initial, 3.)
    df = create_br_df(br)
    if length(br.specialpoint) == 4
        push!(res, @LArray [kr_new[i],Vmax_init_new[i],Km_init_new[i],d_new[i],krep_new[i],ktag_new[i],kma_new[i],kmb_new[i],kdeg_new[i],kin_new[i]] (:kr,:Vmax_init,:Km_init,:d,:krep,:ktag,:km_a,:km_b,:kdeg,:kin))
        p = Plots.plot!(df.kdam, df.rh, size=(900,600), legend=false)
        push!(new_ps, (kr=kr_new[i],Vmax_init = Vmax_init_new[i], Km_init = Km_init_new[i],d = d_new[i], 
        krep = krep_new[i], ktag = ktag_new[i], kdeg = kdeg_new[i], kin = kin_new[i], km_a = kma_new[i], km_b = kmb_new[i]))
    end
end
p










krp = Plots.plot(kr_norm, label="kr", size=(900,600));
scatter!(krp, kr_new, zeros(length(kr_new)));

Vmaxp = Plots.plot(Vmax_init_norm, label="Vmax_init", size=(900,600));
scatter!(Vmaxp, Vmax_init_new, zeros(length(Vmax_init_new)));

Kmp = Plots.plot(Km_init_norm, label="Km_init", size=(900,600));
scatter!(Kmp, Km_init_new, zeros(length(Km_init_new)));

krepp = Plots.plot(krep_norm, label="krep", size=(900,600));
scatter!(krepp, krep_new, zeros(length(krep_new)));

dp = Plots.plot(d_norm, label="d", size=(900,600));
scatter!(dp, d_new, zeros(length(d_new)));

ktagp = Plots.plot(ktag_norm, label="ktag", size=(900,600));
scatter!(ktagp, ktag_new, zeros(length(ktag_new)));

kdegp = Plots.plot(kdeg_norm, label="kdeg", size=(900,600));
scatter!(kdegp, kdeg_new, zeros(length(kdeg_new)));

kinp = Plots.plot(kin_norm, label="kin", size=(900,600));
scatter!(kinp, kin_new, zeros(length(kin_new)));

Plots.plot(krp, Vmaxp, Kmp, krepp, dp, ktagp, kdegp, kinp)

Plots.plot(kinp)




# Define your upper and lower bounds
lower_bound = 0.01
upper_bound = 2

# Define your mean value
mean_value = 0.022

range1 = log(upper_bound) - log(lower_bound)

# Estimate the standard deviation in the log domain
log_std_dev = range1 / (2 * spread_factor)

# Convert the standard deviation back to the original domain
std_dev = exp(log_std_dev)

# Create a truncated log-normal distribution
truncated_lognorm_dist = Truncated(LogNormal(log(mean_value, log_std_dev), std_dev), lower_bound, upper_bound)

# Generate random samples from the truncated log-normal distribution
samples = rand(truncated_lognorm_dist, 10000)

# Plot the histogram of the samples
histogram(samples, bins=50, label="Sampled Distribution")
plot!(title="Truncated Log-Normal Distribution", xlabel="Value", ylabel="Probability Density", legend=true)

plot(truncated_lognorm_dist)