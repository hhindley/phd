using Plots, DifferentialEquations, DataFrames, Statistics, Measures, ProgressBars
include("fullmodel_functions.jl")
include("fullmodel_params.jl")
include("models_new.jl")# plotlyjs()

# tolerance vs. susceptibility
y = [1, 0.1, 0.01, 0.001]
x=-log.(y)

x1=x./2

y1 = [0.01,0.01,0.01,0.01]
x2 = [0,2,4, 7]

y3 = [0.001,10^-2]
x3 = [2.302585,2.302585]

x4 = [4.60517,4.60517]

plot(x,y, yaxis=:log10, yticks=y, legend=false, xlabel="Time (hours)", grid=false,
 ylabel="Fraction of survivors", linecolor=:red, linewidth=3,
    xticks=[0,2,4,6], xlims=(0,6), ylims=(1e-3,1), annotationfontsize=8, annotationcolor=:red)
plot!(x1,y, linecolor=:purple3, linewidth=3)#annotations=([2.4],[0.03],"Susceptible"), annotationfontsize=8, annotationcolor=:purple3)
plot!(x2,y1, linecolor=:grey, linewidth=1)
plot!(x3,y3, linecolor=:grey, linestyle=:dash, linewidth=1)
plot!(x4,y3, linecolor=:grey, linestyle=:dash, linewidth=1)

# dependence on ATP
atp_range_long = range(0,10000, length=10000)
include("fullmodel_params.jl")
rm_a_sols, rtca_sols, rm_b_sols, rtcb_sols, rm_r_sols, rtcr_sols, rdrtca_sols, rtrtcb_sols, ribo_h_sols, ribo_d_sols, ribo_t_sols = @time param_vs_ss(atp_range_long, params_new, 5, RTC_SYSTEM2!, init_rtc, params_new, species_rtc)

function plot_inset()
    p = plot(atp_range_long, [rtca_sols], xlabel="[ATP]", 
    ylabel="steady-state exp level (a.u.)", legend=false, linecolor=:green, 
    annotations=(5.3e3,[2.2e6], "RtcA/B"), annotationcolor=:green, annotationfontsize=14,
    xaxis=(:log, (1, Inf)), yaxis=(:log, (1, Inf)), 
    linewidth=3, grid=false, margin=10mm, yticks=[1,100,10000,1000000]) # previous annotations for a log-log plot annotations=(5.3e3,[2.2e6]
    
    plot!(atp_range_long, rtcr_sols, linecolor=:darkorange2, annotations=(6e3,
    [rtcr_sols[end]+10000], "RtcR"), annotationcolor=:darkorange2, annotationfontsize=14, linewidth=3)
    
    # ylims!(1e3,1.45e6) # previous annotations for a log-log plot annotations=(6e3, [rtcr_sols[end]+10000]
    # yticks!([1e4, 1e6])
    # lens!([1,10],[1,5000], inset = (1, bbox(0.2, 0.1, 0.3, 0.3)), lc=:false)
    # plot!(p[2], atp_range, [rtca_sols], xlims=(0,5), ylims=(0,1000), legend=false, linecolor=:green, xlabel="ATP", ylabel="amount of species", xtickfontsize=7, ytickfontsize=7, xlabelfontsize=7, ylabelfontsize=7, annotations=([atp_range[15]],[rtca_sols[end]-102800], "RtcA/B"), annotationcolor=:green, annotationfontsize=7)
    # plot!(p[2], atp_range, [rtcr_sols], xlims=(0,5), ylims=(0,1000), legend=false, linecolor=:darkorange2, annotations=([atp_range[25]],[rtcr_sols[end]-3400], "RtcR"), annotationcolor=:darkorange2, annotationfontsize=7)
end

p1 = plot_inset()

x = [1,10,100,1000,10000]
y = 10*x 
plot!(x,y)

# model solutions 
sol = simple_solve!(RTC_SYSTEM2!, init_rtc, tspan, params_new)
rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rdrtca, rtrtcb, ribo_h, ribo_d, ribo_t = split_reduced_model(species_rtc, sol)
plot(sol.t, [rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rdrtca, rtrtcb, ribo_d, ribo_h, ribo_t], xaxis=(:log10, (1,Inf)), yaxis=(:log10, (1,Inf)), labels=labels_all)

function plot_sol()
    plot(sol.t, rm_a, linecolor=:green3, annotations=([3.8e2],[rm_a[end]+72], "mRNA-RtcA/B"), annotationcolor=:green3, annotationfontsize=14, xaxis=(:log10, (1,1e3)), 
        margin=5mm, grid=false, yaxis=(:log10, (1,Inf)), linewidth=3, legend=false, xlabel="time (min)", ylabel="exp level (a.u.)")
    plot!(sol.t, rtca, linecolor=:green, linewidth=3, annotations=([6.3e2],[rtca[end]+10000], "RtcA/B"), annotationcolor=:green, annotationfontsize=14)
    plot!(sol.t, rm_r, linecolor=:orange, linewidth=3, annotations=([4.5e2],[rm_r[end]+8], "mRNA-RtcR"), annotationcolor=:orange, annotationfontsize=14)
    plot!(sol.t, rtcr, linecolor=:darkorange2, linewidth=3, annotations=([7.2e2],[rtcr[end]+1000], "RtcR"), annotationcolor=:darkorange2, annotationfontsize=14)
    plot!(sol.t, ribo_h, linecolor=:dodgerblue2, linewidth=3, annotations=([8.3e2],[ribo_h[end]+3], "Rh"), annotationcolor=:dodgerblue2, annotationfontsize=14)
end

p2 = plot_sol()

r_tot = ribo_h + ribo_d + ribo_t
perc_rh = @. ribo_h/r_tot * 100
perc_rd = @. ribo_d/r_tot * 100
perc_rt = @. ribo_t/r_tot * 100 

plot(sol.t, perc_rh, linecolor=:dodgerblue2, linewidth=3, annotations=([8e2],104, "Rh"), annotationcolor=:dodgerblue2, annotationfontsize=14, 
    xaxis=(:log10, (1,1e3)), grid=false, legend=false, xlabel="time (min)", ylabel="% total ribosomes", margin=5mm)
plot!(sol.t, perc_rd, linecolor=:blue, linewidth=3, annotations=([8e2],10, "Rd"), annotationcolor=:blue, annotationfontsize=14)
plot!(sol.t, perc_rt, linecolor=:deepskyblue, linewidth=3, annotations=([8e2],4, "Rt"), annotationcolor=:deepskyblue, annotationfontsize=14)

plot(sol.t, ribo_h, linecolor=:dodgerblue2, linewidth=3, annotations=([8e2],[ribo_h[end]+0.3], "Rh"), annotationcolor=:dodgerblue2, annotationfontsize=10, 
    xaxis=(:log10, (1,1e3)), grid=false, legend=false, xlabel="time (minutes)", ylabel="% total ribosomes")
plot!(sol.t, ribo_d, linecolor=:blue, linewidth=3, annotations=([8e2],[ribo_d[end]+1], "Rd"), annotationcolor=:blue, annotationfontsize=14)
plot!(sol.t, ribo_t, linecolor=:deepskyblue, linewidth=3, annotations=([8e2],[ribo_t[end]+0.3], "Rt"), annotationcolor=:deepskyblue, annotationfontsize=10)





