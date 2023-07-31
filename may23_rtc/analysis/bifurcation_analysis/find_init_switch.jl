using Parameters, CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, Statistics
# using Plots
using PlotlyJS
include("/home/holliehindley/phd/may23_rtc/functions/solving.jl"); include("/home/holliehindley/phd/may23_rtc/functions/set_ups.jl"); include("/home/holliehindley/phd/may23_rtc/functions/plotting.jl"); 
include("/home/holliehindley/phd/may23_rtc/functions/sweep_params.jl"); include("/home/holliehindley/phd/may23_rtc/models/rtc_orig.jl"); include("/home/holliehindley/phd/may23_rtc/models/atp_lam_kin_t.jl"); 
include("/home/holliehindley/phd/may23_rtc/models/single_t.jl"); include("/home/holliehindley/phd/may23_rtc/models/combinations_t.jl");



t, atp_t, lam_t, kin_t = set_time_vars("/home/holliehindley/phd/data/atp_for_rtcmodel.csv")

p = make_subplots(rows=3, cols=1, shared_xaxes=true, vertical_spacing=0.08, subplot_titles=["λ" "ATP" "kin"])
add_trace!(p, (scatter(x=t, y=lam_t)), row=1, col=1)
add_trace!(p, (scatter(x=t, y=atp_t)), row=2, col=1)
add_trace!(p, (scatter(x=t, y=kin_t)), row=3, col=1)
relayout!(p, showlegend=false, xaxis_range=(0,1500))
p

tspan = (0, lam_t[end])
# lam = 0.033; atp = 4000; kin = 0.054;

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
initial = @SVector [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, rh_0, rd_0, rt_0]




# res = change_param_timevars(kdam_range, :kdam, rtc_all_t!, initial, all_species, atp_t, lam_t, kin_t)


tspan = (0,1e6)
params = @LArray [L, c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 1., ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)


initial = @SVector [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, rh_0, rd_0, rt_0]
initial = [0, 0, 0.04, 0, 0, 0, 11.29, 0, 0];
solu = sol(rtc_model, initial, tspan, params);
plotly_plot_sol(solu, "", "", "")

initial1 = [0, 0, 0, 0, 0, 0, 11.29, 0, 0];
solu1 = sol(rtc_model, initial1, tspan, params);
plotly_plot_sol(solu1, "", "", "")


# work out which species affect switch 
ranges = [range(0,1000,length=50),range(0,1000,length=50),range(0,0.1,length=50),range(0,0.02,length=50),range(0,1,length=50),range(0,0.05,length=50),range(0,50,length=50),range(0,200,length=50),range(0,150,length=50)]

function get_init_switch(specie, ranges)
    # ranges = [range(0,1000,length=50),range(0,1000,length=50),range(0,0.1,length=50),range(0,0.1,length=50),range(0,1,length=50),range(0,1,length=50),range(0,50,length=50),range(0,200,length=50),range(0,150,length=50)]
    res_df = DataFrame(rm_a_0=[], rtca_0=[], rm_b_0=[], rtcb_0=[], rm_r_0=[], rtcr_0=[], rh_0=[], rd_0=[], rt_0=[])
    for (i,col,range) in zip(range(1,length(initial)), eachcol(res_df), ranges)
        res=[]
        # print(col)
        initial = [0, 0, 0, 0, 0, 0, 11.29, 0, 0]
        for j in range
            initial[i] = j
            # @show initial
            solu = sol(rtc_model, initial, tspan, params)
            rh_ss = get_ssval(solu,specie)
            push!(res, rh_ss)
        end
        push!(col, ([res...]...))
    end
    return res_df
end

all_df = []
for i in all_species
    push!(all_df, get_init_switch(i,ranges))
end


df_range = DataFrame(rma_range=[], rtca_range=[], rmb_range=[], rtcb_range=[], rmr_range=[], rtcr_range=[], rh_range=[], rd_range=[], rt_range=[])
for (range,col) in zip(ranges,eachcol(df_range))
    # print(([r...]...))
    push!(col, ([range...]...))
    # println(col)

end


rm_a_switch_res = all_df[1]; rtca_switch_res = all_df[2]; rm_b_switch_res = all_df[3]; rtcb_switch_res = all_df[4]; rm_r_switch_res = all_df[5]; rtcr_switch_res = all_df[6]; rh_switch_res = all_df[7]; rd_switch_res = all_df[8] ; rt_switch_res = all_df[9]

function plot_species_switch(xaxis,specie)
    t1 = scatter(x=xaxis,y=rm_a_switch_res[:,specie], name="rm_a")
    t2 = scatter(x=xaxis,y=rtca_switch_res[:,specie], name="rtca")
    t3 = scatter(x=xaxis,y=rm_b_switch_res[:,specie], name="rm_b")
    t4 = scatter(x=xaxis,y=rtcb_switch_res[:,specie], name="rtcb")
    t5 = scatter(x=xaxis,y=rm_r_switch_res[:,specie], name="rm_r")
    t6 = scatter(x=xaxis,y=rtcr_switch_res[:,specie], name="rtcr")
    t7 = scatter(x=xaxis,y=rh_switch_res[:,specie], name="rh")
    t8 = scatter(x=xaxis,y=rd_switch_res[:,specie], name="rd")
    t9 = scatter(x=xaxis,y=rt_switch_res[:,specie], name="rt")
    return plot([t1,t2,t3,t4,t5,t6,t7,t8,t9], Layout(xaxis_title="$specie", yaxis_title="steady-state conc"))#, title="changing initial value of $specie and getting steady state value to see which bs branch it gives"))
end

rma_plot = plot_species_switch(df_range.rma_range,:rm_a_0)
rtca_plot = plot_species_switch(df_range.rtca_range,:rtca_0)
rmb_plot = plot_species_switch(df_range.rmb_range,:rm_b_0)
rtcb_plot = plot_species_switch(df_range.rtcb_range,:rtcb_0)
rmr_plot = plot_species_switch(df_range.rmr_range,:rm_r_0)
rtcr_plot = plot_species_switch(df_range.rtcr_range,:rtcr_0)
rh_plot = plot_species_switch(df_range.rh_range,:rh_0)
rd_plot = plot_species_switch(df_range.rd_range,:rd_0)
rt_plot = plot_species_switch(df_range.rt_range,:rt_0)

p = [rma_plot rtca_plot rmb_plot; rtcb_plot rmr_plot rtcr_plot; rh_plot rd_plot rt_plot]
relayout!(p,legend=false, title_text="changing intitial values to see where and which variables cause the switch")
p
function plot_res(res)
    #plotly plotting
    trace1 = scatter(x=kdam_range, y=res[:rtca], name="RtcA")
    trace2 = scatter(x=kdam_range, y=res[:rtcb], name="RtcB")
    trace3 = scatter(x=kdam_range, y=res[:rm_a], name="mRNA RtcA")
    trace4 = scatter(x=kdam_range, y=res[:rm_b], name="mRNA RtcB")
    trace5 = scatter(x=kdam_range, y=res[:rtcr], name="RtcR")
    trace6 = scatter(x=kdam_range, y=res[:rm_r], name="mRNA RtcR")
    trace7 = scatter(x=kdam_range, y=res[:rt], name="Rt")
    trace8 = scatter(x=kdam_range, y=res[:rd], name="Rd")
    trace9 = scatter(x=kdam_range, y=res[:rh], name="Rh")

    return plot([trace1,trace2,trace3,trace4,trace5,trace6,trace7,trace8,trace9], Layout(title="Steady-state", xaxis_title="kdam", yaxis_title="concentration"))
end

kdam_range = range(0, 3,length=50)
initial = [0, 0, 0.1, 0, 0, 0, 100, 0, 0];
res = change_param(kdam_range, :kdam, rtc_model, initial, get_ssval, all_species, params)

plot_res(res)





# using Plots
#plots plotting
p3 = Plots.plot(kdam_range, res[:rm_a], label="mRNA RtcBA")
p3 = Plots.plot!(kdam_range, res[:rtca], label="RtcA")
p3 = Plots.plot!(kdam_range, res[:rtcb], label="RtcB")
# plot!(kdam_range, res[:rm_b], label="RtcA")
p3 = Plots.plot!(kdam_range, res[:rtcr], label="RtcR")
# p3 = plot!([2.01774465,2.01774465],[0,0.061])
# p3 = plot!([0.63571146,0.63571146],[0,0.061])

p4 = Plots.plot(kdam_range, res[:rh], label="Rh")
p4 = Plots.plot!(kdam_range, res[:rt], label="Rt")
p4 = Plots.plot!(kdam_range, res[:rd], label="Rd")
# p4 = plot!([2.01774465,2.01774465],[0,3.1])
# p4 = plot!([0.63571146,0.63571146],[0,3.1])
