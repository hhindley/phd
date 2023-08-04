using Parameters, CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, Statistics
# using Plots
using PlotlyJS
include("/home/holliehindley/phd/may23_rtc/functions/solving.jl"); include("/home/holliehindley/phd/may23_rtc/functions/set_ups.jl"); include("/home/holliehindley/phd/may23_rtc/functions/plotting.jl"); 
include("/home/holliehindley/phd/may23_rtc/functions/sweep_params.jl"); include("/home/holliehindley/phd/may23_rtc/models/rtc_orig.jl"); include("/home/holliehindley/phd/may23_rtc/models/atp_lam_kin_t.jl"); 
include("/home/holliehindley/phd/may23_rtc/models/single_t.jl"); include("/home/holliehindley/phd/may23_rtc/models/combinations_t.jl");



# t, atp_t, lam_t, kin_t = set_time_vars("/home/holliehindley/phd/data/atp_for_rtcmodel.csv")

# p = make_subplots(rows=3, cols=1, shared_xaxes=true, vertical_spacing=0.08, subplot_titles=["λ" "ATP" "kin"])
# add_trace!(p, (scatter(x=t, y=lam_t)), row=1, col=1)
# add_trace!(p, (scatter(x=t, y=atp_t)), row=2, col=1)
# add_trace!(p, (scatter(x=t, y=kin_t)), row=3, col=1)
# relayout!(p, showlegend=false, xaxis_range=(0,1500))
# p

# tspan = (0, lam_t[end])
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

tspan = (0,1e9)
params1 = @LArray [L, c, kr, Vmax_init, Km_init, 0.05623413251903491, 0.010000000000000002, θtscr, g_max, θtlr, km_a, km_b, d, krep, 1., ktag, kdeg, 0.022222222, 3578.9473684210525, na, nb, nr, 0.014] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
initial =  [rm_a_0, rtca_0, rm_b_0, rtcb_0, rm_r_0, rtcr_0, rh_0, rd_0, rt_0]
initial = [0., 0., 0., 0., 0., 0., 11.29, 0., 0.]

solu = sol(rtc_model, initial, tspan, params1)
plotly_plot_sol(solu, "", "", "")
ss_init = ss_init_vals(solu)
soluss = sol(rtc_model, ss_init, tspan, params1)
plotly_plot_sol(soluss, "", "", "")

initial1 = [0., 0., 1., 0., 0., 0., 11.29, 0., 0.]
solu1 = sol(rtc_model, initial1, tspan, params1)
ss_init_upper = ss_init_vals(solu1)

#define ss for each branch, remember that rt and rd are going to be opposite to the rest of the species
lower_branch = DataFrame(species=["rm_a","rtca","rm_b","rtcb","rm_r","rtcr","rh","rd","rt"],ss_val=[ss_init[1],ss_init[2],ss_init[3],ss_init[4],ss_init[5],ss_init[6],ss_init[7],ss_init[8],ss_init[9]],parameter=["ATP","kin","ω_ab","ω_r","λ","kdam",NaN,NaN,NaN], param_val=[3578.9473684210525,0.022222222,0.05623413251903491,0.010000000000000002,0.014,1.,NaN, NaN, NaN])
upper_branch = DataFrame(species=["rm_a","rtca","rm_b","rtcb","rm_r","rtcr","rh","rd","rt"],ss_val=[ss_init_upper[1],ss_init_upper[2],ss_init_upper[3],ss_init_upper[4],ss_init_upper[5],ss_init_upper[6],ss_init_upper[7],ss_init_upper[8],ss_init_upper[9]],parameter=["ATP","kin","ω_ab","ω_r","λ","kdam",NaN,NaN,NaN], param_val=[3578.9473684210525,0.022222222,0.05623413251903491,0.010000000000000002,0.014,1.,NaN, NaN, NaN])
# CSV.write("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/lower_branch.csv", lower_branch)
# CSV.write("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/upper_branch.csv", upper_branch)


kdam_range = range(0,3,length=20)
res = change_param(kdam_range, :kdam, rtc_model, ss_init, get_ssval, all_species, params1)
plot_change_param_sols(kdam_range, res, :kdam, "","")

res1 = change_param(kdam_range, :kdam, rtc_model, initial, get_ssval, all_species, params1)
plot_change_param_sols(kdam_range, res1, :kdam, "","")

ss_init






lower_branch = DataFrame(CSV.File("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/lower_branch.csv"))
upper_branch = DataFrame(CSV.File("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/upper_branch.csv"))


function set_ss_range(branch_df, specie, n, l)
    a = range((branch_df[branch_df.species .== specie,:].ss_val-(branch_df[branch_df.species .== specie,:].ss_val)),branch_df[branch_df.species .== specie,:].ss_val,length=Int(l/2))
    b = range((branch_df[branch_df.species .== specie,:].ss_val),(branch_df[branch_df.species .== specie,:].ss_val+(n*branch_df[branch_df.species .== specie,:].ss_val)),length=Int(l/2))
    return [(vcat(a,b)...)...]
end

# function set_ss_range(branch_df, specie, n)
#     b = range((branch_df[branch_df.species .== specie,:].ss_val),(branch_df[branch_df.species .== specie,:].ss_val+(n*branch_df[branch_df.species .== specie,:].ss_val)),length=500)
#     return [(b...)...]
# end
function get_all_ranges(branch, n, l)
    all_ranges=[]
    species=["rm_a","rtca","rm_b","rtcb","rm_r","rtcr","rh","rd","rt"]
    for i in species
        push!(all_ranges, set_ss_range(branch, i, n, l))
    end

    return all_ranges
end

function get_ss_change_init(init, range, l)
    initial = deepcopy(init)
    res=[]
    for specie in all_species
        for i in range
            initial[7] = i
            @show initial
            solu = sol(rtc_model, initial, tspan, params1)
            push!(res, get_ssval(solu,specie))
            # push!(res, [get_ssval(solu,specie) for specie in all_species])
        end
    end
    # @show res
    res = Float64.(res)
    res = (reshape(res, (l,9)))
    # @show typeof(res)
    df = DataFrame(res,all_species)
    return df
end


function get_rh_init_switch_all_ranges(ranges, branch, l)
    res=[]
    for (range,i) in zip(ranges,range(1,9))
        initial = deepcopy(branch.ss_val)
        # res=[]
        for j in range
            initial[i] = j
            @show initial
            solu = sol(rtc_model, initial, tspan, params1)
            push!(res, get_ssval(solu,:rh))
        end
        # push!(all_res,res)
    end

    res = Float64.(res)
    res = (reshape(res, (l,9)))
    return DataFrame(res,all_species)
end

# function get_all_init_switch_result(ranges, init)
#     all=[]
#     for range in ranges
#         push!(all, get_ss_change_init(init, range))
#     end
#     return all
# end

function upper_or_lower(df, state,l)
    arr=[]
    state1 = "$state"
    if state1 == "lower"
        for col in eachcol(df)
            for i in col
                for j in i 
                    if round(j;digits=10) == round(col[1];digits=10)
                        push!(arr, 0)
                    else
                        push!(arr, 1)
                    end
                end
            end
        end
    else
        for col in eachcol(df)
            for i in col
                for j in i 
                    if round(j;digits=10) == round(col[1];digits=10)
                        push!(arr, 1)
                    else
                        push!(arr, 0)
                    end
                end
            end
        end
    end
    arr = reshape(arr, (l,9))
    df = DataFrame(arr,all_species)
    return df
end

# function all_upper_lower(all, state)
#     dfs=[]
#     for df in all
#         push!(dfs, upper_or_lower(df, state))
#     end
#     return dfs
# end

function set_shared_range(n,l)
    x1 = range(-100,0, length=Int(l/2))
    x2 = range(0,n*100, length=Int(l/2))
    return vcat(x1,x2)
end
# function set_shared_range(n)
#     # x1 = range(-100,0, length=200)
#     x2 = range(0,n*100, length=500)
#     return x2
# end
n = 620; l = 500;
lower_ranges = get_all_ranges(lower_branch, n, l)
# all = get_all_init_switch_result(lower_ranges, lower_branch.ss_val)
all = get_rh_init_switch_all_ranges(lower_ranges, lower_branch, l)
dfs = upper_or_lower(all, "lower", l)
shared_range = set_shared_range(n,l)

p = plot([scatter(x=shared_range,y=dfs.rm_a, name="rm_a"),scatter(x=shared_range,y=dfs.rtca, name="rtca"),scatter(x=shared_range,y=dfs.rm_b, name="rm_b"),
scatter(x=shared_range,y=dfs.rtcb,name="rtcb"),scatter(x=shared_range,y=dfs.rm_r,name="rm_r"),scatter(x=shared_range,y=dfs.rtcr,name="rtcr"),
scatter(x=shared_range,y=dfs.rh,name="rh"),scatter(x=shared_range,y=dfs.rd,name="rd"),scatter(x=shared_range,y=dfs.rt,name="rt")],
Layout(xaxis_title="% increase from initial ss val", yaxis_title="Branch",
yaxis_tickvals=[0,1], yaxis_ticktext=["lower","upper"], title="switch from lower to upper (rh)"
))


open("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/init_switch.html", "w") do io
    PlotlyBase.to_html(io, p.plot)
end


n1 = 7700; l1=1000;
upper_ranges = get_all_ranges(upper_branch, n1,l1)
all1 = get_rh_init_switch_all_ranges(upper_ranges, upper_branch, l1)
dfs1 = upper_or_lower(all1, "upper",l1)
shared_range1 = set_shared_range(n1,l1)

p = plot([scatter(x=shared_range1,y=dfs1.rm_a, name="rm_a"),scatter(x=shared_range1,y=dfs1.rtca, name="rtca"),scatter(x=shared_range1,y=dfs1.rm_b, name="rm_b"),
scatter(x=shared_range1,y=dfs1.rtcb,name="rtcb"),scatter(x=shared_range1,y=dfs1.rm_r,name="rm_r"),scatter(x=shared_range1,y=dfs1.rtcr,name="rtcr"),
scatter(x=shared_range1,y=dfs1.rh,name="rh"),scatter(x=shared_range1,y=dfs1.rd,name="rd"),scatter(x=shared_range1,y=dfs1.rt,name="rt")],
Layout(xaxis_title="% increase from initial ss val", yaxis_title="Branch",
yaxis_tickvals=[0,1], yaxis_ticktext=["lower","upper"], title="switch from upper to lower (rh)"
))





switch1s=[]
switch2s=[]
for i in range(1,9)
    if length(dfs[i][dfs[i].rh .== 0,:].rh) == 500
        Nothing
    else
        switch1 = lower_ranges[i][length(dfs[i][dfs[i].rh .== 0,:].rh)]
        switch2 = lower_ranges[i][length(dfs[i][dfs[i].rh .== 0,:].rh)+1]
        push!(switch1s, switch1)
        push!(switch2s, switch2)
    end
end
switch1s
switch2s
df = DataFrame(species=["rm_a","rm_b","rm_r","rh","rd","rt"])
insertcols!(df,2, "ss_vals" => [lower_branch.ss_val[1], lower_branch.ss_val[3], lower_branch.ss_val[5], lower_branch.ss_val[7], lower_branch.ss_val[8], lower_branch.ss_val[9]], "switch1" => switch1s, "switch2" => switch2s)

upper_branch.ss_val[9]

length(dfs[2][dfs[2].rh .== 0,:].rh)

lower_ranges[1][length(dfs[1][dfs[1].rh .== 0,:].rh)]
lower_ranges[1][length(dfs[1][dfs[1].rh .== 0,:].rh)+1]

lower_ranges[2][178]
lower_ranges[2][179]



plot([scatter(x=lower_ranges[7],y=dfs[7].rh),scatter(x=lower_ranges[7],y=dfs[7].rd),scatter(x=lower_ranges[7],y=dfs[7].rt),scatter(x=lower_ranges[7],y=dfs[7].rm_a),scatter(x=lower_ranges[7],y=dfs[7].rtca),scatter(x=lower_ranges[7],y=dfs[7].rm_b),scatter(x=lower_ranges[7],y=dfs[7].rtcb),scatter(x=lower_ranges[7],y=dfs[7].rtcr)])














# rh2 = get_init_switch(:rh, lower_ranges, lower_branch.ss_val)
# new_df = DataFrame(rm_a_0=[], rtca_0=[], rm_b_0=[], rtcb_0=[], rm_r_0=[], rtcr_0=[], rh_0=[], rd_0=[], rt_0=[])
# for (col, new_col) in zip(eachcol(rh2),eachcol(new_df))
#     for val in col
#         if val == col[1]
#             push!(new_col, 0)
#         else
#             push!(new_col, 1)
#         end
#     end
# end
# new_df

# x1 = range(-100,0, length=10)
# x2 = range(0,100000, length=10)
# shared_range = vcat(x1,x2)

# insertcols!(new_df,1,"xaxis"=>vcat(x1,x2))
# new_df
# p2 = plot(scatter(x=new_df.xaxis, y=new_df.rt_0, mode="markers"))
# # p2 = plot(scatter(x=lower_ranges[9], y=new_df.rt_0, mode="markers"))


# p1 = (scatter(x=new_df.xaxis, y=new_df.rm_a_0, name="rm_a"))#, mode="markers"))
# p2 = (scatter(x=new_df.xaxis, y=new_df.rtca_0, name="rtca"))#, mode="markers"))
# p3 = (scatter(x=new_df.xaxis, y=new_df.rm_b_0, name="rm_b"))#, mode="markers"))
# p4 = (scatter(x=new_df.xaxis, y=new_df.rtcb_0, name="rtcb"))#, mode="markers"))
# p5 = (scatter(x=new_df.xaxis, y=new_df.rm_r_0, name="rm_r"))#, mode="markers"))
# p6 = (scatter(x=new_df.xaxis, y=new_df.rtcr_0, name="rtcr"))#, mode="markers"))
# p7 = (scatter(x=new_df.xaxis, y=new_df.rh_0, name="rh"))#, mode="markers"))
# p8 = (scatter(x=new_df.xaxis, y=new_df.rd_0, name="rd"))#, mode="markers"))
# p9 = (scatter(x=new_df.xaxis, y=new_df.rt_0, name="rt"))#, mode="markers"))
# plot([p1,p2,p3,p4,p5,p6,p7,p8,p9])



# all_df = []
# for i in all_species
#     push!(all_df, get_init_switch(i,lower_ranges, lower_branch.ss_val))
# end

# all_df = []
# for i in all_species
#     push!(all_df, get_init_switch(i,upper_ranges, upper_branch.ss_val))
# end

# plot(scatter(x=lower_ranges[7],y=all_df_lower[7].rh_0))

# plot(scatter(x=upper_ranges[7],y=all_df_upper[7].rh_0))


















# df_range = DataFrame(rma_range=[], rtca_range=[], rmb_range=[], rtcb_range=[], rmr_range=[], rtcr_range=[], rh_range=[], rd_range=[], rt_range=[])
# for (range,col) in zip(lower_ranges,eachcol(df_range))
#     # print(([r...]...))
#     push!(col, ([range...]...))
#     # println(col)

# end
# df_range

# rm_a_switch_res = all_df[1]; rtca_switch_res = all_df[2]; rm_b_switch_res = all_df[3]; rtcb_switch_res = all_df[4]; rm_r_switch_res = all_df[5]; rtcr_switch_res = all_df[6]; rh_switch_res = all_df[7]; rd_switch_res = all_df[8] ; rt_switch_res = all_df[9]

# function plot_species_switch(xaxis,specie)
#     t1 = scatter(x=xaxis,y=rm_a_switch_res[:,specie], name="rm_a", line_color=:red)#, mode="markers")
#     t2 = scatter(x=xaxis,y=rtca_switch_res[:,specie], name="rtca", line_color=:blue)#, mode="markers")
#     t3 = scatter(x=xaxis,y=rm_b_switch_res[:,specie], name="rm_b", line_color=:green)#, mode="markers")
#     t4 = scatter(x=xaxis,y=rtcb_switch_res[:,specie], name="rtcb", line_color=:orange)#, mode="markers")
#     t5 = scatter(x=xaxis,y=rm_r_switch_res[:,specie], name="rm_r", line_color=:purple)#, mode="markers")
#     t6 = scatter(x=xaxis,y=rtcr_switch_res[:,specie], name="rtcr", line_color=:hotpink)#, mode="markers")
#     t7 = scatter(x=xaxis,y=rh_switch_res[:,specie], name="rh", line_color=:brown)#, mode="markers")
#     t8 = scatter(x=xaxis,y=rd_switch_res[:,specie], name="rd", line_color=:grey)#, mode="markers")
#     t9 = scatter(x=xaxis,y=rt_switch_res[:,specie], name="rt", line_color=:aqua)#, mode="markers")
#     return plot([t1,t2,t3,t4,t5,t6,t7,t8,t9], Layout(xaxis_title="$specie", yaxis_title="steady-state conc"))#, title="changing initial value of $specie and getting steady state value to see which bs branch it gives"))
# end

# rma_plot = plot_species_switch(df_range.rma_range,:rm_a_0);rtca_plot = plot_species_switch(df_range.rtca_range,:rtca_0);rmb_plot = plot_species_switch(df_range.rmb_range,:rm_b_0);rtcb_plot = plot_species_switch(df_range.rtcb_range,:rtcb_0);rmr_plot = plot_species_switch(df_range.rmr_range,:rm_r_0);rtcr_plot = plot_species_switch(df_range.rtcr_range,:rtcr_0);rh_plot = plot_species_switch(df_range.rh_range,:rh_0);rd_plot = plot_species_switch(df_range.rd_range,:rd_0);rt_plot = plot_species_switch(df_range.rt_range,:rt_0);

# p = [rma_plot rtca_plot rmb_plot; rtcb_plot rmr_plot rtcr_plot; rh_plot rd_plot rt_plot]
# relayout!(p,legend=false, title_text="changing intitial values to see where and which variables cause the switch - init as ss vals")
# p

# open("/home/holliehindley/phd/may23_rtc/analysis/bifurcation_analysis/plots/switch1.html", "w") do io
#     PlotlyBase.to_html(io, p.plot)
# end

# function plot_res(res)
#     #plotly plotting
#     trace1 = scatter(x=kdam_range, y=res[:rtca], name="RtcA")
#     trace2 = scatter(x=kdam_range, y=res[:rtcb], name="RtcB")
#     trace3 = scatter(x=kdam_range, y=res[:rm_a], name="mRNA RtcA")
#     trace4 = scatter(x=kdam_range, y=res[:rm_b], name="mRNA RtcB")
#     trace5 = scatter(x=kdam_range, y=res[:rtcr], name="RtcR")
#     trace6 = scatter(x=kdam_range, y=res[:rm_r], name="mRNA RtcR")
#     trace7 = scatter(x=kdam_range, y=res[:rt], name="Rt")
#     trace8 = scatter(x=kdam_range, y=res[:rd], name="Rd")
#     trace9 = scatter(x=kdam_range, y=res[:rh], name="Rh")

#     return plot([trace1,trace2,trace3,trace4,trace5,trace6,trace7,trace8,trace9], Layout(title="Steady-state", xaxis_title="kdam", yaxis_title="concentration"))
# end

# kdam_range = range(0, 3,length=50)
# initial = [0, 0, 0.1, 0, 0, 0, 100, 0, 0];
# res = change_param(kdam_range, :kdam, rtc_model, initial, get_ssval, all_species, params)

# plot_res(res)






