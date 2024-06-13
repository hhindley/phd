using StatsBase, Distributions, Random, DataFrames, CSV, PlotlyJS, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools

include(joinpath(homedir(), "phd/rtc_model/parameters/rtc_params.jl"))
include(joinpath(homedir(), "phd/rtc_model/parameters/rtc_params_molecs.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/indexing.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/hybrid_algo.jl"))
include(joinpath(homedir(), "phd/stochastic_hybrid_code/setup/stoch_model.jl"))


# include("/home/hollie_hindley/Documents/stochastic_hybrid/run_rtc_orig.jl")

n= 100#200 # number of cell cycles
options = Dict(
"threshold"  =>  0.,       # Threshold to decide between determinisitic or stochastic reaction
"FixDetReact"=> [14],# [10,11,12,13,14,15,16,17,18],       # Reactions to be treated determinisitically
    "tspan"     =>   n*log(2)/lam_val,     # Max time for cell cycle
    "samplingFreq"  => 0.1  # for sampling every x mins
)

n*log(2)/lam_val
# X0 = collect(get_X0(indV)')

X0 = collect(get_X0(indV, init_molec)')
par = collect(get_par(indP)')


getssX0 = true
if getssX0
    fout=open("/home/hollie_hindley/Documents/stochastic_hybrid/X0.dat","w")
    propen, S, propList = defineStochModel(par, indV)
    nx = indV.nrOfItems-1
    prop(X) = propen(X[1:nx])
    X0 = hybrid_algo(X0, options, prop, S, out=fout)
    X0[vidx(:V)] = 1
    # run_stoch(X0, 0, 0, "X0")
    # df = DataFrame(CSV.File("/home/hollie_hindley/Documents/stochastic_hybrid/X0.dat", header=["event", "time", "rm_a", "rtca", "rm_b", "rtcb", "rm_r", "rtcr", "rh", "rd", "rt", "volume"]))
    # ss_df = df[1000:end,:]
    # ss = [mean(df[:,col]) for col in names(eachcol(df[:,3:end-2]))]
    # X0 = collect(get_X0(indV, ss)')
    CSV.write("/home/hollie_hindley/Documents/stochastic_hybrid/X0.dat", DataFrame(X0,:auto), header=false)
else
    X0 = CSV.read("/home/hollie_hindley/Documents/stochastic_hybrid/X0.dat", Tables.matrix, header=false)
end


# X0 = collect(get_X0(indV, ssvals_rtc_molec)')
# time_taken = @elapsed run_stoch(0, X0, 20, 0.6, "test_nofloor.dat")
# df_n = DataFrame(CSV.File("/home/hollie_hindley/Documents/stochastic_hybrid/test_nofloor.dat", header=["event", "time", "rm_a", "rtca", "rm_b", "rtcb", "rm_r", "rtcr", "rh", "rd", "rt", "volume", "totprop"]))

# plot(scatter(x=df_n.time, y=df_n.rm_r))

time_taken = @elapsed run_stoch(X0, 20, 0.2, "test_floor.dat")
df = DataFrame(CSV.File("/home/hollie_hindley/Documents/stochastic_hybrid/test_floor.dat", header=["event", "time", "rm_a", "rtca", "rm_b", "rtcb", "rm_r", "rtcr", "rh", "rd", "rt", "volume", "totprop"]))

plot(scatter(x=df.time, y=df.rm_r))#./df.volume))
plot(scatter(x=df.time, y=df.rtca))#./df.volume))
plot(scatter(x=df.time, y=df.rh./df.volume))
plot(scatter(x=df.time, y=df.rtcr))#./df.volume))

df.bins = ceil.(Int, (1:nrow(df))/10)
df_grouped = combine(first, groupby(df, :bins))
plot(histogram(x=df_grouped.rtca, nbinsx=100))#, Layout(yaxis_range=(0,200)))

# # plot([scatter(x=df.time, y=df.rd./df.volume),scatter(x=df.time, y=df.rh./df.volume)])

# alpha = df.rt ./par[pidx(:kr)]
# fa = @. (1+alpha)^6/(par[pidx(:L)]*((1+par[pidx(:c)]*alpha)^6)+(1+alpha)^6)
# ra = @. fa*df.rtcr

# plot(scatter(x=df.time, y=fa))
# plot([scatter(x=df.time,y=df.rtcr),scatter(x=df.time,y=fa)])
# # for i in 1:20
# #     time_taken = @elapsed run_stoch(1, X0, 200, 0.6, "run_multiple/test_floor_$i.dat")
# # end
# # dfs=[]
# # folder = "/home/hollie_hindley/Documents/stochastic_hybrid/run_multiple"
# # files = readdir(folder)
# # for file in files
# #     filepath = joinpath(folder,file)
# #     push!(dfs, DataFrame(CSV.File(filepath, header=["event", "time", "rm_a", "rtca", "rm_b", "rtcb", "rm_r", "rtcr", "rh", "rd", "rt", "volume", "totprop"])))
# # end

# function df_sort(df)
#     df.event = [split(replace(i, r"[\[\]\(Any)]" => ""), ",") for i in df.event]
#     df.event = [parse.(Float64, subarray) for subarray in df.event]
#     df_p = filter(row -> length(row[:event]) > 1, df)
#     df_r = filter(row -> length(row[:event]) == 1 && row[:event][1] != 0, df)
#     df_r.event = first.(df_r.event)

#     props = df_p.event
#     react_names = [:tscr_ab, :tscr_r, :tlr_a, :tlr_b, :tlr_r, :Vinflux, :Vdam, :Vtag, :Vrep, :deg_rd, :deg_rma, :deg_rmb, :deg_rmr, :V]

#     df_props = DataFrame([name => Float64[] for name in react_names])

#     for i in props
#         push!(df_props, [i[j] for j in 1:length(props[end])])
#     end
#     return df_props, df_r, df_p
# end

# plot([scatter(x=dfs[i][1:100:end,:time], y=dfs[i][1:100:end,:rm_a]) for i in 1:length(dfs)])

# df

# df_props, df_r, df_p = df_sort(df)
# df_props1, df_r1, df_p1 = df_sort(df_n)

# react_names = [:tscr_ab, :tscr_r, :tlr_a, :tlr_b, :tlr_r, :Vinflux, :Vdam, :Vtag, :Vrep, :deg_rd, :deg_rma, :deg_rmb, :deg_rmr]

# # s1 = plot(scatter(x=df_p.time, y=df_p.rtca), Layout(title="flooring"));
# # pr1 = plot(scattergl(x=df_r.time, y=df_r.event, mode="markers", marker_color=df_r.event, marker=attr(colorscale="Viridis")), Layout(yaxis=attr(tickvals=range(1,13), ticktext=react_names),));# title="flooring"))
# # p_props = plot([scatter(x=df_p.time, y=df_props[1:100:end,col], name="$col", mode="markers") for col in names(eachcol(df_props[:,1:end-1]))], );#Layout(title="flooring"))
# # t1 = plot(scatter(x=df.time, y=df.totprop));#, Layout(title="flooring"));

# # [s1; pr1; p_props; t1]

# s1 = plot(scatter(x=df_p.time, y=df_p.rtca), Layout(title="flooring"));
# s2 = plot(scatter(x=df_p1.time, y=df_p1.rtca), Layout(title="no flooring"));

# pr1 = plot(scattergl(x=df_r.time, y=df_r.event, mode="markers", marker_color=df_r.event, marker=attr(colorscale="Viridis")), Layout(yaxis=attr(tickvals=range(1,13), ticktext=react_names),));# title="flooring"))
# pr2 = plot(scattergl(x=df_r1.time, y=df_r1.event, mode="markers", marker_color=df_r1.event, marker=attr(colorscale="Viridis")), Layout(yaxis=attr(tickvals=range(1,13), ticktext=react_names),));# title="no flooring"));

# p_props = plot([scatter(x=df_p.time, y=df_props[1:100:end,col], name="$col", mode="markers") for col in names(eachcol(df_props[:,1:end-1]))], );#Layout(title="flooring"))
# p_props1 = plot([scatter(x=df_p1.time, y=df_props1[1:100:end,col], name="$col", mode="markers") for col in names(eachcol(df_props1[:,1:end-1]))],);# Layout(title="no flooring"));

# t1 = plot(scatter(x=df.time, y=df.totprop));#, Layout(title="flooring"));
# t2 = plot(scatter(x=df_n.time, y=df_n.totprop));#, Layout(title="no flooring"));

# [s1 s2; pr1 pr2; p_props p_props1; t1 t2]

# d

# plot(scatter(x=df.time, y=df.rt), Layout(title="flooring"))

# # fl=[]; nf=[];
# # for i in 1:20
# #     time_taken = @elapsed run_stoch(0, X0, 20, 0.6, "test_nofloor")
# #     df_n = DataFrame(CSV.File("/home/hollie_hindley/Documents/stochastic_hybrid/test_nofloor.dat", header=["event", "time", "rm_a", "rtca", "rm_b", "rtcb", "rm_r", "rtcr", "rh", "rd", "rt", "volume", "totprop"]))

# #     time_taken = @elapsed run_stoch(1, X0, 20, 0.6, "test_floor")
# #     df = DataFrame(CSV.File("/home/hollie_hindley/Documents/stochastic_hybrid/test_floor.dat", header=["event", "time", "rm_a", "rtca", "rm_b", "rtcb", "rm_r", "rtcr", "rh", "rd", "rt", "volume", "totprop"]))

# #     df_props, df_r, df_p = df_sort(df)
# #     df_props1, df_r1, df_p1 = df_sort(df_n)

# #     push!(fl, length(df_r.time))
# #     push!(nf, length(df_r1.time))
# # end

# # plot([scatter(x=1:50, y=fl, name="flooring"), scatter(x=1:50, y=nf, name="no flooring")], Layout(yaxis_title="Number of stoch reactions"))







# df_n.event = [split(replace(i, r"[\[\]\(Any)]" => ""), ",") for i in df_n.event]
# df_n.event = [parse.(Float64, subarray) for subarray in df_n.event]
# df_ns = filter(row -> length(row[:event]) > 1, df_n)

# df.event = [split(replace(i, r"[\[\]\(Any)]" => ""), ",") for i in df.event]
# df.event = [parse.(Float64, subarray) for subarray in df.event]
# df_s = filter(row -> length(row[:event]) > 1, df)

# function plotprops(df)
#     props = df.event
#     react_names = [:tscr_ab, :tscr_r, :tlr_a, :tlr_b, :tlr_r, :Vinflux, :Vdam, :Vtag, :Vrep, :deg_rd, :deg_rma, :deg_rmb, :deg_rmr, :V]

#     df_props = DataFrame([name => Float64[] for name in react_names])

#     for i in props
#         push!(df_props, [i[j] for j in 1:length(props[end])])
#     end
#     return plot([scatter(x=df.time, y=df_props[1:end,col], name="$col", mode="markers") for col in names(eachcol(df_props[:,1:end-1]))])

# end
# p_props = plotprops(df_ns);
# p_props1 = plotprops(df_s);
# [plot(scatter(x=df_n.time,y=df_n.rm_a)) plot(scatter(x=df.time,y=df.rm_a))]



# df.event = [eval(Meta.parse(i)) for i in df.event]
# df_r = filter(row -> length(row[:event]) == 1 && row[:event][1] != 0, df)

# df_n.event = [eval(Meta.parse(i)) for i in df_n.event]
# df_rn = filter(row -> length(row[:event]) == 1 && row[:event][1] != 0, df_n)

# react_names = [:tscr_ab, :tscr_r, :tlr_a, :tlr_b, :tlr_r, :Vinflux, :Vdam, :Vtag, :Vrep, :deg_rd, :deg_rma, :deg_rmb, :deg_rmr]
# pr1 = plot(scattergl(x=df_r.time, y=df_r.event, mode="markers", marker_color=df_r.event, marker=attr(colorscale="Viridis")), Layout(yaxis=attr(tickvals=range(1,13), ticktext=react_names)));
# pr2 = plot(scattergl(x=df_rn.time, y=df_rn.event, mode="markers", marker_color=df_rn.event, marker=attr(colorscale="Viridis")), Layout(yaxis=attr(tickvals=range(1,13), ticktext=react_names)));
# [plot(scatter(x=df_n.time, y=df_n.rm_a)) plot(scatter(x=df.time, y=df.rm_a));pr2 pr1]


# df_t = filter(row -> length(row[:event]) == 2, df)

# before = first.(df_t.event)
# after = last.(df_t.event)
# df_t

# df.event = [eval(Meta.parse(i)) for i in df.event]
# dfa = filter(row -> length(row[:event]) > 1, df)

# rma = []; rtca = []; rmb = []; rtcb = []; rmr = []; rtcr = []; rh = []; rd = []; rt = []
# for i in dfa.event
#     push!(rma, i[1])
#     push!(rtca, i[2])
#     push!(rmb, i[3])    
#     push!(rtcb, i[4])
#     push!(rmr, i[5])
#     push!(rtcr, i[6])
#     push!(rh, i[7])
#     push!(rd, i[8])
#     push!(rt, i[9])
# end

# plot([scatter(x=df.time, y=df.rm_a), scatter(x=dfa.time, y=rma)])
# plot([scatter(x=df.time, y=df.rtca), scatter(x=dfa.time, y=rtca)])
# plot([scatter(x=df.time, y=df.rm_b), scatter(x=dfa.time, y=rmb)])
# plot([scatter(x=df.time, y=df.rtcb), scatter(x=dfa.time, y=rtcb)])
# plot([scatter(x=df.time, y=df.rm_r), scatter(x=dfa.time, y=rmr)])
# plot([scatter(x=df.time, y=df.rtcr), scatter(x=dfa.time, y=rtcr)])
# plot([scatter(x=df.time, y=df.rh), scatter(x=dfa.time, y=rh)])
# plot([scatter(x=df.time, y=df.rd), scatter(x=dfa.time, y=rd)])
# plot([scatter(x=df.time, y=df.rt), scatter(x=dfa.time, y=rt)])


# before = first.(df_s.event)
# after = last.(df_s.event)
# plot([scatter(x=df_t.time, y=before./df_t.volume),scatter(x=df_t.time, y=after./df_t.volume)])
# # df_f = DataFrame(CSV.File("/home/hollie_hindley/Documents/stochastic_hybrid/test.dat", header=["event", "time", "rm_a", "rtca", "rm_b", "rtcb", "rm_r", "rtcr", "rh", "rd", "rt", "volume"]))

# # df_rfloor = filter(row -> length(row[:event]) == 1 && row[:event][1] != 0, df_floor)
# # df_r = filter(row -> length(row[:event]) == 1 && row[:event][1] != 0, df)

# [plot(scatter(x=df_n.time,y=df_n.rm_a)) plot(scatter(x=df.time,y=df.rm_a))]
# [plot(scatter(x=df_n.time,y=df_n.event)) plot(scatter(x=df.time,y=df.event))]

# [plot(scatter(x=df.time,y=df.rh)) plot(scatter(x=df_floor.time,y=df_floor.rh))]

# [plot(scatter(x=df_r.time,y=df_r.rm_a)) plot(scatter(x=df_rfloor.time,y=df_rfloor.rm_a))]
# [plot(scatter(x=df.time,y=df.rh)) plot(scatter(x=df_floor.time,y=df_floor.rh))]

# df_n.event = [split(replace(i, r"[\[\]\(Any)]" => ""), ",") for i in df_n.event]
# df_n.event = [parse.(Float64, subarray) for subarray in df_n.event]

# df_t = filter(row -> row[:event] == [0], df)
# df_v = filter(row -> row[:event] == [20], df)
# df_p1 = filter(row -> length(row[:event]) > 1, df_n)
# df_r = filter(row -> length(row[:event]) == 1 && row[:event][1] != 0, df)
# df_r.event = first.(df_r.event)
# # plot([scattergl(x=df.time, y=df[:,col], name="$(names(df)[i])", marker_color=colours[i], legendgroup="$i", showlegend=true) for (col, i) in zip(names(eachcol(df[:,3:end-2])), range(3,length(names(df))-2))])#, title="kdam = $(params_rtc[kdam])"))
# # plot([scattergl(x=df.time, y=df[:,col] ./df.volume, name="$(names(df)[i])", marker_color=colours[i], legendgroup="$i", showlegend=true) for (col, i) in zip(names(eachcol(df[:,3:end-2])), range(3,length(names(df))-2))])#, title="kdam = $(params_rtc[kdam])"))
# df_p.event[120]
# df_p.event = [split(replace(i, r"[\[\]\(LinearAlgebra.Adjoint{Float64})]" => ""), ",") for i in df_p.event]

# # df_p.totprop = map(x -> x[1], df_p.event)
# # df_p.xi = map(x -> length(x) > 1 ? x[2] : NaN, df_p.event)
# # plot([scatter(x=df_p.time, y=df_p.totprop), scatter(x=df_p.time, y=df_p.xi)])

# plot(scatter(x=df_p.time, y=df_p.rd))# ./df_p.volume))
# mean(df.rtca[10000:end])#./df.volume[10000:end])

# X0

# c1 = [split(replace(i, r"" => ""), ",") for i in c]
# df1.event = [parse.(Float64, subarray) for subarray in df1.event]

# df[:,:rm_b]
# any(df_grouped[:,:rm_a] .< 0)


# df.bins = ceil.(Int, (1:nrow(df))/100)
# df_grouped = combine(first, groupby(df, :bins))
# bins = 25#length(df_grouped.bins)
# h_rh = plot(histogram(x=df_grouped.rh ./df_grouped.volume, nbinsx=bins));
# h_rma = plot(histogram(x=df_grouped.rm_a ./df_grouped.volume, nbinsx=bins));
# h_rmb = plot(histogram(x=df_grouped.rm_b ./df_grouped.volume, nbinsx=bins));
# h_rmr = plot(histogram(x=df_grouped.rm_r ./df_grouped.volume, nbinsx=bins));
# h_rtca = plot(histogram(x=df_grouped.rtca ./df_grouped.volume, nbinsx=bins));
# h_rtcb = plot(histogram(x=df_grouped.rtcb ./df_grouped.volume, nbinsx=bins));
# h_rtcr = plot(histogram(x=df_grouped.rtcr ./df_grouped.volume, nbinsx=bins));
# h_rd = plot(histogram(x=df_grouped.rd ./df_grouped.volume, nbinsx=bins));
# h_rt = plot(histogram(x=df_grouped.rt ./df_grouped.volume, nbinsx=bins));


# plot(histogram(x=df_grouped.rh, nbinsx=bins))
# plot(histogram(x=df_grouped.rm_a, nbinsx=bins))
# plot(histogram(x=df_grouped.rm_b, nbinsx=bins))
# plot(histogram(x=df_grouped.rm_r, nbinsx=bins))
# plot(histogram(x=df_grouped.rtca, nbinsx=bins))
# plot(histogram(x=df_grouped.rtcb, nbinsx=bins))
# plot(histogram(x=df_grouped.rtcr, nbinsx=bins))
# plot(histogram(x=df_grouped.rd, nbinsx=bins))
# plot(histogram(x=df_grouped.rt, nbinsx=bins))

# mean(df.rm_a[10000:end])


# time_taken = @elapsed run_stoch(X0, 10, 0.6, "thresh10_kdam06")
# df = DataFrame(CSV.File("/home/hollie_hindley/Documents/stochastic_hybrid/X0_init/thresh10_kdam01.dat", header=["event", "time", "rm_a", "rtca", "rm_b", "rtcb", "rm_r", "rtcr", "rh", "rd", "rt", "volume"]))

# time_taken = @elapsed run_stoch(X0, 15, 0.6, "thresh15_kdam06")
# df1 = DataFrame(CSV.File("/home/hollie_hindley/Documents/stochastic_hybrid/X0_init/thresh15_kdam01.dat", header=["event", "time", "rm_a", "rtca", "rm_b", "rtcb", "rm_r", "rtcr", "rh", "rd", "rt", "volume"]))

# time_taken = @elapsed run_stoch(X0, 20, 0.6, "thresh20_kdam06")
# df2 = DataFrame(CSV.File("/home/hollie_hindley/Documents/stochastic_hybrid/X0_init/thresh20_kdam01.dat", header=["event", "time", "rm_a", "rtca", "rm_b", "rtcb", "rm_r", "rtcr", "rh", "rd", "rt", "volume"]))

# time_taken = @elapsed run_stoch(X0, 100, 0.6, "thresh100_kdam06")
# df3 = DataFrame(CSV.File("/home/hollie_hindley/Documents/stochastic_hybrid/X0_init/thresh50_kdam01.dat", header=["event", "time", "rm_a", "rtca", "rm_b", "rtcb", "rm_r", "rtcr", "rh", "rd", "rt", "volume"]))

# s = :rm_a
# p1 = plot(scatter(x=df.time,y=df[:,s] ./df.volume))
# p2 = plot(scatter(x=df1.time,y=df1[:,s] ./df1.volume));
# p3 = plot(scatter(x=df2.time,y=df2[:,s] ./df2.volume))
# p4 = plot(scatter(x=df3.time,y=df3[:,s] ./df3.volume))

# p5 = plot(scatter(x=df5.time,y=df5[:,s] ./df5.volume));
# p5
# [p1 p2 p3 p4]



# df_props, df_r, df_p = df_sort(df)
# df_props1, df_r1, df_p1 = df_sort(df1)
# df_props2, df_r2, df_p2 = df_sort(df2)
# df_props3, df_r3, df_p3 = df_sort(df3)

# react_names = [:tscr_ab, :tscr_r, :tlr_a, :tlr_b, :tlr_r, :Vinflux, :Vdam, :Vtag, :Vrep, :deg_rd, :deg_rma, :deg_rmb, :deg_rmr]
# pr1 = plot(scattergl(x=df_r.time, y=df_r.event, mode="markers", marker_color=df_r.event, marker=attr(colorscale="Viridis")), Layout(yaxis=attr(tickvals=range(1,13), ticktext=react_names)));
# pr2 = plot(scattergl(x=df_r1.time, y=df_r1.event, mode="markers", marker_color=df_r1.event, marker=attr(colorscale="Viridis")), Layout(yaxis=attr(tickvals=range(1,13), ticktext=react_names)));
# pr3 = plot(scattergl(x=df_r2.time, y=df_r2.event, mode="markers", marker_color=df_r2.event, marker=attr(colorscale="Viridis")), Layout(yaxis=attr(tickvals=range(1,13), ticktext=react_names)));
# pr4 = plot(scattergl(x=df_r3.time, y=df_r3.event, mode="markers", marker_color=df_r3.event, marker=attr(colorscale="Viridis")), Layout(yaxis=attr(tickvals=range(1,13), ticktext=react_names)));

# [pr1 pr2 pr3 pr4]
# length(df_props.tlr_a),length(df_props1.tlr_a),length(df_props2.tlr_a),length(df_props3.tlr_a)

# p_props = plot([scatter(x=df_p.time, y=df_props[1:100:end,col], name="$col", mode="markers") for col in names(eachcol(df_props[:,1:end-1]))])
# p_props1 = plot([scatter(x=df_p1.time, y=df_props1[1:100:end,col], name="$col", mode="markers") for col in names(eachcol(df_props1[:,1:end-1]))]);
# p_props2 = plot([scatter(x=df_p2.time, y=df_props2[1:100:end,col], name="$col", mode="markers") for col in names(eachcol(df_props2[:,1:end-1]))]);
# p_props3 = plot([scatter(x=df_p3.time, y=df_props3[1:100:end,col], name="$col", mode="markers") for col in names(eachcol(df_props3[:,1:end-1]))]);

# [p_props p_props1 p_props2 p_props3]

# function plot_hist(df)
#     df.bins = ceil.(Int, (1:nrow(df))/100)
#     df_grouped = combine(first, groupby(df, :bins))
#     bins = 59#length(df_grouped.bins)
#     h_rh = plot(histogram(x=df_grouped.rh ./df_grouped.volume, nbinsx=bins, name="rh"));
#     h_rma = plot(histogram(x=df_grouped.rm_a ./df_grouped.volume, nbinsx=bins, name="rm_a"));
#     h_rmb = plot(histogram(x=df_grouped.rm_b ./df_grouped.volume, nbinsx=bins, name="rm_b"));
#     h_rmr = plot(histogram(x=df_grouped.rm_r ./df_grouped.volume, nbinsx=bins, name="rm_r"));
#     h_rtca = plot(histogram(x=df_grouped.rtca ./df_grouped.volume, nbinsx=bins, name="rtca"));
#     h_rtcb = plot(histogram(x=df_grouped.rtcb ./df_grouped.volume, nbinsx=bins, name="rtcb"));
#     h_rtcr = plot(histogram(x=df_grouped.rtcr ./df_grouped.volume, nbinsx=bins, name="rtcr"));
#     h_rd = plot(histogram(x=df_grouped.rd ./df_grouped.volume, nbinsx=bins, name="rd"));
#     h_rt = plot(histogram(x=df_grouped.rt ./df_grouped.volume, nbinsx=bins, name="rt"));
#     return h_rh, h_rma, h_rmb, h_rmr, h_rtca, h_rtcb, h_rtcr, h_rd, h_rt
# end
# h_rh, h_rma, h_rmb, h_rmr, h_rtca, h_rtcb, h_rtcr, h_rd, h_rt = plot_hist(df);
# [h_rma h_rmb h_rmr; h_rtca h_rtcb h_rtcr; h_rh h_rd h_rt] 

# h_rh1, h_rma1, h_rmb1, h_rmr1, h_rtca1, h_rtcb1, h_rtcr1, h_rd1, h_rt1 = plot_hist(df1);
# [h_rma1 h_rmb1 h_rmr1; h_rtca1 h_rtcb1 h_rtcr1; h_rh1 h_rd1 h_rt1] 

# h_rh2, h_rma2, h_rmb2, h_rmr2, h_rtca2, h_rtcb2, h_rtcr2, h_rd2, h_rt2 = plot_hist(df2);
# [h_rma2 h_rmb2 h_rmr2; h_rtca2 h_rtcb2 h_rtcr2; h_rh2 h_rd2 h_rt2] 

# h_rh3, h_rma3, h_rmb3, h_rmr3, h_rtca3, h_rtcb3, h_rtcr3, h_rd3, h_rt3 = plot_hist(df3);
# [h_rma3 h_rmb3 h_rmr3; h_rtca3 h_rtcb3 h_rtcr3; h_rh3 h_rd3 h_rt3] 

# h_rh5, h_rma5, h_rmb5, h_rmr5, h_rtca5, h_rtcb5, h_rtcr5, h_rd5, h_rt5 = plot_hist(df5);
# [h_rma5 h_rmb5 h_rmr5; h_rtca5 h_rtcb5 h_rtcr5; h_rh5 h_rd5 h_rt5] 

# plot([scatter(x=df_r.time, y=repeat([1],length(df_r.time)), mode="markers"),
#       scatter(x=df_r1.time, y=repeat([2],length(df_r1.time)), mode="markers"),
#       scatter(x=df_r2.time, y=repeat([3],length(df_r2.time)), mode="markers"),
#       scatter(x=df_r3.time, y=repeat([4],length(df_r3.time)), mode="markers")])





# plot(scatter(x=df.time, y=df.rm_a./df.volume))
# plot(scatter(x=df.time, y=df.rtca))
# plot(scatter(x=df.time,y=df.rtca./df.volume))
# plot(scatter(x=df.time,y=df.rh ./df.volume))
# plot(scatter(x=df.time,y=df.rt ./df.volume))
# plot(scatter(x=df.time,y=df.rd ./df.volume))


# plot([scatter(x=df.time, y=df.rm_a./df.volume),scattergl(x=df.time, y=df.rh./df.volume)])

# filter(row -> row[:event] == 2, df_r).event
# plot([scattergl(x=df.time,y=df.rm_r),scattergl(x=filter(row -> row[:event] == 2, df_r).time,y=filter(row -> row[:event] == 2, df_r).event, mode="markers")])
# plot([scatter(x=df.time,y=df.rm_a),scattergl(x=filter(row -> row[:event] == 1, df_r).time,y=filter(row -> row[:event] == 1, df_r).event, mode="markers")])
# plot([scatter(x=df.time,y=df.rtca),scattergl(x=filter(row -> row[:event] == 3, df_r).time,y=filter(row -> row[:event] == 3, df_r).event, mode="markers")])



# react_names = [:tscr_ab, :tscr_r, :tlr_a, :tlr_b, :tlr_r, :Vinflux, :Vdam, :Vtag, :Vrep, :deg_rd, :deg_rma, :deg_rmb, :deg_rmr]
# plot(scattergl(x=df_r.time, y=df_r.event, mode="markers", marker_color=df_r.event, marker=attr(colorscale="Viridis")), Layout(yaxis=attr(tickvals=range(1,13), ticktext=react_names)))

# function plotprops(df)
#     props = df.event
#     react_names = [:tscr_ab, :tscr_r, :tlr_a, :tlr_b, :tlr_r, :Vinflux, :Vdam, :Vtag, :Vrep, :deg_rd, :deg_rma, :deg_rmb, :deg_rmr, :V]

#     df_props = DataFrame([name => Float64[] for name in react_names])

#     for i in props
#         push!(df_props, [i[j] for j in 1:length(props[end])])
#     end
#     return plot([scatter(x=df.time, y=df_props[1:end,col], name="$col", mode="markers") for col in names(eachcol(df_props[:,1:end-1]))])

# end
# p_props = plotprops(df_ns);
# p_props1 = plotprops(df_s);
# [p_props p_props1]

# open("/home/hollie_hindley/Documents/stochastic_hybrid/p_props.html", "w") do io
#     PlotlyBase.to_html(io, p_props.plot)
# end


# plot([scatter(x=df.time, y=df.rm_r), scattergl(x=df_p.time,y=df_props[:, :tscr_r], mode="markers")])

# time_taken = @elapsed run_stoch(X0, 0, 0.6, "det")

# df = DataFrame(CSV.File("/home/hollie_hindley/Documents/stochastic_hybrid/det.dat", header=["event", "time", "rm_a", "rtca", "rm_b", "rtcb", "rm_r", "rtcr", "rh", "rd", "rt", "volume"]))

# plot([scatter(x=df_rtc.time, y=df_rtc.rh), scatter(x=df.time, y=df.rh ./df.volume)])
# plot([scatter(x=df_rtc.time, y=df_rtc.rm_a), scatter(x=df.time, y=df.rm_a ./df.volume)])
# plot([scatter(x=df_rtc.time, y=df_rtc.rtca), scatter(x=df.time, y=df.rtca ./df.volume)])
# plot([scatter(x=df_rtc.time, y=df_rtc.rm_b), scatter(x=df.time, y=df.rm_b ./df.volume)])
# plot([scatter(x=df_rtc.time, y=df_rtc.rtcb), scatter(x=df.time, y=df.rtcb ./df.volume)])
# plot([scatter(x=df_rtc.time, y=df_rtc.rm_r), scatter(x=df.time, y=df.rm_r ./df.volume)])
# plot([scatter(x=df_rtc.time, y=df_rtc.rtcr), scatter(x=df.time, y=df.rtcr ./df.volume)])
# plot([scatter(x=df_rtc.time, y=df_rtc.rt), scatter(x=df.time, y=df.rt ./df.volume)])
# plot([scatter(x=df_rtc.time, y=df_rtc.rd), scatter(x=df.time, y=df.rd ./df.volume)])



# plot([scatter(x=df_rtc.time, y=df_rtc.rm_a), scatter(x=df.time, y=df.rm_a)])
# plot([scatter(x=df_rtc.time, y=df_rtc.rtca), scatter(x=df.time, y=df.rtca)])
# plot([scatter(x=df_rtc.time, y=df_rtc.rm_b), scatter(x=df.time, y=df.rm_b)])
# plot([scatter(x=df_rtc.time, y=df_rtc.rtcb), scatter(x=df.time, y=df.rtcb)])
# plot([scatter(x=df_rtc.time, y=df_rtc.rm_r), scatter(x=df.time, y=df.rm_r)])
# plot([scatter(x=df_rtc.time, y=df_rtc.rtcr), scatter(x=df.time, y=df.rtcr)])
# plot([scatter(x=df_rtc.time, y=df_rtc.rh), scatter(x=df.time, y=df.rh)])

