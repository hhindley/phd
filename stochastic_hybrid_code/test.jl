using StatsBase, Distributions, Random, DataFrames, CSV, PlotlyJS, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools

PATH = "/home/hollie_hindley/Documents"

include("$PATH/paper/model_params_funcs_2024/params.jl")
include("$PATH/paper/model_params_funcs_2024/rtc_params_molecs.jl")
include("$PATH/stochastic_hybrid/indexing.jl")
include("$PATH/stochastic_hybrid/hybrid_algo.jl")
include("$PATH/stochastic_hybrid/stoch_model.jl")

# include("/home/hollie_hindley/Documents/stochastic_hybrid/run_rtc_orig.jl")

n= 200 # number of cell cycles
options = Dict(
"threshold"  =>  0.,       # Threshold to decide between determinisitic or stochastic reaction
"FixDetReact"=> [14],# [10,11,12,13,14,15,16,17,18],       # Reactions to be treated determinisitically
    "tspan"     =>   n*log(2)/lam_val,     # Max time for cell cycle
    "samplingFreq"  => 1#0.1  # for sampling every x mins
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

hinit=[3.39e-6/sf,90.7e-6/sf,3.39e-6/sf,75.16e-6/sf,4.4e-6/sf,12.57e-6/sf,0.182/sf,5.11/sf,12.58/sf]
X0 = collect(get_X0(indV, hinit)')

X0

include("$PATH/stochastic_hybrid/hybrid_algo.jl")

# X0 = collect(get_X0(indV, ssvals_rtc_molec)')
time_taken = @elapsed run_stoch(X0, 20, 0.6, "test_nofloor")
df_n = DataFrame(CSV.File("/home/hollie_hindley/Documents/stochastic_hybrid/test_nofloor.dat", header=["event", "time", "rm_a", "rtca", "rm_b", "rtcb", "rm_r", "rtcr", "rh", "rd", "rt", "volume", "totprop"]))

time_taken = @elapsed run_stoch(X0, 20, 0.6, "test_floor")
df = DataFrame(CSV.File("/home/hollie_hindley/Documents/stochastic_hybrid/test_floor.dat", header=["event", "time", "rm_a", "rtca", "rm_b", "rtcb", "rm_r", "rtcr", "rh", "rd", "rt", "volume", "totprop"]))

[plot(scatter(x=df_n.time, y=df_n.totprop)) plot(scatter(x=df_n.time, y=df_n.rm_a)) plot(scatter(x=df_n.time, y=df_n.rtca));
plot(scatter(x=df.time, y=df.totprop)) plot(scatter(x=df.time, y=df.rm_a)) plot(scatter(x=df.time, y=df.rtca));]

[plot(scatter(x=df_n.time, y=df_n.rm_a)) plot(scatter(x=df.time, y=df.rm_a))]

df_n.event = [split(replace(i, r"[\[\]\(Any)]" => ""), ",") for i in df_n.event]
df_n.event = [parse.(Float64, subarray) for subarray in df_n.event]
df_ns = filter(row -> length(row[:event]) > 1, df_n)

df.event = [split(replace(i, r"[\[\]\(Any)]" => ""), ",") for i in df.event]
df.event = [parse.(Float64, subarray) for subarray in df.event]
df_s = filter(row -> length(row[:event]) > 1, df)

function plotprops(df)
    props = df.event
    react_names = [:tscr_ab, :tscr_r, :tlr_a, :tlr_b, :tlr_r, :Vinflux, :Vdam, :Vtag, :Vrep, :deg_rd, :deg_rma, :deg_rmb, :deg_rmr, :V]

    df_props = DataFrame([name => Float64[] for name in react_names])

    for i in props
        push!(df_props, [i[j] for j in 1:length(props[end])])
    end
    return plot([scatter(x=df.time, y=df_props[1:end,col], name="$col", mode="markers") for col in names(eachcol(df_props[:,1:end-1]))])

end
p_props = plotprops(df_ns);
p_props1 = plotprops(df_s);
[plot(scatter(x=df_n.time,y=df_n.rm_a)) plot(scatter(x=df.time,y=df.rm_a))]



df.event = [eval(Meta.parse(i)) for i in df.event]
df_r = filter(row -> length(row[:event]) == 1 && row[:event][1] != 0, df)

df_n.event = [eval(Meta.parse(i)) for i in df_n.event]
df_rn = filter(row -> length(row[:event]) == 1 && row[:event][1] != 0, df_n)

react_names = [:tscr_ab, :tscr_r, :tlr_a, :tlr_b, :tlr_r, :Vinflux, :Vdam, :Vtag, :Vrep, :deg_rd, :deg_rma, :deg_rmb, :deg_rmr]
pr1 = plot(scattergl(x=df_r.time, y=df_r.event, mode="markers", marker_color=df_r.event, marker=attr(colorscale="Viridis")), Layout(yaxis=attr(tickvals=range(1,13), ticktext=react_names)));
pr2 = plot(scattergl(x=df_rn.time, y=df_rn.event, mode="markers", marker_color=df_rn.event, marker=attr(colorscale="Viridis")), Layout(yaxis=attr(tickvals=range(1,13), ticktext=react_names)));
[plot(scatter(x=df_n.time, y=df_n.rm_a)) plot(scatter(x=df.time, y=df.rm_a));pr2 pr1]


df_t = filter(row -> length(row[:event]) == 2, df)

before = first.(df_t.event)
after = last.(df_t.event)
df_t

df.event = [eval(Meta.parse(i)) for i in df.event]
dfa = filter(row -> length(row[:event]) > 1, df)

rma = []; rtca = []; rmb = []; rtcb = []; rmr = []; rtcr = []; rh = []; rd = []; rt = []
for i in dfa.event
    push!(rma, i[1])
    push!(rtca, i[2])
    push!(rmb, i[3])    
    push!(rtcb, i[4])
    push!(rmr, i[5])
    push!(rtcr, i[6])
    push!(rh, i[7])
    push!(rd, i[8])
    push!(rt, i[9])
end

plot([scatter(x=df.time, y=df.rm_a), scatter(x=dfa.time, y=rma)])
plot([scatter(x=df.time, y=df.rtca), scatter(x=dfa.time, y=rtca)])
plot([scatter(x=df.time, y=df.rm_b), scatter(x=dfa.time, y=rmb)])
plot([scatter(x=df.time, y=df.rtcb), scatter(x=dfa.time, y=rtcb)])
plot([scatter(x=df.time, y=df.rm_r), scatter(x=dfa.time, y=rmr)])
plot([scatter(x=df.time, y=df.rtcr), scatter(x=dfa.time, y=rtcr)])
plot([scatter(x=df.time, y=df.rh), scatter(x=dfa.time, y=rh)])
plot([scatter(x=df.time, y=df.rd), scatter(x=dfa.time, y=rd)])
plot([scatter(x=df.time, y=df.rt), scatter(x=dfa.time, y=rt)])


before = first.(df_s.event)
after = last.(df_s.event)
plot([scatter(x=df_t.time, y=before./df_t.volume),scatter(x=df_t.time, y=after./df_t.volume)])
# df_f = DataFrame(CSV.File("/home/hollie_hindley/Documents/stochastic_hybrid/test.dat", header=["event", "time", "rm_a", "rtca", "rm_b", "rtcb", "rm_r", "rtcr", "rh", "rd", "rt", "volume"]))

# df_rfloor = filter(row -> length(row[:event]) == 1 && row[:event][1] != 0, df_floor)
# df_r = filter(row -> length(row[:event]) == 1 && row[:event][1] != 0, df)

[plot(scatter(x=df_n.time,y=df_n.rm_a)) plot(scatter(x=df.time,y=df.rm_a))]
[plot(scatter(x=df_n.time,y=df_n.event)) plot(scatter(x=df.time,y=df.event))]

[plot(scatter(x=df.time,y=df.rh)) plot(scatter(x=df_floor.time,y=df_floor.rh))]

[plot(scatter(x=df_r.time,y=df_r.rm_a)) plot(scatter(x=df_rfloor.time,y=df_rfloor.rm_a))]
[plot(scatter(x=df.time,y=df.rh)) plot(scatter(x=df_floor.time,y=df_floor.rh))]

df_n.event = [split(replace(i, r"[\[\]\(Any)]" => ""), ",") for i in df_n.event]
df_n.event = [parse.(Float64, subarray) for subarray in df_n.event]

df_t = filter(row -> row[:event] == [0], df)
df_v = filter(row -> row[:event] == [20], df)
df_p1 = filter(row -> length(row[:event]) > 1, df_n)
df_r = filter(row -> length(row[:event]) == 1 && row[:event][1] != 0, df)
df_r.event = first.(df_r.event)
# plot([scattergl(x=df.time, y=df[:,col], name="$(names(df)[i])", marker_color=colours[i], legendgroup="$i", showlegend=true) for (col, i) in zip(names(eachcol(df[:,3:end-2])), range(3,length(names(df))-2))])#, title="kdam = $(params_rtc[kdam])"))
# plot([scattergl(x=df.time, y=df[:,col] ./df.volume, name="$(names(df)[i])", marker_color=colours[i], legendgroup="$i", showlegend=true) for (col, i) in zip(names(eachcol(df[:,3:end-2])), range(3,length(names(df))-2))])#, title="kdam = $(params_rtc[kdam])"))
df_p.event[120]
df_p.event = [split(replace(i, r"[\[\]\(LinearAlgebra.Adjoint{Float64})]" => ""), ",") for i in df_p.event]

# df_p.totprop = map(x -> x[1], df_p.event)
# df_p.xi = map(x -> length(x) > 1 ? x[2] : NaN, df_p.event)
# plot([scatter(x=df_p.time, y=df_p.totprop), scatter(x=df_p.time, y=df_p.xi)])

plot(scatter(x=df_p.time, y=df_p.rd))# ./df_p.volume))
mean(df.rtca[10000:end])#./df.volume[10000:end])

X0

c1 = [split(replace(i, r"" => ""), ",") for i in c]
df1.event = [parse.(Float64, subarray) for subarray in df1.event]

df[:,:rm_b]
any(df_grouped[:,:rm_a] .< 0)


df.bins = ceil.(Int, (1:nrow(df))/100)
df_grouped = combine(first, groupby(df, :bins))
bins = 25#length(df_grouped.bins)
h_rh = plot(histogram(x=df_grouped.rh ./df_grouped.volume, nbinsx=bins));
h_rma = plot(histogram(x=df_grouped.rm_a ./df_grouped.volume, nbinsx=bins));
h_rmb = plot(histogram(x=df_grouped.rm_b ./df_grouped.volume, nbinsx=bins));
h_rmr = plot(histogram(x=df_grouped.rm_r ./df_grouped.volume, nbinsx=bins));
h_rtca = plot(histogram(x=df_grouped.rtca ./df_grouped.volume, nbinsx=bins));
h_rtcb = plot(histogram(x=df_grouped.rtcb ./df_grouped.volume, nbinsx=bins));
h_rtcr = plot(histogram(x=df_grouped.rtcr ./df_grouped.volume, nbinsx=bins));
h_rd = plot(histogram(x=df_grouped.rd ./df_grouped.volume, nbinsx=bins));
h_rt = plot(histogram(x=df_grouped.rt ./df_grouped.volume, nbinsx=bins));


plot(histogram(x=df_grouped.rh, nbinsx=bins))
plot(histogram(x=df_grouped.rm_a, nbinsx=bins))
plot(histogram(x=df_grouped.rm_b, nbinsx=bins))
plot(histogram(x=df_grouped.rm_r, nbinsx=bins))
plot(histogram(x=df_grouped.rtca, nbinsx=bins))
plot(histogram(x=df_grouped.rtcb, nbinsx=bins))
plot(histogram(x=df_grouped.rtcr, nbinsx=bins))
plot(histogram(x=df_grouped.rd, nbinsx=bins))
plot(histogram(x=df_grouped.rt, nbinsx=bins))

mean(df.rm_a[10000:end])


time_taken = @elapsed run_stoch(X0, 10, 0.6, "thresh10_kdam06")
df = DataFrame(CSV.File("/home/hollie_hindley/Documents/stochastic_hybrid/X0_init/thresh10_kdam01.dat", header=["event", "time", "rm_a", "rtca", "rm_b", "rtcb", "rm_r", "rtcr", "rh", "rd", "rt", "volume"]))

time_taken = @elapsed run_stoch(X0, 15, 0.6, "thresh15_kdam06")
df1 = DataFrame(CSV.File("/home/hollie_hindley/Documents/stochastic_hybrid/X0_init/thresh15_kdam01.dat", header=["event", "time", "rm_a", "rtca", "rm_b", "rtcb", "rm_r", "rtcr", "rh", "rd", "rt", "volume"]))

time_taken = @elapsed run_stoch(X0, 20, 0.6, "thresh20_kdam06")
df2 = DataFrame(CSV.File("/home/hollie_hindley/Documents/stochastic_hybrid/X0_init/thresh20_kdam01.dat", header=["event", "time", "rm_a", "rtca", "rm_b", "rtcb", "rm_r", "rtcr", "rh", "rd", "rt", "volume"]))

time_taken = @elapsed run_stoch(X0, 100, 0.6, "thresh100_kdam06")
df3 = DataFrame(CSV.File("/home/hollie_hindley/Documents/stochastic_hybrid/X0_init/thresh50_kdam01.dat", header=["event", "time", "rm_a", "rtca", "rm_b", "rtcb", "rm_r", "rtcr", "rh", "rd", "rt", "volume"]))

s = :rm_a
p1 = plot(scatter(x=df.time,y=df[:,s] ./df.volume))
p2 = plot(scatter(x=df1.time,y=df1[:,s] ./df1.volume));
p3 = plot(scatter(x=df2.time,y=df2[:,s] ./df2.volume))
p4 = plot(scatter(x=df3.time,y=df3[:,s] ./df3.volume))

p5 = plot(scatter(x=df5.time,y=df5[:,s] ./df5.volume));
p5
[p1 p2 p3 p4]

function df_sort(df)
    df.event = [split(replace(i, r"[\[\]\(Any)]" => ""), ",") for i in df.event]
    df.event = [parse.(Float64, subarray) for subarray in df.event]
    df_p = filter(row -> length(row[:event]) > 1, df)
    df_r = filter(row -> length(row[:event]) == 1 && row[:event][1] != 0, df)
    df_r.event = first.(df_r.event)

    props = df_p.event
    react_names = [:tscr_ab, :tscr_r, :tlr_a, :tlr_b, :tlr_r, :Vinflux, :Vdam, :Vtag, :Vrep, :deg_rd, :deg_rma, :deg_rmb, :deg_rmr, :V]

    df_props = DataFrame([name => Float64[] for name in react_names])

    for i in props
        push!(df_props, [i[j] for j in 1:length(props[end])])
    end
    return df_props, df_r, df_p
end

df_props, df_r, df_p = df_sort(df)
df_props1, df_r1, df_p1 = df_sort(df1)
df_props2, df_r2, df_p2 = df_sort(df2)
df_props3, df_r3, df_p3 = df_sort(df3)

react_names = [:tscr_ab, :tscr_r, :tlr_a, :tlr_b, :tlr_r, :Vinflux, :Vdam, :Vtag, :Vrep, :deg_rd, :deg_rma, :deg_rmb, :deg_rmr]
pr1 = plot(scattergl(x=df_r.time, y=df_r.event, mode="markers", marker_color=df_r.event, marker=attr(colorscale="Viridis")), Layout(yaxis=attr(tickvals=range(1,13), ticktext=react_names)));
pr2 = plot(scattergl(x=df_r1.time, y=df_r1.event, mode="markers", marker_color=df_r1.event, marker=attr(colorscale="Viridis")), Layout(yaxis=attr(tickvals=range(1,13), ticktext=react_names)));
pr3 = plot(scattergl(x=df_r2.time, y=df_r2.event, mode="markers", marker_color=df_r2.event, marker=attr(colorscale="Viridis")), Layout(yaxis=attr(tickvals=range(1,13), ticktext=react_names)));
pr4 = plot(scattergl(x=df_r3.time, y=df_r3.event, mode="markers", marker_color=df_r3.event, marker=attr(colorscale="Viridis")), Layout(yaxis=attr(tickvals=range(1,13), ticktext=react_names)));

[pr1 pr2 pr3 pr4]
length(df_props.tlr_a),length(df_props1.tlr_a),length(df_props2.tlr_a),length(df_props3.tlr_a)

p_props = plot([scatter(x=df_p.time, y=df_props[1:100:end,col], name="$col", mode="markers") for col in names(eachcol(df_props[:,1:end-1]))])
p_props1 = plot([scatter(x=df_p1.time, y=df_props1[1:100:end,col], name="$col", mode="markers") for col in names(eachcol(df_props1[:,1:end-1]))]);
p_props2 = plot([scatter(x=df_p2.time, y=df_props2[1:100:end,col], name="$col", mode="markers") for col in names(eachcol(df_props2[:,1:end-1]))]);
p_props3 = plot([scatter(x=df_p3.time, y=df_props3[1:100:end,col], name="$col", mode="markers") for col in names(eachcol(df_props3[:,1:end-1]))]);

[p_props p_props1 p_props2 p_props3]

function plot_hist(df)
    df.bins = ceil.(Int, (1:nrow(df))/100)
    df_grouped = combine(first, groupby(df, :bins))
    bins = 59#length(df_grouped.bins)
    h_rh = plot(histogram(x=df_grouped.rh ./df_grouped.volume, nbinsx=bins, name="rh"));
    h_rma = plot(histogram(x=df_grouped.rm_a ./df_grouped.volume, nbinsx=bins, name="rm_a"));
    h_rmb = plot(histogram(x=df_grouped.rm_b ./df_grouped.volume, nbinsx=bins, name="rm_b"));
    h_rmr = plot(histogram(x=df_grouped.rm_r ./df_grouped.volume, nbinsx=bins, name="rm_r"));
    h_rtca = plot(histogram(x=df_grouped.rtca ./df_grouped.volume, nbinsx=bins, name="rtca"));
    h_rtcb = plot(histogram(x=df_grouped.rtcb ./df_grouped.volume, nbinsx=bins, name="rtcb"));
    h_rtcr = plot(histogram(x=df_grouped.rtcr ./df_grouped.volume, nbinsx=bins, name="rtcr"));
    h_rd = plot(histogram(x=df_grouped.rd ./df_grouped.volume, nbinsx=bins, name="rd"));
    h_rt = plot(histogram(x=df_grouped.rt ./df_grouped.volume, nbinsx=bins, name="rt"));
    return h_rh, h_rma, h_rmb, h_rmr, h_rtca, h_rtcb, h_rtcr, h_rd, h_rt
end
h_rh, h_rma, h_rmb, h_rmr, h_rtca, h_rtcb, h_rtcr, h_rd, h_rt = plot_hist(df);
[h_rma h_rmb h_rmr; h_rtca h_rtcb h_rtcr; h_rh h_rd h_rt] 

h_rh1, h_rma1, h_rmb1, h_rmr1, h_rtca1, h_rtcb1, h_rtcr1, h_rd1, h_rt1 = plot_hist(df1);
[h_rma1 h_rmb1 h_rmr1; h_rtca1 h_rtcb1 h_rtcr1; h_rh1 h_rd1 h_rt1] 

h_rh2, h_rma2, h_rmb2, h_rmr2, h_rtca2, h_rtcb2, h_rtcr2, h_rd2, h_rt2 = plot_hist(df2);
[h_rma2 h_rmb2 h_rmr2; h_rtca2 h_rtcb2 h_rtcr2; h_rh2 h_rd2 h_rt2] 

h_rh3, h_rma3, h_rmb3, h_rmr3, h_rtca3, h_rtcb3, h_rtcr3, h_rd3, h_rt3 = plot_hist(df3);
[h_rma3 h_rmb3 h_rmr3; h_rtca3 h_rtcb3 h_rtcr3; h_rh3 h_rd3 h_rt3] 

h_rh5, h_rma5, h_rmb5, h_rmr5, h_rtca5, h_rtcb5, h_rtcr5, h_rd5, h_rt5 = plot_hist(df5);
[h_rma5 h_rmb5 h_rmr5; h_rtca5 h_rtcb5 h_rtcr5; h_rh5 h_rd5 h_rt5] 

plot([scatter(x=df_r.time, y=repeat([1],length(df_r.time)), mode="markers"),
      scatter(x=df_r1.time, y=repeat([2],length(df_r1.time)), mode="markers"),
      scatter(x=df_r2.time, y=repeat([3],length(df_r2.time)), mode="markers"),
      scatter(x=df_r3.time, y=repeat([4],length(df_r3.time)), mode="markers")])





plot(scatter(x=df.time, y=df.rm_a./df.volume))
plot(scatter(x=df.time, y=df.rtca))
plot(scatter(x=df.time,y=df.rtca./df.volume))
plot(scatter(x=df.time,y=df.rh ./df.volume))
plot(scatter(x=df.time,y=df.rt ./df.volume))
plot(scatter(x=df.time,y=df.rd ./df.volume))


plot([scatter(x=df.time, y=df.rm_a./df.volume),scattergl(x=df.time, y=df.rh./df.volume)])

filter(row -> row[:event] == 2, df_r).event
plot([scattergl(x=df.time,y=df.rm_r),scattergl(x=filter(row -> row[:event] == 2, df_r).time,y=filter(row -> row[:event] == 2, df_r).event, mode="markers")])
plot([scatter(x=df.time,y=df.rm_a),scattergl(x=filter(row -> row[:event] == 1, df_r).time,y=filter(row -> row[:event] == 1, df_r).event, mode="markers")])
plot([scatter(x=df.time,y=df.rtca),scattergl(x=filter(row -> row[:event] == 3, df_r).time,y=filter(row -> row[:event] == 3, df_r).event, mode="markers")])



react_names = [:tscr_ab, :tscr_r, :tlr_a, :tlr_b, :tlr_r, :Vinflux, :Vdam, :Vtag, :Vrep, :deg_rd, :deg_rma, :deg_rmb, :deg_rmr]
plot(scattergl(x=df_r.time, y=df_r.event, mode="markers", marker_color=df_r.event, marker=attr(colorscale="Viridis")), Layout(yaxis=attr(tickvals=range(1,13), ticktext=react_names)))

function plotprops(df)
    props = df.event
    react_names = [:tscr_ab, :tscr_r, :tlr_a, :tlr_b, :tlr_r, :Vinflux, :Vdam, :Vtag, :Vrep, :deg_rd, :deg_rma, :deg_rmb, :deg_rmr, :V]

    df_props = DataFrame([name => Float64[] for name in react_names])

    for i in props
        push!(df_props, [i[j] for j in 1:length(props[end])])
    end
    return plot([scatter(x=df.time, y=df_props[1:end,col], name="$col", mode="markers") for col in names(eachcol(df_props[:,1:end-1]))])

end
p_props = plotprops(df_ns);
p_props1 = plotprops(df_s);
[p_props p_props1]

open("/home/hollie_hindley/Documents/stochastic_hybrid/p_props.html", "w") do io
    PlotlyBase.to_html(io, p_props.plot)
end


plot([scatter(x=df.time, y=df.rm_r), scattergl(x=df_p.time,y=df_props[:, :tscr_r], mode="markers")])

time_taken = @elapsed run_stoch(X0, 0, 0.6, "det")

df = DataFrame(CSV.File("/home/hollie_hindley/Documents/stochastic_hybrid/det.dat", header=["event", "time", "rm_a", "rtca", "rm_b", "rtcb", "rm_r", "rtcr", "rh", "rd", "rt", "volume"]))

plot([scatter(x=df_rtc.time, y=df_rtc.rh), scatter(x=df.time, y=df.rh ./df.volume)])
plot([scatter(x=df_rtc.time, y=df_rtc.rm_a), scatter(x=df.time, y=df.rm_a ./df.volume)])
plot([scatter(x=df_rtc.time, y=df_rtc.rtca), scatter(x=df.time, y=df.rtca ./df.volume)])
plot([scatter(x=df_rtc.time, y=df_rtc.rm_b), scatter(x=df.time, y=df.rm_b ./df.volume)])
plot([scatter(x=df_rtc.time, y=df_rtc.rtcb), scatter(x=df.time, y=df.rtcb ./df.volume)])
plot([scatter(x=df_rtc.time, y=df_rtc.rm_r), scatter(x=df.time, y=df.rm_r ./df.volume)])
plot([scatter(x=df_rtc.time, y=df_rtc.rtcr), scatter(x=df.time, y=df.rtcr ./df.volume)])
plot([scatter(x=df_rtc.time, y=df_rtc.rt), scatter(x=df.time, y=df.rt ./df.volume)])
plot([scatter(x=df_rtc.time, y=df_rtc.rd), scatter(x=df.time, y=df.rd ./df.volume)])



plot([scatter(x=df_rtc.time, y=df_rtc.rm_a), scatter(x=df.time, y=df.rm_a)])
plot([scatter(x=df_rtc.time, y=df_rtc.rtca), scatter(x=df.time, y=df.rtca)])
plot([scatter(x=df_rtc.time, y=df_rtc.rm_b), scatter(x=df.time, y=df.rm_b)])
plot([scatter(x=df_rtc.time, y=df_rtc.rtcb), scatter(x=df.time, y=df.rtcb)])
plot([scatter(x=df_rtc.time, y=df_rtc.rm_r), scatter(x=df.time, y=df.rm_r)])
plot([scatter(x=df_rtc.time, y=df_rtc.rtcr), scatter(x=df.time, y=df.rtcr)])
plot([scatter(x=df_rtc.time, y=df_rtc.rh), scatter(x=df.time, y=df.rh)])


plot(scatter(x=df.time, y=df.volume))
# sf1 = @. (1e6/(6.022e23*(df.volume*1e-15)))
# atp1 = @. (par[pidx(:atp)]/sf1)#/df.volume
# p=plot(scatter(x=df.time,y=sf1))
# p1=plot(scatter(x=df.time,y=atp))
# [p; p1]

# plot([scatter(x=df_rtc.time, y=repeat([atp_val_molec], length(df_rtc.time))), scatter(x=df.time, y=repeat([mean(atp ./df.volume)],length(df.time)))])

# v_sf = @. 1/df.volume

# alpha1, fa1, ra1 = all_vars(df)
# Voc1 = @. par[pidx(:Vmax_init)]*atp1/((par[pidx(:Km_init)]/v_sf)+atp1) # uM min-1
# sig_o1 = @. ra1*Voc1/par[pidx(:k_diss)] 

# tscr_ab1 = @. sig_o1*par[pidx(:ω_ab)]*atp1/(par[pidx(:θtscr)]/v_sf+atp1) # uM min-1
# tscr_r1 = @. par[pidx(:ω_r)]*atp1/(par[pidx(:θtscr)]+atp1)


# tlr_el1 = @. (par[pidx(:g_max)]*atp1/((par[pidx(:θtlr)]/v_sf)+atp1))/df.volume

# tlr_a1 = @. (1/par[pidx(:na)])*par[pidx(:kc)]*df.rh*df.rm_a*tlr_el1
# tlr_b1 = @. (1/par[pidx(:nb)])*par[pidx(:kc)]*df.rh*df.rm_b*tlr_el1
# tlr_r1 = @. (1/par[pidx(:nr)])*par[pidx(:kc)]*df.rh*df.rm_r*tlr_el1

# Vrep1 = @. df.rtcb*df.rt*par[pidx(:krep)]/(df.rt+par[pidx(:km_b)]/v_sf) # uM min-1 
# Vtag1 = @. df.rtca*df.rd*par[pidx(:ktag)]/(df.rd+par[pidx(:km_a)]/v_sf) # uM min-1 

# Vdam1 = @. par[pidx(:kdam)]*df.rh # uM min-1

# Vinflux = @. par[pidx(:kin)] * par[pidx(:g_max)]*atp1/(par[pidx(:θtlr)]/v_sf+atp1) # uM min-1 


# plot([scatter(x=df_rtc.time, y=alpha), scatter(x=df.time, y=alpha1 ./df.volume)])
# plot([scatter(x=df_rtc.time, y=fa), scatter(x=df.time, y=fa1 ./df.volume)])
# plot([scatter(x=df_rtc.time, y=ra), scatter(x=df.time, y=ra1 ./df.volume)])
# plot([scatter(x=df_rtc.time, y=repeat([Voc], length(df_rtc.time))), scatter(x=df.time, y=Voc1 ./df.volume)])
# plot([scatter(x=df_rtc.time, y=sig_o), scatter(x=df.time, y=sig_o1 ./df.volume)])
# plot([scatter(x=df_rtc.time, y=tscr_ab), scatter(x=df.time, y=tscr_ab1 ./df.volume)])
# plot([scatter(x=df_rtc.time, y=tlr_a), scatter(x=df.time, y=tlr_a1 ./df.volume)])
# plot([scatter(x=df_rtc.time, y=tlr_b), scatter(x=df.time, y=tlr_b1 ./df.volume)])
# plot([scatter(x=df_rtc.time, y=tlr_r), scatter(x=df.time, y=tlr_r1 ./df.volume)])
# plot([scatter(x=df_rtc.time, y=Vrep), scatter(x=df.time, y=Vrep1 ./df.volume)])
# plot([scatter(x=df_rtc.time, y=Vtag), scatter(x=df.time, y=Vtag1 ./df.volume)])
# plot([scatter(x=df_rtc.time, y=Vdam), scatter(x=df.time, y=Vdam1 ./df.volume)])

# plot(scatter(x=df.time, y=df.volume))

# par[pidx(:Vmax_init)]

# params_rtc_molec[Vmax_init]

# par[pidx(:Km_init)]
# params_rtc_molec[Km_init]

# params_rtc_molec[atp]
# atp_val_molec
# atp1[1]


# params_rtc_molec[Vmax_init]*atp_val_molec/(params_rtc_molec[Km_init]+atp_val_molec) # uM min-1 
# Voc = par[pidx(:Vmax_init)]*atp_val_molec/(par[pidx(:Km_init)]+atp_val_molec) # uM min-1 

# par[pidx(:Vmax_init)]*atp1[1]/(par[pidx(:Km_init)]+atp1[1])






# # df.group = ceil.(Int, (1:nrow(df)) / 20)
# # df_new = df[:, Not(:event)]
# # df_avg = combine(groupby(df_new, :group), names(df_new, Not(:group)) .=> mean)
# # plot([scattergl(x=df_avg.time_mean, y=df_avg[:,col], name="$(names(df_avg)[i])", legendgroup="$i", showlegend=true) for (col, i) in zip(names(eachcol(df_avg[:,4:end-2])), range(4,length(names(df_avg))-2))])#, title="kdam = $(params_rtc[kdam])"))

# # avg = mean.(eachcol(df))

# df_v = filter(row -> row[:event] == 2, df)
# df_s = filter(row -> length(row[:event]) == 1, df)

# plot(scatter(x=df.time, y=df.volume))


# df_p = filter(row -> length(row[:event]) > 2, df)

# plot(scatter(x=df.time, y=df.volume))
# plot([scattergl(x=df.time, y=df[:,col], name="$(names(df)[i])", marker_color=colours[i], legendgroup="$i", showlegend=true) for (col, i) in zip(names(eachcol(df[:,3:end-2])), range(3,length(names(df))-2))])#, title="kdam = $(params_rtc[kdam])"))


# plot([scattergl(x=df_v.time, y=df_v[:,col], name="$(names(df_v)[i])", marker_color=colours[i], legendgroup="$i", showlegend=true) for (col, i) in zip(names(eachcol(df_v[:,3:end-2])), range(3,length(names(df_v))-2))])#, title="kdam = $(params_rtc[kdam])"))

# plot([scatter(x=df.time, y=df[:,col]) for col in names(eachcol(df[:,3:end]))])



# props = [split(replace(i, r"[\[\]\(Any)]" => ""), ",") for i in df_p.event]
# props = [parse.(Float64, subarray) for subarray in props]
# props = df_p.event
# react_names = [:tscr_ab, :tscr_r, :tlr_a, :tlr_b, :tlr_r, :Vinflux, :Vdam, :Vtag, :Vrep, :deg_rd, :deg_rma, :deg_rmb, :deg_rmr, :V]

# df_props = DataFrame([name => Float64[] for name in react_names])

# for i in props
#     push!(df_props, [i[j] for j in 1:length(props[end])])
# end

# plot([scattergl(x=df_p.time, y=df_props[:,col], name="$col", mode="markers") for col in names(eachcol(df_props))])
