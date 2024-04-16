using StatsBase, Distributions, Random, DataFrames, CSV, PlotlyJS, DifferentialEquations, OrderedCollections, ProgressBars

PATH = "/home/holliehindley/phd"

include("$PATH/paper/model_params_funcs_2024/params.jl")
include("$PATH/paper/model_params_funcs_2024/rtc_params_molecs.jl")
include("$PATH/stochastic_hybrid/indexing.jl")
include("$PATH/stochastic_hybrid/hybrid_algo.jl")
include("$PATH/stochastic_hybrid/stoch_model.jl")
include("$PATH/stochastic_hybrid/indexing.jl")

options = Dict(
"threshold"  =>  0.,       # Threshold to decide between determinisitic or stochastic reaction
"FixDetReact"=>  [],       # Reactions to be treated determinisitically
    "tspan"     =>   1e6,     # Max time for cell cycle
    "samplingFreq"  => 10     # for sampling every x mins
)


X0 = collect(get_X0(indV)')
par = collect(get_par(indP)')


getssX0 = false

if getssX0
    fout=open("$PATH/stochastic_hybrid/all_X0.dat","w")
    propen, S, propList = defineStochModel(par, indV)
    nx = indV.nrOfItems-1
    prop(X) = propen(X[1:nx])
    X0 = hybrid_algo(X0, options, prop, S, out=fout)
    CSV.write("$PATH/stochastic_hybrid/X0.dat", DataFrame(X0,:auto), header=false)
else
    X0 = CSV.read("$PATH/stochastic_hybrid/X0.dat", Tables.matrix, header=false)
end

 
function prop(X)
    nx = indV.nrOfItems - 1
    propen(X[1:nx]) 
end
function run_stoch(X0, thresh, kdam)
    par[pidx(:kdam)] = kdam
    threshold = thresh # set to zero for a deterministic result
    options["threshold"] = threshold
    fout=open("$PATH/stochastic_hybrid/test1.dat","w")
    global propen, S, propList  = defineStochModel(par, indV)    
    solu = hybrid_algo(X0, options, prop, S, out=fout)
    close(fout)
    df = DataFrame(CSV.File("$PATH/stochastic_hybrid/test1.dat", header=["event", "time", "rm_a", "rtca", "rm_b", "rtcb", "rm_r", "rtcr", "rh", "rd", "rt", "TotProp"]))
return df
end
# plot(scatter(x=df.time, y=df.event), Layout(xaxis_type="log"))

options["threshold"] = 0
df = @time run_stoch(X0, options["threshold"], 0.1)

high_thresh = 5
df3 = @time run_stoch(X0, high_thresh, 0.1)

low_thresh = 0.01
df4 = @time run_stoch(X0, low_thresh, 0.1)


p_det = plot([scattergl(x=df.time, y=df[:,col], name="$(names(df)[i])", marker_color=colours[i], legendgroup="$i", showlegend=true) for (col, i) in zip(names(eachcol(df[:,3:end-1])), range(3,length(names(df))-1))], Layout(yaxis_tickformat=".2e", title="flooring stoch"))#, title="kdam = $(params_rtc[kdam])"))
p_species = plot([scattergl(x=df3.time, y=df3[:,col], name="$(names(df3)[i])", marker_color=colours[i], legendgroup="$i", showlegend=true) for (col, i) in zip(names(eachcol(df3[:,3:end-1])), range(3,length(names(df3))-1))], Layout(yaxis_tickformat=".2e", title="high threshold"))#, title="kdam = $(params_rtc[kdam])"))
p_species1 = plot([scattergl(x=df4.time, y=df4[:,col], name="$(names(df4)[i])", marker_color=colours[i], legendgroup="$i", showlegend=true) for (col, i) in zip(names(eachcol(df4[:,3:end-1])), range(3,length(names(df4))-1))], Layout(yaxis_tickformat=".2e", title="low threshold"))#, title="kdam = $(params_rtc[kdam])"))


open("/home/hollie_hindley/Documents/stochastic_hybrid/higher_threshold/result.html", "w") do io
    PlotlyBase.to_html(io, p_species.plot)
end
open("/home/hollie_hindley/Documents/stochastic_hybrid/low_threshold/result.html", "w") do io
    PlotlyBase.to_html(io, p_species1.plot)
end
react_names = [:tscr_a, :tscr_b, :tscr_r, :tlr_a, :tlr_b, :tlr_r, :Vinflux, :Vdam, :Vtag, :Vrep, :dil_rh, :dil_rd, :dil_rt, :dil_rma, :dil_rmb, :dil_rmr, :dil_rtca, :dil_rtcb, :dil_rtcr, :deg_rd, :deg_rma, :deg_rmb, :deg_rmr]

function df_rearrange(df3)
    df3.event = [split(replace(i, r"[\[\]\(Any)]" => ""), ",") for i in df3.event]
    df3.event = [parse.(Float64, subarray) for subarray in df3.event]

    df_stoch = filter(row -> length(row.event) == 2, df3)
    df_props = filter(row -> length(row.event) > 2, df3)
    df_stoch.stochcum, df_stoch.stochreact = broadcast(x -> x[1], df_stoch.event), broadcast(x -> x[2], df_stoch.event)

    select!(df_props, Not(:TotProp))
    df_props.totprop = map(x -> x[1], df_props.event)
    df_props.xi = map(x -> x[2], df_props.event)
    map!(x -> x[3:end], df_props.event, df_props.event)

    df_react = DataFrame([name => Float64[] for name in react_names])
    for vec in df_props.event
        push!(df_react, vec)
    end
    df_react_red = select(df_react, Not([:tscr_r, :Vinflux]))
    return df_stoch, df_props, df_react_red
end 

df_stoch_high, df_props_high, df_react_high = df_rearrange(df3)
df_stoch_low, df_props_low, df_react_low = df_rearrange(df4)

p_species = plot([scattergl(x=df_props_high.time, y=df_props_high[:,col], name="$(names(df_props_high)[i])", marker_color=colours[i], legendgroup="$i", showlegend=true) for (col, i) in zip(names(eachcol(df_props_high[:,3:end-2])), range(3,length(names(df_props_high))-2))], Layout(yaxis_tickformat=".2e", title="high threshold"))#, title="kdam = $(params_rtc[kdam])"))
p_species1 = plot([scattergl(x=df_props_low.time, y=df_props_low[:,col], name="$(names(df_props_low)[i])", marker_color=colours[i], legendgroup="$i", showlegend=true) for (col, i) in zip(names(eachcol(df_props_low[:,3:end-2])), range(3,length(names(df_props_low))-2))], Layout(yaxis_tickformat=".2e", title="low threshold"))#, title="kdam = $(params_rtc[kdam])"))


open("/home/hollie_hindley/Documents/stochastic_hybrid/higher_threshold/result.html", "w") do io
    PlotlyBase.to_html(io, p_species.plot)
end
open("/home/hollie_hindley/Documents/stochastic_hybrid/low_threshold/result.html", "w") do io
    PlotlyBase.to_html(io, p_species1.plot)
end

totprop_xi_high = plot([scattergl(x=df_props_high.time, y=df_props_high.totprop,name="totprop"), scattergl(x=df_props_high.time, y=df_props_high.xi, name="xi")], Layout(xaxis_title="time", yaxis_title="propensity&xi",title="high stoch"));
totprop_xi_low = plot([scattergl(x=df_props_low.time, y=df_props_low.totprop,name="totprop"), scattergl(x=df_props_low.time, y=df_props_low.xi, name="xi")], Layout(xaxis_title="time", yaxis_title="propensity&xi",title="low stoch"));

open("/home/hollie_hindley/Documents/stochastic_hybrid/high_threshold/totprop_xi_high.html", "w") do io
    PlotlyBase.to_html(io, totprop_xi_high.plot)
end
open("/home/hollie_hindley/Documents/stochastic_hybrid/low_threshold/totprop_xi_low.html", "w") do io
    PlotlyBase.to_html(io, totprop_xi_low.plot)
end



p_stochcum_high = plot(scattergl(x=df_stoch_high.time, y=df_stoch_high.stochcum), Layout(xaxis_title="time", yaxis_title="cumulative stochastic events", title="high stoch"));
p_stochcum_low = plot(scattergl(x=df_stoch_low.time, y=df_stoch_low.stochcum), Layout(xaxis_title="time", yaxis_title="cumulative stochastic events", title="low stoch"));

open("/home/hollie_hindley/Documents/stochastic_hybrid/higher_threshold/p_stochcum_high.html", "w") do io
    PlotlyBase.to_html(io, p_stochcum_high.plot)
end
open("/home/hollie_hindley/Documents/stochastic_hybrid/low_threshold/p_stochcum_low.html", "w") do io
    PlotlyBase.to_html(io, p_stochcum_low.plot)
end



p_stochreact_high = plot(scattergl(x=df_stoch_high.time, y=df_stoch_high.stochreact, mode="markers", marker_color=df_stoch_high.stochreact, marker=attr(colorscale="Viridis")), Layout(xaxis_title="time",yaxis=attr(tickvals=range(1,23), ticktext=react_names, title="high stoch")))
p_stochreact_low = plot(scattergl(x=df_stoch_low.time, y=df_stoch_low.stochreact, mode="markers", marker_color=df_stoch_low.stochreact, marker=attr(colorscale="Viridis")), Layout(xaxis_title="time",yaxis=attr(tickvals=range(1,23), ticktext=react_names, title="low stoch"));)

open("/home/hollie_hindley/Documents/stochastic_hybrid/higher_threshold/p_stochreact_high.html", "w") do io
    PlotlyBase.to_html(io, p_stochreact_high.plot)
end
open("/home/hollie_hindley/Documents/stochastic_hybrid/low_threshold/p_stochreact_low.html", "w") do io
    PlotlyBase.to_html(io, p_stochreact_low.plot)
end



p_props_high = plot([scattergl(x=df_props_high.time, y=df_react_high[:,col], mode="markers", name="$(names(df_react_high)[i])", legendgroup="$i", showlegend=true) for (col, i) in zip(names(eachcol(df_react_high)), range(1,length(names(df_react_high))))], Layout(xaxis_title="reaction propensities", yaxis_title="time", yaxis_tickformat=".2e", title="Reaction propensities - high stoch", yaxis_type="log"));#, title="kdam = $(params_rtc[kdam])"))
add_trace!(p_props_high, scattergl(x=[df_props_high.time[2],df_props_high.time[end]], y=[high_thresh,high_thresh], mode="lines", marker_color="black", name="threshold", marker_width=4))
p_props_high

p_props_low = plot([scattergl(x=df_props_low.time, y=df_react_low[:,col], mode="markers", name="$(names(df_react_low)[i])", legendgroup="$i", showlegend=true) for (col, i) in zip(names(eachcol(df_react_low)), range(1,length(names(df_react_low))))], Layout(xaxis_title="reaction propensities", yaxis_title="time", yaxis_tickformat=".2e", title="Reaction propensities - low stoch", yaxis_type="log"));#, title="kdam = $(params_rtc[kdam])"))
add_trace!(p_props_low, scattergl(x=[df_props_low.time[2],df_props_low.time[end]], y=[low_thresh,low_thresh], mode="lines", marker_color="black", name="threshold", marker_width=4))

open("/home/hollie_hindley/Documents/stochastic_hybrid/higher_threshold/p_props_high.html", "w") do io
    PlotlyBase.to_html(io, p_props_high.plot)
end
open("/home/hollie_hindley/Documents/stochastic_hybrid/low_threshold/p_props_low.html", "w") do io
    PlotlyBase.to_html(io, p_props_low.plot)
end




# looking at xi and totProp

df3.totprop, df3.xi = map(x -> x[1], df3.event), map(x -> x[end], df3.event)
plot([scatter(x=df3.time, y=df3.totprop, name="totprop"), scatter(x=df3.time, y=df3.xi, name="xi")], Layout(xaxis_type="log"))


# looking at saving interval 

freq(t) = t <= 2000 ? t/100 : log(1 + ((t-2000)^40))
plot(scatter(x=df3.time, y=freq.(df3.time)), Layout(xaxis_type="log"))


open("/home/hollie_hindley/Documents/stochastic_hybrid/p_props.html", "w") do io
    PlotlyBase.to_html(io, p_props.plot)
end



prop_arr1 = df.event
values1 = [split(replace(i, r"[\[\]]" => ""), ",") for i in prop_arr1]
arr1 = [parse.(Float64, i) for i in values1]

react_names = [:tscr_a, :tscr_b, :tscr_r, :tlr_a, :tlr_b, :tlr_r, :Vinflux, :Vdam, :Vtag, :Vrep, :dil_rh, :dil_rd, :dil_rt, :dil_rma, :dil_rmb, :dil_rmr, :dil_rtca, :dil_rtcb, :dil_rtcr, :deg_rd, :deg_rma, :deg_rmb, :deg_rmr]
df_react1 = DataFrame([name => Float64[] for name in react_names])

for vec in arr1[2:end]
    push!(df_react1, vec)
end

p2 = plot([scattergl(x=df.time, y=df_react1[:,col], name="$(names(df_react1)[i])", legendgroup="$i", showlegend=true) for (col, i) in zip(names(eachcol(df_react1)), range(1,length(names(df_react1))))], Layout(xaxis_type="log", yaxis_tickformat=".2e", title="Reaction propensities"));#, title="kdam = $(params_rtc[kdam])"))
add_trace!(p2, scatter(x=[df.time[2],df.time[end]], y=[0,0], mode="lines", marker_color="black", name="threshold"))
p2



plot([scatter(x=df3.time, y=df_react.tscr_a)], Layout(xaxis_type="log"))



p_rtc1 = plot([scattergl(x=df.time, y=col, name="$(names(df)[i])", legendgroup="$i") for (col, i) in zip(eachcol(df[:,3:end-1]), range(3,length(names(df))-1))], Layout(xaxis_type="log", yaxis_tickformat=".2e"))#, title="kdam = $(params_rtc[kdam])"))

colours =["#636EFA", "#EF553B", "#00CC96", "#AB63FA", "#FFA15A", "#19D3F3", "#FF6692", "#B6E880", "#FF97FF", "#FECB52", :blue]

p = plot([scattergl(x=df[df.event .== 0, :time], y=df[df.event .== 0, col], name="$(names(df)[i])", marker_color=colours[i], legendgroup="$i") for (col, i) in zip(names(eachcol(df[:,3:end-1])), range(3,length(names(df))-1))], Layout(xaxis_type="log", yaxis_tickformat=".2e", title="deterministic"))#, title="kdam = $(params_rtc[kdam])"))
p1 = plot([scattergl(x=df3[df3.event .== 0, :time], y=df3[df3.event .== 0, col], name="$(names(df3)[i])", marker_color=colours[i], legendgroup="$i", showlegend=false) for (col, i) in zip(names(eachcol(df3[:,3:end-1])), range(3,length(names(df3))-1))], Layout(xaxis_type="log", yaxis_tickformat=".2e", title="flooring stoch"))#, title="kdam = $(params_rtc[kdam])"))

pl = [p p1]

plot(scattergl(x=df[df.event .!= 0, :time], y=df[df.event .!= 0, :event], mode="markers", marker_color=df[df.event .!= 0, :event], marker=attr(colorscale="Viridis")), Layout(xaxis_type="log"))

p = plot(scattergl(x=df1[df1.event .!= 0, :time], y=df1[df1.event .!= 0, :event], mode="markers", marker_color=df1[df1.event .!= 0, :event], marker=attr(colorscale="Viridis")), Layout(xaxis_type="log"));
p1 = plot(scattergl(x=df2[df2.event .!= 0, :time], y=df2[df2.event .!= 0, :event], mode="markers", marker_color=df2[df2.event .!= 0, :event], marker=attr(colorscale="Viridis")), Layout(xaxis_type="log"));
p3 = plot(scattergl(x=df3[df3.event .!= 0, :time], y=df3[df3.event .!= 0, :event], mode="markers", marker_color=df3[df3.event .!= 0, :event], marker=attr(colorscale="Viridis")), Layout(xaxis_type="log"))

[p p1 p3]



react_names = [:thres_val, :tscr_a, :tscr_b, :tscr_r, :tlr_a, :tlr_b, :tlr_r, :Vinflux, :Vdam, :Vtag, :Vrep, :dil_rh, :dil_rd, :dil_rt, :dil_rma, :dil_rmb, :dil_rmr, :dil_rtca, :dil_rtcb, :dil_rtcr, :deg_rd, :deg_rma, :deg_rmb, :deg_rmr]
df_react = DataFrame([name => Float64[] for name in react_names])
thresh_range = 10 .^ range(log10(0.1),log10(50), length=5)
for i in ProgressBar(thresh_range)
    X0 = CSV.read("$PATH/stochastic_hybrid/X0.dat", Tables.matrix, header=false)
    df = run_stoch(X0, i, 0.1)
    stochReact = df[df.event .!="0",:event]
    stochReact = [parse.(Bool, split(replace(replace(i, "Bool[" => ""), "]" => ""), ", ")) for i in stochReact]
    # contain_detReact = [count(x -> x == 0, stochReact[i]) for i in range(1,length(stochReact))]
    det_react = [findall(x->x==0, i) for i in stochReact]
    det_react = filter(!isempty, det_react)
    react_num_det=[]
    for num in range(1,23)
        push!(react_num_det, sum([count(x -> x == num, i) for i in det_react]))
    end
    pushfirst!(react_num_det, i)
    push!(df_react, react_num_det)
end

df_react
tot_det = sum.(eachrow(df_react[:,2:end]))
df_react.tot_det = tot_det
df_react

react_names = [:tscr_a, :tscr_b, :tscr_r, :tlr_a, :tlr_b, :tlr_r, :Vinflux, :Vdam, :Vtag, :Vrep, :dil_rh, :dil_rd, :dil_rt, :dil_rma, :dil_rmb, :dil_rmr, :dil_rtca, :dil_rtcb, :dil_rtcr, :deg_rd, :deg_rma, :deg_rmb, :deg_rmr]
df_react = DataFrame([name => Float64[] for name in react_names])
print("running")
for i in range(1,100)
    df = run_stoch(X0, 5, 0.1)
    arr = [count(x -> x == i, df[!,:event]) for i in range(1,23)]
    push!(df_react, arr)
end
print("finished")
CSV.write("/home/hollie_hindley/Documents/stochastic_hybrid/reaction_freq.csv",df_react)
print("file saved")

df_react = DataFrame(CSV.File("/home/hollie_hindley/Documents/stochastic_hybrid/reaction_freq.csv"))
plot([scatter(x=range(1,100), y=col, name="$(names(df_react)[i])") for (col, i) in zip(eachcol(df_react),range(1,length(react_names)))])#, Layout(yaxis_type="log"))





# plot(scatter(x=df.time, y=df.rh), Layout(xaxis_type="log", yaxis_tickformat=".2e"))
# plot(scatter(x=df[df.event .== 0, :time], y=df[df.event .== 0, :rh]), Layout(xaxis_type="log", yaxis_tickformat=".2e"))



