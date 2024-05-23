using StatsBase, Distributions, Random, DataFrames, CSV, PlotlyJS, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools

PATH = "/home/holliehindley/phd/rtc_model/stochastic_hybrid1"

include("/home/holliehindley/phd/rtc_model/parameters/rtc_params.jl")
include("/home/holliehindley/phd/rtc_model/parameters/rtc_params_molecs.jl")
include("$PATH/indexing.jl")
include("$PATH/hybrid_algo.jl")
include("$PATH/stoch_model.jl")

include("/home/holliehindley/phd/rtc_model/stochastic_hybrid1/run_rtc_orig.jl")


options = Dict(
"threshold"  =>  0.,       # Threshold to decide between determinisitic or stochastic reaction
"FixDetReact"=> [14],# [10,11,12,13,14,15,16,17,18],       # Reactions to be treated determinisitically
    "tspan"     =>   200*log(2)/lam_val,     # Max time for cell cycle
    "samplingFreq"  => 2  # for sampling every x mins
)

X0 = collect(get_X0(indV)')
par = collect(get_par(indP)')

getssX0 = true
if getssX0
    fout=open("/home/holliehindley/phd/rtc_model/stochastic_hybrid1/X0.dat","w")
    propen, S, propList = defineStochModel(par, indV)
    nx = indV.nrOfItems-1
    prop(X) = propen(X[1:nx])
    X0 = hybrid_algo(X0, options, prop, S, out=fout)
    X0[vidx(:V)] = 1
    CSV.write("/home/holliehindley/phd/rtc_model/stochastic_hybrid1/X0.dat", DataFrame(X0,:auto), header=false)
else
    X0 = CSV.read("/home/holliehindley/phd/rtc_model/stochastic_hybrid1/X0.dat", Tables.matrix, header=false)
end
X0


function prop(X)
    nx = indV.nrOfItems - 1
    propen(X[1:nx]) 
end
function run_stoch(X0, thresh, kdam)
    par[pidx(:kdam)] = kdam
    threshold = thresh # set to zero for a deterministic result
    options["threshold"] = threshold
    fout=open("/home/holliehindley/phd/rtc_model/stochastic_hybrid1/test1.dat","w")
    global propen, S, propList  = defineStochModel(par, indV)    
    solu = hybrid_algo(X0, options, prop, S, out=fout)
    close(fout)
end
X0
include("/home/holliehindley/phd/rtc_model/stochastic_hybrid1/stoch_model.jl")
time_taken = @elapsed run_stoch(X0, 0, 0.1)

df = DataFrame(CSV.File("/home/holliehindley/phd/rtc_model/stochastic_hybrid1/test1.dat", header=["event", "time", "rm_a", "rtca", "rm_b", "rtcb", "rm_r", "rtcr", "rh", "rd", "rt", "volume"]))
# plot([scattergl(x=df.time, y=df[:,col], name="$(names(df)[i])", marker_color=colours[i], legendgroup="$i", showlegend=true) for (col, i) in zip(names(eachcol(df[:,3:end-2])), range(3,length(names(df))-2))])#, title="kdam = $(params_rtc[kdam])"))
# plot([scattergl(x=df.time, y=df[:,col] ./df.volume, name="$(names(df)[i])", marker_color=colours[i], legendgroup="$i", showlegend=true) for (col, i) in zip(names(eachcol(df[:,3:end-2])), range(3,length(names(df))-2))])#, title="kdam = $(params_rtc[kdam])"))

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

sf1 = @. (1e6/(6.022e23*(df.volume*1e-15)))
atp1 = @. (par[pidx(:atp)]/sf1)#/df.volume
p=plot(scatter(x=df.time,y=sf1))
p1=plot(scatter(x=df.time,y=atp))
[p; p1]

plot([scatter(x=df_rtc.time, y=repeat([atp_val_molec], length(df_rtc.time))), scatter(x=df.time, y=repeat([mean(atp ./df.volume)],length(df.time)))])

v_sf = @. 1/df.volume

alpha1, fa1, ra1 = all_vars(df)
Voc1 = @. par[pidx(:Vmax_init)]*atp1/((par[pidx(:Km_init)]/v_sf)+atp1) # uM min-1
sig_o1 = @. ra1*Voc1/par[pidx(:k_diss)] 

tscr_ab1 = @. sig_o1*par[pidx(:ω_ab)]*atp1/(par[pidx(:θtscr)]/v_sf+atp1) # uM min-1
tscr_r1 = @. par[pidx(:ω_r)]*atp1/(par[pidx(:θtscr)]+atp1)


tlr_el1 = @. (par[pidx(:g_max)]*atp1/((par[pidx(:θtlr)]/v_sf)+atp1))/df.volume

tlr_a1 = @. (1/par[pidx(:na)])*par[pidx(:kc)]*df.rh*df.rm_a*tlr_el1
tlr_b1 = @. (1/par[pidx(:nb)])*par[pidx(:kc)]*df.rh*df.rm_b*tlr_el1
tlr_r1 = @. (1/par[pidx(:nr)])*par[pidx(:kc)]*df.rh*df.rm_r*tlr_el1

Vrep1 = @. df.rtcb*df.rt*par[pidx(:krep)]/(df.rt+par[pidx(:km_b)]/v_sf) # uM min-1 
Vtag1 = @. df.rtca*df.rd*par[pidx(:ktag)]/(df.rd+par[pidx(:km_a)]/v_sf) # uM min-1 

Vdam1 = @. par[pidx(:kdam)]*df.rh # uM min-1

Vinflux = @. par[pidx(:kin)] * par[pidx(:g_max)]*atp1/(par[pidx(:θtlr)]/v_sf+atp1) # uM min-1 


plot([scatter(x=df_rtc.time, y=alpha), scatter(x=df.time, y=alpha1 ./df.volume)])
plot([scatter(x=df_rtc.time, y=fa), scatter(x=df.time, y=fa1 ./df.volume)])
plot([scatter(x=df_rtc.time, y=ra), scatter(x=df.time, y=ra1 ./df.volume)])
plot([scatter(x=df_rtc.time, y=repeat([Voc], length(df_rtc.time))), scatter(x=df.time, y=Voc1 ./df.volume)])
plot([scatter(x=df_rtc.time, y=sig_o), scatter(x=df.time, y=sig_o1 ./df.volume)])
plot([scatter(x=df_rtc.time, y=tscr_ab), scatter(x=df.time, y=tscr_ab1 ./df.volume)])
plot([scatter(x=df_rtc.time, y=tlr_a), scatter(x=df.time, y=tlr_a1 ./df.volume)])
plot([scatter(x=df_rtc.time, y=tlr_b), scatter(x=df.time, y=tlr_b1 ./df.volume)])
plot([scatter(x=df_rtc.time, y=tlr_r), scatter(x=df.time, y=tlr_r1 ./df.volume)])
plot([scatter(x=df_rtc.time, y=Vrep), scatter(x=df.time, y=Vrep1 ./df.volume)])
plot([scatter(x=df_rtc.time, y=Vtag), scatter(x=df.time, y=Vtag1 ./df.volume)])
plot([scatter(x=df_rtc.time, y=Vdam), scatter(x=df.time, y=Vdam1 ./df.volume)])

plot(scatter(x=df.time, y=df.volume))

par[pidx(:Vmax_init)]

params_rtc_molec[Vmax_init]

par[pidx(:Km_init)]
params_rtc_molec[Km_init]

params_rtc_molec[atp]
atp_val_molec
atp1[1]


params_rtc_molec[Vmax_init]*atp_val_molec/(params_rtc_molec[Km_init]+atp_val_molec) # uM min-1 
Voc = par[pidx(:Vmax_init)]*atp_val_molec/(par[pidx(:Km_init)]+atp_val_molec) # uM min-1 

par[pidx(:Vmax_init)]*atp1[1]/(par[pidx(:Km_init)]+atp1[1])






# df.group = ceil.(Int, (1:nrow(df)) / 20)
# df_new = df[:, Not(:event)]
# df_avg = combine(groupby(df_new, :group), names(df_new, Not(:group)) .=> mean)
# plot([scattergl(x=df_avg.time_mean, y=df_avg[:,col], name="$(names(df_avg)[i])", legendgroup="$i", showlegend=true) for (col, i) in zip(names(eachcol(df_avg[:,4:end-2])), range(4,length(names(df_avg))-2))])#, title="kdam = $(params_rtc[kdam])"))

# avg = mean.(eachcol(df))

df_v = filter(row -> row[:event] == 2, df)
# df_s = filter(row -> row[:event] == 1, df)

plot(scatter(x=df.time, y=df.volume))


df_p = filter(row -> length(row[:event]) > 1, df)

plot(scatter(x=df.time, y=df.volume))
plot([scattergl(x=df.time, y=df[:,col], name="$(names(df)[i])", marker_color=colours[i], legendgroup="$i", showlegend=true) for (col, i) in zip(names(eachcol(df[:,3:end-2])), range(3,length(names(df))-2))])#, title="kdam = $(params_rtc[kdam])"))


plot([scattergl(x=df_v.time, y=df_v[:,col], name="$(names(df_v)[i])", marker_color=colours[i], legendgroup="$i", showlegend=true) for (col, i) in zip(names(eachcol(df_v[:,3:end-2])), range(3,length(names(df_v))-2))])#, title="kdam = $(params_rtc[kdam])"))

plot([scatter(x=df.time, y=df[:,col]) for col in names(eachcol(df[:,3:end]))])



props = [split(replace(i, r"[\[\]\(Any)]" => ""), ",") for i in df_p.event]
props = [parse.(Float64, subarray) for subarray in props]

react_names = [:tscr_ab, :tscr_r, :tlr_a, :tlr_b, :tlr_r, :Vinflux, :Vdam, :Vtag, :Vrep, :deg_rd, :deg_rma, :deg_rmb, :deg_rmr, :V]

df_props = DataFrame([name => Float64[] for name in react_names])

for i in props
    push!(df_props, [i[j] for j in 1:length(props[end])])
end

plot([scattergl(x=df_p.time, y=df_props[:,col], name="$col", mode="markers") for col in names(eachcol(df_props))])

