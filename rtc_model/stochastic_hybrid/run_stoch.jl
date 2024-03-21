using StatsBase, Distributions, Random, DataFrames, CSV, PlotlyJS, DifferentialEquations

PATH = "/home/holliehindley/phd"

include("$PATH/rtc_model/parameters/rtc_params.jl")
include("$PATH/rtc_model/stochastic_hybrid/indexing.jl")
include("$PATH/rtc_model/stochastic_hybrid/hybrid_algo.jl")
include("$PATH/rtc_model/stochastic_hybrid/stoch_model.jl")
include("$PATH/rtc_model/stochastic_hybrid/indexing.jl")

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
    fout=open("$PATH/rtc_model/stochastic_hybrid/all_X0.dat","w")
    propen, S, propList = defineStochModel(par, indV)
    nx = indV.nrOfItems-1
    prop(X) = propen(X[1:nx])
    X0 = hybrid_algo(X0, options, prop, S, out=fout)
    CSV.write("$PATH/rtc_model/stochastic_hybrid/X0.dat", DataFrame(X0,:auto), header=false)
else
    X0 = CSV.read("$PATH/rtc_model/stochastic_hybrid/X0.dat", Tables.matrix, header=false)
end


par[pidx(:kdam)] = 0.1 
threshold = 5 # set to zero for a deterministic result
options["threshold"] = threshold

fout=open("$PATH/rtc_model/stochastic_hybrid/test1.dat","w")
global propen, S, propList  = defineStochModel(par, indV)
nx = indV.nrOfItems - 1
prop(X) = propen(X[1:nx]) 
@time solu = hybrid_algo(X0, options, prop, S, out=fout)
close(fout)


StochSim = DataFrame(CSV.File("$PATH/rtc_model/stochastic_hybrid/test1.dat", header=["event", "time", "rm_a", "rtca", "rm_b", "rtcb", "rm_r", "rtcr", "rh", "rd", "rt", "TotProp"]))

plot(scatter(x=StochSim.time, y=StochSim.event), Layout(xaxis_type="log",yaxis_type="log"))

p_rtc1 = plot([scatter(x=StochSim.time, y=col, name="$(names(StochSim)[i])", legendgroup="$i") for (col, i) in zip(eachcol(StochSim[:,3:end-1]), range(3,length(names(StochSim))-1))], Layout(xaxis_type="log", yaxis_tickformat=".2e"))#, title="kdam = $(params_rtc[kdam])"))

0 in StochSim[30:300,:event] 

count(x -> x == 1, StochSim[!,:event])
count(x -> x == 0, StochSim[!,:event])
StochSim[30:300,:event] 

p_rtc1 = plot([scatter(x=StochSim.time, y=StochSim[StochSim.event .== 0, col], name="$(names(StochSim)[i])", legendgroup="$i") for (col, i) in zip(eachcol(StochSim[:,3:end-1]), range(3,length(names(StochSim))-1))], Layout(xaxis_type="log", yaxis_tickformat=".2e"))#, title="kdam = $(params_rtc[kdam])"))

plot(scatter(x=StochSim.time, y=StochSim.rh), Layout(xaxis_type="log", yaxis_tickformat=".2e"))
plot(scatter(x=StochSim[StochSim.event .== 0, :time], y=StochSim[StochSim.event .== 0, :rh]), Layout(xaxis_type="log", yaxis_tickformat=".2e"))
# isStochReact = determine_partitioning(X0, prop, 1e-8, [])

# a = prop(X0)

# a[isStochReact]


# a0_s = sum(a[isStochReact])

# a_d = a .* .~isStochReact

# dySpecies = S*a_d

# dfposdt = []
# if length(X0) > indV.nrOfItems # why would u be bigger than number of species? 
#     fpos = X0[(indV.nrOfItems + 1):end]
#     dfposdt = ones(length(fpos))
# end
# dySpecies

# ODEProblem(odefunction!, X0, tspan, isStochReact)
