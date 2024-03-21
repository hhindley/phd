using StatsBase, Distributions, Random, DataFrames, CSV, Plots

PATH = "/home/holliehindley/phd"

include("$PATH/growth_model/elena_stochastic/Index.jl") # add indexing functions
include("$PATH/growth_model/elena_stochastic/InitialConditions.jl") # add initial conditions functions 
include("$PATH/growth_model/elena_stochastic/Parameters.jl") # add parameter functions
include("$PATH/growth_model/elena_stochastic/ModelDefinition.jl")
include("$PATH/growth_model/elena_stochastic/ODEint.jl")
include("$PATH/growth_model/elena_stochastic/HybridAlgo.jl")


# Set options 4 simulation

    options = Dict(
    "threshold"  =>  0.,       # Threshold to decide between determinisitic or stochastic reaction
    "FixDetReact"=>  [],       # Reactions to be treated determinisitically
     "tspan"     =>   1000, #1e6,     # Max time for cell cycle
     "samplingFreq"  => 10     # for sampling every x mins
    )

    getssX0 = false 

    X0 = collect(get_X0(I)') # get initial conditions (has to be a row vector)
    par = get_parameters(J) # get parameters 
    propen, S, propList  = defineStochModel(par, I)
        
    nx = I.nrOfItems - 1
    prop(X) = propen(X[1:nx]) 
    
    HybridAlgo(X0, par, options, prop, S, out=stdout)

    # get X0s --------------------------------
    
    ns = 7
    par[pidx("ns")] = ns

    if getssX0
    # Get deterministic steady state conditions for X0
    fout=open("$PATH/growth_model/elena_stochastic/all_X0ns$ns.dat","w")
    par[pidx("ns")] = ns
    global propen, S, propList  = defineStochModel(par, I)
    nx = I.nrOfItems - 1
    prop(X) = propen(X[1:nx]) 
    X0 = HybridAlgo(X0, par, options, prop, S,out=fout)
    CSV.write("$PATH/growth_model/elena_stochastic/X0ns$ns.dat", DataFrame(X0,:auto), header=false)
    else 
    X0 = CSV.read("$PATH/growth_model/elena_stochastic/X0ns$ns.dat", Tables.matrix, header=false) 
    end 

    # Simulate model -------------------------
    ## Redefine simulation options
    threshold = 1 # Normally at 50 for my model but changed it for testing purposes
    options["threshold"] = threshold
 

    fout=open("$PATH/growth_model/elena_stochastic/test1.dat","w")
    global propen, S, propList  = defineStochModel(par, I)
    nx = I.nrOfItems - 1
    prop(X) = propen(X[1:nx]) 
    sol = HybridAlgo(X0, par, options, prop, S,out=fout)
    close(fout)

   ## Plot results ------------------------------------
   GetX0ss4plotting = DataFrame(CSV.File("HybridStochAlgo/StochAlgo_Isolated/all_X0ns$ns.dat",  header=["event", "time", "lambda","NUTR", "ATP", "mENZ", "mTRP", "mHKP", "mRIB", "RIB_mENZ", "RIB_mTRP", "RIB_mHKP", "RIB_mRIB", "AB_RIB_mENZ", "AB_RIB_mTRP", "AB_RIB_mHKP", "AB_RIB_mRIB", "ENZ", "TRP", "HKP", "RIB", "AB", "TotProp"]) )
   plot(GetX0ss4plotting.time/60, GetX0ss4plotting.lambda*60)

   StochSim = DataFrame(CSV.File("HybridStochAlgo/StochAlgo_Isolated/test1.dat", header=["event", "time", "lambda", "NUTR", "ATP", "mENZ", "mTRP", "mHKP", "mRIB", "RIB_mENZ", "RIB_mTRP", "RIB_mHKP", "RIB_mRIB", "AB_RIB_mENZ", "AB_RIB_mTRP", "AB_RIB_mHKP", "AB_RIB_mRIB", "ENZ", "TRP", "HKP", "RIB", "AB", "TotProp"]))
   adapttime = StochSim.time .+ GetX0ss4plotting.time[end]
   plot!(adapttime/60, StochSim.lambda*60, label="threshold = $threshold", lw=3)

