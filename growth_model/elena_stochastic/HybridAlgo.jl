using StatsBase
using Distributions
using Random
using DataFrames


include("Index.jl") # add indexing functions
include("InitialConditions.jl") # add initial conditions functions 
include("Parameters.jl") # add parameter functions
include("ModelDefinition.jl")
include("ODEint.jl")

# Create a struct for saving variables
struct Output
    t::Matrix{Float64}
    X::Matrix{Float64}
    X0::Matrix{Float64}
    lam::Vector{Float64}
end


function Lambda(par,X)
    Kgamma =  par[pidx("Kp")]
    gamma      =  par[pidx("gmax")] * X[sidx("ATP")] /( Kgamma +  X[sidx("ATP")]);
    lam        = (((X[sidx("RIB_mENZ")] + X[sidx("RIB_mTRP")] + X[sidx("RIB_mHKP")] + X[sidx("RIB_mRIB")]).*gamma)./(par[pidx("Mref")]));
    return lam
end

function dump(event,time,state;output=stdout)

    lambda = Lambda(par,state)

    print(output,"$event\t$time\t$lambda\t")
    
    states = state[1:I.nrOfItems-1]
    for s in states
        print(output,"$s\t")
    end
    print(output,"\n")
    
    flush(output)
end


function HybridAlgo(X0, par, options, prop, S; out=stdout)

# propen, S, propList  = defineStochModel(par, I)

 # Assign variables
    t0          = 0                       # normally 0

     # Hyrid Algo
    thresh      = options["threshold"]          # Threshold to decide between determinisitic or stochastic reaction
    FixDetReact = options["FixDetReact"]        # Reactions to be treated determinisitically
    tint        = options["tspan"]              # Max time for cell cycle (imp to define in main file)
    tf          = tint                          # final simulation time (it is an adaptable variable)
    sampleFreq  = options["samplingFreq"]       # imp to define so that you don't create too big of a matrix
    
    global t_save =  sampleFreq



     # Vectors for storage 
    s0          = t0        # current system time
    Xs          = X0        # matrix storing all simulated system states (one row per time point)


    getpartitioning(X) = determine_partitioning(X, prop, thresh, FixDetReact)
     # safeguard against (mis-)evaluating propensities with the expanded state 
     # vector which contains the cumulative propensities. 
     nx = I.nrOfItems - 1

     # determine initial partitioning into deterministic and stochastic reactions
    isStochReact = getpartitioning(X0)  

    global stochReact = false
    global savedat    = false

    dump(0,t0,X0;output=out) # store system state at the start of the simulation


        while (tf > s0)  #   while final simulation time not reached

        #   If final time is not reached
        #   draw random number to determine the time of the next reaction
        #   xi is used as a threshold to stop integration with the event function

        exp1 = Exponential(1)
        global xi = rand(exp1,1) 
    
        #   the ODE solver uses the event function to stop the integration when a
        #   stochastic reaction should fire: when total stochastic propensity is
        #   < xi
        #   the event function can also detect when a new replication fork
        #   needs to be initiated or when the cell has to divide 

        #    dump(1,tinit[end],X0;output=out) # If you want to store all stoch. events, uncomment this line
       
        X0[sidx("TotProp")] = 0;

        tspan = [s0,tf]
        prob = ODEProblem(odefunction!, X0, tspan, isStochReact)

        stochreact_cb = DiscreteCallback(conditionStoch, affect_term!, save_positions=(false, false))
        sampling_cb   = DiscreteCallback(conditionSampling, affect_term!, save_positions=(false, false))
        # callbacks = [stochreact_cb, sampling_cb]

        # Solve the ODE with the conditions
        sol = solve(prob, Rodas5(),callback=CallbackSet(stochreact_cb, sampling_cb), save_everystep = false, abstol=1e-12,reltol=1e-9)#, dtmax = 0.05)
        ## Changed callback specification from original code

        # Access the solution
        ss = sol.t
        Y = vcat(sol.u...)'

        #   store state update 
        #   only store the values related to the simulated species, not
        #   integrated propensities (these are re-set to 0 at each iteration)
        X0       = Y[:,end]'
        # Species  = vcat(Species, Y[1:nx,end]')
        # tvec     = hcat(tvec, ss[end])
        s0       = ss[end]  

        ##  Saving
        if savedat  
        dump(0,ss[end],X0;output=out)
        global t_save = ss[end] + sampleFreq
        savedat = false
        end 

        # @show X0[sidx("TotProp")]


        if stochReact  
        #   re-calculate all current reaction propensities
        a  = prop(X0) 
        
        #   set all reaction propensities of deterministic reactions to 0 as
        #   these cannot fire now
        a[.!isStochReact] .= 0 

        #   calculate the total stochastic reaction propensity
        a0 = sum(a) 
        
        #   determine which stochastic reaction fires
        j = sample(1:length(a), Weights(a ./ a0), 1, replace=false)[1]
        # global whatstochreact = cat(whatstochreact,j, dims = 1 )

        # println(j)
                
        #   update the system time (creating a duplicate time)
        # ts = [ts  s0] 

        Update = S[:,j] 
        # Update = vcat(Update, zeros(length(X0) - length(Update)))  #To allow the sum when new replication forks have happened

        #   update the system state
        X0  =  X0 + Update' 

        #   update the partitioning
        isStochReact = getpartitioning(X0) 

        stochReact = false
        
    end #end stoch            

end  # end while loop

return X0
end  # end function






