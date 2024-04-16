
struct Output
    t::Matrix{Float64}
    X::Matrix{Float64}
    X0::Matrix{Float64}
    lam::Vector{Float64}
end

function hybrid_algo(X0, options, prop, S; out=stdout)
    # save_when = 10 .^ range(log10(1e-6), log10(1e6), length=1000)

    t0 = 0

    thresh = options["threshold"]
    FixDetReact = options["FixDetReact"]
    tf = options["tspan"]
    tint = options["tspan"]
    # sampleFreq = options["samplingFreq"]

    sampleFreq = 1e-2
    freq(t) = t <= 2000 ? t/100 : log(1 + ((t-2000)^40))

    # freq(t) = log(1 + (t^2))#/10)

    global t_save = sampleFreq

    s0 = t0 # system time
    Xs = X0
    isStochReact = determine_partitioning(floor.(X0), prop, thresh, FixDetReact)

    global stochReact = false
    global savedat = false

    dump(0,t0,X0;output=out)

    # Random.seed!(123)
    nu = 1 

    while (tf > s0)
        exp1 = Exponential(1)
        global xi = rand(exp1, 1) # get random number from exp distribution
        # @show xi
        X0[vidx(:totProp)] = 0 # set total propensities to zero 
        tspan = [s0,tf]

        prob = ODEProblem(odefunction!, X0, tspan, isStochReact) # sets ode prob for deterministic reactions

        stochreact_cb = DiscreteCallback(conditionStoch, affect_term!, save_positions=(false,false)) # is called at every time step during integration, if total propensities are greater than the random number then terminate integration and set stochReact to true
        sampling_cb = DiscreteCallback(conditionSampling, affect_term!, save_positions=(false,false)) # is called at every time step during integration, if t is greater than sampling frequency then make savedat = true and terminate integration


        solu = solve(prob, Rodas4(), callback=CallbackSet(stochreact_cb, sampling_cb), save_everystep=false) # solve for deterministic reactions and stop when one of the callbacks is reached


        ss = solu.t # time solution
        Y = vcat(solu.u...)' # species solution

        X0 = Y[:,end]' # new initials

        s0 = ss[end] # new time point 
        
        # @show s0, freq(s0)
        if savedat # if savedat is true 
            # println(ss[end])
            # dump(prop(X0),ss[end],X0;output=out)
            dump(ss[end],ss[end],X0;output=out)

            # dump(stochReact,ss[end],X0;output=out)
            # dump([X0[vidx(:totProp)],xi,prop(X0)],ss[end],X0;output=out)

            # @show X0[vidx(:totProp)]
            global t_save = ss[end] +  freq(ss[end]) # update t_save to be a distance from the current timestep 
            savedat = false
        end
        
        

        if stochReact
            # X0_round = round.(X0)
            # X0_round = ceil.(X0)
            X0_floor = floor.(X0)
            a = prop(X0_floor)
            # a = prop(X0)
            a[.!isStochReact] .= 0 # zeros out propensities of deterministic reactions

            a0 = sum(a)

            j = sample(1:length(a), Weights(a./a0), 1, replace=false)[1] # randomly selects a reaction index based on the probabilities given by the propensities 

            Update = S[:,j] # selects the chosen reaction 

            X0 = X0 + Update' # updates the state variables by adding stoichiometric changes from the selected reaction

            # dump(prop(X0),ss[end],X0;output=out)

            nu += 1
            # dump([nu,j],ss[end],X0;output=out)

            isStochReact = determine_partitioning(X0, prop, thresh, FixDetReact) # updates status of stochastic reactions based on new X0
            

            stochReact = false
        end
    end
    
    
    return X0
end


function determine_partitioning(u, prop, thresh, FixDetReact)

    isStochReact = prop(u) .< thresh # if the propensity is less than the threshold then the reaction will be stochastic 

    if !isempty(FixDetReact)
        isStochReact[FixDetReact] .= 0
    end

    return isStochReact
end
# dump(0,t0,X0;output=out)
function dump(event, time, state; output=stdout)
    # rh = state[vidx(:rh)]

    print(output, "$time\t$(state[vidx(:rm_a)])\t$(state[vidx(:rh)])")
#     print(output, "$event\t$time\t")
#     states = state[1:indV.nrOfItems-1]
# # 
#     for s in states 
#         print(output, "$s\t")
#     end
    print(output, "\n")
    flush(output)
end

function odefunction!(du, u, isStochReact, t) # dont need to pass parameters here because that was done when building the stochastic model so this information is already stored within the propensities 
    u[u .< 0] .= 0 # zero any negatives (non-negativity)

    a = prop(u) # calculate reaction propensities

    a0_s = sum(a[isStochReact]) # total reaction propensities for stochastic reactions

    a_d = a .* .~isStochReact # .~ flips the boolean array so true becomes false and vice versa and then sets the stochastic reactions to zero so they are not counted as deterministic 

    dySpecies = S*a_d # calculates the derivatives of the deterministic reactions by doing the dot product (for each species (i), multiply each stoichiometric coeff in the ith row of S with the corresponding propensity in a_d and then sum the products)
    # dySpecies is a vector of the rate of change of all species for deterministic reactions (du in normal ode models)
    dySpecies[end] = a0_s
    du .= vcat(dySpecies)'

    return du
end

function conditionStoch(u, t, integrator) 
    if u[vidx(:totProp)] >= xi[1]
        global stochReact = true
    else
        global stochReact = false
    end

    return stochReact

end

function conditionSampling(u, t, integrator)

    if t >= t_save
        global savedat = true
    else
        global savedat = false
    end

    return savedat
end


function affect_term!(integrator) 
    terminate!(integrator)
end