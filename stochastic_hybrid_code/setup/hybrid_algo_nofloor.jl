
struct Output
    t::Matrix{Float64}
    X::Matrix{Float64}
    X0::Matrix{Float64}
    lam::Vector{Float64}
end



function prop(X)
    nx = indV.nrOfItems - 1
    propen(X[1:nx]) 
end
function run_stoch(X0, thresh, kdam, file)
    par[pidx(:kdam)] = kdam
    threshold = thresh # set to zero for a deterministic result
    options["threshold"] = threshold
    fout=open("$file","w")
    # fout=open("/home/hollie_hindley/Documents/stochastic_hybrid/$file","w")
    global propen, S, propList  = defineStochModel(par, indV)    
    solu = hybrid_algo(X0, options, prop, S, out=fout)
    close(fout)
end

function flooring(X0)
    X0_new = deepcopy(X0)
    for i in range(1,9)
        if X0_new[i] < floor(X0_new[i]) + 0.98
            X0_new[i] = floor(X0_new[i])
            # X0_new[end] = floor(X0_new[end])
        else
            X0_new[i] = ceil(X0_new[i])
            # X0_new[end] = ceil(X0_new[end])
        end
    end
    
    return X0_new
end


function hybrid_algo(X0, options, prop, S; out=stdout)
    # save_when = 10 .^ range(log10(1e-6), log10(1e6), length=1000)

    t0 = 0

    thresh = options["threshold"]
    FixDetReact = options["FixDetReact"]
    tf = options["tspan"]
    tint = options["tspan"]
    sampleFreq = options["samplingFreq"]

    # sampleFreq = 1e-2
    # freq(t) = t <= 2000 ? t/100 : log(1 + ((t-2000)^40))

    # freq(t) = log(1 + (t^2))#/10)

    global t_save = sampleFreq

    s0 = t0 # system time
    Xs = X0
    # isStochReact = determine_partitioning(floor.(X0), prop, thresh, FixDetReact)
    partition(X0) = determine_partitioning(X0, prop, thresh, FixDetReact)
    isStochReact = partition((X0))

    div_times = [i*(log(2)/par[pidx(:lam)]) for i in 1:floor(tf/(log(2)/par[pidx(:lam)]))]

    global stochReact = false
    global savedat = false
    global division = false

    dump([0],t0,X0;output=out)

    # Random.seed!(111)
    nu = 0

    while ~division && (tf > s0)
        exp1 = Exponential(1)
        global xi = rand(exp1, 1) # get random number from exp distribution
        # @show xi, X0[vidx(:totProp)]

        X0[vidx(:totProp)] = 0 # set total propensities to zero 
        tspan = [s0,tf]

        prob = ODEProblem(odefunction!, X0, tspan, isStochReact) # sets ode prob for deterministic reactions

        stochreact_cb = DiscreteCallback(conditionStoch, affect_term!, save_positions=(false,false)) # is called at every time step during integration, if total propensities are greater than the random number then terminate integration and set stochReact to true
        sampling_cb = DiscreteCallback(conditionSampling, affect_term!, save_positions=(false,false)) # is called at every time step during integration, if t is greater than sampling frequency then make savedat = true and terminate integration
        
        condition_div(u,t,integrator) = t in div_times

        division_cb = DiscreteCallback(condition_div, affect_division!, save_positions=(false,false)) # is called at every time step during integration, if t is greater than division time then make division = true and terminate integration
        # division_cb = PresetTimeCallback(div_times, affect_division!) # is called at every time step during integration, if t is greater than division time then make division = true and terminate integration

        solu = solve(prob, Rodas4(), callback=CallbackSet(division_cb, stochreact_cb, sampling_cb), tstops = div_times, save_everystep=false) # solve for deterministic reactions and stop when one of the callbacks is reached

        ss = solu.t # time solution
        Y = vcat(solu.u...)' # species solution

        X0 = Y[:,end]' # new initials

        s0 = ss[end] # new time point 
        # @show division
        if division
            # dump(20,ss[end],X0;output=out)
            X0[vidx(:V)] = X0[vidx(:V)]/2
            for j in 1:length(X0)-2
                # Partition the molecules stochastically using a binomial distribution
                # pDiv=rand(Beta(200.,200.))
                pDiv=0.5
                # X0[j] = X0[j] > 0. ? float(rand(Binomial(round(Int, X0[j]), pDiv))) : 0                
                X0[j] = X0[j] / 2
            end
            # dump(20,ss[end],X0;output=out)
            division = false
        end

        # @show s0, t_save
        if savedat # if savedat is true 
            # println(ss[end])
            # nu += 1

            # dump(X0[vidx(:totProp)],ss[end],X0;output=out)
            dump(prop(X0),ss[end],X0;output=out)

            # dump(stochReact,ss[end],X0;output=out)
            # dump([X0[vidx(:totProp)][1],xi[1]],ss[end],X0;output=out)

            # @show X0[vidx(:totProp)], xi

            global t_save = ss[end] +  sampleFreq #freq(ss[end]) # update t_save to be a distance from the current timestep 

            savedat = false
        end


        # @show stochReact
        if stochReact

            a = prop((X0))
            a[.!isStochReact] .= 0 # zeros out propensities of deterministic reactions

            a0 = sum(a)

            j = sample(1:length(a), Weights(a./a0), 1, replace=false)[1] # randomly selects a reaction index based on the probabilities given by the propensities 

            Update = S[:,j] # selects the chosen reaction 

            X0 = X0 + Update' # updates the state variables by adding stoichiometric changes from the selected reaction
            clamp!(X0, 0, Inf) # clamps the state variables to be non-negative
            # dump([a1,a2],ss[end],X0;output=out)

            # nu += 1
            dump([j], ss[end],X0;output=out)

            isStochReact = partition((X0)) # updates status of stochastic reactions based on new X0

            stochReact = false
        end

    end
    
    
    return X0
end

function determine_partitioning(u, prop, thresh, FixDetReact)

    isStochReact = prop((u)) .< thresh # if the propensity is less than the threshold then the reaction will be stochastic 

    if !isempty(FixDetReact)
        for i in FixDetReact
            isStochReact[i] = 0
        end
    end

    return isStochReact
end
# dump(0,t0,X0;output=out)
function dump(event, time, state; output=stdout)
    # rh = state[vidx(:rh)]

    # print(output, "$time\t$(state[vidx(:rm_a)])\t$(state[vidx(:rh)])")
    print(output, "$event\t$time\t")
    states = state[1:indV.nrOfItems]
# 
    for s in states 
        print(output, "$s\t")
    end
    print(output, "\n")
    flush(output)
end

function odefunction!(du, u, isStochReact, t) # dont need to pass parameters here because that was done when building the stochastic model so this information is already stored within the propensities 
    u[u .< 0] .= 0 # zero any negatives (non-negativity)

    a = prop((u)) # calculate reaction propensities
    a2 = prop(u)

    a0_s = sum(a[isStochReact]) # total reaction propensities for stochastic reactions

    a_d = a2 .* .~isStochReact # .~ flips the boolean array so true becomes false and vice versa and then sets the stochastic reactions to zero so they are not counted as deterministic 

    dySpecies = S*a_d # calculates the derivatives of the deterministic reactions by doing the dot product (for each species (i), multiply each stoichiometric coeff in the ith row of S with the corresponding propensity in a_d and then sum the products)
    # dySpecies is a vector of the rate of change of all species for deterministic reactions (du in normal ode models)
    dySpecies[end] = a0_s # updates totprop

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

function affect_division!(integrator)
    terminate!(integrator)

    global division = true 

    return division 
end