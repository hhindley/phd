

function defineStochModel(par, indV)
    n_react = 23
    nu = zeros(indV.nrOfItems, n_react)
    propList = Vector{Function}(undef, n_react)

    alpha(X) = X[vidx(:rt)]/par[pidx(:kr)] # unitless
    fa(X) = (1+alpha(X))^6/(par[pidx(:L)]*((1+par[pidx(:c)]*alpha(X))^6)+(1+alpha(X))^6) # unitless 
    ra(X) = fa(X)*X[vidx(:rtcr)] # uM 
    
    # transcription
    Voc = par[pidx(:Vmax_init)]*par[pidx(:atp)]/(par[pidx(:Km_init)]+par[pidx(:atp)]) # uM min-1 
    sig_o(X) = ra(X)*Voc/par[pidx(:k_diss)] # uM

    tscr_ab(X) = sig_o(X)*par[pidx(:ω_ab)]*par[pidx(:atp)]/(par[pidx(:θtscr)]+par[pidx(:atp)]) # uM min-1
    tscr_r = par[pidx(:ω_r)]*par[pidx(:atp)]/(par[pidx(:θtscr)]+par[pidx(:atp)]) # uM min-1

    tlr_el = par[pidx(:g_max)]*par[pidx(:atp)]/(par[pidx(:θtlr)]+par[pidx(:atp)])
    tlr_a(X) = (1/par[pidx(:na)])*par[pidx(:kc)]*X[vidx(:rh)]*X[vidx(:rm_a)]*tlr_el
    tlr_b(X) = (1/par[pidx(:nb)])*par[pidx(:kc)]*X[vidx(:rh)]*X[vidx(:rm_b)]*tlr_el
    tlr_r(X) = (1/par[pidx(:nr)])*par[pidx(:kc)]*X[vidx(:rh)]*X[vidx(:rm_r)]*tlr_el

    # # ribosomes
    Vrep(X) = X[vidx(:rtcb)]*X[vidx(:rt)]*par[pidx(:krep)]/(X[vidx(:rt)]+par[pidx(:km_b)]) # uM min-1 
    Vdam(X) = par[pidx(:kdam)]*X[vidx(:rh)] # uM min-1
    Vinflux = par[pidx(:kin)] * par[pidx(:g_max)]*par[pidx(:atp)]/(par[pidx(:θtlr)]+par[pidx(:atp)]) # uM min-1 
    Vtag(X) = X[vidx(:rtca)]*X[vidx(:rd)]*par[pidx(:ktag)]/(X[vidx(:rd)]+par[pidx(:km_a)]) # uM min-1 

    # stoichiometric vectors and reaction propensities
    
    # transcription  
    # reaction 1: * → rm_a
    mu = 1
    nu[vidx(:rm_a), mu] = 1
    propList[mu] = X -> tscr_ab(X)

    # reaction 2: * → rm_b
    mu += 1
    nu[vidx(:rm_b), mu] = 1
    propList[mu] = X -> tscr_ab(X)

    # reaction 3: * → rm_r
    mu += 1
    nu[vidx(:rm_r), mu] = 1
    propList[mu] = X -> tscr_r

    # translation
    # reaction 4: rm_a → rtca
    mu += 1
    nu[vidx(:rtca), mu] = 1
    propList[mu] = X -> tlr_a(X)

    # reaction 5: rm_b → rtcb
    mu += 1
    nu[vidx(:rtcb), mu] = 1
    propList[mu] = X -> tlr_b(X)

    # reaction 6: rm_r → rtcr
    mu += 1
    nu[vidx(:rtcr), mu] = 1
    propList[mu] = X -> tlr_r(X)

    # influx of rh 
    # reaction 7: * → rh
    mu += 1
    nu[vidx(:rh), mu] = 1
    propList[mu] = X -> Vinflux

    # damage
    # reaction 8: rh → rd
    mu += 1
    nu[vidx(:rh), mu] = -1
    nu[vidx(:rd), mu] = 1
    propList[mu] = X -> Vdam(X)

    # tagging
    # reaction 9: rd → rt
    mu += 1
    nu[vidx(:rd), mu] = -1
    nu[vidx(:rt), mu] = 1
    propList[mu] = X -> Vtag(X)

    # repair
    # reaction 10: rt → rh
    mu += 1
    nu[vidx(:rt), mu] = -1
    nu[vidx(:rh), mu] = 1
    propList[mu] = X -> Vrep(X)

    # dilution
    # reaction 11: rh → *
    mu += 1
    nu[vidx(:rh), mu] = -1
    propList[mu] = X -> par[pidx(:lam)] * X[vidx(:rh)]

    # reaction 12: rd → *
    mu += 1
    nu[vidx(:rd), mu] = -1
    propList[mu] = X -> par[pidx(:lam)] * X[vidx(:rd)]

    # reaction 13: rt → *
    mu += 1
    nu[vidx(:rt), mu] = -1
    propList[mu] = X -> par[pidx(:lam)] * X[vidx(:rt)]

    # reaction 14: rm_a → *
    mu += 1
    nu[vidx(:rm_a), mu] = -1
    propList[mu] = X -> par[pidx(:lam)] * X[vidx(:rm_a)]

    # reaction 15: rm_b → *
    mu += 1
    nu[vidx(:rm_b), mu] = -1
    propList[mu] = X -> par[pidx(:lam)] * X[vidx(:rm_b)]

    # reaction 16: rm_r → *
    mu += 1
    nu[vidx(:rm_r), mu] = -1
    propList[mu] = X -> par[pidx(:lam)] * X[vidx(:rm_r)]

    # reaction 17: rtca → *
    mu += 1
    nu[vidx(:rtca), mu] = -1
    propList[mu] = X -> par[pidx(:lam)] * X[vidx(:rtca)]

    # reaction 18: rtcb → *
    mu += 1
    nu[vidx(:rtcb), mu] = -1
    propList[mu] = X -> par[pidx(:lam)] * X[vidx(:rtcb)]

    # reaction 19: rtcr → *
    mu += 1
    nu[vidx(:rtcr), mu] = -1
    propList[mu] = X -> par[pidx(:lam)] * X[vidx(:rtcr)]

    # degradation
    # reaction 20: rd → *
    mu += 1
    nu[vidx(:rd), mu] = -1
    propList[mu] = X -> par[pidx(:kdeg)] * X[vidx(:rd)]

    # reaction 21: rm_a → *
    mu += 1
    nu[vidx(:rm_a), mu] = -1
    propList[mu] = X -> par[pidx(:d)] * X[vidx(:rm_a)]

    # reaction 22: rm_b → *
    mu += 1
    nu[vidx(:rm_b), mu] = -1
    propList[mu] = X -> par[pidx(:d)] * X[vidx(:rm_b)]

    # reaction 23: rm_r → *
    mu += 1
    nu[vidx(:rm_r), mu] = -1
    propList[mu] = X -> par[pidx(:d)] * X[vidx(:rm_r)]

    propFun(X) = [max(f(X),0) for f in propList]
    return propFun, nu, propList
end





