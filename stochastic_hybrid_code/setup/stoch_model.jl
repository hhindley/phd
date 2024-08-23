function defineStochModel(par, indV)
    n_react = 14
    nu = zeros(indV.nrOfItems, n_react)
    propList = Vector{Function}(undef, n_react)

    Vref = 1
    # v_sf(X) = Vref/X[vidx(:V)]
    v_sf(X) = X[vidx(:V)]

    sf(X) = (1e6/(6.022e23*(X[vidx(:V)]*1e-15)))
    atp(X) = par[pidx(:atp)]/sf(X)

    alpha(X) = X[vidx(:rt)]/par[pidx(:kr)] * v_sf(X) # unitless
    fa(X) = (1+alpha(X))^6/(par[pidx(:L)]*((1+par[pidx(:c)]*alpha(X))^6)+(1+alpha(X))^6)  # unitless 
    ra(X) = fa(X)*X[vidx(:rtcr)] # uM 
    
    # transcription
    Voc(X) = (par[pidx(:Vmax_init)]*atp(X)/((par[pidx(:Km_init)]/v_sf(X))+atp(X)))
    sig_o(X) = ra(X)*Voc(X)/par[pidx(:k_diss)] # uM

    tscr_ab(X) = (sig_o(X)*par[pidx(:ω_ab)]*atp(X)/(par[pidx(:θtscr)]/v_sf(X)+atp(X)) )# uM min-1
    tscr_r(X) = (par[pidx(:ω_r)]/v_sf(X)*atp(X)/((par[pidx(:θtscr)]/v_sf(X))+atp(X))) # uM min-1

    tlr_el(X) = (par[pidx(:g_max)]*atp(X)/((par[pidx(:θtlr)]/v_sf(X))+atp(X)))
    tlr_a(X) = (1/par[pidx(:na)])*par[pidx(:kc)]*X[vidx(:rh)]*X[vidx(:rm_a)]*tlr_el(X)
    tlr_b(X) = (1/par[pidx(:nb)])*par[pidx(:kc)]*X[vidx(:rh)]*X[vidx(:rm_b)]*tlr_el(X)
    tlr_r(X) = (1/par[pidx(:nr)])*par[pidx(:kc)]*X[vidx(:rh)]*X[vidx(:rm_r)]*tlr_el(X)

    # # ribosomes
    Vrep(X) = X[vidx(:rtcb)]*X[vidx(:rt)]*par[pidx(:krep)]/(X[vidx(:rt)]+par[pidx(:km_b)]/v_sf(X)) # uM min-1 
    Vdam(X) = par[pidx(:kdam)]*X[vidx(:rh)] # uM min-1
    Vinflux(X) = par[pidx(:kin)] * par[pidx(:g_max)]*atp(X)/(par[pidx(:θtlr)]/v_sf(X)+atp(X)) # uM min-1 
    Vtag(X) = X[vidx(:rtca)]*X[vidx(:rd)]*par[pidx(:ktag)]/(X[vidx(:rd)]+par[pidx(:km_a)]/v_sf(X)) # uM min-1 


    # 1st order - do nothing
    # 2nd order - divide by v_sf
    # zero order - multiply by v_sf # we do the opposite because we are using v_sf(X) which already divides by V

    # stoichiometric vectors and reaction propensities
    # transcription  
    # reaction 1: * → rm_a - what order?? 1st order with respect to sig_o?
    mu = 1
    nu[vidx(:rm_a), mu] = 1
    nu[vidx(:rm_b), mu] = 1
    propList[mu] = X -> tscr_ab(X) 


    # reaction 2: * → rm_r - zero order 
    mu += 1
    nu[vidx(:rm_r), mu] = 1
    propList[mu] = X -> tscr_r(X) 

    # translation
    # reaction 3: rm_a + rh → rtca - 2nd order 
    mu += 1
    nu[vidx(:rtca), mu] = 1
    propList[mu] = X -> tlr_a(X) * v_sf(X)

    # reaction 4: rm_b + rh → rtcb - 2nd order 
    mu += 1
    nu[vidx(:rtcb), mu] = 1
    propList[mu] = X -> tlr_b(X) * v_sf(X)

    # reaction 5: rm_r + rh → rtcr - 2nd order 
    mu += 1
    nu[vidx(:rtcr), mu] = 1
    propList[mu] = X -> tlr_r(X) * v_sf(X)

    # influx of rh 
    # reaction 6: * → rh - zero order
    mu += 1
    nu[vidx(:rh), mu] = 1
    propList[mu] = X -> Vinflux(X) / v_sf(X)

    # damage
    # reaction 7: rh → rd - 1st order 
    mu += 1
    nu[vidx(:rh), mu] = -1
    nu[vidx(:rd), mu] = 1
    propList[mu] = X -> Vdam(X) 

    # tagging
    # reaction 8: rd + rtca → rt - 2nd order 
    mu += 1
    nu[vidx(:rd), mu] = -1
    nu[vidx(:rt), mu] = 1
    propList[mu] = X -> Vtag(X) 

    # repair
    # reaction 9: rt + rtcb → rh - 2nd order 
    mu += 1
    nu[vidx(:rt), mu] = -1
    nu[vidx(:rh), mu] = 1
    propList[mu] = X -> Vrep(X) 

    # degradation
    # reaction 10: rd → * - 1st order 
    mu += 1
    nu[vidx(:rd), mu] = -1
    propList[mu] = X -> par[pidx(:kdeg)] * X[vidx(:rd)]

    # reaction 11: rm_a → * - 1st order 
    mu += 1
    nu[vidx(:rm_a), mu] = -1
    propList[mu] = X -> par[pidx(:d)] * X[vidx(:rm_a)]

    # reaction 12: rm_b → * - 1st order
    mu += 1
    nu[vidx(:rm_b), mu] = -1
    propList[mu] = X -> par[pidx(:d)] * X[vidx(:rm_b)]

    # reaction 13: rm_r → * - 1st order 
    mu += 1
    nu[vidx(:rm_r), mu] = -1
    propList[mu] = X -> par[pidx(:d)] * X[vidx(:rm_r)]

    # reaction 14: V growth
    mu += 1
    nu[vidx(:V), mu] = 1
    propList[mu] = X -> X[vidx(:V)] * par[pidx(:lam)]   

    propFun(X) = [max(f(X),0) for f in propList]
    return propFun, nu, propList
end

