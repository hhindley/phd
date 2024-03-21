#### Author: Elena Pascual Garcia 
#### Date: 14th Dec 2023
#### Project: Hybrid Stochastic Bacterial Growth & Division Model + Antibiotics

### Definition ODE model in matrix form. Input: parameter vector + indexing structure
### Output: propFun = propensities of each reaction
###         nu = stochiometrix matrix
###         propList = reaction list
### Usage instructions can be found at the end of the script

function defineStochModel(par, I)

    nrOfRct = 51
    nu = zeros(I.nrOfItems, nrOfRct)
    propList = Vector{Function}(undef, nrOfRct)
    

    NA = par[pidx("NA")]  # Avogadro number
    SF  = NA * 1e-6 * 1e-15
    SF2 = NA * 1e-6 * 1e-15

    # Parameter names
    NUTR_ext = par[pidx("NUTR_ext")]
    ns       = par[pidx("ns")]
    vt       = par[pidx("vt")]
    vm       = par[pidx("vm")]
    Km       = par[pidx("Km")]
    Kt       = par[pidx("Kt")]
    thetar   = par[pidx("thetar")]
    thetax   = par[pidx("thetax")]
    we       = par[pidx("we")]
    wr       = par[pidx("wr")]
    wq       = par[pidx("wq")]
    Kq       = par[pidx("Kq")]
    nq       = par[pidx("nq")]
    gmax     = par[pidx("gmax")]
    Kp       = par[pidx("Kp")]
    nx       = par[pidx("nx")]
    nr       = par[pidx("nr")]
    kb       = par[pidx("kb")]
    ku       = par[pidx("ku")]
    dm       = par[pidx("dm")]
    Mref     = par[pidx("Mref")]
    AB_ext   = par[pidx("AB_ext")] * SF2
    kon      = par[pidx("kon")]
    KD       = par[pidx("KD")]
    pin      = par[pidx("pin")]
    pout     = par[pidx("pout")]


    # Pre-calculations

    # Define Koff
    koff = kon * KD
    kon /= SF

    # Transport and metabolism rates 
    vimp(X) =  X[sidx("TRP")] * vt * NUTR_ext / (Kt + NUTR_ext)
    vcat(X) =  X[sidx("ENZ")] * vm * X[sidx("NUTR")] / (Km + X[sidx("NUTR")])

 
    # Transcription rates 
    Wm(X) = we  * X[sidx("ATP")] / (thetax  + X[sidx("ATP")])
    Wt(X) = we  * X[sidx("ATP")] / (thetax  + X[sidx("ATP")])
    Wq(X) = wq  * X[sidx("ATP")] / (thetax  + X[sidx("ATP")]) / (1 + (X[sidx("HKP")] / (Kq )) ^ nq)
    Wr(X) = wr  * X[sidx("ATP")] / (thetar  + X[sidx("ATP")])

    # Translation rates
    Kgamma(X) = Kp 
    gamma(X)  = gmax * X[sidx("ATP")] / (Kgamma(X) + X[sidx("ATP")])

    Vm(X) = gamma(X) / nx * X[sidx("RIB_mENZ")]
    Vt(X) = gamma(X) / nx * X[sidx("RIB_mTRP")]
    Vq(X) = gamma(X) / nx * X[sidx("RIB_mHKP")]
    Vr(X) = gamma(X) / nr * X[sidx("RIB_mRIB")]



    ttrate(X)   = (X[sidx("RIB_mENZ")] + X[sidx("RIB_mTRP")] + X[sidx("RIB_mHKP")] + X[sidx("RIB_mRIB")])*gamma(X); 
    lam(X)      = ttrate(X)/Mref;
    # Stochiometric vectors and reaction propensities, listed by reactions ----
   
    ### NUTR -------------------------------------------------------    
   
    mu = 1 # vimp: * -> NUTR
    nu[sidx("NUTR"), mu] = 1
    propList[mu] = X -> vimp(X)

    mu += 1 # vcat: NUTR -> ATP
    nu[sidx("NUTR"), mu] = -1
    nu[sidx("ATP"), mu] = ns
    propList[mu] = X -> vcat(X)

    mu += 1 # NUTR -> *
    nu[sidx("NUTR"), mu] = -1
    propList[mu] = X -> lam(X)*X[sidx("NUTR")]


    ### ATP -------------------------------------------------------------
    # ns*vcat: NUTR -> ATP (described)
    # - ttrate --> described later on

    mu += 1 # ATP -> *
    nu[sidx("ATP"), mu] = -1
    propList[mu] = X -> lam(X)*X[sidx("ATP")]


    ### mENZ ------------------------------------------------------------

    mu += 1 # + Wm : * -> mENZ 
    nu[sidx("mENZ"), mu] = 1
    propList[mu] = X -> Wm(X)

    mu += 1 # # + ku*RIB_mENZ: RIB_mENZ --> RIB+mENZ  
    nu[sidx("mENZ"), mu] = 1
    nu[sidx("RIB"), mu] = 1
    nu[sidx("RIB_mENZ"), mu] = -1
    propList[mu] = X -> ku * X[sidx("RIB_mENZ")]

    mu += 1 # - kb*mENZ*RIB: RIB+mENZ --> RIB_mENZ
    nu[sidx("mENZ"), mu] = -1
    nu[sidx("RIB"), mu] = -1
    nu[sidx("RIB_mENZ"), mu] = 1
    propList[mu] = X -> kb * X[sidx("mENZ")] * X[sidx("RIB")] 

    mu += 1 # +Vm:  RIB_mENZ+ATP --> RIB+mENZ+ENZ 
    nu[sidx("mENZ"), mu] = 1
    nu[sidx("RIB"), mu] = 1
    nu[sidx("ENZ"), mu] = 1
    nu[sidx("RIB_mENZ"), mu] = -1
    nu[sidx("ATP"), mu] = -nx
    propList[mu] = X -> Vm(X)

    mu += 1 # degradation: mENZ: mENZ -> * 
    nu[sidx("mENZ"), mu] = -1
    propList[mu] = X -> dm * X[sidx("mENZ")]

    mu += 1 # dilution: mENZ -> *
    nu[sidx("mENZ"), mu] = -1
    propList[mu] = X -> lam(X)*X[sidx("mENZ")]


    ### mTRP -------------------------------------------------------------
    mu += 1 # + Wt : * -> mTRP 
    nu[sidx("mTRP"), mu] = 1
    propList[mu] = X -> Wt(X)

    mu += 1 # + ku*RIB_mTRP: RIB_mTRP --> RIB+mTRP
    nu[sidx("mTRP"), mu] = 1
    nu[sidx("RIB"), mu] = 1
    nu[sidx("RIB_mTRP"), mu] = -1
    propList[mu] = X -> ku * X[sidx("RIB_mTRP")]

    mu += 1 # - kb*mTRP*RIB: RIB+mTRP --> RIB_mTRP
    nu[sidx("mTRP"), mu] = -1
    nu[sidx("RIB"), mu] = -1
    nu[sidx("RIB_mTRP"), mu] = 1
    propList[mu] = X -> kb * X[sidx("mTRP")] * X[sidx("RIB")] 

    mu += 1 # +Vt:  RIB_mTRP+ATP --> RIB+TRP+mTRP 
    nu[sidx("mTRP"), mu] = 1
    nu[sidx("RIB"), mu] = 1
    nu[sidx("TRP"), mu] = 1
    nu[sidx("RIB_mTRP"), mu] = -1
    nu[sidx("ATP"), mu] = -nx
    propList[mu] = X -> Vt(X)

    mu += 1 # degradation mTRP: mTRP -> *  
    nu[sidx("mTRP"), mu] = -1
    propList[mu] = X -> dm * X[sidx("mTRP")]

    mu += 1 # dilution: mTRP -> *
    nu[sidx("mTRP"), mu] = -1
    propList[mu] = X -> lam(X)*X[sidx("mTRP")]


    ### mHKP --------------------------------------------------------------
    mu += 1 # + Wq : * -> mHKP  
    nu[sidx("mHKP"), mu] = 1
    propList[mu] = X -> Wq(X)

    mu += 1 # + ku*RIB_mHKP: RIB_mHKP --> RIB+mHKP 
    nu[sidx("mHKP"), mu] = 1
    nu[sidx("RIB"), mu] = 1
    nu[sidx("RIB_mHKP"), mu] = -1
    propList[mu] = X -> ku * X[sidx("RIB_mHKP")]

    mu += 1 # - kb*mHKP*RIB: RIB+mHKP --> RIB_mHKP 
    nu[sidx("mHKP"), mu] = -1
    nu[sidx("RIB"), mu] = -1
    nu[sidx("RIB_mHKP"), mu] = 1
    propList[mu] = X -> kb * X[sidx("mHKP")] * X[sidx("RIB")] 

    mu += 1 # + Vq:  RIB_mHKP+ATP --> RIB+mHKP+HKP  
    nu[sidx("mHKP"), mu] = 1
    nu[sidx("RIB"), mu] = 1
    nu[sidx("HKP"), mu] = 1
    nu[sidx("RIB_mHKP"), mu] = -1
    nu[sidx("ATP"), mu] = -nx
    propList[mu] = X -> Vq(X)

    mu += 1 # degradation mHKP: mHKP -> *
    nu[sidx("mHKP"), mu] = -1
    propList[mu] = X -> dm * X[sidx("mHKP")]

    mu += 1 # dilution: mHKP -> *
    nu[sidx("mHKP"), mu] = -1
    propList[mu] = X -> lam(X)*X[sidx("mHKP")]

    
    ### mRIB ---------------------------------------------------------------
    mu += 1 # + Wr : * -> mRIB 
    nu[sidx("mRIB"), mu] = 1
    propList[mu] = X -> Wr(X)

    mu += 1 # + ku*RIB_mRIB: RIB_mRIB --> RIB+mRIB 
    nu[sidx("mRIB"), mu] = 1
    nu[sidx("RIB"), mu] = 1
    nu[sidx("RIB_mRIB"), mu] = -1
    propList[mu] = X -> ku * X[sidx("RIB_mRIB")]

    mu += 1 # - kb*mRIB*RIB: RIB+mRIB --> RIB_mRIB 
    nu[sidx("mRIB"), mu] = -1
    nu[sidx("RIB"), mu] = -1
    nu[sidx("RIB_mRIB"), mu] = 1
    propList[mu] = X -> kb * X[sidx("mRIB")] * X[sidx("RIB")] 

    mu += 1 # +Vr:  RIB_mRIB+ATP --> RIB+mRIB+RIB
    nu[sidx("mRIB"), mu] = 1
    nu[sidx("RIB"), mu] = 2
    nu[sidx("RIB_mRIB"), mu] = -1
    nu[sidx("ATP"), mu] = -nr
    propList[mu] = X -> Vr(X)

    mu += 1 # degradation mTRP: mRIB -> * 
    nu[sidx("mRIB"), mu] = -1
    propList[mu] = X -> dm * X[sidx("mRIB")]

    mu += 1 # dilution: mRIB -> *
    nu[sidx("mRIB"), mu] = -1
    propList[mu] = X -> lam(X)*X[sidx("mRIB")]


    ### RIB_mENZ ---------------------------------------------------------------

    # kb*mENZ*RIB: RIB+mENZ --> RIB_mENZ (described)
    # - ku*RIB_mENZ: RIB_mENZ --> RIB+mENZ (described)
    # -Vm:  RIB_mENZ+ATP --> RIB+mENZ+ENZ (described)

    mu += 1 # + koff*AB_RIB_mENZ:  AB_RIB_mENZ -> AB+RIB_mENZ  
    nu[sidx("AB_RIB_mENZ"), mu] = -1
    nu[sidx("AB"), mu] = 1
    nu[sidx("RIB_mENZ"), mu] = 1
    propList[mu] = X -> koff * X[sidx("AB_RIB_mENZ")]

    mu += 1 # + kon*AB*RIB_mENZ:  AB+RIB_mENZ --> AB_RIB_mENZ
    nu[sidx("AB_RIB_mENZ"), mu] = 1
    nu[sidx("AB"), mu] = -1
    nu[sidx("RIB_mENZ"), mu] = -1
    propList[mu] = X -> kon * X[sidx("AB")] * X[sidx("RIB_mENZ")] 

    mu += 1 # dilution: RIB_mENZ -> *
    nu[sidx("RIB_mENZ"), mu] = -1
    propList[mu] = X -> lam(X)*X[sidx("RIB_mENZ")]


    ### RIB_mTRP ---------------------------------------------------------------

    # kb*mTRP*RIB: RIB+mTRP --> RIB_mTRP (described)
    # - ku*RIB_mTRP: RIB_mTRP --> RIB+mTRP (described)
    # -Vt:  RIB_mTRP+ATP --> RIB+mTRP+TRP (described)
    
    mu += 1 # + koff*AB_RIB_mTRP:  AB_RIB_mTRP -> AB+RIB_mTRP
    nu[sidx("AB_RIB_mTRP"), mu] = -1
    nu[sidx("AB"), mu] = 1
    nu[sidx("RIB_mTRP"), mu] = 1
    propList[mu] = X -> koff * X[sidx("AB_RIB_mTRP")]

    mu += 1 # + kon*AB*RIB_mTRP:  AB+RIB_mTRP --> AB_RIB_mTRP
    nu[sidx("AB_RIB_mTRP"), mu] = 1
    nu[sidx("AB"), mu] = -1
    nu[sidx("RIB_mTRP"), mu] = -1
    propList[mu] = X -> kon * X[sidx("AB")] * X[sidx("RIB_mTRP")] 

    mu += 1 # dilution: RIB_mTRP -> *
    nu[sidx("RIB_mTRP"), mu] = -1
    propList[mu] = X -> lam(X)*X[sidx("RIB_mTRP")]

    ### RIB_mHKP ---------------------------------------------------------------

    # kb*mHKP*RIB: RIB+mHKP --> RIB_mHKP (described)
    # - ku*RIB_mHKP: RIB_mHKP --> RIB+mHKP (described)
    # -Vq:  RIB_mHKP+ATP --> RIB+mHKP+HKP (described)    
    
    mu += 1 # + koff*AB_RIB_mHKP:  AB_RIB_mHKP -> AB+RIB_mHKP
    nu[sidx("AB_RIB_mHKP"), mu] = -1
    nu[sidx("AB"), mu] = 1
    nu[sidx("RIB_mHKP"), mu] = 1
    propList[mu] = X -> koff * X[sidx("AB_RIB_mHKP")]

    mu += 1 # + kon*AB*RIB_mHKP:  AB+RIB_mHKP --> AB_RIB_mHKP
    nu[sidx("AB_RIB_mHKP"), mu] = 1
    nu[sidx("AB"), mu] = -1
    nu[sidx("RIB_mHKP"), mu] = -1
    propList[mu] = X -> kon * X[sidx("AB")] * X[sidx("RIB_mHKP")] 

    mu += 1 # dilution: RIB_mHKP -> *
    nu[sidx("RIB_mHKP"), mu] = -1
    propList[mu] = X -> lam(X)*X[sidx("RIB_mHKP")]


    ### RIB_mRIB ---------------------------------------------------------------

    # kb*mRIB*RIB: RIB+mRIB --> RIB_mRIB (described)
    # - ku*RIB_mRIB: RIB_mRIB --> RIB+mRIB (described)
    # -Vr:  RIB_mRIB+ATP --> RIB+mRIB+RIB (described)
    
    mu += 1 # + koff*AB_RIB_mRIB:  AB_RIB_mRIB -> AB+RIB_mRIB 
    nu[sidx("AB_RIB_mRIB"), mu] = -1
    nu[sidx("AB"), mu] = 1
    nu[sidx("RIB_mRIB"), mu] = 1
    propList[mu] = X -> koff * X[sidx("AB_RIB_mRIB")]

    mu += 1 # + kon*AB*RIB_mRIB:  AB+RIB_mRIB --> AB_RIB_mRIB 
    nu[sidx("AB_RIB_mRIB"), mu] = 1
    nu[sidx("AB"), mu] = -1
    nu[sidx("RIB_mRIB"), mu] = -1
    propList[mu] = X -> kon * X[sidx("AB")] * X[sidx("RIB_mRIB")] 

    mu += 1 # dilution: RIB_mRIB -> *
    nu[sidx("RIB_mRIB"), mu] = -1
    propList[mu] = X -> lam(X)*X[sidx("RIB_mRIB")]


    ### AB ---------------------------------------------------------------

    # kon*AB*(RIB_mENZ+RIB_mTRP+RIB_mHKP+RIB_mRIB) (described)
    # koff*(str_rmm+str_rmt+str_rmq+str_rmr) (described)
    
    mu = mu+1;  # pin * AB: * --> AB  
    nu[sidx("AB"),mu]     = 1;
    propList[mu] = X -> pin * AB_ext;
    
    mu = mu+1;  # - pout * AB: AB --> *  
    nu[sidx("AB"),mu]    = -1;
    propList[mu] = X -> pout * X[sidx("AB")];
    
    mu += 1 # dilution: AB -> *
    nu[sidx("AB"), mu] = -1
    propList[mu] = X -> lam(X)*X[sidx("AB")]

    ### AB_RIB_mENZ ---------------------------------------------------------------
    
    # + kon*AB*RIB_mENZ (described)
    # - koff*str_rmm  (described)

    mu += 1 # dilution: AB_RIB_mENZ -> *
    nu[sidx("AB_RIB_mENZ"), mu] = -1
    propList[mu] = X -> lam(X)*X[sidx("AB_RIB_mENZ")]

    ### AB_RIB_mTRP ---------------------------------------------------------------
    
    # + kon*AB*RIB_mTRP (described)
    # - koff*AB_RIB_mTRP  (described)

    mu += 1 # dilution: AB_RIB_mTRP -> *
    nu[sidx("AB_RIB_mTRP"), mu] = -1
    propList[mu] = X -> lam(X)*X[sidx("AB_RIB_mTRP")]
    
    ### AB_RIB_mHKP ---------------------------------------------------------------
    
    # + kon*AB*RIB_mHKP (described)
    # - koff*AB_RIB_mHKP  (described)

    mu += 1 # dilution: AB_RIB_mHKP -> *
    nu[sidx("AB_RIB_mHKP"), mu] = -1
    propList[mu] = X -> lam(X)*X[sidx("AB_RIB_mHKP")]
    
    ### AB_RIB_mRIB ---------------------------------------------------------------
    
    # + kon*AB*RIB_mRIB (described)
    # - koff*AB_RIB_mRIB  (described)

    mu += 1 # dilution: AB_RIB_mRIB -> *
    nu[sidx("AB_RIB_mRIB"), mu] = -1
    propList[mu] = X -> lam(X)*X[sidx("AB_RIB_mRIB")]
    
    
    ### ENZ ---------------------------------------------------------------
    
    # + Vm:  RIB_mENZ+ATP --> RIB+mENZ+ENZ (described)
    
    mu += 1 # dilution: ENZ -> *
    nu[sidx("ENZ"), mu] = -1
    propList[mu] = X -> lam(X)*X[sidx("ENZ")]
    
    ### TRP ---------------------------------------------------------------
    
    # + Vt:  RIB_mTRP+ATP --> RIB+mTRP+TRP (described)
    
    mu += 1 # dilution: TRP -> *
    nu[sidx("TRP"), mu] = -1
    propList[mu] = X -> lam(X)*X[sidx("TRP")]
    
    ### HKP ---------------------------------------------------------------
    
    # + Vq:  RIB_mHKP+ATP --> RIB+mHKP+HKP (described)
    
    mu += 1 # dilution: HKP -> *
    nu[sidx("HKP"), mu] = -1
    propList[mu] = X -> lam(X)*X[sidx("HKP")]
    
    
    ### RIB ---------------------------------------------------------------
    
    # Vr:  RIB_mRIB+ATP --> RIB+mRIB+RIB (described)
    # Vm + Vt + Vq + Vr (described)
    # ku * (RIB_mENZ+RIB_mTRP+RIB_mHKP+RIB_mRIB) (described)
    # kb * RIB * (mENZ+mTRP+mHKP+mRIB) (described)
    
    mu += 1 # dilution: RIB -> *
    nu[sidx("RIB"), mu] = -1
    propList[mu] = X -> lam(X)*X[sidx("RIB")]
    
    # combine all function handles and enforce positivity of propensities
    propFun(X) = [max(f(X), 0) for f in propList]
    return propFun, nu, propList # HH dont need to return propList? 

 end


#  ##### ============ Usage =================
#  (propFun, nu, propList) = defineStochModel(par, I)
# propFun(X0)