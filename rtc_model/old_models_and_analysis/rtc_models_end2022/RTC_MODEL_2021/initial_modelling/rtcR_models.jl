function full_model!(dydt, initial, params, t) 
    k1, k2, k3, k4, katp = params 

    a, b, c, d, e, p = initial 

    dadt, dbdt, dcdt, dddt, dedt, dpdt = zeros(length(dydt))

    dydt[1] = -k1*a*b^6 + k2*c 
    dydt[2] = -6k1*a*b^6 + 6*k2*c
    dydt[3] = k1*a*b^6 - k2*c - k3*c*d + k4*e + katp*e 
    dydt[4] = -k3*c*d + k4*e 
    dydt[5] = k3*c*d - k4*e - katp*e 
    dydt[6] = katp*e    
end


function hill!(dydt, initial, params, t) 
    k1, k2, k3, k4, katp = params 

    a, b, c, d, e, p = initial 

    dadt, dbdt, dcdt, dddt, dedt, dpdt = zeros(length(dydt))
    
    ctot =  c+e_0#(a0 + b0)/7 #a0+b0+7d_0 #(a0 + b0)/7 #a0+c_0+e_0 
    atot = a_0+c_0+e_0
    hmax = k1*atot
    kb =  k1/k2 #c/(a*b^6)
    hill = hmax*kb*b^6/(1+kb*b^6)

    dydt[1] = -k1*a*b^6 + k2*c 
    dydt[2] = -6k1*a*b^6 + 6*k2*c
    dydt[3] = hill - k2*c - k3*c*d + k4*e + katp*e 
    dydt[4] = -k3*c*d + k4*e 
    dydt[5] = k3*c*d - k4*e - katp*e 
    dydt[6] = katp*e    
end

function mass_cons_a!(dydt, initial, params, t)
    
    k1, k2, k3, k4, katp, ctot, a_0, b_0 = params 

    b, c, d, e, p = initial 

    dbdt, dcdt, dddt, dedt, dpdt = zeros(length(dydt))

    a = (6*a_0-b_0+b)/6

    dydt[1] = -6k1*a*b^6 + 6*k2*c
    dydt[2] = k1*a*b^6 - k2*c - k3*c*d + k4*e + katp*e 
    dydt[3] = -k3*c*d + k4*e 
    dydt[4] = k3*c*d - k4*e - katp*e 
    dydt[5] = katp*e    
end

function mass_cons_a_with_hill!(dydt, initial, params, t)
    
    k1, k2, k3, k4, katp, ctot, a_0, b_0 = params 

    b, c, d, e, p = initial 

    dbdt, dcdt, dddt, dedt, dpdt = zeros(length(dydt))

    a = (6*a_0-b_0+b)/6
    kb =  c/(a*b^6)
    ctot =  a_0+c_0+e_0#(a0 + b0)/7 #a0+b0+7d_0 #(a0 + b0)/7 #a0+c_0+e_0 
    atot = a_0+c_0+e_0
    hmax = k1*ctot
    hill = hmax*kb*b^6/(1+kb*b^6)

    dydt[1] = -6k1*a*b^6 + 6*k2*c
    dydt[2] = hill - k2*c - k3*c*d + k4*e + katp*e 
    dydt[3] = -k3*c*d + k4*e 
    dydt[4] = k3*c*d - k4*e - katp*e 
    dydt[5] = katp*e    
end

function mass_cons_hill_mm!(dydt, initial, params, t) # BEST 
    
    k1, k2, k3, k4, katp, ctot, a_0, b_0 = params 

    b, c, d, p = initial 

    dbdt, dcdt, dddt, dpdt = zeros(length(dydt))

    ctot = c+e_0 #a0+b0+7d_0 #(a0 + b0)/7 #a0+c_0+e_0 
    a = (6*a_0-b_0+b)/6
    e = (k3*d*ctot/(k4+katp+(k3*d)))
    atot = a_0+c_0+e_0
    hmax = k1*atot
    kb = k1/k2 # c/(a*b^6)
    hill = hmax*kb*b^6/(1+kb*b^6)

    dydt[1] = -6k1*a*b^6 + 6*k2*c
    dydt[2] = hill - k2*c - k3*c*d + k4*e + katp*e 
    dydt[3] = -k3*c*d + k4*e  
    dydt[4] = katp*e    
end

function simple_odes_mm_e_mass_cons_a!(dydt, initial, params, t) 
    k1, k2, k3, k4, katp = params 

    b, c, d, p = initial 

    dbdt, dcdt, dddt, dpdt = zeros(length(dydt))

    a = (6*a_0-b_0+b)/6
    e = (k3*d*c #=c here should be ctot=#/(k4+katp+(k3*d))) # ctot = c because ctot = c0+e0 but that would be 0 so ends up being just c for initial conc of c because in this scenario e isnt there yet 

    dydt[1] = -6k1*a*b^6 + 6*k2*c
    dydt[2] = k1*a*b^6 - k2*c - k3*c*d + k4*e + katp*e 
    dydt[3] = -k3*c*d + k4*e 
    dydt[4] = katp*e    
end



function comparison_rtc!(dydt, initial, params, t) # BEST 
    
    k1, k2, k3, k4, katp, a_0, b_0 = params 

    b, c, d, p = initial 

    dbdt, dcdt, dddt, dpdt = zeros(length(dydt))

    ctot = c+e_0 #a0+b0+7d_0 #(a0 + b0)/7 #a0+c_0+e_0 
    a =  (6*a_0-b_0+b)/6
    e = (k3*d*ctot/(k4+katp+(k3*d)))
    atot = a_0+c_0+e_0
    hmax = k1*atot
    kb = k1/k2 # c/(a*b^6)
    hill = hmax*kb*b^6/(1+kb*b^6)

    dydt[1] = -6*k1*a*b^6 + 6*k2*c
    dydt[2] = hill - k2*c - k3*c*d + k4*e + katp*e 
    dydt[3] = -k3*c*d + k4*e  
    dydt[4] = katp*e    
end