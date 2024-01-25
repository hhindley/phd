function growth_model(dydt, initial, params, t)
    dm, kb, ku, f, thetar, s0, gmax, thetax, Kt, M, we, Km, vm, nx, Kq, vt, wr, wq, nq, nr, ns, Kgamma = params
    cr, em, cq, ct, et, cm, mt, mm, q, si, mq, mr, r, a = initial

    dcr, dem, dcq, dct, det, dcm, dmt, dmm, dq, dsi, dmq, dmr, dr, da = zeros(length(dydt))
    
    # Kgamma = 7*sf #gmax/Kp
    gamma = gmax*a/(Kgamma+a) # aa min-1 molec-1 (aa min-1 molec-1(in terms of uM so this does cancel out below))
    ttrate = (cq+cr+ct+cm)*gamma # aa min-1 (aa min-1)
    lam = ttrate/M # min-1 (min-1)
    vimp = (et*vt*s0/(Kt+s0)) # molec min-1 (uM min-1)
    nucat = (em*vm*si/(Km+si)) # molec min-1 (uM min-1)

    dydt[1] =  +r*mr*kb - cr*ku - cr*lam - gamma/nr*cr - f*cr # molec min-1
    dydt[2] = - lam*em + cm*gamma/nx # molec min-1
    dydt[3] =  +r*mq*kb - cq*ku - cq*lam - gamma/nx*cq - f*cq # molec min-1
    dydt[4] =  +r*mt*kb - ct*ku - ct*lam - gamma/nx*ct - f*ct # molec min-1
    dydt[5] =  - lam*et + gamma/nx*ct # molec min-1
    dydt[6] =   +r*mm*kb - cm*ku - cm*lam - gamma/nx*cm - f*cm # molec min-1
    dydt[7] =   +(we*a/(thetax+a)) - mt*dm - mt*lam - r*mt*kb + ct*ku + gamma/nx*ct # molec min-1
    dydt[8] =   +(we*a/(thetax+a)) - mm*dm - mm*lam - r*mm*kb + cm*ku + gamma/nx*cm # molec min-1
    dydt[9] = - lam*q + gamma/nx*cq # molec min-1
    dydt[10] = - lam*si + vimp - nucat # molec min-1
    dydt[11] =   +(wq*a/(thetax+a)/(1+(q/Kq)^nq)) - mq*dm - mq*lam - r*mq*kb + cq*ku + gamma/nx*cq # molec min-1
    dydt[12] =   +(wr*a/(thetar+a)) - mr*dm - mr*lam - r*mr*kb + cr*ku + gamma/nr*cr # molec min-1
    dydt[13] =  - lam*r - r*mr*kb - r*mt*kb - r*mm*kb - r*mq*kb + cr*ku + ct*ku + cm*ku + cq*ku + gamma/nr*cr + gamma/nr*cr + gamma/nx*ct + gamma/nx*cm + gamma/nx*cq # molec min-1
    dydt[14] =  +ns*nucat - ttrate - lam*a ### aa min-1
end


function odemodel!(dydt, initial, params, t)
    b, dm, kb, ku, f, thetar, k_cm, s0, gmax, cl, thetax, Kt, M, we, Km, vm, nx, Kq, Kp, vt, wr, wq, wp, nq, nr, ns, Kgamma = params
    cr, em, cp, cq, ct, et, cm, mt, mm, q, p, si, mq, mp, mr, r, a = initial

    dcr, dem, dcp, dcq, dct, det, dcm, dmt, dmm, dq, dp, dsi, dmq, dmp, dmr, dr, da = zeros(length(dydt))
    
    # Kgamma = 7*sf #gmax/Kp
    gamma = gmax*a/(Kgamma+a)
    ttrate = (cq+cr+cp+ct+cm)*gamma
    lam = ttrate/M
    vimp = (et*vt*s0/(Kt+s0))
    nucat = (em*vm*si/(Km+si))

    dydt[1] =  +r*mr*kb - cr*ku - cr*lam - gamma/nr*cr - f*cr
    dydt[2] = - lam*em + cm*gamma/nx
    dydt[3] =   +r*mp*kb - cp*ku - cp*lam - gamma/nx*cp - f*cp
    dydt[4] =  +r*mq*kb - cq*ku - cq*lam - gamma/nx*cq - f*cq
    dydt[5] =  +r*mt*kb - ct*ku - ct*lam - gamma/nx*ct - f*ct
    dydt[6] =  - lam*et + gamma/nx*ct
    dydt[7] =   +r*mm*kb - cm*ku - cm*lam - gamma/nx*cm - f*cm
    dydt[8] =   +(we*a/(thetax+a)) - mt*dm - mt*lam - r*mt*kb + ct*ku + gamma/nx*ct
    dydt[9] =   +(we*a/(thetax+a)) - mm*dm - mm*lam - r*mm*kb + cm*ku + gamma/nx*cm
    dydt[10] = - lam*q + gamma/nx*cq
    dydt[11] =  - lam*p + gamma/nx*cp
    dydt[12] = - lam*si + vimp - nucat
    dydt[13] =   +(wq*a/(thetax+a)/(1+(q/Kq)^nq)) - mq*dm - mq*lam - r*mq*kb + cq*ku + gamma/nx*cq
    dydt[14] = +(wp*a/(thetax+a)) - mp*dm - mp*lam - r*mp*kb + cp*ku + gamma/nx*cp
    dydt[15] =   +(wr*a/(thetar+a)) - mr*dm - mr*lam - r*mr*kb + cr*ku + gamma/nr*cr
    dydt[16] =  - lam*r - r*mr*kb - r*mt*kb - r*mm*kb - r*mq*kb - r*mp*kb + cr*ku + ct*ku + cm*ku + cq*ku + cp*ku + gamma/nr*cr + gamma/nr*cr + gamma/nx*ct + gamma/nx*cm + gamma/nx*cq + gamma/nx*cp
    dydt[17] =  +ns*nucat - ttrate - lam*a
end

function odemodelfull!(dydt, initial, params, t)
    b, dm, kb, ku, f, thetar, k_cm, s0, gmax, cl, thetax, Kt, M, we, Km, vm, nx, Kq, Kp, vt, wr, wq, wp, hq, nr, ns = params
    cr, em, cp, cq, ct, et, cm, mt, mm, q, p, si, mq, mp, mr, r, a, zmr, zmp, zmq, zmt, zmm = initial

    dcr, dem, dcp, dcq, dct, det, dcm, dmt, dmm, dq, dp, dsi, dmq, dmp, dmr, dr, da, dzmr, dzmp, dzmq, dzmt, dzmm = zeros(length(dydt))

    Kgamma = gmax/Kp
    gamma = gmax*a/(Kgamma+a)
    ttrate = (cq+cr+cp+ct+cm)*gamma
    lam = ttrate/M
    vimp = (et*vt*s0)/(Kt+s0)
    nucat = (em*vm*si)/(Km+si)

    dydt[1] =  +r*mr*kb - cr*ku - cr*lam - cr*gamma/nr - f*cr + b*zmr
    dydt[2] = - lam*em + cm*gamma/nx
    dydt[3] =   +r*mp*kb - cp*ku - cp*lam - cp*gamma/nx - f*cp + b*zmp
    dydt[4] =  +r*mq*kb - cq*ku - cq*lam - cq*gamma/nx - f*cq + b*zmq
    dydt[5] =  +r*mt*kb - ct*ku - ct*lam - ct*gamma/nx - f*ct + b*zmt
    dydt[6] =  - lam*et + ct*gamma/nx
    dydt[7] =   +r*mm*kb - cm*ku - cm*lam - cm*gamma/nx - f*cm + b*zmm
    dydt[8] =   +we*a/(thetax+a) - mt*lam - mt*dm - r*mt*kb + ct*ku + ct*gamma/nx
    dydt[9] =   +we*a/(thetax+a) - mm*lam - mt*dm - r*mm*kb + cm*ku + cm*gamma/nx
    dydt[10] = - lam*q + (cq*gamma)/nx
    dydt[11] =  - lam*p + (cp*gamma)/nx
    dydt[12] = - lam*si + vimp - nucat
    dydt[13] =   +wq*a/(thetax+a)/1+(q/Kq)^hq - mq*lam- mt*dm - r*mq*kb + cq*ku + cq*gamma/nx
    dydt[14] = +wp*a/(thetax+a) - mp*lam- mt*dm - r*mp*kb + cp*ku + cp*gamma/nx
    dydt[15] =   +wr*a/(thetar+a) - mr*lam- mt*dm- r*mr*kb + cr*ku + (cr*gamma)/nr
    dydt[16] =  - lam*r - r*mr*kb - r*mt*kb - r*mm*kb - r*mq*kb - r*mp*kb + cr*ku + ct*ku + cm*ku + cq*ku + cp*ku + cr*gamma/nr + cr*gamma/nr + ct*gamma/nx + cm*gamma/nx + cq*gamma/nx + cp*gamma/nx
    dydt[17] =  +ns*nucat - ttrate - lam*a
    dydt[18] = +f*cr-b*zmr-lam*zmr
    dydt[19] = +f*cp-b*zmp-lam*zmp
    dydt[20] = +f*cq-b*zmq-lam*zmq
    dydt[21] = +f*ct-b*zmt-lam*zmt
    dydt[22] = +f*cm-b*zmm-lam*zmm
end

function pop_model!(dydt, initial, params, t)
    dm, kb, ku, f, thetar, gmax, thetax, Kt, M, we, Km, vm, nx, Kq, vt, wr, wq, wp, nq, nr, ns, kin, d_s, d_n, Kgamma = params
    cr, em, cp, cq, ct, et, cm, mt, mm, q, p, si, mq, mp, mr, r, a, s, N = initial

    dcr, dem, dcp, dcq, dct, det, dcm, dmt, dmm, dq, dp, dsi, dmq, dmp, dmr, dr, da, ds, dN = zeros(length(dydt))
    
    # Kgamma = 7*sf;#gmax/Kp
    gamma = gmax*a/(Kgamma+a)
    ttrate = (cq+cr+cp+ct+cm)*gamma
    lam = ttrate/M
    vimp = (et*vt*s/(Kt+s))
    nucat = (em*vm*si/(Km+si))

    dydt[1] =  +r*mr*kb - cr*ku - cr*lam - gamma/nr*cr - f*cr
    dydt[2] = - lam*em + cm*gamma/nx
    dydt[3] =   +r*mp*kb - cp*ku - cp*lam - gamma/nx*cp - f*cp
    dydt[4] =  +r*mq*kb - cq*ku - cq*lam - gamma/nx*cq - f*cq
    dydt[5] =  +r*mt*kb - ct*ku - ct*lam - gamma/nx*ct - f*ct
    dydt[6] =  - lam*et + gamma/nx*ct
    dydt[7] =   +r*mm*kb - cm*ku - cm*lam - gamma/nx*cm - f*cm
    dydt[8] =   +(we*a/(thetax+a)) - mt*dm - mt*lam - r*mt*kb + ct*ku + gamma/nx*ct
    dydt[9] =   +(we*a/(thetax+a)) - mm*dm - mm*lam - r*mm*kb + cm*ku + gamma/nx*cm
    dydt[10] = - lam*q + gamma/nx*cq
    dydt[11] =  - lam*p + gamma/nx*cp
    dydt[12] = - lam*si + vimp - nucat
    dydt[13] =   +(wq*a/(thetax+a)/(1+(q/Kq)^nq)) - mq*dm - mq*lam - r*mq*kb + cq*ku + gamma/nx*cq
    dydt[14] = +(wp*a/(thetax+a)) - mp*dm - mp*lam - r*mp*kb + cp*ku + gamma/nx*cp
    dydt[15] =   +(wr*a/(thetar+a)) - mr*dm - mr*lam - r*mr*kb + cr*ku + gamma/nr*cr
    dydt[16] =  - lam*r - r*mr*kb - r*mt*kb - r*mm*kb - r*mq*kb - r*mp*kb + cr*ku + ct*ku + cm*ku + cq*ku + cp*ku + gamma/nr*cr + gamma/nr*cr + gamma/nx*ct + gamma/nx*cm + gamma/nx*cq + gamma/nx*cp
    dydt[17] =  +ns*nucat - ttrate - lam*a
    dydt[18] = kin - vimp*N - d_s*s
    dydt[19] = lam*N - d_n*N
end




