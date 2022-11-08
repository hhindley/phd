function odemodel!(dydt, initial, params, t)
    b, dm, kb, ku, f, thetar, k_cm, s0, gmax, cl, thetax, Kt, M, we, Km, vm, nx, Kq, Kp, vt, wr, wq, wp, hq, nr, ns = params
    cr, em, cp, cq, ct, et, cm, mt, mm, q, p, si, mq, mp, mr, r, a = initial

    dcr, dem, dcp, dcq, dct, det, dcm, dmt, dmm, dq, dp, dsi, dmq, dmp, dmr, dr, da = zeros(length(dydt))
    
    Kgamma = gmax/Kp
    gamma = gmax*a/(Kgamma+a)
    ttrate = (cq+cr+cp+ct+cm)gamma
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
    dydt[13] =   +(wq*a/(thetax+a)/(1+(q/Kq)^hq)) - mq*dm - mq*lam - r*mq*kb + cq*ku + gamma/nx*cq
    dydt[14] = +(wp*a/(thetax+a)) - mp*dm - mp*lam - r*mp*kb + cp*ku + gamma/nx*cp
    dydt[15] =   +(wr*a/(thetar+a)) - mr*dm - mr*lam - r*mr*kb + cr*ku + gamma/nr*cr
    dydt[16] =  - lam*r - r*mr*kb - r*mt*kb - r*mm*kb - r*mq*kb - r*mp*kb + cr*ku + ct*ku + cm*ku + cq*ku + cp*ku + gamma/nr*cr + gamma/nr*cr + gamma/nx*ct + gamma/nx*cm + gamma/nx*cq + gamma/nx*cp
    dydt[17] =  +ns*nucat - ttrate - lam*a
end

function odemodelfull!(dydt, initial, params, t)
    b, dm, kb, ku, cl, k_cm, f, thetar, s0, gmax, thetax, Kt, M, we, Km, vm, nx, Kq, Kp, vt, wr, wq, wp, hq, nr, ns = params
    cr, em, cp, cq, ct, et, cm, zmr, zmp, zmq, zmt, zmm, mt, mm, q, p, si, mq, mp, mr, r, a = initial

    dcr, dem, dcp, dcq, dct, det, dcm, dzmr, dzmp, dzmq, dzmt, dzmm, dmt, dmm, dq, dp, dsi, dmq, dmp, dmr, dr, da = zeros(length(dydt))

    Kgamma = gmax/Kp
    gamma = gmax*a/(Kgamma+a)
    ttrate = (cq+cr+cp+ct+cm)gamma
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
    dydt[8] = +f*cr - b*zmr - lam*zmr
    dydt[9] = +f*cp - b*zmp - lam*zmp
    dydt[10] = +f*cq - b*zmq - lam*zmq
    dydt[11] = +f*ct - b*zmt - lam*zmt
    dydt[12] = +f*cm - b*zmm - lam*zmm
    dydt[13] =   +(we*a/(thetax+a)) - mt*dm - mt*lam - r*mt*kb + ct*ku + gamma/nx*ct
    dydt[14] =   +(we*a/(thetax+a)) - mm*dm - mm*lam - r*mm*kb + cm*ku + gamma/nx*cm
    dydt[15] = - lam*q + gamma/nx*cq
    dydt[16] =  - lam*p + gamma/nx*cp
    dydt[17] = - lam*si + vimp - nucat
    dydt[18] =   +(wq*a/(thetax+a)/(1+(q/Kq)^hq)) - mq*dm - mq*lam - r*mq*kb + cq*ku + gamma/nx*cq
    dydt[19] = +(wp*a/(thetax+a)) - mp*dm - mp*lam - r*mp*kb + cp*ku + gamma/nx*cp
    dydt[20] =   +(wr*a/(thetar+a)) - mr*dm - mr*lam - r*mr*kb + cr*ku + gamma/nr*cr
    dydt[21] =  - lam*r - r*mr*kb - r*mt*kb - r*mm*kb - r*mq*kb - r*mp*kb + cr*ku + ct*ku + cm*ku + cq*ku + cp*ku + gamma/nr*cr + gamma/nr*cr + gamma/nx*ct + gamma/nx*cm + gamma/nx*cq + gamma/nx*cp
    dydt[22] =  +ns*nucat - ttrate - lam*a

end

function pop_odes!(dYdt, initial, params, t)
    vimp, lam = params
    N, s = initial
 
    dN, ds = zeros(length(dYdt))

    dYdt[1] = + lam*N
    dYdt[2] = - vimp*N #- s

end


function popodemodeldiff!(dydt, initial, params, t)
    b, dm, kb, ku, cl, k_cm, f, thetar, gmax, thetax, Kt, M, we, Km, vm, nx, Kq, Kp, vt, wr, wq, wp, hq, nr, ns = params
    N, s0, cr, em, cp, cq, ct, et, cm, zmr, zmp, zmq, zmt, zmm, mt, mm, q, p, si, mq, mp, mr, r, a = initial

    dN, ds0, dcr, dem, dcp, dcq, dct, det, dcm, dzmr, dzmp, dzmq, dzmt, dzmm, dmt, dmm, dq, dp, dsi, dmq, dmp, dmr, dr, da = zeros(length(dydt))

    Kgamma = 3e8
    gamma = gmax*a/(Kgamma+a)
    ttrate = (cq+cr+cp+ct+cm)gamma
    lam = ttrate/M
    vimp = (et*vt*s0/(Kt+s0))
    nucat = (em*vm*si/(Km+si))

    dydt[1] = +lam*N
    dydt[2] = -vimp*N - s0
    dydt[3] =  +r*mr*kb - cr*ku - cr*lam - gamma/nr*cr - f*cr
    dydt[4] = - lam*em + cm*gamma/nx
    dydt[5] =   +r*mp*kb - cp*ku - cp*lam - gamma/nx*cp - f*cp
    dydt[6] =  +r*mq*kb - cq*ku - cq*lam - gamma/nx*cq - f*cq
    dydt[7] =  +r*mt*kb - ct*ku - ct*lam - gamma/nx*ct - f*ct
    dydt[8] =  - lam*et + gamma/nx*ct
    dydt[9] =   +r*mm*kb - cm*ku - cm*lam - gamma/nx*cm - f*cm
    dydt[10] = +f*cr - b*zmr - lam*zmr
    dydt[11] = +f*cp - b*zmp - lam*zmp
    dydt[12] = +f*cq - b*zmq - lam*zmq
    dydt[13] = +f*ct - b*zmt - lam*zmt
    dydt[14] = +f*cm - b*zmm - lam*zmm
    dydt[15] =   +(we*a/(thetax+a)) - mt*dm - mt*lam - r*mt*kb + ct*ku + gamma/nx*ct
    dydt[16] =   +(we*a/(thetax+a)) - mm*dm - mm*lam - r*mm*kb + cm*ku + gamma/nx*cm
    dydt[17] = - lam*q + gamma/nx*cq
    dydt[18] =  - lam*p + gamma/nx*cp
    dydt[19] = - lam*si + vimp - nucat
    dydt[20] =   +(wq*a/(thetax+a)/(1+(q/Kq)^hq)) - mq*dm - mq*lam - r*mq*kb + cq*ku + gamma/nx*cq
    dydt[21] = +(wp*a/(thetax+a)) - mp*dm - mp*lam - r*mp*kb + cp*ku + gamma/nx*cp
    dydt[22] =   +(wr*a/(thetar+a)) - mr*dm - mr*lam - r*mr*kb + cr*ku + gamma/nr*cr
    dydt[23] =  - lam*r - r*mr*kb - r* mt*kb - r*mm*kb - r*mq*kb - r*mp*kb + cr*ku + ct*ku + cm*ku + cq*ku + cp*ku + gamma/nr*cr + gamma/nr*cr + gamma/nx*ct + gamma/nx*cm + gamma/nx*cq + gamma/nx*cp
    dydt[24] =  +ns*nucat - ttrate - lam*a


end

function popodemodelfull!(dydt, initial, params, t)
    b, dm, kb, ku, cl, k_cm, f, thetar, gmax, thetax, Kt, M, we, Km, vm, nx, Kq, Kp, vt, wr, wq, wp, hq, nr, ns = params
    cr, em, cp, cq, ct, et, cm, zmr, zmp, zmq, zmt, zmm, mt, mm, q, p, si, mq, mp, mr, r, a, N, s0 = initial

    dcr, dem, dcp, dcq, dct, det, dcm, dzmr, dzmp, dzmq, dzmt, dzmm, dmt, dmm, dq, dp, dsi, dmq, dmp, dmr, dr, da, dN, ds0 = zeros(length(dydt))

    Kgamma = 3e8
    gamma = gmax*a/(Kgamma+a)
    ttrate = (cq+cr+cp+ct+cm)gamma
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
    dydt[8] = +f*cr - b*zmr - lam*zmr
    dydt[9] = +f*cp - b*zmp - lam*zmp
    dydt[10] = +f*cq - b*zmq - lam*zmq
    dydt[11] = +f*ct - b*zmt - lam*zmt
    dydt[12] = +f*cm - b*zmm - lam*zmm
    dydt[13] =   +(we*a/(thetax+a)) - mt*dm - mt*lam - r*mt*kb + ct*ku + gamma/nx*ct
    dydt[14] =   +(we*a/(thetax+a)) - mm*dm - mm*lam - r*mm*kb + cm*ku + gamma/nx*cm
    dydt[15] = - lam*q + gamma/nx*cq
    dydt[16] =  - lam*p + gamma/nx*cp
    dydt[17] = - lam*si + vimp - nucat
    dydt[18] =   +(wq*a/(thetax+a)/(1+(q/Kq)^hq)) - mq*dm - mq*lam - r*mq*kb + cq*ku + gamma/nx*cq
    dydt[19] = +(wp*a/(thetax+a)) - mp*dm - mp*lam - r*mp*kb + cp*ku + gamma/nx*cp
    dydt[20] =   +(wr*a/(thetar+a)) - mr*dm - mr*lam - r*mr*kb + cr*ku + gamma/nr*cr
    dydt[21] =  - lam*r - r*mr*kb - r* mt*kb - r*mm*kb - r*mq*kb - r*mp*kb + cr*ku + ct*ku + cm*ku + cq*ku + cp*ku + gamma/nr*cr + gamma/nr*cr + gamma/nx*ct + gamma/nx*cm + gamma/nx*cq + gamma/nx*cp
    dydt[22] =  +ns*nucat - ttrate - lam*a
    dydt[23] = +lam*N
    dydt[24] = - vimp*N - s0

end


function odeHeteroExpress!(dydt, initial, params, t)
    b, dm, kb, ku, cl, k_cm, f, thetar, s0, gmax, thetax, Kt, M, we, Km, vm, nx, Kq, Kp, vt, wr, wq, wp, hq, nr, ns = params
    cr, em, cp, cq, ct, et, cm, zmr, zmp, zmq, zmt, zmm, mt, mm, q, p, si, mq, mp, mr, r, a = initial

    dcr, dem, dcp, dcq, dct, det, dcm, dzmr, dzmp, dzmq, dzmt, dzmm, dmt, dmm, dq, dp, dsi, dmq, dmp, dmr, dr, da, dmg1, dmg2, dmg3, dcg1, dcg2, dcg3, dg1, dg2, dg3 = zeros(length(dydt))

    Kgamma = gmax/Kp
    gamma = gmax*a/(Kgamma+a)
    ttrate = (cq+cr+cp+ct+cm)gamma
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
    dydt[8] = +f*cr - b*zmr - lam*zmr
    dydt[9] = +f*cp - b*zmp - lam*zmp
    dydt[10] = +f*cq - b*zmq - lam*zmq
    dydt[11] = +f*ct - b*zmt - lam*zmt
    dydt[12] = +f*cm - b*zmm - lam*zmm
    dydt[13] =   +(we*a/(thetax+a)) - mt*dm - mt*lam - r*mt*kb + ct*ku + gamma/nx*ct
    dydt[14] =   +(we*a/(thetax+a)) - mm*dm - mm*lam - r*mm*kb + cm*ku + gamma/nx*cm
    dydt[15] = - lam*q + gamma/nx*cq
    dydt[16] =  - lam*p + gamma/nx*cp
    dydt[17] = - lam*si + vimp - nucat
    dydt[18] =   +(wq*a/(thetax+a)/(1+(q/Kq)^hq)) - mq*dm - mq*lam - r*mq*kb + cq*ku + gamma/nx*cq
    dydt[19] = +(wp*a/(thetax+a)) - mp*dm - mp*lam - r*mp*kb + cp*ku + gamma/nx*cp
    dydt[20] =   +(wr*a/(thetar+a)) - mr*dm - mr*lam - r*mr*kb + cr*ku + gamma/nr*cr
    dydt[21] =  - lam*r - r*mr*kb - r*mt*kb - r*mm*kb - r*mq*kb - r*mp*kb + cr*ku + ct*ku + cm*ku + cq*ku + cp*ku + gamma/nr*cr + gamma/nr*cr + gamma/nx*ct + gamma/nx*cm + gamma/nx*cq + gamma/nx*cp
    dydt[22] =  +ns*nucat - ttrate - lam*a
    dydt[23] = +(wg*a/(thetax+a))/(1+(g1/Kq)^h) + ku*cg1 + (gamma*a/ng)*cg1 - kb*r*mg1 - dmg*mg1 - lam*mg1 
    dydt[24] = +(wg*a/(thetax+a))/(1+(g2/Kq)^h) + ku*cg2 + (gamma*a/ng)*cg2 - kb*r*mg2 - dmg*mg2 - lam*mg2
    dydt[25] = +(wg*a/(thetax+a))/(1+(g3/Kq)^h) + ku*cg3 + (gamma*a/ng)*cg3 - kb*r*mg3 - dmg*mg3 - lam*mg3 
    dydt[26] = +kb*r*mg1 - (gamma*a/ng)*cg3 - ku*cg1 - lam*cg1
    dydt[27] = +kb*r*mg2 - (gamma*a/ng)*cg3 - ku*cg2 - lam*cg2
    dydt[28] = +kb*r*mg3 - (gamma*a/ng)*cg3 - ku*cg3 - lam*cg3
    dydt[29] = +(gamma*a/ng)*cg1 - lam*g1 - dg*g1
    dydt[30] = +(gamma*a/ng)*cg2 - lam*g2 - dg*g2
    dydt[31] = +(gamma*a/ng)*cg3 - lam*g3 - dg*g3  
end