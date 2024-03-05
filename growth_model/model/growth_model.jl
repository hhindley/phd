include("/home/holliehindley/phd/general_funcs/all_model_funcs.jl")

indexof(sym, syms) = findfirst(isequal(sym),syms)

@variables t 
@parameters d kb ku thetar s0 gmax thetax Kt M we Km vm nx Kq vt wr wq nq nr ns Kgamma abx kon
species_gm1 = @syms mr(t) mt(t) mm(t) mq(t) cr(t) ct(t) cm(t) cq(t) et(t) em(t) q(t) r(t) zmr(t) zmt(t) zmm(t) zmq(t) si(t) a(t) 

species_gm = [Symbol(i) for i in species_gm1]


D = Differential(t)

@mtkmodel GROWTH_MODEL begin
    @parameters begin
        d
        kb
        ku
        thetar
        s0
        gmax
        thetax
        Kt
        M
        we
        Km
        vm
        nx
        Kq
        vt
        wr
        wq
        nq
        nr
        ns
        Kgamma
        abx
        kon
    end
    @variables begin
        mr(t) 
        mt(t) 
        mm(t) 
        mq(t) 
        cr(t) 
        ct(t) 
        cm(t) 
        cq(t) 
        et(t) 
        em(t) 
        q(t) 
        r(t) 
        zmr(t) 
        zmt(t) 
        zmm(t) 
        zmq(t) 
        si(t) 
        a(t) 

        rhs_mr(t) 
        rhs_mt(t) 
        rhs_mm(t) 
        rhs_mq(t) 
        rhs_cr(t) 
        rhs_ct(t) 
        rhs_cm(t) 
        rhs_cq(t) 
        rhs_et(t) 
        rhs_em(t) 
        rhs_q(t) 
        rhs_r(t) 
        rhs_zmr(t) 
        rhs_zmt(t) 
        rhs_zmm(t) 
        rhs_zmq(t) 
        rhs_si(t) 
        rhs_a(t) 
        
        gamma(t)
        ttrate(t)
        lam(t)
        vimp(t)
        nucat(t)

    end

    @equations begin
        gamma ~ gmax*a/(Kgamma+a)
        ttrate ~ (cr + ct + cm + cq)*gamma
        lam ~ ttrate/M
        vimp ~ (et*vt*s0/(Kt+s0))
        nucat ~ (em*vm*si/(Km+si))

        rhs_mt ~   +ω_p(we, thetax, a) - deg(mt) - dil(mt,lam) - rh_bind(mt,r) + rh_unbind(ct) + v_x(ct,nx,gamma)
        rhs_mm ~   +ω_p(we, thetax, a) - deg(mm) - dil(mm,lam) - rh_bind(mm,r) + rh_unbind(cm) + v_x(cm,nx,gamma)
        rhs_mq ~   +ω_q(wq, thetax, a, q) - deg(mq) - dil(mq,lam) - rh_bind(mq,r) + rh_unbind(cq) + v_x(cq,nx,gamma)
        rhs_mr ~   +ω_p(wr, thetar, a) - deg(mr) - dil(mr,lam) - rh_bind(mr,r) + rh_unbind(cr) + v_x(cr,nr,gamma)

        rhs_cr ~  +rh_bind(mr,r) - rh_unbind(cr) - dil(cr,lam) - v_x(cr,nr,gamma) - zm(cr)
        rhs_cq ~  +rh_bind(mq,r) - rh_unbind(cq) - dil(cq,lam) - v_x(cq,nx,gamma) - zm(cq)
        rhs_ct ~  +rh_bind(mt,r) - rh_unbind(ct) - dil(ct,lam) - v_x(ct,nx,gamma) - zm(ct)
        rhs_cm ~   +rh_bind(mm,r) - rh_unbind(cm) - dil(cm,lam) - v_x(cm,nx,gamma) - zm(cm)


        rhs_em ~ - dil(em,lam) + v_x(cm,nx,gamma)
        rhs_et ~  - dil(et,lam) + v_x(ct,nx,gamma)
        rhs_q ~ - dil(q,lam) + v_x(cq,nx,gamma)
        rhs_r ~  - dil(r,lam) - rh_bind(mr,r) - rh_bind(mt,r) - rh_bind(mm,r) - rh_bind(mq,r) + rh_unbind(cr) + rh_unbind(ct) + rh_unbind(cm) + rh_unbind(cq) + v_x(cr,nr,gamma) + v_x(cr,nr,gamma) + v_x(ct,nx,gamma) + v_x(cm,nx,gamma) + v_x(cq,nx,gamma)

        rhs_si ~ - dil(si,lam) + vimp - nucat
        rhs_a ~  +ns*nucat - ttrate - dil(a,lam)

        rhs_zmr ~  zm(cr) - dil(zmr,lam)
        rhs_zmq ~ zm(cq) - dil(zmq,lam)
        rhs_zmt ~ zm(ct) - dil(zmt,lam)
        rhs_zmm ~ zm(cm) - dil(zmm,lam)
        
        D(mr) ~ rhs_mr
        D(mt) ~ rhs_mt
        D(mm) ~ rhs_mm
        D(mq) ~ rhs_mq

        D(cr) ~ rhs_cr
        D(ct) ~ rhs_ct
        D(cm) ~ rhs_cm
        D(cq) ~ rhs_cq

        D(et) ~ rhs_et
        D(em) ~ rhs_em
        D(q) ~ rhs_q
        D(r) ~ rhs_r

        D(zmr) ~ rhs_zmr
        D(zmt) ~ rhs_zmt
        D(zmm) ~ rhs_zmm
        D(zmq) ~ rhs_zmq

        D(si) ~ rhs_si
        D(a) ~ rhs_a


    end
end

@mtkbuild growth_model = GROWTH_MODEL()

init_gm = [growth_model.mr=>0.0, growth_model.mt=>0.0, growth_model.mm=>0.0, growth_model.mq=>0.0, 
            growth_model.cr=>0.0, growth_model.ct=>0.0, growth_model.cm=>0.0, growth_model.cq=>0.0, 
            growth_model.et=>0.0, growth_model.em=>0.0, growth_model.q=>0.0, growth_model.r=>10.0,
            growth_model.zmr=>0.0, growth_model.zmt=>0.0, growth_model.zmm=>0.0, growth_model.zmq=>0.0,
            growth_model.si=>0.0, growth_model.a=>1000.0]

init_gm_uM = [growth_model.mr=>0.0, growth_model.mt=>0.0, growth_model.mm=>0.0, growth_model.mq=>0.0, 
growth_model.cr=>0.0, growth_model.ct=>0.0, growth_model.cm=>0.0, growth_model.cq=>0.0, 
growth_model.et=>0.0, growth_model.em=>0.0, growth_model.q=>0.0, growth_model.r=>0.0166,
growth_model.zmr=>0.0, growth_model.zmt=>0.0, growth_model.zmm=>0.0, growth_model.zmq=>0.0,
growth_model.si=>0.0, growth_model.a=>1.66]


params_gm = Dict(d=>d_val, kb=>kb_val, ku=>ku_val, thetar=>thetar_val, s0=>s0_val, gmax=>gmax_val, thetax=>thetax_val, Kt=>Kt_val, M=>M_val, we=>we_val, Km=>Km_val, vm=>vm_val, nx=>nx_val, Kq=>Kq_val, vt=>vt_val, wr=>wr_val, wq=>wq_val, nq=>nq_val, nr=>nr_val, ns=>ns_val, Kgamma=>Kgamma_val, abx=>abx_val, kon=>kon_val)
params_gm_uM = Dict(d=>d_val, kb=>kb_uM_val, ku=>ku_val, thetar=>thetar_uM_val, s0=>s0_uM_val, gmax=>gmax_val, thetax=>thetax_uM_val, Kt=>Kt_uM_val, M=>M_uM_val, we=>we_uM_val, Km=>Km_uM_val, vm=>vm_val, nx=>nx_val, Kq=>Kq_uM_val, vt=>vt_val, wr=>wr_uM_val, wq=>wq_uM_val, nq=>nq_val, nr=>nr_val, ns=>ns_val, Kgamma=>Kgamma_uM_val, abx=>abx_val, kon=>kon_val)

ssvals_gm = steady_states(growth_model, init_gm, params_gm)


# prob2 = ODEProblem(growth_model, init_gm, tspan, params_gm; jac=true);
# solu2 = solve(prob2, Rodas4(), abstol=1e-12, reltol=1e-9);


# function growth_model(dydt, initial, params, t)
#     dm, kb, ku, f, thetar, s0, gmax, thetax, Kt, M, we, Km, vm, nx, Kq, vt, wr, wq, nq, nr, ns, Kgamma = params
#     cr, em, cq, ct, et, cm, mt, mm, q, si, mq, mr, r, a = initial

#     dcr, dem, dcq, dct, det, dcm, dmt, dmm, dq, dsi, dmq, dmr, dr, da = zeros(length(dydt))
    
#     # Kgamma = 7*sf #gmax/Kp
#     gamma = gmax*a/(Kgamma+a) # aa min-1 molec-1 (aa min-1 molec-1(in terms of uM so this does cancel out below))
#     ttrate = (cq+cr+ct+cm)*gamma # aa min-1 (aa min-1)
#     lam = ttrate/M # min-1 (min-1)
#     vimp = (et*vt*s0/(Kt+s0)) # molec min-1 (uM min-1)
#     nucat = (em*vm*si/(Km+si)) # molec min-1 (uM min-1)

#     dydt[1] =  +r*mr*kb - cr*ku - cr*lam - gamma/nr*cr - f*cr # molec min-1
#     dydt[2] = - lam*em + cm*gamma/nx # molec min-1
#     dydt[3] =  +r*mq*kb - cq*ku - cq*lam - gamma/nx*cq - f*cq # molec min-1
#     dydt[4] =  +r*mt*kb - ct*ku - ct*lam - gamma/nx*ct - f*ct # molec min-1
#     dydt[5] =  - lam*et + gamma/nx*ct # molec min-1
#     dydt[6] =   +r*mm*kb - cm*ku - cm*lam - gamma/nx*cm - f*cm # molec min-1
#     dydt[7] =   +(we*a/(thetax+a)) - mt*dm - mt*lam - r*mt*kb + ct*ku + gamma/nx*ct # molec min-1
#     dydt[8] =   +(we*a/(thetax+a)) - mm*dm - mm*lam - r*mm*kb + cm*ku + gamma/nx*cm # molec min-1
#     dydt[9] = - lam*q + gamma/nx*cq # molec min-1
#     dydt[10] = - lam*si + vimp - nucat # molec min-1
#     dydt[11] =   +(wq*a/(thetax+a)/(1+(q/Kq)^nq)) - mq*dm - mq*lam - r*mq*kb + cq*ku + gamma/nx*cq # molec min-1
#     dydt[12] =   +(wr*a/(thetar+a)) - mr*dm - mr*lam - r*mr*kb + cr*ku + gamma/nr*cr # molec min-1
#     dydt[13] =  - lam*r - r*mr*kb - r*mt*kb - r*mm*kb - r*mq*kb + cr*ku + ct*ku + cm*ku + cq*ku + gamma/nr*cr + gamma/nr*cr + gamma/nx*ct + gamma/nx*cm + gamma/nx*cq # molec min-1
#     dydt[14] =  +ns*nucat - ttrate - lam*a ### aa min-1
# end

function growth_model_incl_abx(dydt, initial, params, t)
    # params
    dm, kb, ku, thetar, s0, gmax, thetax, Kt, M, we, Km, vm, nx, Kq, vt, wr, wq, nq, nr, ns, Kgamma, Cm, k_cm = params

    cr, em, cq, ct, et, cm, mt, mm, q, si, mq, mr, r, a, zmr, zmq, zmt, zmm = initial

    # equations
    gamma = (gmax*a)/(Kgamma+a)
    ttrate = (cr + ct + cm + cq)*gamma
    lam = ttrate/M
    vimp = (et*vt*s0/(Kt+s0))
    nucat = (em*vm*si/(Km+si))

    dydt[1] =  +r*mr*kb - cr*ku - cr*lam - gamma/nr*cr - cr*Cm*k_cm
    dydt[2] = - lam*em + cm*gamma/nx
    dydt[3] =  +r*mq*kb - cq*ku - cq*lam - gamma/nx*cq - cq*Cm*k_cm
    dydt[4] =  +r*mt*kb - ct*ku - ct*lam - gamma/nx*ct - ct*Cm*k_cm
    dydt[5] =  - lam*et + gamma/nx*ct
    dydt[6] =   +r*mm*kb - cm*ku - cm*lam - gamma/nx*cm - cm*Cm*k_cm
    dydt[7] =   +(we*a/(thetax+a)) - mt*dm - mt*lam - r*mt*kb + ct*ku + gamma/nx*ct
    dydt[8] =   +(we*a/(thetax+a)) - mm*dm - mm*lam - r*mm*kb + cm*ku + gamma/nx*cm
    dydt[9] = - lam*q + gamma/nx*cq
    dydt[10] = - lam*si + vimp - nucat
    dydt[11] =   +(wq*a/(thetax+a)/(1+(q/Kq)^nq)) - mq*dm - mq*lam - r*mq*kb + cq*ku + gamma/nx*cq
    dydt[12] =   +(wr*a/(thetar+a)) - mr*dm - mr*lam - r*mr*kb + cr*ku + gamma/nr*cr
    dydt[13] =  - lam*r - r*mr*kb - r*mt*kb - r*mm*kb - r*mq*kb + cr*ku + ct*ku + cm*ku + cq*ku + gamma/nr*cr + gamma/nr*cr + gamma/nx*ct + gamma/nx*cm + gamma/nx*cq
    dydt[14] =  +ns*nucat - ttrate - lam*a
    dydt[15] =  cr*Cm*k_cm - lam*zmr
    dydt[16] = cq*Cm*k_cm - lam*zmq
    dydt[17] = ct*Cm*k_cm - lam*zmt
    dydt[18] = cm*Cm*k_cm - lam*zmm
end

# function odemodel!(dydt, initial, params, t)
#     b, dm, kb, ku, f, thetar, k_cm, s0, gmax, cl, thetax, Kt, M, we, Km, vm, nx, Kq, Kp, vt, wr, wq, wp, nq, nr, ns, Kgamma = params
#     cr, em, cp, cq, ct, et, cm, mt, mm, q, p, si, mq, mp, mr, r, a = initial

#     dcr, dem, dcp, dcq, dct, det, dcm, dmt, dmm, dq, dp, dsi, dmq, dmp, dmr, dr, da = zeros(length(dydt))
    
#     # Kgamma = 7*sf #gmax/Kp
#     gamma = gmax*a/(Kgamma+a)
#     ttrate = (cq+cr+cp+ct+cm)*gamma
#     lam = ttrate/M
#     vimp = (et*vt*s0/(Kt+s0))
#     nucat = (em*vm*si/(Km+si))

#     dydt[1] =  +r*mr*kb - cr*ku - cr*lam - gamma/nr*cr - f*cr
#     dydt[2] = - lam*em + cm*gamma/nx
#     dydt[3] =   +r*mp*kb - cp*ku - cp*lam - gamma/nx*cp - f*cp
#     dydt[4] =  +r*mq*kb - cq*ku - cq*lam - gamma/nx*cq - f*cq
#     dydt[5] =  +r*mt*kb - ct*ku - ct*lam - gamma/nx*ct - f*ct
#     dydt[6] =  - lam*et + gamma/nx*ct
#     dydt[7] =   +r*mm*kb - cm*ku - cm*lam - gamma/nx*cm - f*cm
#     dydt[8] =   +(we*a/(thetax+a)) - mt*dm - mt*lam - r*mt*kb + ct*ku + gamma/nx*ct
#     dydt[9] =   +(we*a/(thetax+a)) - mm*dm - mm*lam - r*mm*kb + cm*ku + gamma/nx*cm
#     dydt[10] = - lam*q + gamma/nx*cq
#     dydt[11] =  - lam*p + gamma/nx*cp
#     dydt[12] = - lam*si + vimp - nucat
#     dydt[13] =   +(wq*a/(thetax+a)/(1+(q/Kq)^nq)) - mq*dm - mq*lam - r*mq*kb + cq*ku + gamma/nx*cq
#     dydt[14] = +(wp*a/(thetax+a)) - mp*dm - mp*lam - r*mp*kb + cp*ku + gamma/nx*cp
#     dydt[15] =   +(wr*a/(thetar+a)) - mr*dm - mr*lam - r*mr*kb + cr*ku + gamma/nr*cr
#     dydt[16] =  - lam*r - r*mr*kb - r*mt*kb - r*mm*kb - r*mq*kb - r*mp*kb + cr*ku + ct*ku + cm*ku + cq*ku + cp*ku + gamma/nr*cr + gamma/nr*cr + gamma/nx*ct + gamma/nx*cm + gamma/nx*cq + gamma/nx*cp
#     dydt[17] =  +ns*nucat - ttrate - lam*a
# end

# function odemodelfull!(dydt, initial, params, t)
#     b, dm, kb, ku, f, thetar, k_cm, s0, gmax, cl, thetax, Kt, M, we, Km, vm, nx, Kq, Kp, vt, wr, wq, wp, hq, nr, ns = params
#     cr, em, cp, cq, ct, et, cm, mt, mm, q, p, si, mq, mp, mr, r, a, zmr, zmp, zmq, zmt, zmm = initial

#     dcr, dem, dcp, dcq, dct, det, dcm, dmt, dmm, dq, dp, dsi, dmq, dmp, dmr, dr, da, dzmr, dzmp, dzmq, dzmt, dzmm = zeros(length(dydt))

#     Kgamma = gmax/Kp
#     gamma = gmax*a/(Kgamma+a)
#     ttrate = (cq+cr+cp+ct+cm)*gamma
#     lam = ttrate/M
#     vimp = (et*vt*s0)/(Kt+s0)
#     nucat = (em*vm*si)/(Km+si)

#     dydt[1] =  +r*mr*kb - cr*ku - cr*lam - cr*gamma/nr - f*cr + b*zmr
#     dydt[2] = - lam*em + cm*gamma/nx
#     dydt[3] =   +r*mp*kb - cp*ku - cp*lam - cp*gamma/nx - f*cp + b*zmp
#     dydt[4] =  +r*mq*kb - cq*ku - cq*lam - cq*gamma/nx - f*cq + b*zmq
#     dydt[5] =  +r*mt*kb - ct*ku - ct*lam - ct*gamma/nx - f*ct + b*zmt
#     dydt[6] =  - lam*et + ct*gamma/nx
#     dydt[7] =   +r*mm*kb - cm*ku - cm*lam - cm*gamma/nx - f*cm + b*zmm
#     dydt[8] =   +we*a/(thetax+a) - mt*lam - mt*dm - r*mt*kb + ct*ku + ct*gamma/nx
#     dydt[9] =   +we*a/(thetax+a) - mm*lam - mt*dm - r*mm*kb + cm*ku + cm*gamma/nx
#     dydt[10] = - lam*q + (cq*gamma)/nx
#     dydt[11] =  - lam*p + (cp*gamma)/nx
#     dydt[12] = - lam*si + vimp - nucat
#     dydt[13] =   +wq*a/(thetax+a)/1+(q/Kq)^hq - mq*lam- mt*dm - r*mq*kb + cq*ku + cq*gamma/nx
#     dydt[14] = +wp*a/(thetax+a) - mp*lam- mt*dm - r*mp*kb + cp*ku + cp*gamma/nx
#     dydt[15] =   +wr*a/(thetar+a) - mr*lam- mt*dm- r*mr*kb + cr*ku + (cr*gamma)/nr
#     dydt[16] =  - lam*r - r*mr*kb - r*mt*kb - r*mm*kb - r*mq*kb - r*mp*kb + cr*ku + ct*ku + cm*ku + cq*ku + cp*ku + cr*gamma/nr + cr*gamma/nr + ct*gamma/nx + cm*gamma/nx + cq*gamma/nx + cp*gamma/nx
#     dydt[17] =  +ns*nucat - ttrate - lam*a
#     dydt[18] = +f*cr-b*zmr-lam*zmr
#     dydt[19] = +f*cp-b*zmp-lam*zmp
#     dydt[20] = +f*cq-b*zmq-lam*zmq
#     dydt[21] = +f*ct-b*zmt-lam*zmt
#     dydt[22] = +f*cm-b*zmm-lam*zmm
# end

# function pop_model!(dydt, initial, params, t)
#     dm, kb, ku, f, thetar, gmax, thetax, Kt, M, we, Km, vm, nx, Kq, vt, wr, wq, wp, nq, nr, ns, kin, d_s, d_n, Kgamma = params
#     cr, em, cp, cq, ct, et, cm, mt, mm, q, p, si, mq, mp, mr, r, a, s, N = initial

#     dcr, dem, dcp, dcq, dct, det, dcm, dmt, dmm, dq, dp, dsi, dmq, dmp, dmr, dr, da, ds, dN = zeros(length(dydt))
    
#     # Kgamma = 7*sf;#gmax/Kp
#     gamma = gmax*a/(Kgamma+a)
#     ttrate = (cq+cr+cp+ct+cm)*gamma
#     lam = ttrate/M
#     vimp = (et*vt*s/(Kt+s))
#     nucat = (em*vm*si/(Km+si))

#     dydt[1] =  +r*mr*kb - cr*ku - cr*lam - gamma/nr*cr - f*cr
#     dydt[2] = - lam*em + cm*gamma/nx
#     dydt[3] =   +r*mp*kb - cp*ku - cp*lam - gamma/nx*cp - f*cp
#     dydt[4] =  +r*mq*kb - cq*ku - cq*lam - gamma/nx*cq - f*cq
#     dydt[5] =  +r*mt*kb - ct*ku - ct*lam - gamma/nx*ct - f*ct
#     dydt[6] =  - lam*et + gamma/nx*ct
#     dydt[7] =   +r*mm*kb - cm*ku - cm*lam - gamma/nx*cm - f*cm
#     dydt[8] =   +(we*a/(thetax+a)) - mt*dm - mt*lam - r*mt*kb + ct*ku + gamma/nx*ct
#     dydt[9] =   +(we*a/(thetax+a)) - mm*dm - mm*lam - r*mm*kb + cm*ku + gamma/nx*cm
#     dydt[10] = - lam*q + gamma/nx*cq
#     dydt[11] =  - lam*p + gamma/nx*cp
#     dydt[12] = - lam*si + vimp - nucat
#     dydt[13] =   +(wq*a/(thetax+a)/(1+(q/Kq)^nq)) - mq*dm - mq*lam - r*mq*kb + cq*ku + gamma/nx*cq
#     dydt[14] = +(wp*a/(thetax+a)) - mp*dm - mp*lam - r*mp*kb + cp*ku + gamma/nx*cp
#     dydt[15] =   +(wr*a/(thetar+a)) - mr*dm - mr*lam - r*mr*kb + cr*ku + gamma/nr*cr
#     dydt[16] =  - lam*r - r*mr*kb - r*mt*kb - r*mm*kb - r*mq*kb - r*mp*kb + cr*ku + ct*ku + cm*ku + cq*ku + cp*ku + gamma/nr*cr + gamma/nr*cr + gamma/nx*ct + gamma/nx*cm + gamma/nx*cq + gamma/nx*cp
#     dydt[17] =  +ns*nucat - ttrate - lam*a
#     dydt[18] = kin - vimp*N - d_s*s
#     dydt[19] = lam*N - d_n*N
# end




