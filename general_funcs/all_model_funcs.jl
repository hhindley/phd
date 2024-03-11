ω_rtcBA(θ_x, a, sig_o) = sig_o * w_BA*a/(θ_x + a) # uM min-1

ω_p(w_x, θ_x, a) = w_x*a/(θ_x + a) # uM min-1 

ω_q(w_q, θ_x, a, q) = (w_q*a/(θ_x + a))/(1+(q/Kq)^nq) # uM min-1

v_x(c_x, nx, gamma) = c_x*gamma/nx # uM min-1 

dil(x,lam) = lam*x # uM min-1

deg(x) = d*x # uM min-1

rh_bind(m_x,rh) = kb*rh*m_x # uM min-1

rh_unbind(c_x) = ku*c_x # uM min-1

zm(c_x) = c_x*abx*kon # uM min-1 

zm_diss(z_x) = koff*z_x # uM min-1

tlr(rm_x, nx, rh, tlr_el) = (1/nx)*kc*rh*rm_x*tlr_el # uM min-1 