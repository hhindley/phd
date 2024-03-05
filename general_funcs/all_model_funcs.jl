ω_rtcBA(θ_x, a, sig_o) = sig_o * w_BA*a/(θ_x + a)

ω_p(w_x, θ_x, a) = w_x*a/(θ_x + a)

ω_q(w_q, θ_x, a, q) = (w_q*a/(θ_x + a))/(1+(q/Kq)^nq)

v_x(c_x, nx, gamma) = c_x*gamma/nx

dil(x,lam) = lam*x

deg(x) = d*x

rh_bind(m_x,rh) = kb*rh*m_x

rh_unbind(c_x) = ku*c_x

zm(c_x) = c_x*abx*kon

zm_diss(z_x) = koff*z_x

tlr(rm_x, nx, rh, tlr_el) = (1/nx)*kc*rh*rm_x*tlr_el 