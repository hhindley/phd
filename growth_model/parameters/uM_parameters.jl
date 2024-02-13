include("/home/holliehindley/phd/growth_model/parameters/parameters.jl")

# parameters
sf = 1e6/(6.022e23*1e-15)

s0_uM= s0*sf; # uM
Kgamma_uM = Kgamma*sf # uM
Kt_uM= Kt*sf; # uM
Km_uM= Km*sf; # uM
wr_uM= wr*sf; # uM min-1
we_uM= we*sf # uM min-1
wq_uM= wq*sf; # uM min-1 
thetar_uM= thetar*sf # uM
thetax_uM= thetax*sf # uM
Kq_uM= Kq*sf; # uM
kb_uM= kb/sf; # min-1 uM-1
M_uM= M*sf; # (uM of) aa 

Kp_uM= Kp*sf;

# params_init= [b, dm, kb, ku, f_init, thetar, k_cm, s0, gmax, cl, thetax, Kt, M, we, Km, vm, nx, Kq, Kp, vt, wr, wq, wp, hq, nr, ns]
# params= [b, dm, kb, ku, f, thetar, k_cm, s0, gmax, cl, thetax, Kt, M, we, Km, vm, nx, Kq, Kp, vt, wr, wq, wp, hq, nr, ns]
params_uM = @LArray [dm, kb_uM, ku, f, thetar_uM, s0_uM, gmax, thetax_uM, Kt_uM, M_uM, we_uM, Km_uM, vm, nx, Kq_uM, vt, wr_uM, wq_uM, nq, nr, ns, Kgamma_uM] (:dm, :kb, :ku, :f, :thetar, :s0, :gmax, :thetax, :Kt, :M, :we, :Km, :vm, :nx, :Kq, :vt, :wr, :wq, :nq, :nr, :ns, :Kgamma)

# params_uM = [0, 0.1, 0.00166, 0.00166, 0, 15594.7, 0, 1661, 2.0923, 0, 160.01, 1.66, 1e8, 0.00687, 1.66, 5800, 300, 252.77, 0.299, 726, 1.54, 1.576, 0.00687, 4, 7459, 0.5, 255.7]

r_0_uM= 0.0166
a_0_uM= 1.66

init_gm_uM= [cr_0, em_0, cq_0, ct_0, et_0, cm_0, mt_0, mm_0, q_0, si_0, mq_0, mr_0, r_0_uM, a_0_uM]
