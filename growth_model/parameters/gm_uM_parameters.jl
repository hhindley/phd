include("$PATH/growth_model/parameters/growth_model_params.jl")

# parameters
sf = 1e6/(6.022e23*1e-15)

s0_uM_val= s0_val*sf; # uM
Kgamma_uM_val = Kgamma_val*sf # uM
Kt_uM_val= Kt_val*sf; # uM
Km_uM_val= Km_val*sf; # uM
wr_uM_val= wr_val*sf; # uM min-1
we_uM_val= we_val*sf # uM min-1
wq_uM_val= wq_val*sf; # uM min-1 
thetar_uM_val= thetar_val*sf # uM
thetax_uM_val= thetax_val*sf # uM
Kq_uM_val= Kq_val*sf; # uM
kb_uM_val= kb_val/sf; # min-1 uM-1
M_uM_val= M_val*sf; # (uM of) aa 

# Kp_uM= Kp*sf;

# params_init= [b, dm, kb, ku, f_init, thetar, k_cm, s0, gmax, cl, thetax, Kt, M, we, Km, vm, nx, Kq, Kp, vt, wr, wq, wp, hq, nr, ns]
# params= [b, dm, kb, ku, f, thetar, k_cm, s0, gmax, cl, thetax, Kt, M, we, Km, vm, nx, Kq, Kp, vt, wr, wq, wp, hq, nr, ns]


# params_uM = [0, 0.1, 0.00166, 0.00166, 0, 15594.7, 0, 1661, 2.0923, 0, 160.01, 1.66, 1e8, 0.00687, 1.66, 5800, 300, 252.77, 0.299, 726, 1.54, 1.576, 0.00687, 4, 7459, 0.5, 255.7]

# r_0_uM= 0.0166
# a_0_uM= 1.66

# init_gm_uM= [cr_0, em_0, cq_0, ct_0, et_0, cm_0, mt_0, mm_0, q_0, si_0, mq_0, mr_0, r_0_uM, a_0_uM]

# solu_gm_uM = sol(growth_model, init_gm_uM, tspan, params_gm_uM)
# ssvals_gm_uM = get_all_ssvals(solu_gm_uM, species_gm)