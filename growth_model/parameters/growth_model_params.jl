# parameters
sf_gm = 22000; # molec aa-1

s0_val= 1e4; # molecs
d_val= 0.1; # min-1
ns_val= 0.5; #0.01
nr_val= 7549.0; # aa molecs-1
nx_val= 300.0; # aa molecs-1
gmax_val= 1260.0; # aa min-1 molecs-1
Kgamma_val = 7*sf_gm; # molecs 
vt_val= 726.0; # min-1 
Kt_val= 1.0e3; # molecs
vm_val= 5800.0; # min-1
Km_val= 1.0e3; # molecs
wr_val= 929.9678874564831; # molecs min-1
we_val= 4.139172187824451; # molecs min-1
wq_val= 948.9349882947897; # molecs min-1
thetar_val= 426.8693338968694*sf_gm; # molecs 
thetax_val= 4.379733394834643*sf_gm; # molecs
Kq_val= 1.522190403737490e+05; # molecs
nq_val= 4;
kb_val= 1; # min-1 molecs-1
ku_val= 1.0; # min-1
M_val= 1.0e8; # (molecs of) aa
kon_val= 0.005990373118888; # min-1 uM-1

abx_val=0 # uM
# f= 0#1;
# b= 0;
# cl= 0;
# Kp= 180.1378030928276; # molecs
# wp= 0.0; # molecs min-1
# kin = 0;
# d_s = 0;
# d_n = 0;
# params_init= [b, dm, kb, ku, f_init, thetar, k_cm, s0, gmax, cl, thetax, Kt, M, we, Km, vm, nx, Kq, Kp, vt, wr, wq, wp, hq, nr, ns]
tspan = (0, 1e9);



# species_gm = [:cr, :em, :cq, :ct, :et, :cm, :mt, :mm, :q, :si, :mq, :mr, :r, :a, :zmr, :zmq, :zmt, :zmm]


# params= [dm, kb, ku, f, thetar, s0, gmax, thetax, Kt, M, we, Km, vm, nx, Kq, vt, wr, wq, nq, nr, ns, Kgamma]
# params_abx = [dm, kb, ku, thetar, s0, gmax, thetax, Kt, M, we, Km, vm, nx, Kq, vt, wr, wq, nq, nr, ns, Kgamma, Cm, k_cm]
# params_ext = [b, dm, kb, ku, f, thetar, k_cm, s0, gmax, cl, thetax, Kt, M, we, Km, vm, nx, Kq, Kp, vt, wr, wq, wp, nq, nr, ns, Kgamma]
# pop_params = [dm, kb, ku, f, thetar, gmax, thetax, Kt, M, we, Km, vm, nx, Kq, vt, wr, wq, wp, nq, nr, ns, kin, d_s, d_n, Kgamma]

# cr_0= 0.
# em_0= 0.
# cp_0= 0.
# cq_0= 0.
# ct_0= 0.
# et_0= 0. #0
# cm_0= 0.
# mt_0= 0.
# mm_0= 0.
# q_0= 0.
# p_0= 0.
# si_0= 0.
# mq_0= 0.
# mp_0= 0.
# mr_0= 0.
# r_0= 10.0
# a_0= 1000.0
# s_0 = 1e4#1e6
# N_0 = 1
# init1= [cr_0, em_0, cp_0, cq_0, ct_0, et_0, cm_0, mt_0, mm_0, q_0, p_0, si_0, mq_0, mp_0, mr_0, r_0, a_0]
# init_gm= [cr_0, em_0, cq_0, ct_0, et_0, cm_0, mt_0, mm_0, q_0, si_0, mq_0, mr_0, r_0, a_0]
# pop_init= [cr_0, em_0, cp_0, cq_0, ct_0, et_0, cm_0, mt_0, mm_0, q_0, p_0, si_0, mq_0, mp_0, mr_0, r_0, a_0, s_0, N_0]


# cr_0= 0.
# em_0= 0.
# cp_0= 0.
# cq_0= 0.
# ct_0= 0.
# et_0= 0.
# cm_0= 0.
# mt_0= 0.
# mm_0= 0.
# q_0= 0.
# p_0= 0.
# si_0= 0.
# mq_0= 0.
# mp_0= 0.
# mr_0= 0.
# r_0= 10.0
# a_0= 1000.0
# zmr_0= 0.
# zmp_0= 0. 
# zmq_0= 0.
# zmt_0= 0.
# zmm_0= 0.

# # initfull= [cr_0, em_0, cp_0, cq_0, ct_0, et_0, cm_0, mt_0, mm_0, q_0, p_0, si_0, mq_0, mp_0, mr_0, r_0, a_0, zmr_0, zmp_0, zmq_0, zmt_0, zmm_0]

# init_abx=[cr_0, em_0, cq_0, ct_0, et_0, cm_0, mt_0, mm_0, q_0, si_0, mq_0, mr_0, r_0, a_0, zmr_0, zmq_0, zmt_0, zmm_0]