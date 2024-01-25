# parameters
sf_gm = 22000; # molec aa-1

s0= 1e4; # molecs
dm= 0.1; # min-1
ns= 0.5; #0.01
nr= 7549.0; # aa molecs-1
nx= 300.0; # aa molecs-1
gmax= 1260.0; # aa min-1 molecs-1
Kgamma = 7*sf_gm; # molecs 
vt= 726.0; # min-1 
Kt= 1.0e3; # molecs
vm= 5800.0; # min-1
Km= 1.0e3; # molecs
wr= 929.9678874564831; # molecs min-1
we= 4.139172187824451; # molecs min-1
wq= 948.9349882947897; # molecs min-1
thetar= 426.8693338968694*sf_gm; # molecs 
thetax= 4.379733394834643*sf_gm; # molecs
Kq= 1.522190403737490e+05; # molecs
nq= 4;
kb= 1; # min-1 molecs-1
ku= 1.0; # min-1
M= 1.0e8; # (molecs of) aa
k_cm= 0.005990373118888; # min-1 uM-1


f= 0#1;
b= 0;
cl= 0;
Kp= 180.1378030928276; # molecs
wp= 0.0; # molecs min-1
kin = 0;
d_s = 0;
d_n = 0;
# params_init= [b, dm, kb, ku, f_init, thetar, k_cm, s0, gmax, cl, thetax, Kt, M, we, Km, vm, nx, Kq, Kp, vt, wr, wq, wp, hq, nr, ns]

params= [dm, kb, ku, f, thetar, s0, gmax, thetax, Kt, M, we, Km, vm, nx, Kq, vt, wr, wq, nq, nr, ns, Kgamma]
params_ext = [b, dm, kb, ku, f, thetar, k_cm, s0, gmax, cl, thetax, Kt, M, we, Km, vm, nx, Kq, Kp, vt, wr, wq, wp, nq, nr, ns, Kgamma]
pop_params = [dm, kb, ku, f, thetar, gmax, thetax, Kt, M, we, Km, vm, nx, Kq, vt, wr, wq, wp, nq, nr, ns, kin, d_s, d_n, Kgamma]

cr_0= 0.
em_0= 0.
cp_0= 0.
cq_0= 0.
ct_0= 0.
et_0= 0. #0
cm_0= 0.
mt_0= 0.
mm_0= 0.
q_0= 0.
p_0= 0.
si_0= 0.
mq_0= 0.
mp_0= 0.
mr_0= 0.
r_0= 10.0
a_0= 1000.0
s_0 = 1e4#1e6
N_0 = 1
init1= [cr_0, em_0, cp_0, cq_0, ct_0, et_0, cm_0, mt_0, mm_0, q_0, p_0, si_0, mq_0, mp_0, mr_0, r_0, a_0]
init_gm= [cr_0, em_0, cq_0, ct_0, et_0, cm_0, mt_0, mm_0, q_0, si_0, mq_0, mr_0, r_0, a_0]
pop_init= [cr_0, em_0, cp_0, cq_0, ct_0, et_0, cm_0, mt_0, mm_0, q_0, p_0, si_0, mq_0, mp_0, mr_0, r_0, a_0, s_0, N_0]


cr_0= 0.
em_0= 0.
cp_0= 0.
cq_0= 0.
ct_0= 0.
et_0= 0.
cm_0= 0.
mt_0= 0.
mm_0= 0.
q_0= 0.
p_0= 0.
si_0= 0.
mq_0= 0.
mp_0= 0.
mr_0= 0.
r_0= 10.0
a_0= 1000.0
zmr_0= 0.
zmp_0= 0. 
zmq_0= 0.
zmt_0= 0.
zmm_0= 0.

initfull= [cr_0, em_0, cp_0, cq_0, ct_0, et_0, cm_0, mt_0, mm_0, q_0, p_0, si_0, mq_0, mp_0, mr_0, r_0, a_0, zmr_0, zmp_0, zmq_0, zmt_0, zmm_0]
