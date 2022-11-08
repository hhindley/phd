# params

s = 1e4
dm = 0.1
ns = 0.5
nr = 7459 
nx = 300
gmax = 1260
Kgamma = 7
vt = 726 
Kt = 1000
vm = 5800
Km = 1000 
wr = 930
we = 4.14
wq = 948.93
thetar = 426.87
thetanr = 4.38
Kq = 152219
hq = 4
kb = 1
ku = 1
M = 1e8
kcm = 0.00599
params = [s, dm, ns, nr, nx, gmax, Kgamma, vt, Kt, vm, Km, wr, we, wq, thetar, thetanr, Kq, hq, kb, ku, M, kcm]



# rates
lam = 2 
vimp = 2 
dm = 1
kb = 1
ku = 1
wx = 1
vx = 10
rates = [lam, vimp, dm, kb, ku, wx, vx]

# initial values
r_0 = 0
et_0 = 0
em_0 = 0
q_0 = 0
si_0 = 0
a_0 = 0
mr_0 = 0
mt_0 = 0
mm_0 = 0
mq_0 = 0
cr_0 = 0
ct_0 = 0
cm_0 = 0
cq_0 = 0
init = [r_0, et_0, em_0, q_0, si_0, a_0, mr_0, mt_0, mm_0, mq_0, cr_0, ct_0, cm_0, cq_0]

















from cmath import e
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


# define the function
def model(t, variables, s, d_m, n_s, n_r, n_t, n_m, n_q, gamma_max, K_gamma, v_t, K_t, v_m, K_m, w_r, w_e, w_t, w_m, \
    w_q, theta_r, theta_nr, K_q, h_q, k_b, k_u, M, k_cm):

    s_i, a, r, e_t, e_m, q, m_r, m_t, m_m, m_q, c_r, c_t, c_m, c_q =  variables

    # helper function
    def v_imp(e_t, s):
        return e_t * v_t * s / (K_t + s)


    def v_cat(e_m, s_i):
        return e_m * v_m * s_i / (K_m + s_i)


    def gamma(a):
        return gamma_max * a / (K_gamma + a)


    def v_x(x, c_x, a):
        if x == 'r':
            return c_x * gamma(a) / n_r
        elif x == 't':
            return c_x * gamma(a) / n_t
        elif x == 'm':
            return c_x * gamma(a) / n_m
        else:
            return c_x * gamma(a) / n_q

    def omega(x, a):
        if x == 'r':
            return w_r * a / (theta_r + a)
        elif x == 't':
            return w_t * a / (theta_nr + a)
        else:
            return w_m * a / (theta_nr + a)


    def iota(q):
        return 1 / (1 + (q + K_q) ** h_q)


    def omega_q(q, a):
        return w_q * a * iota(q) / (theta_nr + a)

    R_t = c_r + c_t + c_m + c_q
    growth_rate = gamma(a) * R_t / M

    dsidt = v_imp(e_t, s) - v_cat(e_m, s_i) - growth_rate * s_i
    dadt = n_s * v_cat(e_m, s_i) - n_r * v_x('r', c_r, a) - n_t * v_x('t', c_t, a) - n_m * v_x('m', c_m, a) - n_q * v_x('q', c_m, a) - growth_rate * a
    drdt = v_x('r', c_r, a) - growth_rate * r + v_x('r', c_r, a) - k_b * r * m_r + k_u * c_r + v_x('t', c_t, a) - k_b * r * m_t + k_u * c_t + v_x('m', c_m, a) - k_b * r * m_m + k_u * c_m + v_x('q', c_q, a) - k_b * r * m_q + k_u * c_q
    detdt = v_x('t', c_t, a) - growth_rate * e_t
    demdt = v_x('m', c_m, a) - growth_rate * e_m
    dqdt = v_x('q', c_q, a) - growth_rate * q
    dmrdt = omega('r', a) - (growth_rate + d_m) * m_r - k_b * r * m_r + k_u * c_r + v_x('r', c_r, a)
    dmtdt = omega('t', a) - (growth_rate + d_m) * m_t - k_b * r * m_t + k_u * c_t + v_x('t', c_t, a)
    dmmdt = omega('m', a) - (growth_rate + d_m) * m_m - k_b * r * m_m + k_u * c_m + v_x('m', c_t, a)
    dmqdt = omega_q(q, a) - (growth_rate + d_m) * m_q - k_b * r * m_q + k_u * c_q + v_x('q', c_t, a)
    dcrdt = k_b * m_r * r - k_u * c_r - v_x('r', c_r, a) - growth_rate * c_r
    dctdt = k_b * m_t * r - k_u * c_t - v_x('t', c_t, a) - growth_rate * c_t
    dcmdt = k_b * m_m * r - k_u * c_m - v_x('m', c_m, a) - growth_rate * c_m
    dcqdt = k_b * m_q * r - k_u * c_q - v_x('q', c_q, a) - growth_rate * c_q

    return [dsidt, dadt, drdt, detdt, demdt, dqdt, dmrdt, dmtdt, dmmdt, dmqdt, dcrdt, dctdt, dcmdt, dcqdt]

# params
s = 10 ** 4  # external nutrient
d_m = 0.1  # mRNA-degradation rate
n_s = 0.5  # nutrient efficiency
n_r = 7459  # ribosome length
n_t, n_m, n_q = 300, 300, 300  # length of non-ribosomal proteins
gamma_max = 1260  # max. transl. elongation rate
K_gamma = 7  # transl. elongation threshold
v_t = 726  # max. nutrient import rate
K_t = 1000  # nutrient import threshold
v_m = 5800  # max. enzymatic rate
K_m = 1000  # enzymatic threshold
w_r = 930  # max. ribosome transcription rate
w_e, w_t, w_m = 4.14, 4.14, 4.14  # max. enzyme transcription rate
w_q = 948.93  # max. q-transcription rate
theta_r = 426.87  # ribosome transcription threshold
theta_nr = 4.38  # non-ribosomal transcription threshold
K_q = 152219  # q-autoinhibition threshold
h_q = 4  # q-autoinhibition Hill coeff.
k_b = 1  # mRNA-ribosome binding rate
k_u = 1  # mRNA-ribosome unbinding rate
M = 10 ** 8  # total cell mass
k_cm = 0.00599  # chloramphenicol-binding rate

p = (s, d_m, n_s, n_r, n_t, n_m, n_q, gamma_max, K_gamma, v_t, K_t, v_m, K_m, w_r, w_e, w_t, w_m, \
    w_q, theta_r, theta_nr, K_q, h_q, k_b, k_u, M, k_cm)

y0 = [0 for i in range(14)]
y0[1] = 1000 # initial value of energy
y0[2] = 10 # initial value of ribosomes

t_span = (0.0, 900.0)
t = np.arange(0.0, 900.0, 0.1)

result_solve_ivp = solve_ivp(model, t_span, y0, args=p, t_eval=t)
t = result_solve_ivp.t

s_i = result_solve_ivp.y[0]
a = result_solve_ivp.y[1]
r = result_solve_ivp.y[2]
e_t = result_solve_ivp.y[3]
e_m = result_solve_ivp.y[4]
q = result_solve_ivp.y[5]
m_r = result_solve_ivp.y[6]
m_t = result_solve_ivp.y[7]
m_m = result_solve_ivp.y[8]
m_q = result_solve_ivp.y[9]
c_r = result_solve_ivp.y[10]
c_t = result_solve_ivp.y[11]
c_m = result_solve_ivp.y[12]
c_q = result_solve_ivp.y[13]

print("finish")


fig, ax = plt.subplots(14, 1, figsize=(15, 35), sharex=True)
ax[0].plot(t, s_i, 'b-', label='s_i')
ax[0].set_ylabel('growth_rate')
ax[0].legend()

ax[1].plot(t, a, 'g-', label = 'a')
ax[1].set_ylabel('growth_rate')
ax[1].legend()

ax[2].plot(t, r, 'r-', label='r')
ax[2].set_ylabel('growth_rate')
ax[2].legend()

ax[3].plot(t, e_t, 'b-', label='e_t')
ax[3].set_ylabel('growth_rate')
ax[3].legend()

ax[4].plot(t, e_m, 'b-', label='e_m')
ax[4].set_ylabel('growth_rate')
ax[4].legend()

ax[5].plot(t, q, 'b-', label='q')
ax[5].set_ylabel('growth_rate')
ax[5].legend()

ax[6].plot(t, m_r, 'b-', label='m_r')
ax[6].set_ylabel('growth_rate')
ax[6].legend()

ax[7].plot(t, m_t, 'b-', label='m_t')
ax[7].set_ylabel('growth_rate')
ax[7].legend()

ax[8].plot(t, m_m, 'b-', label='m_m')
ax[8].set_ylabel('growth_rate')
ax[8].legend()

ax[9].plot(t, m_q, 'b-', label='m_q')
ax[9].set_ylabel('growth_rate')
ax[9].legend()

ax[10].plot(t, c_r, 'b-', label='c_r')
ax[10].set_ylabel('growth_rate')
ax[10].legend()

ax[11].plot(t, c_t, 'b-', label='c_t')
ax[11].set_ylabel('growth_rate')
ax[11].legend()

ax[12].plot(t, c_m, 'b-', label='c_m')
ax[12].set_ylabel('growth_rate')
ax[12].legend()

ax[13].plot(t, c_q, 'b-', label='c_q')
ax[13].set_xlabel('time(s)')
ax[13].set_ylabel('growth_rate')
ax[13].legend()

plt.savefig('model_15min.png')
plt.show()
