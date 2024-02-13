sf = 1e6/(6.022e23*1e-15) # multiply to go from molecs/cell to uM 

1e6* sf
4.38*sf*22000

w_rh = 930 # molecs/min*cell
w_t = 4.14 # molecs/min*cell, same as w_m
w_q = 948.93 # molecs/min*cell

θ_rh = 426.87 # molecs/cell
θ_nr = 4.38 # molecs/cell
Kq = 152219 # molecs/cell

gmax = 1260 # aa/min*molecs
Kgamma = 7 # molecs/cell  

kb = 1 # cell/min*molecs 

s0 = 1e4 # molecs 

Kt = 1000 # molecs 

Km = 1000 # molecs/cell 

# NEW units
θ_rh_NEW = θ_rh*sf*22000 # uM
θ_nr_NEW = θ_nr*sf*22000 # uM
Kq_NEW = Kq*sf # uM
Kgamma_NEW = Kgamma*sf*22000 # uM  
Km_NEW = Km*sf # uM


w_rh_NEW = w_rh*sf # uM/min
w_t_NEW = w_t*sf # uM/min
w_q_NEW = w_q*sf # uM/min


gmax_NEW = gmax*sf # aa/min*molecs

kb_NEW = kb*sf # 1/min*molecs 


s0_NEW = s0*sf # uM 

Kt_NEW = Kt*sf # uM 

