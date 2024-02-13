a_0 = 5 # 10 
b_0 = 10 # 10
c_0 = 0 # 0
d_0 = 5 # 5
p_0 = 0 # 0
e_0 = 0 #0
init = [b_0, c_0, d_0, p_0];
init2 = [b_0, c_0, d_0, e_0, p_0];
init3 = [a_0, b_0, c_0, d_0, p_0];
init_full = [a_0, b_0, c_0, d_0, e_0, p_0]

k1 = 0.1 # 0.1
k2 = 0.1 # 0.1
k3 = 0.1 # 0.1
k4 = 0.1 # 0.1
katp = 10 # 10
# a0 = 5 # 10
# b0 = 10 # 10
params = [k1, k2, k3, k4, katp, a_0, b_0];
params_full = [k1, k2, k3, k4, katp];
comp_params = [k1, k2, k3, k4, katp, a_0, b_0];


tspan = (0.0, 100);

species = ["b-Rt" "c-Ra" "d-\\sigma_pc" "p-\\sigma*"];
species3 = ["b-Rt" "c-Ra" "d-\\sigma_pc" "e-Ra\\sigma_pc" "p-\\sigma*"];
species2 = [:b, :c, :d, :p];
species_full = ["a-Ri" "b-Rt" "c-Ra" "d-\\sigma_pc" "e-Ra\\sigma_pc" "p-\\sigma*"] 
species_full2 = [:a, :b, :c, :d, :e, :p]
species4 = [:b, :c, :d, :e, :p]
species5 = ["a-Ri" "b-Rt" "c-Ra" "d-\\sigma_pc" "p-\\sigma*"]