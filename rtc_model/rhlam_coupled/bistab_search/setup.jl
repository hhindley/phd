include(joinpath(homedir(), "phd/rtc_model/rhlam_coupled/rhlam_model.jl"))

lam_c_vals = 10 .^range(log10(1e-8),log10(1), length=10) #8e-7
kin_c_vals = 10 .^range(log10(1e-6),log10(1), length=10) #1.5e-5
wab_vals = 10 .^range(log10(1e-6),log10(1), length=10) #1.3e-5
kdam_range = range(0,100,length=20)
