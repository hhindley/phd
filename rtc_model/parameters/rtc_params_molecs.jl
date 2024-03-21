sf = 1e6/(6.022e23*1e-15) # multiply to go from molecs/cell to uM # uM molec-1

kr_val_molec = kr_val/sf; # uM

Km_init_val_molec = Km_init_val/sf; # uM
θtscr_val_molec = θtscr_val/sf; # uM
θtlr_val_molec = θtlr_val/sf; # uM
# k_b = 17.7; 

atp_val_molec = atp_val/sf;#4000;#2500; # uM
km_a_val_molec = km_a_val/sf; # uM
km_b_val_molec = km_b_val/sf; # uM


kin_val_molec = kin_val/sf #2.381 # uM molec aa-1 

ω_r_val_molec = ω_r_val/sf; #0.01/1e4 #0.0019*6 #70.53; #0.0019*6#79.43865871861044; #0.0019*6;  # uM min-1

# kb_gm_val = 1/sf # min-1 uM-1
# ku_gm_val = 1. # min-1

kc_val_molec = kc_val*sf # kc (uM-1) is the binding constant (kb/ku) (uM-1) for ribosomes to mRNAs - kc/sf = uM-1 


