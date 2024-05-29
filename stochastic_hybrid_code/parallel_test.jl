using Distributed
nprocs()
addprocs(10)

@everywhere include("/home/hollie_hindley/Documents/phd/stochastic_hybrid_code/setup/file_funcs.jl")

@elapsed df_props, df_reacts, df_res = loadsort_all_arrow_files("/home/hollie_hindley/Documents/stochastic_hybrid/thresh_test_arrow")

df_props
