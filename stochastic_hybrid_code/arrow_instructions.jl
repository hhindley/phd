using StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Arrow, FilePathsBase

# this line needs to be added in at the end of whenever you run a simulation that creates many .dat files 
arrow_conv("/home/hollie_hindley/Documents/stochastic_hybrid/thresh_test", "/home/hollie_hindley/Documents/stochastic_hybrid/thresh_test_arrow")


dfs = load_arrow_files("/home/hollie_hindley/Documents/phd/stochastic_hybrid_code/testing_conv_arrow")

# now need to convert all .dat files to .arrow files that I plan to analyse