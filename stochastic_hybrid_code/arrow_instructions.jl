using StatsBase, Distributions, Random, DataFrames, CSV, DifferentialEquations, OrderedCollections, ProgressBars, BenchmarkTools, Arrow, FilePathsBase

# this line needs to be added in at the end of whenever you run a simulation that creates many .dat files 
arrow_conv("/home/hollie_hindley/Documents/phd/stochastic_hybrid_code/testing_conv", "/home/hollie_hindley/Documents/phd/stochastic_hybrid_code/testing_conv_arrow")

# this is how you load in the arrow files as dataframes - should be really quick now 
function load_arrow_files(folder_path)
    tabs = [Arrow.Table(joinpath(folder_path, file)) for file in readdir(folder_path)]
    dfs = [DataFrame(tab) for tab in tabs]
    return dfs
end

dfs = load_arrow_files("/home/hollie_hindley/Documents/phd/stochastic_hybrid_code/testing_conv_arrow")

# now need to convert all .dat files to .arrow files that I plan to analyse