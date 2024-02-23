using ModelingToolkit, DifferentialEquations, PlotlyJS, LinearAlgebra, DataFrames, LabelledArrays, Printf, BifurcationKit

include("/home/holliehindley/phd/rtc_model/models/rtc_orig.jl")
include("/home/holliehindley/phd/general_funcs/solving.jl")
include("/home/holliehindley/phd/rtc_model/parameters/params.jl")
include("/home/holliehindley/phd/rtc_model/functions/bf_funcs/bf_funcs.jl")


br = get_br(rtc_model, ssvals_rtc, params_rtc, 1.5)

