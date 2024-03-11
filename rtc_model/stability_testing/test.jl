using ModelingToolkit, DifferentialEquations, PlotlyJS, LinearAlgebra, DataFrames, LabelledArrays, Printf, BifurcationKit

include("$PATHrtc_model/models/rtc_orig.jl")
include("$PATHgeneral_funcs/solving.jl")
include("$PATHrtc_model/parameters/params.jl")
include("$PATHrtc_model/functions/bf_funcs/bf_funcs.jl")


br = get_br(rtc_model, ssvals_rtc, params_rtc, 1.5)

