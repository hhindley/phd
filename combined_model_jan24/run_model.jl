using Parameters, CSV, DataFrames, DifferentialEquations, StaticArrays, LabelledArrays, BenchmarkTools, OrderedCollections, DataInterpolations, Statistics
using Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra, Printf
using PlotlyJS, ProgressBars

include("/home/holliehindley/phd/may23_rtc/functions/solving.jl")
include("/home/holliehindley/phd/combined_model_jan24/params.jl")
include("/home/holliehindley/phd/combined_model_jan24/combined_model.jl")

comb_species = [:m_rh, :m_t, :m_m, :m_q, :m_R, :m_A, :m_B, :c_rh, :c_t, :c_m, :c_q, :c_R, :c_A, :c_B, :et, :em, :q, :R, :A, :B, :rh, :rt, :rd, :z_rh, :z_t, :z_m, :z_q, :z_R, :z_A, :z_B, :si, :a]
tspan = (0,1e9)

solu = sol(combined_model, comb_init, tspan, comb_params)
# solu = simple_solve!(combined_model, comb_init, tspan, comb_params)

df = create_solu_df(solu, comb_species)

plot([scatter(x=df.time, y=col, name="$(names(df)[i])") for (col, i) in zip(eachcol(df[:,2:end]), range(2,length(names(df))))], Layout(xaxis_type="log"))

alpha = df.rt/kr 
fa = @. (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6)
ra = @. fa*df.R

ω_p(w_x, θ_x) = @. w_x*df.a/(θ_x + df.a)
plot(scatter(x=df.time, y=ω_p(w_rh, θ_nr)), Layout(xaxis_type="log"))

ω_q(θ_x) = @. (w_q*df.a/(θ_x + df.a))/(1+(df.q/Kq)^hq) # transcription rate for q # uM min-1
plot(scatter(x=df.time, y=ω_q(θ_nr)), Layout(xaxis_type="log"))

ω_rtcBA(θ_x) = @. (w_BA*df.a/(θ_x + df.a)) * ra*Vmax_init*df.a/(Km_init+df.a) # transcription rate for RtcBA # uM min-1 ???
plot(scatter(x=df.time, y=ω_rtcBA(θ_nr)), Layout(xaxis_type="log"))

gamma = @. gmax*df.a/(Kgamma+df.a) # aa min-1
v_x(c_x, nx) = @. df[:,c_x]*gamma/nx # uM min-1
ttrate = @. (df.c_q+df.c_rh+df.c_t+df.c_m+df.c_R+df.c_A+df.c_B)*gamma
plot(scatter(x=df.time, y=v_x(:c_rh, nrh)), Layout(xaxis_type="log"))
plot(scatter(x=df.time, y=v_x(:c_t, nx)), Layout(xaxis_type="log"))
plot(scatter(x=df.time, y=v_x(:c_B, nx)), Layout(xaxis_type="log"))
plot(scatter(x=df.time, y=ttrate), Layout(xaxis_type="log"))

λ = @. ttrate/M # min-1
plot(scatter(x=df.time, y=λ), Layout(xaxis_type="log"))

dil(x) = @. λ*x # uM min-1 
plot(scatter(x=df.time, y=dil(df.a)), Layout(xaxis_type="log"))
plot(scatter(x=df.time, y=dil(df.rh)), Layout(xaxis_type="log"))
plot(scatter(x=df.time, y=dil(df.B)), Layout(xaxis_type="log"))

deg(x) = @. d*x # uM min-1 
plot(scatter(x=df.time, y=deg(df.m_q)), Layout(xaxis_type="log"))
plot(scatter(x=df.time, y=deg(df.m_rh)), Layout(xaxis_type="log"))
plot(scatter(x=df.time, y=deg(df.m_B)), Layout(xaxis_type="log"))

rh_bind(m_x) = @. kb*df.rh*m_x # uM min-1 
plot(scatter(x=df.time, y=rh_bind(df.m_rh)), Layout(xaxis_type="log"))
plot(scatter(x=df.time, y=rh_bind(df.m_B)), Layout(xaxis_type="log"))

rh_unbind(c_x) = @. ku*c_x # uM min-1 
plot(scatter(x=df.time, y=rh_unbind(df.c_B)), Layout(xaxis_type="log"))
plot(scatter(x=df.time, y=rh_unbind(df.c_rh)), Layout(xaxis_type="log"))

zm(c_x) = @. c_x*abx*kon # uM min-1
plot(scatter(x=df.time, y=zm(df.c_B)), Layout(xaxis_type="log"))
plot(scatter(x=df.time, y=zm(df.c_A)), Layout(xaxis_type="log"))

zm_diss(z_x) = @. koff*z_x # uM min-1 
plot(scatter(x=df.time, y=zm_diss(df.z_B)), Layout(xaxis_type="log"))

vimp = @. (df.et*vt*s0/(Kt+s0)) # uM min-1 
vcat = @. (df.em*vm*df.si/(Km+df.si)) # uM min-1 
plot(scatter(x=df.time, y=vcat), Layout(xaxis_type="log"))

A1 = @. (df.a*df.A)/(df.a+(km_a*df.rd)) # uM
B1 = @. (df.a*df.B)/(df.a+(km_b*df.rt)) # uM
plot(scatter(x=df.time, y=A1), Layout(xaxis_type="log"))
plot(scatter(x=df.time, y=B1), Layout(xaxis_type="log"))

Vrep = @. krep*B1*df.rt # uM min-1 technically should be uM^2 min-1
Vdam = @. (df.z_rh + df.z_t + df.z_m + df.z_q + df.z_R + df.z_A + df.z_B)*koff*kdam # uM min-1
Vtag = @. ktag*A1*df.rd # uM min-1 
plot(scatter(x=df.time, y=Vrep), Layout(xaxis_type="log"))
plot(scatter(x=df.time, y=Vdam), Layout(xaxis_type="log"))
plot(scatter(x=df.time, y=Vtag), Layout(xaxis_type="log"))

