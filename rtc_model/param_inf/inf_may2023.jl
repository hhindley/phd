using CSV, DifferentialEquations, Statistics, BayesOpt, BlackBoxOptim, PyCall, StaticArrays, LabelledArrays, DataInterpolations, PlotlyJS
include("$PATHrtc_models/Oct2022_model/rtc_model.jl")
include("$PATHrtc_models/sol_species_funcs.jl")
include("$PATHrtc_models/params_init_tspan.jl")
include("$PATHParam_inf/inf_setup.jl")


csv_a_conc = DataFrame(CSV.File("$PATHdata/df_final_conc_a.csv"))
csv_b_conc = DataFrame(CSV.File("$PATHdata/df_final_conc_b.csv"))
t = [4,6,8,10,12,24]
csv_a_std_conc = DataFrame(CSV.File("$PATHdata/df_final_conc_a_sd.csv"))  
csv_b_std_conc = DataFrame(CSV.File("$PATHdata/df_final_conc_b_sd.csv"))  

csv_a_copy = DataFrame(CSV.File("$PATHdata/df_final_copy_a.csv"))
csv_b_copy = DataFrame(CSV.File("$PATHdata/df_final_copy_b.csv"))
csv_a_std_copy = DataFrame(CSV.File("$PATHdata/df_final_copy_a_sd.csv"))  
csv_b_std_copy = DataFrame(CSV.File("$PATHdata/df_final_copy_b_sd.csv"))  


rtca_a = csv_a_conc.RtcA
rtca_b = csv_b_conc.RtcA
rtca_a_sd = csv_a_std_conc.RtcA
rtca_b_sd = csv_b_std_conc.RtcA
rtcb_a = csv_a_conc.RtcB
rtcb_b = csv_b_conc.RtcB
rtcb_a_sd = csv_a_std_conc.RtcB
rtcb_b_sd = csv_b_std_conc.RtcB
rtcr_a = csv_a_conc.RtcR
rtcr_b = csv_b_conc.RtcR
rtcr_a_sd = csv_a_std_conc.RtcR
rtcr_b_sd = csv_b_std_conc.RtcR

# set atp and lam curves from files 
csv_atp = DataFrame(CSV.File("$PATHdata/atp_for_rtcmodel.csv"))

csv_atp.atp = csv_atp.atp/5
rh = 11.29
lam_t = QuadraticInterpolation(csv_atp."gr",csv_atp."t")
atp_t = QuadraticInterpolation(csv_atp."atp",csv_atp."t")
kin_model = @. lam_t*rh/g_max
kin_t = QuadraticInterpolation(kin_model, csv_atp."t") 


tspan2 = (0, 1440)
t_2 = [4,6,8,10,12,24]*60

function plot_all_data()
    rtca_conc_a = scatter(x=t*60, y=csv_a_conc.RtcA, error_y=attr(type="data", array=csv_a_std_conc.RtcA), name="RtcA")
    rtcb_conc_a = scatter(x=t*60, y=csv_a_conc.RtcB, error_y=attr(type="data", array=csv_a_std_conc.RtcB), name="RtcB")
    rtcr_conc_a = scatter(x=t*60, y=csv_a_conc.RtcR, error_y=attr(type="data", array=csv_a_std_conc.RtcR), name="RtcR")
    p_conc_a = plot([rtca_conc_a, rtcb_conc_a, rtcr_conc_a], Layout(title="Concentration (μM) samples A", xaxis_title="Time (hours)", yaxis_title="Concentration (μM)"))

    rtca_conc_b = scatter(x=t*60, y=csv_b_conc.RtcA, error_y=attr(type="data", array=csv_b_std_conc.RtcA), name="RtcA")
    rtcb_conc_b = scatter(x=t*60, y=csv_b_conc.RtcB, error_y=attr(type="data", array=csv_b_std_conc.RtcB), name="RtcB")
    rtcr_conc_b = scatter(x=t*60, y=csv_b_conc.RtcR, error_y=attr(type="data", array=csv_b_std_conc.RtcR), name="RtcR")
    p_conc_b = plot([rtca_conc_b, rtcb_conc_b, rtcr_conc_b], Layout(title="Concentration (μM) samples A", xaxis_title="Time (hours)", yaxis_title="Concentration (μM)"))

    return [p_conc_a; p_conc_b]
end

inf_data = plot_all_data()
open("./data_figs/inf_data.html", "w") do io
    PlotlyBase.to_html(io, inf_data.plot)
end

ω_ab = 0.01; ω_r = 0.01;
params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin_t, atp_t, na, nb, nr, lam_t] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
solu = sol_with_t(rtc_all_t!, initial, params, tspan2, t_2)
plotly_plot_sol(solu, "", "", "sol")
rm_a = get_curve(solu, :rm_a)
rm_b = get_curve(solu, :rm_b)
rm_r = get_curve(solu, :rm_r)


dataset_a = [rtca_a, rtcb_a, rtcr_a]
dataset_a_var = [rtca_a_sd.^2, rtcb_a_sd.^2, rtcr_a_sd.^2]
dataset_b = [rtca_b[2:6], rtcb_b[2:6], rtcr_b[2:6]]
dataset_b_var = [(rtca_b_sd.^2)[2:6], (rtcb_b_sd.^2)[2:6], (rtcr_b_sd.^2)[2:6]]

# model_data_a = [rm_a, rm_b, rm_r]
# model_data_b = [rm_a[2:6], rm_b[2:6], rm_r[2:6]]
# dataset_a = [rtca_a, rtcb_a]#, rtcr_a]
# dataset_a_var = [rtca_a_sd.^2, rtcb_a_sd.^2]#, rtcr_a_sd.^2]

# error = 0
# for rtc in collect(1:3)
#     error_a = sum([(k*abs2(i-j))/k*6 for (i,j,k) in zip(dataset_a[rtc], model_data_a[rtc], dataset_a_var[rtc])])
#     error_b = sum([(k*abs2(i-j))/k*6 for (i,j,k) in zip(dataset_b[rtc], model_data_b[rtc], dataset_b_var[rtc])])

#     # error_a = sum([abs2((i-j)/k) for (i,j,k) in zip(dataset_a[rtc], model_data[rtc], dataset_a_var[rtc])])
#     # error_b = sum([abs2((i-j)/k) for (i,j,k) in zip(dataset_b[rtc], model_data[rtc], dataset_b_var[rtc])])
#     tot_error = error_a + error_b
#     error += tot_error
# end
# error

# param_range = collect(0:0.01:0.3)
# sweep_paramx2_few_t(rtc_all_t!, tspan2, t_2, lam_t, atp_t, kin_t, :rm_a, get_ssval, :ω_ab, :ω_r, param_range, param_range)
# sweep_paramx2_few_t(rtc_all_t!, tspan2, t_2, lam_t, atp_t, kin_t, :rm_r, get_ssval, :ω_ab, :ω_r, param_range, param_range)

# dataset_a = [rtca_a, rtcb_a]#, rtcr_a]
# dataset_a_var = [rtca_a_sd.^2, rtcb_a_sd.^2]#, rtcr_a_sd.^2]
# dataset_b = [rtca_b, rtcb_b]#, rtcr_b]
# dataset_b_var = [rtca_b_sd.^2, rtcb_b_sd.^2]#, rtcr_b_sd.^2]

# collect(1:3)

function sse(model, initial, params, tspan2, t_2)
    solu = sol_with_t(model, initial, params, tspan2, t_2)

    rm_a = get_curve(solu, :rm_a)
    rm_b = get_curve(solu, :rm_b)
    rm_r = get_curve(solu, :rm_r)
    
    model_data_a = [rm_a, rm_b, rm_r]
    model_data_b = [rm_a[2:6], rm_b[2:6], rm_r[2:6]]
    
    error = 0
    for rtc in collect(1:2)
        error_a = sum([(abs2(i-j))*k for (i,j,k) in zip(dataset_a[rtc], model_data_a[rtc], dataset_a_var[rtc])])
        error_b = sum([(abs2(i-j))*k for (i,j,k) in zip(dataset_b[rtc], model_data_b[rtc], dataset_b_var[rtc])])
        tot_error = error_a + error_b
        error += tot_error
    end
    return error
end

function wmse(model, initial, params, tspan2, t_2)
    solu = sol_with_t(model, initial, params, tspan2, t_2)

    rm_a = get_curve(solu, :rm_a)
    rm_b = get_curve(solu, :rm_b)
    rm_r = get_curve(solu, :rm_r)
    
    model_data_a = [rm_a, rm_b, rm_r]
    model_data_b = [rm_a[2:6], rm_b[2:6], rm_r[2:6]]

    error = 0
    for rtc in collect(1:3)
        error_a = sum([(k*abs2(i-j))/k*6 for (i,j,k) in zip(dataset_a[rtc], model_data_a[rtc], dataset_a_var[rtc])])
        error_b = sum([(k*abs2(i-j))/k*6 for (i,j,k) in zip(dataset_b[rtc], model_data_b[rtc], dataset_b_var[rtc])])
        tot_error = error_a + error_b
        error += tot_error
    end
    return error
end


function bo_sse(;ω_ab, ω_r)
    obj = sse(rtc_all_t!, initial, [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin_t, atp, na, nb, nr, lam_t], tspan2, t_2)
    return -obj/216
end

function bo_wmse(;ω_ab, ω_r)
    obj = wmse(rtc_all_t!, initial, [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin_t, atp_t, na, nb, nr, lam_t], tspan2, t_2)
    return -obj
end

py"""
param_range_ω = {'ω_ab': (0, 1), 'ω_r': (0, 1)}
"""


# import bayes_opt package from python
bayes_opt = pyimport("bayes_opt")

# setting the optimizer 
optimizer_sse = bayes_opt.BayesianOptimization(f=bo_sse, pbounds=py"param_range_ω")#, random_state=27, verbose=2) # verbose = 1 prints only when a maximum is observed (pink)
optimizer_wmse = bayes_opt.BayesianOptimization(f=bo_wmse, pbounds=py"param_range_ω")#, random_state=27, verbose=2) # verbose = 1 prints only when a maximum is observed (pink)

# timing the process and maximising the optimizer 
function timer()
    optimizer_sse.maximize(init_points=2, n_iter=100, acq="ei", xi=0.01) #kappa=2, xi = 0.0 (prefer exploitation), xi = 0.1 (prefer exploration)
end

@time timer()
print(optimizer_sse.max)



using BlackBoxOptim
function rtc_bo_ω(x)
    # obj = wmse(rtc_all_t!, initial, [L, c, kr, Vmax_init, Km_init, x[1], x[2], θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin_t, atp_t, na, nb, nr, lam_t], tspan2, t_2)

    obj = sse(rtc_model1!, initial, [L, c, kr, Vmax_init, Km_init, x[1], x[2], θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam_t], tspan2, t_2)

    return (obj)
end

res = bboptimize(rtc_bo_ω; SearchRange = [(0, 1), (0, 1)], Method = :adaptive_de_rand_1_bin, MaxSteps = 1000)


print("ω_ab = ", best_candidate(res)[1])
print("ω_r = ", best_candidate(res)[2])
print("error = ", best_fitness(res))



using BayesOpt

config = ConfigParameters()         # initiates parameters 

config.n_iterations = 300
config.n_iter_relearn = 50
config.n_init_samples = 500     # how to change parameter values from original
config.noise = 1e-14

optimizer, optimum = bayes_optimization(rtc_bo_ω, [0,0], [1,1], config)

optimum

using BayesianOptimization, GaussianProcesses, Distributions

# Choose as a model an elastic GP with input dimensions 2.
# The GP is called elastic, because data can be appended efficiently.
model = ElasticGPE(2,                            # 2 input dimensions
                   mean = MeanConst(0.),         
                   kernel = SEArd([0., 0.], 5.),
                   logNoise = 0.,
                   capacity = 3000)              # the initial capacity of the GP is 3000 samples.
set_priors!(model.mean, [Normal(1, 2)])

# Optimize the hyperparameters of the GP using maximum a posteriori (MAP) estimates every 50 steps
modeloptimizer = NoModelOptimizer()#every = 50, noisebounds = [-4, 3],       # bounds of the logNoise
                                # kernbounds = [[-1, -1, 0], [4, 4, 10]],  # bounds of the 3 parameters GaussianProcesses.get_param_names(model.kernel)
                                # maxeval = 40)
opt = BOpt(rtc_bo_ω,
           model,
           ExpectedImprovement(),                   # type of acquisition
           modeloptimizer,                        
           [0,0], [1,1],                     # lowerbounds, upperbounds         
           repetitions = 1,                          # evaluate the function for each input 5 times
           maxiterations = 300,                      # evaluate at 100 input positions
           sense = Min,                              # minimize the function
           acquisitionoptions = (method = :LD_LBFGS, # run optimization of acquisition function with NLopts :LD_LBFGS method
                                 restarts = 5,       # run the NLopt method from 5 random initial conditions each time.
                                 maxtime = 0.1,      # run the NLopt method for at most 0.1 second each time
                                 maxeval = 1000),    # run the NLopt methods for at most 1000 iterations (for other options see https://github.com/JuliaOpt/NLopt.jl)
            verbosity = Progress)

result = boptimize!(opt)









function plot_data_vs_model(solu, log, vary)
    rm_a = get_curve(solu, :rm_a)
    rm_b = get_curve(solu, :rm_b)
    rm_r = get_curve(solu, :rm_r)
    model_rtca = scatter(x=solu.t, y=rm_a, name="Model RtcA")
    model_rtcb = scatter(x=solu.t, y=rm_b, name="Model RtcB")
    model_rtcr = scatter(x=solu.t, y=rm_r, name="Model RtcR")

    rtca_conc_a = scatter(x=t*60, y=csv_a_conc.RtcA, error_y=attr(type="data", array=csv_a_std_conc.RtcA), name="RtcA_a")
    rtcb_conc_a = scatter(x=t*60, y=csv_a_conc.RtcB, error_y=attr(type="data", array=csv_a_std_conc.RtcB), name="RtcB_a")
    rtcr_conc_a = scatter(x=t*60, y=csv_a_conc.RtcR, error_y=attr(type="data", array=csv_a_std_conc.RtcR), name="RtcR_a")
    rtca_conc_b = scatter(x=t[2:6]*60, y=csv_b_conc.RtcA[2:6], error_y=attr(type="data", array=csv_b_std_conc.RtcA[2:6]), name="RtcA_b")
    rtcb_conc_b = scatter(x=t[2:6]*60, y=csv_b_conc.RtcB[2:6], error_y=attr(type="data", array=csv_b_std_conc.RtcB[2:6]), name="RtcB_b")
    rtcr_conc_b = scatter(x=t[2:6]*60, y=csv_b_conc.RtcR[2:6], error_y=attr(type="data", array=csv_b_std_conc.RtcR[2:6]), name="RtcR_b")


    return plot([model_rtca, model_rtcb, model_rtcr, rtca_conc_b, rtcb_conc_b, rtcr_conc_b, rtca_conc_a, rtcb_conc_a, rtcr_conc_a], Layout(title="Data vs. model - $vary", xaxis_title="Time (mins)", yaxis_title="Concentration (μM)", yaxis_type=log))
end



function rtc_bo_ω(x)
    # obj = wmse(rtc_all_t!, initial, [L, c, kr, Vmax_init, Km_init, x[1], x[2], θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin_t, atp_t, na, nb, nr, lam_t], tspan2, t_2)

    obj = sse(rtc_model1!, initial, [L, c, kr, Vmax_init, Km_init, x[1], x[2], θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam_t], tspan2, t_2)

    return (obj)
end

res = bboptimize(rtc_bo_ω; SearchRange = [(0, 1), (0, 1)], Method = :adaptive_de_rand_1_bin, MaxSteps = 1000)


print("ω_ab = ", best_candidate(res)[1])
print("ω_r = ", best_candidate(res)[2])
print("error = ", best_fitness(res))

ω_ab = best_candidate(res)[1]; ω_r = best_candidate(res)[2];

ω_ab = 0.00032; ω_r = 0.00012;
lam = 0.033; kin = 0.054; atp = 4000;

include("$PATHrtc_models/params_init_tspan.jl")

param = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin_t, atp, na, nb, nr, lam_t] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
solu = sol_with_t(rtc_all_t!, initial, param, tspan2, t_2)

tspan=(0,1e9)
solu = sol(lam_kin_t, initial, tspan, param)

plotly_plot_sol(solu, "log", "", "sol")

data_vs_model = plot_data_vs_model(solu, "log", "all_t")

savefig(data_vs_model, "data_vs_model.svg")


sopen("./data_figs/data_vs_model.html", "w") do io
    PlotlyBase.to_html(io, data_vs_model.plot)
end




