using CSV, TotalLeastSquares, Statistics, BayesOpt, BlackBoxOptim, PyCall, StaticArrays
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

function plot_individual_data(csv, csv_std, title, type)
    p = make_subplots(rows=4, cols=1, shared_xaxes=true, vertical_spacing=0.06, y_title="$type", x_title="Time (hours)", subplot_titles=["RtcA" "RtcB" "RtcR" "OmpA"])
    add_trace!(p, (scatter(x=t, y=csv.RtcA, error_y=attr(type="data", array=csv_std.RtcA))), row=1, col=1)
    add_trace!(p, (scatter(x=t, y=csv.RtcB, error_y=attr(type="data", array=csv_std.RtcB))), row=2, col=1)
    add_trace!(p, (scatter(x=t, y=csv.RtcR, error_y=attr(type="data", array=csv_std.RtcR))), row=3, col=1)
    add_trace!(p, (scatter(x=t, y=csv.OmpA, error_y=attr(type="data", array=csv_std.OmpA))), row=4, col=1)
    relayout!(p, showlegend=false, title_text="$title")
    return p
end

p_conc_a = plot_individual_data(csv_a_conc, csv_a_std_conc, "Dataset A conc", "Concentration (μM)")
p_conc_b = plot_individual_data(csv_b_conc, csv_b_std_conc, "Dataset B conc", "Concentration (μM)")
p_copy_a = plot_individual_data(csv_a_copy, csv_a_std_copy, "Dataset A copy number", "Copy number")
p_copy_b = plot_individual_data(csv_b_copy, csv_b_std_copy, "Dataset B copy number", "Copy number")


function plot_all_data()
    rtca_copy_a = scatter(x=t, y=csv_a_copy.RtcA, error_y=attr(type="data", array=csv_a_std_copy.RtcA), name="RtcA")
    rtcb_copy_a = scatter(x=t, y=csv_a_copy.RtcB, error_y=attr(type="data", array=csv_a_std_copy.RtcB), name="RtcB")
    rtcr_copy_a = scatter(x=t, y=csv_a_copy.RtcR, error_y=attr(type="data", array=csv_a_std_copy.RtcR), name="RtcR")
    ompa_copy_a = scatter(x=t, y=csv_a_copy.OmpA, error_y=attr(type="data", array=csv_a_std_copy.OmpA), name="OmpA")
    p_copy_a = plot([rtca_copy_a, rtcb_copy_a, rtcr_copy_a, ompa_copy_a], Layout(title="Copy number samples A", yaxis_type="log", xaxis_title="Time (hours)", yaxis_title="Copy number"))

    rtca_copy_b = scatter(x=t, y=csv_b_copy.RtcA, error_y=attr(type="data", array=csv_b_std_copy.RtcA), name="RtcA")
    rtcb_copy_b = scatter(x=t, y=csv_b_copy.RtcB, error_y=attr(type="data", array=csv_b_std_copy.RtcB), name="RtcB")
    rtcr_copy_b = scatter(x=t, y=csv_b_copy.RtcR, error_y=attr(type="data", array=csv_b_std_copy.RtcR), name="RtcR")
    ompa_copy_b = scatter(x=t, y=csv_b_copy.OmpA, error_y=attr(type="data", array=csv_b_std_copy.OmpA), name="OmpA")
    p_copy_b = plot([rtca_copy_b, rtcb_copy_b, rtcr_copy_b, ompa_copy_b], Layout(title="Copy number samples B", yaxis_type="log", xaxis_title="Time (hours)", yaxis_title="Copy number"))

    rtca_conc_a = scatter(x=t, y=csv_a_conc.RtcA, error_y=attr(type="data", array=csv_a_std_conc.RtcA), name="RtcA")
    rtcb_conc_a = scatter(x=t, y=csv_a_conc.RtcB, error_y=attr(type="data", array=csv_a_std_conc.RtcB), name="RtcB")
    rtcr_conc_a = scatter(x=t, y=csv_a_conc.RtcR, error_y=attr(type="data", array=csv_a_std_conc.RtcR), name="RtcR")
    ompa_conc_a = scatter(x=t, y=csv_a_conc.OmpA, error_y=attr(type="data", array=csv_a_std_conc.OmpA), name="OmpA")
    p_conc_a = plot([rtca_conc_a, rtcb_conc_a, rtcr_conc_a, ompa_conc_a], Layout(title="Concentration (μM) samples A", yaxis_type="log", xaxis_title="Time (hours)", yaxis_title="Concentration (μM)"))

    rtca_conc_b = scatter(x=t, y=csv_b_conc.RtcA, error_y=attr(type="data", array=csv_b_std_conc.RtcA), name="RtcA")
    rtcb_conc_b = scatter(x=t, y=csv_b_conc.RtcB, error_y=attr(type="data", array=csv_b_std_conc.RtcB), name="RtcB")
    rtcr_conc_b = scatter(x=t, y=csv_b_conc.RtcR, error_y=attr(type="data", array=csv_b_std_conc.RtcR), name="RtcR")
    ompa_conc_b = scatter(x=t, y=csv_b_conc.OmpA, error_y=attr(type="data", array=csv_b_std_conc.OmpA), name="OmpA")
    p_conc_b = plot([rtca_conc_b, rtcb_conc_b, rtcr_conc_b, ompa_conc_b], Layout(title="Concentration (μM) samples A", yaxis_type="log", xaxis_title="Time (hours)", yaxis_title="Concentration (μM)"))

    return [p_copy_a p_conc_a; p_copy_b p_conc_b]
end
plot_all_data()

rtca_a = csv_a.RtcA
rtca_b = csv_b.RtcA
rtca_a_sd = csv_a_std.RtcA
rtca_b_sd = csv_b_std.RtcA
rtcb_a = csv_a.RtcB
rtcb_b = csv_b.RtcB
rtcb_a_sd = csv_a_std.RtcB
rtcb_b_sd = csv_b_std.RtcB
rtcr_a = csv_a.Rtc,
rtcr_b = csv_b.RtcR
rtcr_a_sd = csv_a_std.RtcR
rtcr_b_sd = csv_b_std.RtcR

csv_gr = DataFrame(CSV.File("$PATHdata/results_colD_grfit.csv")) # read csv to a datafram
csv_gr = select!(csv_gr, Not(["log(OD)", "log(OD) error", "gr error", "od"]))

# set atp and lam curves from files 
lam_t, new_df = extend_gr_curve(csv_gr)
lam_t[lam_t.< 0] .= 0 #zero(eltype(lam_colD))

tspan2 = (0, 1440)
t_2 = [4,6,8,10,12,24]*60





data = vcat(rtca_a, rtcb_a, rtcr_a, rtca_b, rtcb_b, rtcr_b)
time = hcat(t,t,t,t,t,t)


atp = 4000; kin = 0.054;
params = @LArray [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam_t] (:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam)
solu = sol(rtc_model1!, initial, tspan2, params)

rm_a = get_curve(solu, :rm_a)
rm_b = get_curve(solu, :rm_b)
rm_r = get_curve(solu, :rm_r)

# compute the residuals between the model and the data using TLS
b = vcat(rm_a, rm_b, rm_r)
@show length(b)
B = [time ones(size(time,1))] \ [b data]




# Define the TLS objective function
function TLS_objective(time, data, rtc_model1!, initial, params, tspan2, t_2)

    # solve the ODE system using the given parameter values
    solu = sol_with_t(rtc_model1!, initial, params, tspan2, t_2)

    rm_a = get_curve(solu, :rm_a)
    rm_b = get_curve(solu, :rm_b)
    rm_r = get_curve(solu, :rm_r)
    
    # compute the residuals between the model and the data using TLS
    b = vcat(rm_a, rm_b, rm_r)
    @show length(b)
    B = [time ones(size(time,1))] \ [b data]
    return norm([b data] - [time ones(size(time,1))] * B)
end

dataset_a = [rtca_a, rtcb_a, rtcr_a]
dataset_a_var = [rtca_a_sd.^2, rtcb_a_sd.^2, rtcr_a_sd.^2]
dataset_b = [rtca_b, rtcb_b, rtcr_b]
dataset_b_var = [rtca_b_sd.^2, rtcb_b_sd.^2, rtcr_b_sd.^2]

collect(1:3)

function sse(rtc_model1!, initial, params, tspan2, t_2)
    solu = sol_with_t(rtc_model1!, initial, params, tspan2, t_2)

    rm_a = get_curve(solu, :rm_a)
    rm_b = get_curve(solu, :rm_b)
    rm_r = get_curve(solu, :rm_r)
    
    model_data = [rm_a, rm_b, rm_r]
    error = 0
    for rtc in collect(1:3)
        error_a = sum([abs2((i-j)/k) for (i,j,k) in zip(dataset_a[rtc], model_data[rtc], dataset_a_var[rtc])])
        error_b = sum([abs2((i-j)/k) for (i,j,k) in zip(dataset_b[rtc], model_data[rtc], dataset_b_var[rtc])])
        tot_error = error_a + error_b
        error += tot_error
    end
    @show error
end


function bo(;ω_ab, ω_r)
    obj = sse(rtc_model1!, initial, [L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam_t], tspan2, t_2)
    return -obj/216
end
py"""
param_range_ω = {'ω_ab': (0, 100), 'ω_r': (0, 100)}
"""


# import bayes_opt package from python
bayes_opt = pyimport("bayes_opt")

# setting the optimizer 
optimizer = bayes_opt.BayesianOptimization(f=bo, pbounds=py"param_range_ω", random_state=27, verbose=2) # verbose = 1 prints only when a maximum is observed (pink)

# timing the process and maximising the optimizer 
function timer()
    optimizer.maximize(init_points=2, n_iter=150, acq="ei", xi=0.1) #kappa=2, xi = 0.0 (prefer exploitation), xi = 0.1 (prefer exploration)
end

@time timer()
print(optimizer.max)













function tls_objective(x::Vector, X::Matrix, sigma::Vector)
    # Compute the TLS error for a given set of parameters x
    # X is your independent variable(s)
    # sigma is the standard deviation values for each data point
    # add bias term, if necessary
    X = hcat(X, ones(size(X, 1)))
    A, residuals, rank, s = ldiv!(qr(X / diagm(sqrt.(sigma))), y / sqrt.(sigma))
    error = sqrt(mean(residuals.^2))
    return error
end

csv_gr = DataFrame(CSV.File("$PATHdata/results_colD_grfit.csv")) # read csv to a datafram
csv_gr = select!(csv_gr, Not(["log(OD)", "log(OD) error", "gr error", "od"]))

# set atp and lam curves from files 
lam_t, new_df = extend_gr_curve(csv_gr)
lam_t[lam_t.< 0] .= 0 #zero(eltype(lam_colD))

tspan2 = (0, 1440)
t_2 = [4,6,8,10,12,24]*60

function rtc_bo_ω(;ω_ab, ω_r)
    # @show ω_ab, ω_r
    obj_wt_ab = obj(rtc_model1!, initial, ([L, c, kr, Vmax_init, Km_init, ω_ab, ω_r, θtscr, g_max, θtlr, km_a, km_b, d, krep, kdam, ktag, kdeg, kin, atp, na, nb, nr, lam_t]), tspan2, t_2, "mrna", rtca_a, rtca_a_sd)
    return -(obj_wt_ab/36)
end


# in python writing the ranges to search for parameter value
py"""
param_range_ω = {'ω_ab': (0, 50), 'ω_r': (0, 50)}
"""


# import bayes_opt package from python
bayes_opt = pyimport("bayes_opt")

# setting the optimizer 
optimizer = bayes_opt.BayesianOptimization(f=rtc_bo_ω, pbounds=py"param_range_ω", random_state=27, verbose=2) # verbose = 1 prints only when a maximum is observed (pink)

# timing the process and maximising the optimizer 
function timer()
    optimizer.maximize(init_points=2, n_iter=150, acq="ei", xi=0.1) #kappa=2, xi = 0.0 (prefer exploitation), xi = 0.1 (prefer exploration)
end

@time timer()
print(optimizer.max)






