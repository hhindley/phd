using DifferentialEquations, PyCall, DataFrames, Measures, CSV, DiffEqCallbacks, DataInterpolations, LabelledArrays, Statistics, OrderedCollections
using PlotlyJS


include("model.jl")
include("parameters.jl")
include("functions.jl")

function gr_and_atp(solu_ss)
    M = 10^8
    cr = get_curve(solu_ss, :cr); cp = get_curve(solu_ss, :cp); cq = get_curve(solu_ss, :cq); ct = get_curve(solu_ss, :ct); cm = get_curve(solu_ss, :cm);
    Rt = @. cr+cp+cq+ct+cm
    g_max = 1260
    k_g = 7*sf
    a1 = get_curve(solu_ss, :a)
    SF = 1e6/(6.022e23*1e-15)
    a = a1.*SF
    lam = @. (Rt*(g_max*a/(k_g+a)))/(M)
    return lam, a
end
function max_vals_param_sweep(parameter, param_range, ssinit)
    params = @LArray [dm, kb, ku, f, thetar, gmax, thetax, Kt, M, we, Km, vm, nx, Kq, vt, wr, wq, wp, nq, nr, ns, kin, d_s, d_n, Kgamma] (:dm, :kb, :ku, :f, :thetar, :gmax, :thetax, :Kt, :M, :we, :Km, :vm, :nx, :Kq, :vt, :wr, :wq, :wp, :nq, :nr, :ns, :kin, :d_s, :d_n, :Kgamma)
    max_vals = OrderedDict(val => [Dict("λ"=>[],"ATP"=>[])] for val in param_range)
    for (val, i) in zip(param_range, values(max_vals))
        params[parameter] = val
        solu_ss = simple_solve!(pop_model!, ssinit, tspan, params)
        lam, a = gr_and_atp(solu_ss)
        # for i in values(max_vals)
        push!(i[1]["λ"], maximum(lam))
        push!(i[1]["ATP"], maximum(a))

        # end 
    end
    return max_vals
end
function max_vals_init_sweep(param_range)
    params = @LArray [dm, kb, ku, f, thetar, gmax, thetax, Kt, M, we, Km, vm, nx, Kq, vt, wr, wq, wp, nq, nr, ns, kin, d_s, d_n, Kgamma] (:dm, :kb, :ku, :f, :thetar, :gmax, :thetax, :Kt, :M, :we, :Km, :vm, :nx, :Kq, :vt, :wr, :wq, :wp, :nq, :nr, :ns, :kin, :d_s, :d_n, :Kgamma)
    max_vals = OrderedDict(val => [Dict("λ"=>[],"ATP"=>[])] for val in param_range)
    for (s0, i) in zip(param_range, values(max_vals))
        ssinit = get_ss_init(solu, s0, N_0)
        solu_ss = simple_solve!(pop_model!, ssinit, tspan, params)
        lam, a = gr_and_atp(solu_ss)
        # for i in values(max_vals)
        push!(i[1]["λ"], maximum(lam))
        push!(i[1]["ATP"], maximum(a))

        # end 
    end
    return max_vals
end

tspan = (0,1e9)

solu = simple_solve!(odemodel!, init, tspan, params)
plotly_plot_sol(solu, "", "")

s_0 = 8e7; N_0 = 1;
ssinit = get_ss_init(solu, s_0, N_0);
# print(ssinit)

tspan = (0,1e9);
# better to have smaller sf for growth rate 
ns = 3.5; sf = 10000; Kgamma = 7*sf; # change kgamma to change onset of atp/growth rate
pop_params = [dm, kb, ku, f, thetar, gmax, thetax, Kt, M, we, Km, vm, nx, Kq, vt, wr, wq, wp, nq, nr, ns, kin, d_s, d_n, Kgamma ];

solu_ss = simple_solve!(pop_model!, ssinit, tspan, pop_params);
# plotly_plot_sol_pop(solu_ss, "log", "")

lam, a = gr_and_atp(solu_ss);
# print("λ = ", maximum(lam), " and ATP = ", maximum(a))

# s0_range = collect(5e7:500000:8e7)
# max_vals = max_vals_init_sweep(s0_range)

# println(max_vals)




# function sweep_paramx3(param1, param2, param_range1, param_range2, param_range3)
#     all_res = []
#     params = @LArray [dm, kb, ku, f, thetar, gmax, thetax, Kt, M, we, Km, vm, nx, Kq, vt, wr, wq, wp, nq, nr, ns, kin, d_s, d_n, Kgamma] (:dm, :kb, :ku, :f, :thetar, :gmax, :thetax, :Kt, :M, :we, :Km, :vm, :nx, :Kq, :vt, :wr, :wq, :wp, :nq, :nr, :ns, :kin, :d_s, :d_n, :Kgamma)
#     for v in param_range1
#         params[param1] = v
#         res1 = []
#         for i in param_range2
#             params[param2] = i
#             max_vals = OrderedDict(val => [Dict("λ"=>[],"ATP"=>[])] for val in param_range)
#             for (s0, i) in zip(param_range3, values(max_vals))
#                 ssinit = get_ss_init(solu, s0, N_0)
#                 solu_ss = simple_solve!(pop_model!, ssinit, tspan, params)
#                 lam, a = gr_and_atp(solu_ss)
#                 push!(i[1]["λ"], maximum(lam))
#                 push!(i[1]["ATP"], maximum(a))            
#             end
#             push!(res1, max_vals)
#         end
#         push!(all_res, res1)
#     end
#     return all_res
# end

# function find_good_param(max_vals, param_range, parameter)
#     for i in param_range
#         if 0.01 < max_vals[i][1]["λ"][1] && max_vals[i][1]["ATP"][1] < 40000
#             println("$parameter = ", i, " and ", max_vals[i][1])
#         end
#     end
# end

# sf_range = collect(1:1000:10000)
# Kgamma_range = [7*i for i in sf_range]
# param_range = collect(0:1.5:14)
# max_vals = max_vals_param_sweep(:ns, param_range, ssinit)
# find_good_param(max_vals, param_range, "ns")
# println(max_vals)

# max_vals_kgamma = max_vals_param_sweep(:Kgamma, Kgamma_range, ssinit)
# find_good_param(max_vals_kgamma, Kgamma_range, "Kgamma")

# res = sweep_paramx3(:ns, :Kgamma, param_range, Kgamma_range, s0_range)
# res[8]



# using FlexiMaps; 
# s0_range = maprange(log, 1, 1e10, length=10)
# max_vals_s0 = max_vals_init_sweep(s0_range)
# find_good_param(max_vals_s0, s0_range, "s0")

# p_atp1 = plot(scatter(x=solu_ss.t, y=a1), Layout(xaxis_type="log"))
p_atp = plot(scatter(x=solu_ss.t, y=a), Layout(xaxis_range=(0,200), xaxis_type="", title="ATP"))
p_lam = plot(scatter(x=solu_ss.t, y=lam), Layout(xaxis_range=(0,200), xaxis_type="", title="λ"))





# split data for interpolation 
first_a = a[1:argmax(a)]
last_a = a[argmax(a):end]
first_lam = lam[1:argmax(a)]
last_lam = lam[argmax(a):end]

real = (scatter(x=lam, y=a))
first = (scatter(x=first_lam, y=first_a))#, Layout(xaxis_range=(0,200), xaxis_type=""))
last = (scatter(x=last_lam, y=last_a))#, Layout(xaxis_range=(0,200), xaxis_type=""))

plot([first, last])

# load scipy.interpolate
np = pyimport("numpy")
interp = pyimport("scipy.interpolate")

#create interpolation objects 
int_first = interp.interp1d(first_lam, first_a)
int_last = interp.interp1d(last_lam, last_a)

#test to see if it returns the same curve as created with 
test_first = int_first(first_lam)
test_last = int_last(last_lam)

p1 = scatter(x=first_lam, y=test_first)
p2 = scatter(x=last_lam, y=test_last)

plot([p1,p2])

# load data growth rate 
csv_lam = DataFrame(CSV.File("/home/holliehindley/phd/data/results_colD_grfit.csv"))
csv_lam = select!(csv_lam, Not(["log(OD)", "log(OD) error", "gr error", "od"]))
# csv_lam."gr"[csv_lam."gr".< 0] .= 0 #zero(eltype(lam_colD))
plot(scatter(x=csv_lam."t", y=csv_lam."gr"))

function extend_gr_curve(csv)
    mean_gr = mean((csv."gr"[Int64((length(csv."t")*2/3)+1):end]))
    df = DataFrame(t=Float64[], gr=Float64[])
    for t in collect(csv."t"[end]+10:5000:1e9)
        push!(df, [t, csv_lam.gr[28]])
    end    
    new_df = vcat(csv, df)
    lam = QuadraticInterpolation(new_df."gr",new_df."t")
    return lam, new_df
end

lam_ext, new_df = extend_gr_curve(csv_lam)

lam_ext[lam_ext.< 0] .= csv_lam.gr[28]
plot(scatter(x=new_df."t", y=lam_ext), Layout(xaxis_range=(0,10000)))

# split data growth rate 
dlam = new_df."gr"
first_dlam = dlam[1:argmax(dlam)]
last_dlam = dlam[argmax(dlam):end]
lam1 = scatter(x=new_df.t[1:argmax(dlam)], y=first_dlam)
lam2 = scatter(x=new_df.t[argmax(dlam):end], y=last_dlam)
plot([lam1,lam2], Layout(xaxis_range=(0,1300)))


# interpolate to get new atp values 
atp_rtc_first = int_first(first_dlam)
atp_rtc_last = int_last(last_dlam)

atp1 = scatter(x=new_df.t[1:argmax(dlam)], y=atp_rtc_first)
atp2 = scatter(x=new_df.t[argmax(dlam)+1:end], y=atp_rtc_last)
plot([atp1,atp2], Layout(xaxis_range=(0,1000)))

# rejoin sections for full atp array 
new_atp = [atp_rtc_first; atp_rtc_last]
plot(scatter(x=new_df.t, y=new_atp), Layout(xaxis_range=(0,10000)))

# CSV.write("colD_gr_data.csv", new_df)

p = make_subplots(rows=2, cols=2, vertical_spacing=0.08, subplot_titles=["Model λ" "Model ATP"; "Data λ" "Inferred ATP"])
add_trace!(p, (scatter(x=solu_ss.t, y=lam)), row=1, col=1)
add_trace!(p, (scatter(x=solu_ss.t, y=a)), row=2, col=1)
add_trace!(p, (scatter(x=new_df.t, y=new_df.gr)), row=1, col=2)
add_trace!(p, (scatter(x=new_df.t, y=new_atp)), row=2, col=2)
relayout!(p, showlegend=false, xaxis_range=(0,1300), xaxis2_range=(0,1300), xaxis3_range=(0,1300), xaxis4_range=(0,1300))#xaxis_type="log", xaxis2_type="log", xaxis3_type="log", xaxis4_type="log")
p


atp_lam_data = DataFrame(t = new_df.t, atp = new_atp[1:end-1], gr = new_df.gr)

CSV.write("/home/holliehindley/phd/data/atp_for_rtcmodel.csv", atp_lam_data)

