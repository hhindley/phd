using DataFrames, CSV, PyCall

# data
dfc = DataFrame(CSV.File("/home/holliehindley/phd/data/fig4c_bo.csv"))
WT1 = dfc[!,2]; hpx = dfc[!,3]; WT1_std = dfc[!,4]; hpx_std = dfc[!,5]; # promoter expression

dfe = DataFrame(CSV.File("/home/holliehindley/phd/data/fig4e_rtcoff_bo.csv"))
WT2 = dfe[!,2]; hpx_rtcoff = dfe[!,3]; WT2_std = dfe[!,4]; hpx_rtcoff_std = dfe[!,5]; # growth

dff = DataFrame(CSV.File("/home/holliehindley/phd/data/fig4f_rtcon_bo.csv"))
WT3 = dff[!,2]; hpx_rtcon = dff[!,3]; WT3_std = dff[!,4]; hpx_rtcon_std = dff[!,5]; # growth

df2 = DataFrame(CSV.File("/home/holliehindley/phd/data/colD_supf2_bo.csv"))
WT4 = df2[!,2]; WT_colD = df2[!,3]; WT4_std = df2[!,18]; WT_colD_std = df2[!,19]; # growth
nA = df2[!,4]; nA_colD = df2[!,5];  nA_std = df2[!,20]; nA_colD_std = df2[!,21];
nB = df2[!,8]; nB_colD = df2[!,9]; nB_B = df2[!,10]; nB_B_colD = df2[!,11]; nB_Bmut = df2[!,12]; nB_Bmut_colD = df2[!,13];  
nB_std = df2[!,24]; nB_colD_std = df2[!,25]; nB_B_std = df2[!,26]; nB_B_colD_std = df2[!,27]; nB_Bmut_std = df2[!,28]; nB_Bmut_colD_std = df2[!,29];
nR = df2[!,14]; nR_colD = df2[!,15]; nR_std = df2[!,30]; nR_colD_std = df2[!,31]; 



# set time span and how many time points to solve at 
tspan2 = (0, 2880)
t_2 = dfc[!,1]*60

tspan4 = (0, 1297)
t_4 = df2[!,1]*60


function compare_data_and_sol(solu, curve, data, std)#(model, params, init, tspan, t, curve, data, std)
    obj = [] #, obj1 = [], []
    # params_dict, init = choose_param_vector(model)
    # params = set_N(data)
    # print(params)
    # solu = sol_with_t(model, init, params, tspan, t)
    if curve == "mrna"
        rm_a = get_curve(solu, :rm_a)
        # rm_r = get_curve(solu, :rm_r)
        append!(obj, [abs2((i-j)/k) for (i,j,k) in zip(rm_a, data, std)])
        # append!(obj1, [abs2((i-j)/k) for (i,j,k) in zip(rm_r, data, std)])
    return sum(obj)#, sum(obj1)
    elseif curve == "OD"
        OD = get_curve(solu, :OD)
        append!(obj, [abs2((i-j)/k) for (i,j,k) in zip(OD, data, std)])
    return sum(obj)
    end
end
 
function obj(model, init, params, tspan, t, curve, data, std)
    # @show ω_ab
    solu = sol_with_t(model, init, params, tspan, t)
    obj_wt = compare_data_and_sol(solu, curve, data, std)
    return obj_wt
end

function obj_OD(model, init, params, tspan, t, curve, data, std)
    solu = sol_with_t(model, init, params, tspan, t)
    obj_wt = compare_data_and_sol(solu, curve, data, std)
    return obj_wt
end

function set_k(data)
    if data == WT2
        k = findmax(WT2)[1]
    elseif data == WT3 
        k = findmax(WT3)[1]
    else 
        k = findmax(WT4)[1]
    end
    return k
end

function set_OD0(data)
    if data == WT2
        OD_0 = OD_0_wt2
    elseif data == WT3
        OD_0 = OD_0_wt3
    else 
        OD_0 = OD_0_wt4
    end
    return OD_0
end


# analysis 
function results_one_param(optimizer)
    vals, errors = [], []
    for i in collect(1:length(optimizer.space.target))
        a = collect(values(optimizer.res)[i])
        append!(vals, collect(values(a[2][2])))
        append!(errors, -(collect(values(a[1][2]))))
    end

    # reversing sum of squares error calculation to get errors close to zero 
    errors_ori = sqrt.(errors/15)

    # best values 
    best_param = collect(values(collect(values(optimizer.max))[2]))
    best_error = -[collect(values(optimizer.max))[1]]

    # best error with sum of squares reversed 
    best_error_ori = sqrt.(best_error/15)
    return vals, errors, errors_ori, best_param, best_error, best_error_ori
end

function results_two_param(optimizer)
    val_ab, val_r, errors = Array{Float64}(undef, 0), Array{Float64}(undef, 0), Array{Float64}(undef, 0)
    for i in collect(1:length(optimizer.space.target))
        a = collect(values(optimizer.res)[i])
        append!(val_ab, collect(values(a[2][2]))[2])
        append!(val_r, collect(values(a[2][2]))[1])
        append!(errors, -(collect(values(a[1][2]))))
    end

    # best values 
    best_param_ab = collect(values(collect(values(optimizer.max))[2]))[2]
    best_param_r = collect(values(collect(values(optimizer.max))[2]))[1]
    best_error = -[collect(values(optimizer.max))[1]]

    return val_ab, val_r, errors, best_param_ab, best_param_r, best_error
end

function plot_3D(vals_ab, vals_r, errors)
    points = hcat(vals_ab, vals_r)'
    nx = 101
    ny = 200
    x = LinRange((py"param_range_ω"["ω_ab"])[1],(py"param_range_ω"["ω_ab"])[2], nx)
    y = LinRange((py"param_range_ω"["ω_r"])[1],(py"param_range_ω"["ω_r"])[2], ny)
    xy = Iterators.product(x,y)
    xx = getindex.(xy,1)
    yy = getindex.(xy,2)
    points2 = hcat(vec(xx), vec(yy))'

    itp = interpolate(Multiquadratic(), points, errors)
    interpolated = ScatteredInterpolation.evaluate(itp, points2) # evaluate it for how ever many coordinates that you want
    z = reshape(interpolated, nx,ny) # make the interpolated into a grid for plotting 
    return x, y, z
end
