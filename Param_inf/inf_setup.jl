using DataFrames, CSV

# data
dfc = DataFrame(CSV.File("/home/holliehindley/phd/data/fig4c_bo.csv"))
WT1 = dfc[!,2]; hpx = dfc[!,3]; WT1_std = dfc[!,4]; hpx_std = dfc[!,5]; # promoter expression

dfe = DataFrame(CSV.File("/home/holliehindley/phd/data/fig4e_rtcoff_bo.csv"))
WT2 = dfc[!,2]; hpx_rtcoff = dfc[!,3]; WT2_std = dfc[!,4]; hpx_rtcoff_std = dfc[!,5]; # growth

dff = DataFrame(CSV.File("/home/holliehindley/phd/data/fig4f_rtcon_bo.csv"))
WT3 = dfc[!,2]; hpx_rtcon = dfc[!,3]; WT3_std = dfc[!,4]; hpx_rtcon_std = dfc[!,5]; # growth

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


function compare_data_and_sol(model, init, tspan, dict, t, curve, data, std)
    obj = []
    param_vector = @SVector [values(dict)]
    solu = sol_with_t(model, init, tspan, param_vector[1], t)
    if curve == "mrnas"
        rm_a = get_curve(solu, :rm_a)
        rm_b = get_curve(solu, :rm_b)
        append!(obj, [abs2((i-j)/k) for (i,j,k) in zip(rm_a, data, std)])
        append!(obj, [abs2((i-j)/k) for (i,j,k) in zip(rm_b, data, std)])
    elseif curve == "den"
        den = get_curve(solu, :den)
        append!(obj, [abs2((i-j)/k) for (i,j,k) in zip(den, data, std)])
    end
    return obj
end