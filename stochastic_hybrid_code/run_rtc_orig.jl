using ModelingToolkit, DifferentialEquations, PlotlyJS, LinearAlgebra, DataFrames, LabelledArrays, Printf, BifurcationKit, OrderedCollections

include("/home/hollie_hindley/Documents/paper/model_params_funcs_2024/solving.jl")
include("/home/hollie_hindley/Documents/paper/model_params_funcs_2024/rtc_orig.jl")
include("/home/hollie_hindley/Documents/paper/model_params_funcs_2024/params.jl")
include("/home/hollie_hindley/Documents/paper/model_params_funcs_2024/rtc_params_molecs.jl")

# include("/home/hollie_hindley/Documents/paper/model_params_funcs_2024/rtc_with_vol.jl")

solu_rtc_molec = sol(rtc_model, init_rtc_molec, (0,200*log(2)/lam_val), params_rtc_molec)

df_rtc = create_solu_df(solu_rtc_molec, species_rtc)

params_kdam = deepcopy(params_rtc_molec)
params_kdam[kdam] = 0.6
solu_rtc_dam = sol(rtc_model, ssvals_rtc_molec, (0, 1000*log(2)/lam_val), params_kdam)

df_rtc = create_solu_df(solu_rtc_dam, species_rtc)

println(ssvals_rtc_molec)

function all_vars(df)
    alpha = @. df.rt/par[pidx(:kr)] #/ v_sf(X)# unitless
    fa = @. (1+alpha)^6/(par[pidx(:L)]*((1+par[pidx(:c)]*alpha)^6)+(1+alpha)^6) # unitless 
    ra = @. fa*df.rtcr # uM 

    # Voc = par[pidx(:Vmax_init)]*atp/(par[pidx(:Km_init)]/v_sf(X)+atp(X)) # uM min-1 
    # sig_o = ra(X)*Voc(X)/par[pidx(:k_diss)] # uM

    # tscr_ab(X) = sig_o(X)*par[pidx(:ω_ab)]*atp(X)/(par[pidx(:θtscr)]/v_sf(X)+atp(X)) # uM min-1
    # tscr_r(X) = par[pidx(:ω_r)]*atp(X)/(par[pidx(:θtscr)]/v_sf(X)+atp(X))
    return alpha, fa, ra
end

alpha, fa, ra = all_vars(df_rtc)

Voc = par[pidx(:Vmax_init)]*atp_val_molec/(par[pidx(:Km_init)]+atp_val_molec) # uM min-1 
sig_o = @. ra*Voc/par[pidx(:k_diss)] 

tscr_ab = @. sig_o*par[pidx(:ω_ab)]*atp_val_molec/(par[pidx(:θtscr)]+atp_val_molec) # uM min-1
tscr_r = par[pidx(:ω_r)]*atp_val_molec/(par[pidx(:θtscr)]+atp_val_molec) # uM min-1

# colours =["#EF553B",  "#00CC96", "#AB63FA", "#FFA15A", "#19D3F3", "#FF6692", "#B6E880", "#FF97FF", "#FECB52", :blue]

# p_rtc1_molec = plot([scatter(x=df_molec.time, y=col, name="$(names(df_molec)[i])", legendgroup="$i", marker_color=colours[i]) for (col, i) in zip(eachcol(df_molec[:,2:end]), range(2,length(names(df_molec))))], Layout(title="kdam = $(params_rtc_molec[kdam])"))


# solu_rtc_molec = sol(rtc_model_vol, init_rtc_molec_vol, (0,2000), params_rtc_molec)

# df_rtc = create_solu_df(solu_rtc_molec, species_rtc2)

# plot(scatter(x=df_rtc.time, y=df_rtc.v))


function get_X0(indV, ssvals)
    X0 = zeros(indV.nrOfItems)

    X0[vidx(:rm_a)] = ssvals[1]
    X0[vidx(:rtca)] = ssvals[2]
    X0[vidx(:rm_b)] = ssvals[3]
    X0[vidx(:rtcb)] = ssvals[4]
    X0[vidx(:rm_r)] = ssvals[5]
    X0[vidx(:rtcr)] = ssvals[6]
    X0[vidx(:rh)] = ssvals[7]
    X0[vidx(:rd)] = ssvals[8]
    X0[vidx(:rt)] = ssvals[9]
    X0[vidx(:V)] = 1 # added in volume/dilution corrections
    X0[vidx(:totProp)] = 0

    return X0
end

X0 = collect(get_X0(indV, ssvals_rtc_molec)')



