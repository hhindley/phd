
struct MyIndexDict
    indices::Dict{Symbol, Int}
    nrOfItems::Int
    allItems::UnitRange{Int}
end

# Function to create an IndexDict with given names
function create_index_dict(names)
    indices = Dict{Symbol, Int}()
    for (i, name) in enumerate(names)
        indices[Symbol(name)] = i
    end
    return MyIndexDict(indices, length(names), 1:length(names))
end

# Helper function for compact indexing
function index(dict::MyIndexDict, item::Symbol)
    sym_item = Symbol(item)
    haskey(dict.indices, sym_item) || throw(ArgumentError("Item $item not found in the index dictionary"))
    return dict.indices[sym_item]
end

species_rtc = [:rm_a, :rtca, :rm_b, :rtcb, :rm_r, :rtcr, :rh, :rd, :rt]
spec = deepcopy(species_rtc)
params_names = [:L, :c, :kr, :Vmax_init, :Km_init, :ω_ab, :ω_r, :θtscr, :g_max, :θtlr, :km_a, :km_b, :d, :krep, :kdam, :ktag, :kdeg, :kin, :atp, :na, :nb, :nr, :lam, :kc, :k_diss]
push!(spec, :V)
push!(spec, :totProp)
indV = create_index_dict(spec)
indP = create_index_dict(params_names)


# Compact indexing function
vidx(item) = index(indV, item)
pidx(item) = index(indP, item)

sf = 1e6/(6.022e23*1e-15) # multiply to go from molecs/cell to uM # uM molec-1

# function get_X0(indV)
#     X0 = zeros(indV.nrOfItems)

#     X0[vidx(:rm_a)] = 0
#     X0[vidx(:rtca)] = 0
#     X0[vidx(:rm_b)] = 0
#     X0[vidx(:rtcb)] = 0
#     X0[vidx(:rm_r)] = 0
#     X0[vidx(:rtcr)] = 0
#     X0[vidx(:rh)] = 11.29
#     X0[vidx(:rd)] = 0
#     X0[vidx(:rt)] = 0
#     X0[vidx(:totProp)] = 0

#     return X0
# end

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
    X0[vidx(:V)] = 1 #1e-15 # volume of ecoli cell so now volume should be correct in simulations #1 # added in volume/dilution corrections
    X0[vidx(:totProp)] = 0

    return X0
end



# function get_par(indP)
#     par = zeros(indP.nrOfItems)

#     par[pidx(:L)] = L_val
#     par[pidx(:c)] = c_val
#     par[pidx(:kr)] = kr_val
#     par[pidx(:Vmax_init)] = Vmax_init_val
#     par[pidx(:Km_init)] = Km_init_val
#     par[pidx(:ω_ab)] = ω_ab_val
#     par[pidx(:ω_r)] = ω_r_val
#     par[pidx(:θtscr)] = θtscr_val
#     par[pidx(:g_max)] = g_max_val
#     par[pidx(:θtlr)] = θtlr_val
#     par[pidx(:km_a)] = km_a_val
#     par[pidx(:km_b)] = km_b_val
#     par[pidx(:d)] = d_val
#     par[pidx(:krep)] = krep_val
#     par[pidx(:kdam)] = kdam_val
#     par[pidx(:ktag)] = ktag_val
#     par[pidx(:kdeg)] = kdeg_val
#     par[pidx(:kin)] = kin_val
#     par[pidx(:atp)] = atp_val
#     par[pidx(:na)] = nA_val
#     par[pidx(:nb)] = nB_val
#     par[pidx(:nr)] = nR_val
#     par[pidx(:lam)] = lam_val
#     par[pidx(:kc)] = kc_val
#     par[pidx(:k_diss)] = k_diss_val

#     return par

# end

function get_par(indP)
    par = zeros(indP.nrOfItems)

    par[pidx(:L)] = L_val
    par[pidx(:c)] = c_val
    par[pidx(:kr)] = kr_val_molec
    par[pidx(:Vmax_init)] = Vmax_init_val
    par[pidx(:Km_init)] = Km_init_val_molec
    par[pidx(:ω_ab)] = ω_ab_val
    par[pidx(:ω_r)] = ω_r_val_molec
    par[pidx(:θtscr)] = θtscr_val_molec
    par[pidx(:g_max)] = g_max_val
    par[pidx(:θtlr)] = θtlr_val_molec
    par[pidx(:km_a)] = km_a_val_molec
    par[pidx(:km_b)] = km_b_val_molec
    par[pidx(:d)] = d_val
    par[pidx(:krep)] = krep_val
    par[pidx(:kdam)] = kdam_val
    par[pidx(:ktag)] = ktag_val
    par[pidx(:kdeg)] = kdeg_val
    par[pidx(:kin)] = kin_val_molec
    par[pidx(:atp)] = atp_val
    par[pidx(:na)] = nA_val
    par[pidx(:nb)] = nB_val
    par[pidx(:nr)] = nR_val
    par[pidx(:lam)] = lam_val
    par[pidx(:kc)] = kc_val_molec
    par[pidx(:k_diss)] = k_diss_val
    # par[pidx(:vol_max)] = vol_max_val

    return par

end