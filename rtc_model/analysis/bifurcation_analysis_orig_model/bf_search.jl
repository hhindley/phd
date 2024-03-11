using Plots
using Revise, ForwardDiff, Parameters, Setfield, LinearAlgebra, DataFrames

include("$PATHmay23_rtc/analysis/bifurcation_analysis/bf_funcs.jl")

wab_range = 10 .^range(-5,stop=0,length=10)
wr_range = 10 .^(range(-7,stop=0,length=10))

atp_range = range(500,stop=5000,length=10)
kin_range = range(0,stop=0.2,length=10)
lam_range = range(0.001,stop=0.04,length=10)

wr_range = (range(0.0001,stop=0.001,length=10)) # use this range with the above 3 ranges and get ~3000 combos that give bistability
wab_range = range(0.01, stop=1., length=10)

function run_param_search(atp_range, kin_range, lam_range, wr_range, wab_range, params)
    bistab1_df=DataFrame(wab = Float64[], wr = Float64[], kin = Float64[], lam = Float64[], bs_type = Symbol[])
    bistab1 = []
    params = deepcopy(params1)
    for i in wab_range
        params = merge(params, (ω_ab=i,))
        for j in wr_range
            params = merge(params, (ω_r=j,))
            # for k in atp_range 
            #     params = merge(params, (atp=k,))
                for l in kin_range
                    params = merge(params, (kin=l,))

                    for m in lam_range
                        params = merge(params, (lam=m,))
                        # println("lam = $m")
                        # println("kin = $l, lam = $m")
                        # println("atp = $k, kin = $l, lam = $m")
                        # println("wr = $j, atp = $k, kin = $l, lam = $m")
                        # println("wab = $i, wr = $j, atp = $k, kin = $l, lam = $m")

                        br = get_br(params, initial)
                        for b in range(1,length(br.specialpoint))
                            if br.specialpoint[b].type == :endpoint
                                Nothing
                            else
                                # push!(bistab, ("lam = $m", br.specialpoint[b].type))
                                # push!(bistab, ("kin = $l, lam = $m", br.specialpoint[b].type))
                                # push!(bistab, ("atp = $k, kin = $l, lam = $m", br.specialpoint[b].type))
                                # push!(bistab, ("wr = $j, atp = $k, kin = $l, lam = $m", br.specialpoint[b].type))
                                # push!(bistab1, ("wab = $i, wr = $j, kin = $l, lam = $m", br.specialpoint[b].type))
                                # push!(bistab1, ("wab = $i, wr = $j, atp = $k, kin = $l, lam = $m", br.specialpoint[b].type))

                                push!(bistab1_df, (i, j, l, m, br.specialpoint[b].type))
                            end
                        end                     
                    end
                end
            # end
        end
    end
    return bistab1_df
end
bistab1_df = @time run_param_search(atp_range, kin_range, lam_range, wr_range, wab_range, params)


bistab1_df = bistab1_df[1:2:end,:]

CSV.write("$PATHmay23_rtc/analysis/bifurcation_analysis/param_sweep_bs.csv", bistab1_df)
