@everywhere begin
    using Arrow, DataFrames, TableOperations

    function loadandsort_arrow_file(file)
        df = Arrow.Table(file) |> DataFrame

        df_p = filter(row -> length(row[:event]) > 3, df)
        df_r = filter(row -> length(row[:event]) < 3 && row[:event][1] != 0, df)

        react_names = [:tscr_ab, :tscr_r, :tlr_a, :tlr_b, :tlr_r, :Vinflux, :Vdam, :Vtag, :Vrep, :deg_rd, :deg_rma, :deg_rmb, :deg_rmr, :V]

        df_prop = split.(replace.(df_p.event, r"[\[\]\(Any)]" => ""), ",")
        df_prop = permutedims(mapcols(x -> parse.(Float64, x), DataFrame(df_prop, :auto)))
        rename!(df_prop, react_names)
        
        df_results = select(df_p, Not(:event))
        df_r.event = [parse.(Float64, subarray) for subarray in df_r.event]
        
        df_react = combine(groupby(df_r, :event), nrow => :count)[2:12,:]
        insertcols!(df_react, :reaction => react_names[Int64.(df_react.event)])

        return df_prop, df_react, df_results
    end
end

function loadsort_all_arrow_files(folder_path)
    files = readdir(folder_path)[1:5]
    n = length(files)
    df_props = Vector{DataFrame}(undef, n)
    df_reacts = Vector{DataFrame}(undef, n)
    df_res = Vector{DataFrame}(undef, n)

    @sync @distributed for (i, file) in collect(enumerate(files))
        file_name = joinpath(folder_path, file)
        df_props[i], df_reacts[i], df_res[i] = loadandsort_arrow_file(file_name)
    end

    return df_props, df_reacts, df_res
end 

using Distributed
nprocs()
addprocs(10)

@everywhere include("/home/hollie_hindley/Documents/phd/stochastic_hybrid_code/setup/file_funcs.jl")

@elapsed df_props, df_reacts, df_res = loadsort_all_arrow_files("/home/hollie_hindley/Documents/stochastic_hybrid/thresh_test_arrow")

df_props