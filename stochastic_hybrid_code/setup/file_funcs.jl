using CSV, Arrow, FilePathsBase, Distributed, JSON, DataFrames

# function arrow_conv(folder_path, arrow_folder_path)
#     folder_path = folder_path
#     files = readdir(folder_path)
#     arrow_folder_path = arrow_folder_path
#     if !isdir(arrow_folder_path)
#         mkdir(arrow_folder_path)
#     end

#     for file in files
#         # print(joinpath(folder_path, file))
#         open(Arrow.Writer, joinpath(arrow_folder_path, splitext(basename(file))[1] * ".arrow")) do writer
#             for chunk in CSV.Chunks(joinpath(folder_path, file); header=["event", "time", "rm_a", "rtca", "rm_b", "rtcb", "rm_r", "rtcr", "rh", "rd", "rt", "volume", "totprop", "missing"], drop=["missing"])
#                 df_chunk = DataFrame(chunk) 
#                 df_chunk.event = Array{Float64}.(JSON.parse.(df_chunk.event))
#                 Arrow.write(writer, df_chunk)
#             end
#         end
#     end

#     # for file in files
#     #     df = DataFrame(CSV.File(joinpath(folder_path, file), header=["event", "time", "rm_a", "rtca", "rm_b", "rtcb", "rm_r", "rtcr", "rh", "rd", "rt", "volume", "totprop"]))[:,1:end-1]
#     #     # print("$file df complete \n")
#     #     df.event = Array{Float64}.(JSON.parse.(df.event))
#     #     # print("$file df event converted \n")
#     #     Arrow.write(joinpath(arrow_folder_path, splitext(basename(file))[1] * ".arrow"), df)
#     #     # print("$file arrow file created \n")
#     #     # rm(joinpath(folder_path, file), force=true)
#     #     # print("$file removed \n")
#     # end
#     # rm(folder_path, recursive=true, force=true)
# end

function arrow_conv(folder_path, arrow_folder_path)
    folder_path = folder_path
    files = readdir(folder_path)
    arrow_folder_path = arrow_folder_path
    if !isdir(arrow_folder_path)
        mkdir(arrow_folder_path)
    end
    for i in ["results", "props", "reacts"]
        if !isdir(joinpath(arrow_folder_path, i))
            # print(joinpath(arrow_folder_path, i))
            mkdir(joinpath(arrow_folder_path, i))
        end
    end
    react_names=[:tscr_ab, :tscr_r, :tlr_a, :tlr_b, :tlr_r, :Vinflux, :Vdam, :Vtag, :Vrep, :deg_rd, :deg_rma, :deg_rmb, :deg_rmr, :V]

    for file in files
        # print(joinpath(folder_path, file))
        
        results_writer = open(Arrow.Writer, joinpath(joinpath(arrow_folder_path, "results"), splitext(basename(file))[1] * ".arrow"))
        props_writer = open(Arrow.Writer, joinpath(joinpath(arrow_folder_path, "props"), splitext(basename(file))[1] * ".arrow"))
        # reacts_writer = open(Arrow.Writer, joinpath(joinpath(arrow_folder_path, "reacts"), splitext(basename(file))[1] * ".arrow"))
        println("opened writes for $(basename(file))")
    
        df_reacts = DataFrame[]
        for chunk in CSV.Chunks(joinpath(folder_path, file); header=["event", "time", "rm_a", "rtca", "rm_b", "rtcb", "rm_r", "rtcr", "rh", "rd", "rt", "volume", "totprop", "missing"], drop=["missing"])
            df_chunk = DataFrame(chunk) 
            df_chunk.event = Array{Float64}.(JSON.parse.(df_chunk.event))

            dfp_chunk = filter(row -> length(row.event) > 3, df_chunk)
            Arrow.write(results_writer, dfp_chunk)
            # println("written chunks to results file")

            df_prop_chunk = DataFrame(transpose(hcat(dfp_chunk.event...)), react_names)
            Arrow.write(props_writer, df_prop_chunk)
            # println("written chunks to props file")

            dfr_chunk = filter(row -> length(row.event) < 3 && row.event[1] != 0, df_chunk)
            df_react_chunk = combine(groupby(dfr_chunk, :event), nrow => :count)
            if !isempty(df_react_chunk)
                df_react_chunk.event = reduce(vcat, map(v -> round.(Int64, v), collect.(df_react_chunk.event)))
                insertcols!(df_react_chunk, :reaction => react_names[df_react_chunk.event])
                push!(df_reacts, df_react_chunk)
            end
            # println("finished chunk")
        end
        df_react = combine(groupby(vcat(df_reacts...), [:reaction, :event]), :count => sum => :count)
        Arrow.write(joinpath(joinpath(arrow_folder_path, "reacts"), splitext(basename(file))[1] * ".arrow"), df_react)
        # println("written to reacts file")

        close(results_writer)
        close(props_writer)
        println("finished $(basename(file))")
        # close(reacts_writer)
    end

    rm(folder_path, recursive=true, force=true)
end

function load_files(folder_path; dataframe=true)
    files = readdir(folder_path)
    df_list = []
    if dataframe
        for file in files
            df = Arrow.Table(joinpath(folder_path, file)) |> DataFrame
            push!(df_list, df)
        end
    else
        for file in files
            df = Arrow.Table(joinpath(folder_path, file))
            push!(df_list, df)
        end
    end
    return df_list
end

function create_df_reacts(df_reacts)
    df_reacts = [DataFrame(df_reacts[i]) for i in eachindex(df_reacts)]
    for i in eachindex(df_reacts)
        for num in 1:13
            if num ∉ df_reacts[i].event
                push!(df_reacts[i], (event=num, count=0, reaction=react_names[num]))
            end
        end
    end
    return df_reacts
end

function load_hist_files(mainfolder)
    foldernames = readdir(mainfolder)
    dfs = Dict{String, Vector{DataFrame}}()
    for folder in foldernames
        folderpath = joinpath(mainfolder, folder)
        filenames = readdir(folderpath)
        for file in filenames
            prefix = split(file, ".")[1]
            filepath = joinpath(folderpath, file)
            df = DataFrame(CSV.File(filepath))
            df.bin = Array{Float64}.(JSON.parse.(df.bin))
            if haskey(dfs, prefix)
                push!(dfs[prefix], df)
            else
                dfs[prefix] = [df]
            end
        end
    end
    return dfs
end