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

function load_files(folder_path)
    files = readdir(folder_path)
    df_list = []
    for file in files
        df = Arrow.Table(joinpath(folder_path, file)) #|> DataFrame
        push!(df_list, df)
    end
    return df_list
end
# function loadandsort_arrow_file(file)
#     atab = Arrow.Table(file)
#     df_p = atab |> TableOperations.filter(x -> length(x.event) > 3) |> DataFrame
#     df_r = atab |> TableOperations.filter(x -> length(x.event) < 3 && x.event[1] !=0) |> DataFrame
#     df_prop = DataFrame(transpose(hcat(df_p.event...)), [:tscr_ab, :tscr_r, :tlr_a, :tlr_b, :tlr_r, :Vinflux, :Vdam, :Vtag, :Vrep, :deg_rd, :deg_rma, :deg_rmb, :deg_rmr, :V])
#     df_react = combine(groupby(df_r, :event), nrow => :count)#[2:11,:]
#     df_react.event = reduce(vcat, map(v -> round.(Int64, v), collect.(df_react.event)))
#     insertcols!(df_react, :reaction => react_names[df_react.event])

#     return df_prop, df_react, df_p

# end

# function loadandsort_arrow_file(file)
#     df = Arrow.Table(file) |> DataFrame

#     # takes a view of the data and keeps those rows that are less than 3 long and not zero 
#     df_r = @view df[findall(row -> length(row[:event]) < 3 && row[:event][1] != 0, eachrow(df)), :]
    
#     # takes a view of the data and keeps those rows that are more than 3 long
#     df_p = @view df[findall(row -> length(row[:event]) > 3, eachrow(df)), :]
    
#     react_names = [:tscr_ab, :tscr_r, :tlr_a, :tlr_b, :tlr_r, :Vinflux, :Vdam, :Vtag, :Vrep, :deg_rd, :deg_rma, :deg_rmb, :deg_rmr, :V]
    
#     # converts the event column to an array 
#     props = str_to_arr.(df_p.event)
    
#     # creates a new dataframe of the propensities next to the reaction names 
#     df_prop = DataFrame(transpose(hcat(props...)), react_names)
        
#     # look into this? 
#     df_r.event = [parse.(Float64, subarray) for subarray in df_r.event]
#     # new dataframe for reactions 
#     df_react = combine(groupby(df_r, :event), nrow => :count)[2:12,:]
#     insertcols!(df_react, :reaction => react_names[Int64.(df_react.event)])
    
#     # df_results = select!(df_p, Not(:event))

# return df_prop, df_react, df_p
# end

# function loadsort_all_arrow_files(folder_path)
#     # folder_path = folder_path
#     for i in ["results", "props", "reacts"]
#         if !isdir(joinpath(folder_path, i))
#             print(joinpath(folder_path, i))
#             mkdir(joinpath(folder_path, i))
#         end
#     end
#     files = filter(isfile, readdir(folder_path, join=true))#[1:2]
#     # print(files)
#     for file in files
#         # need to add here - check if file in each folder exists, if it does then skip that file in the for loop and if it doesnt then do convert that file
#         if basename(file) ∉ readdir(joinpath(folder_path, "results")) && file ∉ readdir(joinpath(folder_path, "props")) && file ∉ readdir(joinpath(folder_path, "reacts"))
#             print("$(basename(file)) does not yet exist in subfolders, performing task \n")
#             df_prop, df_react, df_res = loadandsort_arrow_file(file)
#             Arrow.write(joinpath(folder_path, "results/$(basename(file))"), df_res::DataFrame)
#             Arrow.write(joinpath(folder_path, "props/$(basename(file))"), df_prop::DataFrame)
#             Arrow.write(joinpath(folder_path, "reacts/$(basename(file))"), df_react::DataFrame)
#             print("finished $(basename(file)) \n")

#         else
#             print("$(basename(file)) already exists in subfolders \n")
#         end
#     end

# end 