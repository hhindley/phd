using CSV, Arrow, FilePathsBase, Distributed, JSON, DataFrames

function arrow_conv(folder_path, arrow_folder_path)
    folder_path = folder_path
    files = readdir(folder_path)
    arrow_folder_path = arrow_folder_path
    if !isdir(arrow_folder_path)
        mkdir(arrow_folder_path)
    end

    for file in files
        df = DataFrame(CSV.File(joinpath(folder_path, file), header=["event", "time", "rm_a", "rtca", "rm_b", "rtcb", "rm_r", "rtcr", "rh", "rd", "rt", "volume", "totprop"]))[:,1:end-1]
        # print("$file df complete \n")
        df.event = Array{Float64}.(JSON.parse.(df.event))
        # print("$file df event converted \n")
        Arrow.write(joinpath(arrow_folder_path, splitext(basename(file))[1] * ".arrow"), df)
        # print("$file arrow file created \n")
        # rm(joinpath(folder_path, file), force=true)
        # print("$file removed \n")
    end
    rm(folder_path, recursive=true, force=true)
end


function loadandsort_arrow_file(file)
    atab = Arrow.Table(file)
    df_p = atab |> TableOperations.filter(x -> length(x.event) > 3) |> DataFrame
    df_r = atab |> TableOperations.filter(x -> length(x.event) < 3 && x.event[1] !=0) |> DataFrame
    df_prop = DataFrame(transpose(hcat(df_p.event...)), [:tscr_ab, :tscr_r, :tlr_a, :tlr_b, :tlr_r, :Vinflux, :Vdam, :Vtag, :Vrep, :deg_rd, :deg_rma, :deg_rmb, :deg_rmr, :V])
    df_react = combine(groupby(df_r, :event), nrow => :count)[2:11,:]
    df_react.event = reduce(vcat, map(v -> round.(Int64, v), collect.(df_react.event)))
    insertcols!(df_react, :reaction => react_names[df_react.event])

    return df_prop, df_react, df_p

end

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

function loadsort_all_arrow_files(folder_path)
    # folder_path = folder_path
    for i in ["results", "props", "reacts"]
        if !isdir(joinpath(folder_path, i))
            print(joinpath(folder_path, i))
            mkdir(joinpath(folder_path, i))
        end
    end
    files = filter(isfile, readdir(folder_path, join=true))#[1:2]
    print(files)
    for file in files
        # file_name = joinpath(folder_path, file)
        # print(file)
        df_prop, df_react, df_res = loadandsort_arrow_file(file)
        Arrow.write(joinpath(folder_path, "results/$(basename(file))"), df_res::DataFrame)
        Arrow.write(joinpath(folder_path, "props/$(basename(file))"), df_prop::DataFrame)
        Arrow.write(joinpath(folder_path, "reacts/$(basename(file))"), df_react::DataFrame)
    
    end

end 