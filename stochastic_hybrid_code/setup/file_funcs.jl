using CSV, Arrow, FilePathsBase, Distributed, JSON, DataFrames

function arrow_conv(folder_path, arrow_folder_path)
    folder_path = folder_path
    files = readdir(folder_path)
    arrow_folder_path = arrow_folder_path
    if !isdir(arrow_folder_path)
        mkdir(arrow_folder_path)
    end
    for i in ["results", "props", "reacts"]
        if !isdir(joinpath(arrow_folder_path, i))
            mkdir(joinpath(arrow_folder_path, i))
        end
    end
    react_names=[:tscr_ab, :tscr_r, :tlr_a, :tlr_b, :tlr_r, :Vinflux, :Vdam, :Vtag, :Vrep, :deg_rd, :deg_rma, :deg_rmb, :deg_rmr, :V]

    for file in files
        results_file_path = joinpath(joinpath(arrow_folder_path, "results"), splitext(basename(file))[1] * ".arrow")
        reacts_file_path = joinpath(joinpath(arrow_folder_path, "reacts"), splitext(basename(file))[1] * ".arrow")
        props_file_path = joinpath(joinpath(arrow_folder_path, "props"), splitext(basename(file))[1] * ".arrow")

        if isfile(results_file_path) && isfile(reacts_file_path) && isfile(props_file_path)
            println("File $(basename(file)) already processed. Skipping.")
            continue
        end

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

        if !isempty(df_reacts)
            df_react = combine(groupby(vcat(df_reacts...), [:reaction, :event]), :count => sum => :count)
            Arrow.write(joinpath(joinpath(arrow_folder_path, "reacts"), splitext(basename(file))[1] * ".arrow"), df_react)
        else
            println("No stochastic reactions happened in $(basename(file))")
        end
        # println("written to reacts file")

        close(results_writer)
        close(props_writer)
        println("finished $(basename(file))")
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
    all_items = readdir(mainfolder)
    foldernames = [item for item in all_items if !occursin("DS", item)]
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

function prod_tot_count(df_reacts)
    tot_counts = Int64[]
    for i in eachindex(df_reacts)
        push!(tot_counts,sum(df_reacts[i].count))
    end
    return tot_counts
end


function LoadDataVars(folder; reacts=true, results=true, props=true)
    folder = folder
    filepath = joinpath(mount_path, folder)
    if isfile(joinpath(filepath, replace(folder, "final_files" => "") * "times.csv"))
        df_times = CSV.File(joinpath(filepath, replace(folder, "final_files" => "") * "times.csv")) |> DataFrame
    else
        df_times = []#CSV.File(timefilepath) |> DataFrame
    end

    if "kdam" in names(df_times) && "thresh" in names(df_times)
        threshold_vals = df_times.thresh
        kdam_vals = df_times.kdam
        titles = ["kdam: $(kdam_vals[i]), thresh: $(threshold_vals[i])" for i in eachindex(threshold_vals)]
        tested_vals = Dict(:threshold=>threshold_vals, :kdam=>kdam_vals)
    elseif "kdam" in names(df_times)
        kdam_vals = df_times.kdam
        titles = ["kdam: $(kdam_vals[i])" for i in eachindex(kdam_vals)]
        tested_vals = Dict(:kdam=>kdam_vals)
    elseif "thresh" in names(df_times)
        threshold_vals = df_times.thresh
        titles = ["threshold: $(round(threshold_vals[i], digits=2))" for i in eachindex(threshold_vals)]
        tested_vals = Dict(:threshold=>threshold_vals)
    end

    if reacts && results && props
        df_reacts = load_files(joinpath(filepath, "reacts"))
        df_reacts = create_df_reacts(df_reacts)
        df_results = load_files(joinpath(filepath, "results"))
        df_props = load_files(joinpath(filepath, "props"), dataframe=false)
    elseif reacts && results
        df_reacts = load_files(joinpath(filepath, "reacts"))
        df_reacts = create_df_reacts(df_reacts)
        df_results = load_files(joinpath(filepath, "results"))
        df_props = []
    elseif reacts && props
        df_reacts = load_files(joinpath(filepath, "reacts"))
        df_reacts = create_df_reacts(df_reacts)
        df_props = load_files(joinpath(filepath, "props"), dataframe=false)
        df_results = []
    elseif results && props
        df_results = load_files(joinpath(filepath, "results"))
        df_props = load_files(joinpath(filepath, "props"), dataframe=false)
        df_reacts = []
    elseif results
        df_results = load_files(joinpath(filepath, "results"))
        df_props = []
        df_reacts = []
    elseif props
        df_props = load_files(joinpath(filepath, "props"), dataframe=false)
        df_results = []
        df_reacts = []
    elseif reacts
        df_reacts = load_files(joinpath(filepath, "reacts"))
        df_reacts = create_df_reacts(df_reacts)
        df_results = []
        df_props = []
    end
    return df_times, tested_vals, titles, df_results, df_reacts, df_props
end

function setup_dicts(folders)
    dict_times = Dict(i => DataFrame() for i in 1:length(folders))
    dict_threshvals = Dict(i => Dict() for i in 1:length(folders))
    dict_titles = Dict(i => String[] for i in 1:length(folders))
    dict_results = Dict(i => [] for i in 1:length(folders))
    dict_reacts = Dict(i => [] for i in 1:length(folders))
    dict_props = Dict(i => [] for i in 1:length(folders))
    dict_counts = Dict(i => [] for i in 1:length(folders))
    dict_hists = Dict{Int, Any}()
    return dict_times, dict_threshvals, dict_titles, dict_results, dict_reacts, dict_props, dict_counts, dict_hists
end

