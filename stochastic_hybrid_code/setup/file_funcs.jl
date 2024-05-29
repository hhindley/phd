using CSV, Arrow, FilePathsBase

function arrow_conv(folder_path, arrow_folder_path)
    folder_path = folder_path
    files = readdir(folder_path)
    arrow_folder_path = arrow_folder_path
    if !isdir(arrow_folder_path)
        mkdir(arrow_folder_path)
    end
    for file in files
        df = DataFrame(CSV.File(joinpath(folder_path, file), header=["event", "time", "rm_a", "rtca", "rm_b", "rtcb", "rm_r", "rtcr", "rh", "rd", "rt", "volume", "totprop"]))
        Arrow.write(joinpath(arrow_folder_path, splitext(basename(file))[1] * ".arrow"), df)
    end
    rm(folder_path, recursive=true, force=true)
end
# this is how you load in the arrow files as dataframes - should be really quick now 
function load_arrow_files(folder_path)
    tabs = [Arrow.Table(joinpath(folder_path, file)) for file in readdir(folder_path)]
    dfs = [DataFrame(tab) for tab in tabs]
    return dfs
end