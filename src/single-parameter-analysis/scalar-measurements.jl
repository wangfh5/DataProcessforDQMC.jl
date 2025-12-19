#= 
单参数分析模块 (Single Parameter Analysis) 之标量量分析

此模块提供了一组函数，用于分析单个参数文件夹中的DQMC模拟数据。
通常在模拟输出目录中直接运行，处理该目录下的 .bin 文件。

主要功能:
- 能量分析 (EnergyAnalysis)
- Renyi负值计算 (RenyiNegativity, RenyiNegativity_all)

=#

export EnergyAnalysis
export RenyiNegativity, RenyiNegativity_all

## -------------------------------------------------------------------------- ##
##                              Renyi Negativity                              ##
## -------------------------------------------------------------------------- ##

"""
    RenyiNegativity(filename, filedir=pwd(); printLA=nothing, startbin=3, endbin=nothing,
                    outlier_mode=:dropmaxmin, outlier_param=1)

Calculate Renyi negativity from expRenyiN*.bin files.

# Arguments
- `filename`: Name of the Renyi negativity file (e.g., "expRenyiN3.bin")
- `filedir`: Directory containing the file (default: current directory)
- `printLA`: Which LA values to print (default: all)
- `startbin`: First bin to include (default: 3)
- `endbin`: Last bin to include, inclusive (default: all bins)
- `outlier_mode`: `:dropmaxmin` (trim extremes) or `:iqrfence` (Tukey IQR fence).
- `outlier_param`: For `:dropmaxmin`, number of max/min values to trim (integer, default 1).
                  For `:iqrfence`, positive Real k multiplier for IQR.

# Returns
- `(expRenyiN, RenyiN)`: Arrays of shape (L+1, 2) with [mean, error] for each LA.
"""
function RenyiNegativity(filename::String, filedir::String=pwd();
    printLA=nothing, startbin::Int=3, endbin::Union{Int,Nothing}=nothing,
    outlier_mode::Symbol=:dropmaxmin, outlier_param::Real=1)

    # Add .bin extension if not present
    if !endswith(filename, ".bin")
        filename = filename * ".bin"
    end

    # 打开文件
    filepath = joinpath(filedir, filename)

    # 使用 readdlm 读取整个文件
    data = readdlm(filepath, Float64)
    nbin = size(data, 1)
    start_idx = max(startbin, 1)
    end_idx = isnothing(endbin) ? nbin : min(endbin, nbin)
    data = data[start_idx:end_idx, :]

    # deduce L and rank from the data
    L = size(data,2) ÷ 2 - 1
    # example filename: expRenyiN3.bin, expRenyiN4_TW.bin
    rank = parse(Int, split(filename, "expRenyiN")[2][1])
    quantity_name = split(filename, ".bin")[1]

    # 输出结果
    expRenyiN = zeros(L+1,2)
    RenyiN = zeros(L+1,2)
    for i in 1:L+1
        vectmp = data[:,2i-1] # take real part
        vectmp = remove_outliers(vectmp, outlier_mode, outlier_param; min_n=5)
        if mean(vectmp) > 0.0
            mn = mean(vectmp)

            # Calculate error with auto_digits
            err = error(vectmp, sigma=1, bessel=true, auto_digits=true)

            expRenyiN[i,:] = [mn, err]
            RenyiN[i,:] = [log(mn)/(1-rank), abs(err/mn/(1-rank))]

            if (printLA == nothing) || ((i-1) ∈ printLA)
                # Format values for display
                formatted_mn, formatted_err = format_value_error(mn, err)
                renyi_val = RenyiN[i,1]
                renyi_err = RenyiN[i,2]
                formatted_renyi, formatted_renyi_err = format_value_error(renyi_val, renyi_err)

                println("$(quantity_name) $(i-1) $formatted_mn ± $formatted_err $formatted_renyi ± $formatted_renyi_err")
            end
        end
    end
    # @show expRenyiN
    # @show RenyiN
    return expRenyiN, RenyiN
end

function RenyiNegativity_all(filedir::String=pwd();maxrank::Int=4, kwargs...)
    # find all the files with name expRenyiN*.bin or expRenyiN*_TW.bin, where * is an integer
    filenames = Base.filter(x->occursin(r"expRenyiN\d+\.bin", x) || occursin(r"expRenyiN\d+_TW\.bin", x), readdir(filedir))
    # sort the filenames
    filenames = sort(filenames)
    # get the number of files
    nfiles = length(filenames)
    # for each file, create an array to store the average and error of Renyi negativity
    # and then append the result to an opened JLD2 file
    file = jldopen("RenyiNall.jld2", "w") do file
        for i in 1:nfiles
            filename = filenames[i]
            suffix = split(filename, "expRenyiN")[2][1:end-4]
            rank = parse(Int, suffix[1])
            if rank <= maxrank
                expRenyiN, RenyiN = RenyiNegativity(filename, filedir;kwargs...)
                # save the result to a JLD2 file
                file["expRenyiN$suffix"] = expRenyiN
                file["RenyiN$suffix"] = RenyiN
            end
        end
    end
end

## -------------------------------------------------------------------------- ##
##                                   Energy                                   ##
## -------------------------------------------------------------------------- ##

"""
    EnergyAnalysis(filename="energy.bin", filedir=pwd();
                   startbin=2, endbin=nothing, outlier_mode=:dropmaxmin, outlier_param=0,
                   columns=[5,7,9], labels=["E_kin", "d_occ", "E_tot"],
                   verbose=true)

Analyze energy.bin file and calculate Monte Carlo averages and errors for specified columns.

# Arguments
- `filename`: Name of the energy file (default: "energy.bin")
- `filedir`: Directory containing the file (default: current directory)
- `startbin`: First bin to include in analysis (default: 2)
- `endbin`: Last bin to include in analysis, inclusive (default: all bins)
- `outlier_mode`: `:dropmaxmin` (trim extremes) or `:iqrfence` (Tukey IQR fence).
- `outlier_param`: For `:dropmaxmin`, number of max/min values to trim (integer, default 0).
                  For `:iqrfence`, positive Real k multiplier for IQR.
- `columns`: Array of column indices to analyze (default: [5,7,9])
- `labels`: Array of labels for the columns (default: ["E_kin", "d_occ", "E_tot"])
- `verbose`: Whether to print results to console (default: true)

# Returns
- Dictionary with column labels as keys and [mean, error] arrays as values

# Example
```julia
# Analyze energy.bin in current directory
results = EnergyAnalysis()

# Analyze specific columns with custom labels
results = EnergyAnalysis(columns=[5,7], labels=["Kinetic", "Potential"])

# Use IQR fence instead of trim
results = EnergyAnalysis(outlier_mode=:iqrfence, outlier_param=10.0)

# Access results
kinetic_energy = results["Kinetic"][1]  # mean value
kinetic_error = results["Kinetic"][2]   # error
```
"""
function EnergyAnalysis(filename::String="energy.bin", filedir::String=pwd();
                       startbin::Int=2, endbin::Union{Int,Nothing}=nothing,
                       outlier_mode::Symbol=:dropmaxmin, outlier_param::Real=0,
                       columns::Vector{Int}=[5,7,9],
                       labels::Vector{String}=["E_kin", "d_occ", "E_tot"],
                       verbose::Bool=true)
    # Validate inputs
    if length(columns) != length(labels)
        error("Number of columns ($(length(columns))) must match number of labels ($(length(labels)))")
    end

    # Add .bin extension if not present
    if !endswith(filename, ".bin")
        filename = filename * ".bin"
    end

    # Open and read the file
    filepath = joinpath(filedir, filename)
    if !isfile(filepath)
        error("File not found: $filepath")
    end

    # Read data
    data = readdlm(filepath, Float64)
    nbin = size(data, 1)

    # Check if file has enough columns
    if size(data, 2) < maximum(columns)
        error("File has only $(size(data, 2)) columns, but analysis requested column $(maximum(columns))")
    end

    # Apply start and end bin selection (inclusive bin indices)
    start_idx = max(startbin, 1)
    end_idx = isnothing(endbin) ? nbin : min(endbin, nbin)
    data = data[start_idx:end_idx, :]

    # Calculate statistics for each column
    results = Dict{String, Vector{Float64}}()

    if verbose
        println("\nMonte Carlo analysis of $filename:")
        println("----------------------------------------")
        println("Using $(size(data, 1)) measurements after filtering")
        println("----------------------------------------")
    end

    for (i, col) in enumerate(columns)
        label = labels[i]
        values = data[:, col]
        # Apply outlier removal via unified interface
        values = remove_outliers(values, outlier_mode, outlier_param; min_n=5)

        # Calculate mean and error
        mean_val = mean(values)

        # Calculate error with appropriate precision
        err_val = error(values, sigma=1, bessel=true, auto_digits=true)

        # Store results
        results[label] = [mean_val, err_val]

        if verbose
            # Format value and error for display
            formatted_val, formatted_err = format_value_error(mean_val, err_val)
            println("$label = $formatted_val ± $formatted_err")
        end
    end

    if verbose
        println("----------------------------------------\n")
    end

    return results
end