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

function RenyiNegativity(filename::String, filedir::String=pwd();
    printLA=nothing, startbin::Int=3, endbin::Union{Int,Nothing}=nothing, dropmaxmin::Int=1)
    # Add .bin extension if not present
    if !endswith(filename, ".bin")
        filename = filename * ".bin"
    end

    # 打开文件
    filepath = joinpath(filedir, filename)

    # 使用 readdlm 读取整个文件
    data = readdlm(filepath, Float64)
    # Apply start and end bin selection
    data = data[startbin:end,:]
    if !isnothing(endbin)
        data = data[1:endbin,:]
    end
    # remove the max and min
    data = filter(data, dropmaxmin)

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
                   startbin=2, endbin=nothing, dropmaxmin=0,
                   columns=[5,7,9,11], labels=["E_kin", "E_pot", "E_tot", "E_tot^2"],
                   verbose=true)

Analyze energy.bin file and calculate Monte Carlo averages and errors for specified columns.

Arguments:
- `filename`: Name of the energy file (default: "energy.bin")
- `filedir`: Directory containing the file (default: current directory)
- `startbin`: First bin to include in analysis (default: 2)
- `endbin`: Last bin to include in analysis (default: all bins)
- `dropmaxmin`: Number of maximum and minimum values to drop (default: 0)
- `columns`: Array of column indices to analyze (default: [5,7,9,11])
- `labels`: Array of labels for the columns (default: ["E_kin", "E_pot", "E_tot", "E_tot^2"])
- `verbose`: Whether to print results to console (default: true)

Returns:
- Dictionary with column labels as keys and [mean, error] arrays as values

Example:
```julia
# Analyze energy.bin in current directory
results = EnergyAnalysis()

# Analyze specific columns with custom labels
results = EnergyAnalysis(columns=[5,7], labels=["Kinetic", "Potential"])

# Access results
kinetic_energy = results["Kinetic"][1]  # mean value
kinetic_error = results["Kinetic"][2]   # error
```
"""
function EnergyAnalysis(filename::String="energy.bin", filedir::String=pwd();
                       startbin::Int=2, endbin::Union{Int,Nothing}=nothing, dropmaxmin::Int=0,
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

    # Check if file has enough columns
    if size(data, 2) < maximum(columns)
        error("File has only $(size(data, 2)) columns, but analysis requested column $(maximum(columns))")
    end

    # Apply start and end bin selection
    data = data[startbin:end,:]
    if !isnothing(endbin)
        data = data[1:endbin,:]
    end

    # Remove the max and min values
    data = filter(data, dropmaxmin)

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