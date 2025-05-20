#= 
单参数分析模块 (Single Parameter Analysis)

此模块提供了一组函数，用于分析单个参数文件夹中的DQMC模拟数据。
通常在模拟输出目录中直接运行，处理该目录下的 .bin 文件。

主要功能:
- 能量分析 (EnergyAnalysis)
- 单轨道关联函数分析 (CorrelationAnalysis)
- 多轨道关联函数分析 (MultiOrbitalCorrelationAnalysis)
- Renyi负值计算 (RenyiNegativity, RenyiNegativity_all)

使用示例:
1. 在Julia REPL中:
   using DataProcessforDQMC
   EnergyAnalysis("energy.bin", ".")

2. 作为脚本运行:
   julia -e 'using DataProcessforDQMC; EnergyAnalysis("energy.bin", ".")'
=#

using Printf
using Statistics
using DelimitedFiles
using DataFrames

export EnergyAnalysis, CorrelationAnalysis, MultiOrbitalCorrelationAnalysis
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

## -------------------------------------------------------------------------- ##
##                      Real-space Correlation functions                      ##
## -------------------------------------------------------------------------- ##

## ---------------- Helper functions for correlation analysis --------------- ##

"""Calculate Euclidean distance between two points"""
function euclidean_distance(x, y)
    return sqrt(x^2 + y^2)
end

"""Calculate i coordinate from imj when j is fixed at (1,1)"""
function calculate_i_coord(imj_x, imj_y, Lx, Ly)
    i_x = mod1(1 + imj_x, Lx)  # mod1 ensures result is in range 1:Lx
    i_y = mod1(1 + imj_y, Ly)  # mod1 ensures result is in range 1:Ly
    return (i_x, i_y)
end

"""
    filter_bins(df, startbin, endbin, dropmaxmin, n_bins, verbose)

Filter bins based on startbin, endbin, and dropmaxmin.

Arguments:
- `df`: DataFrame containing all bin data
- `startbin`: First bin to include in analysis
- `endbin`: Last bin to include in analysis (default: all bins)
- `dropmaxmin`: Number of maximum and minimum values to drop (default: 0)
- `n_bins`: Total number of bins
- `verbose`: Whether to print results to console (default: true)

Returns:
- Filtered DataFrame
"""
function filter_bins(df, startbin, endbin, dropmaxmin, n_bins, verbose; value_columns=[:real_val])
    # 选择bin范围
    if !isnothing(endbin)
        filtered_df = df[(df.bin .>= startbin) .& (df.bin .<= endbin), :]
    else
        filtered_df = df[df.bin .>= startbin, :]
    end
    
    selected_bins = unique(filtered_df.bin)
    selected_count = length(selected_bins)
    
    # 应用max/min过滤（如果需要）
    if dropmaxmin > 0 && selected_count > 2*dropmaxmin
        # 检查所有请求的列是否存在
        valid_columns = [col for col in value_columns if col in propertynames(filtered_df)]
        
        if isempty(valid_columns)
            @warn "None of the specified value_columns exist in the DataFrame. Using first numeric column for filtering."
            # 找到第一个数值列
            numeric_cols = [col for col in propertynames(filtered_df) if eltype(filtered_df[!, col]) <: Number && col != :bin]
            if !isempty(numeric_cols)
                valid_columns = [first(numeric_cols)]
            else
                @error "No numeric columns found for bin filtering"
                return filtered_df
            end
        end
        
        # 对每个bin计算所有列的平均值总和
        bin_avg_df = combine(groupby(filtered_df, :bin)) do group
            total = 0.0
            for col in valid_columns
                total += mean(group[!, col])
            end
            (bin=first(group.bin), bin_avg=total)
        end
        
        # 对bin排序
        sorted_bins = sort(bin_avg_df, :bin_avg).bin
        valid_bins = sorted_bins[(dropmaxmin+1):(end-dropmaxmin)]
        
        # 只保留有效bin的数据
        filtered_df = Base.filter(row -> row.bin in valid_bins, filtered_df)
    end
    
    if verbose
        filtered_count = length(unique(filtered_df.bin))
        
        # 详细的bin选择和过滤信息
        if n_bins > 1
            if startbin > 1 || !isnothing(endbin)
                bin_range = "$startbin-$(isnothing(endbin) ? n_bins : min(endbin, n_bins))"
                println("  Selected bins: $selected_count of $n_bins (range: $bin_range)")
            else
                println("  Selected bins: $selected_count of $n_bins")
            end
            
            if dropmaxmin > 0
                println("  Filtered bins: $filtered_count (removed $dropmaxmin max/min values)")
            else
                println("  No max/min filtering applied")
            end
        end
    end
    
    return filtered_df
end

"""
    calculate_statistics(df, auto_digits)

Calculate statistics for each coordinate in the DataFrame.

Arguments:
- `df`: DataFrame containing bin information
- `auto_digits`: Whether to use automatic precision for error calculation

Returns:
- DataFrame containing statistics for each coordinate
"""
function calculate_statistics(df, auto_digits)
    # 对每个坐标分组
    grouped_df = groupby(df, :coord)
    
    # 对每个组计算统计量
    results_df = combine(grouped_df) do group_df
        values = group_df.complex_val
        real_values = real.(values)
        imag_values = imag.(values)
        
        # 计算均值
        mean_real = mean(real_values)
        mean_imag = mean(imag_values)
        
        # 计算全精度误差
        err_real = error(real_values, sigma=1, bessel=true, auto_digits=false)
        err_imag = error(imag_values, sigma=1, bessel=true, auto_digits=false)
        
        # 计算格式化误差
        err_real_fmt = error(real_values, sigma=1, bessel=true, auto_digits=auto_digits)
        err_imag_fmt = error(imag_values, sigma=1, bessel=true, auto_digits=auto_digits)
        
        # 格式化显示值
        formatted_real, formatted_real_err = format_value_error(mean_real, err_real_fmt)
        formatted_imag, formatted_imag_err = format_value_error(mean_imag, err_imag_fmt)
        
        return (
            # 分开保存实部和虚部用于后续计算
            mean_real = mean_real,
            mean_imag = mean_imag,
            err_real = err_real,
            err_imag = err_imag,
            
            # 格式化值用于显示
            formatted_real = "$(formatted_real) ± $(formatted_real_err)",
            formatted_imag = "$(formatted_imag) ± $(formatted_imag_err)",
            
            # 坐标信息
            imj_x = first(group_df.imj_x),
            imj_y = first(group_df.imj_y)
        )
    end
    
    # 按imj坐标排序
    sort!(results_df, [:imj_x, :imj_y])
    
    return results_df
end


"""
    calculate_multi_orbital_statistics(filtered_df, orbital_labels, auto_digits)

Calculate statistics for each coordinate and each orbital in the DataFrame.

Arguments:
- `filtered_df`: DataFrame containing filtered bin data
- `orbital_labels`: Array of orbital labels
- `auto_digits`: Whether to use automatic precision for error calculation

Returns:
- DataFrame containing statistics for each coordinate and orbital
"""
function calculate_statistics_multi_orbital(filtered_df, orbital_labels, auto_digits)
    # Group by coordinates
    return combine(groupby(filtered_df, :coord)) do group
        result = NamedTuple()
        
        # Basic coordinate information
        result = merge(result, (
            imj_x = first(group.imj_x),
            imj_y = first(group.imj_y)
        ))
        
        # Process each orbital
        for orbital in orbital_labels
            real_col = Symbol("$(orbital)_real")
            imag_col = Symbol("$(orbital)_imag")
            
            real_values = group[!, real_col]
            imag_values = group[!, imag_col]
            
            # Calculate statistics
            mean_real = mean(real_values)
            mean_imag = mean(imag_values)
            
            # Calculate errors (full precision)
            err_real = error(real_values, sigma=1, bessel=true, auto_digits=false)
            err_imag = error(imag_values, sigma=1, bessel=true, auto_digits=false)
            
            # Calculate formatted errors for display
            err_real_fmt = error(real_values, sigma=1, bessel=true, auto_digits=auto_digits)
            err_imag_fmt = error(imag_values, sigma=1, bessel=true, auto_digits=auto_digits)
            
            # Format for display
            formatted_real, formatted_real_err = format_value_error(mean_real, err_real_fmt)
            formatted_imag, formatted_imag_err = format_value_error(mean_imag, err_imag_fmt)
            
            # Add to result
            result = merge(result, NamedTuple{
                (Symbol("$(orbital)_mean_real"),
                 Symbol("$(orbital)_mean_imag"),
                 Symbol("$(orbital)_err_real"),
                 Symbol("$(orbital)_err_imag"),
                 Symbol("$(orbital)_formatted_real"),
                 Symbol("$(orbital)_formatted_imag"))
            }((
                mean_real,
                mean_imag,
                err_real,
                err_imag,
                "$(formatted_real) ± $(formatted_real_err)",
                "$(formatted_imag) ± $(formatted_imag_err)"
            )))
        end
        
        return result
    end
end

"""
    print_correlation_results(results_df, filename, Lx, Ly, dropmaxmin, startbin, endbin, n_bins)

Print single-orbital correlation analysis results.

Arguments:
- `results_df`: DataFrame containing analysis results
- `filename`: Name of the correlation file
- `Lx`, `Ly`: Lattice dimensions
- `dropmaxmin`: Number of maximum and minimum values dropped
- `startbin`, `endbin`: Bin range used
- `n_bins`: Total number of bins
"""
function print_correlation_results(results_df, filename, Lx, Ly, dropmaxmin, startbin, endbin, n_bins)
    println("\nCorrelation analysis of $filename:")
    println("----------------------------------------")
    println("Lattice size: $(Lx)×$(Ly)")
    
    # 打印bin信息
    if n_bins > 1
        if startbin > 1 || !isnothing(endbin)
            bin_range = "$startbin-$(isnothing(endbin) ? n_bins : min(endbin, n_bins))"
            println("Bins: using range $bin_range from $n_bins total, dropped $dropmaxmin max/min values")
        else
            println("Bins: dropped $dropmaxmin max/min values from $n_bins total")
        end
    else
        println("Only one bin found - no filtering applied")
    end
    
    # 创建显示用DataFrame
    display_df = DataFrame(
        "imj" => [string("(", r.imj_x, ",", r.imj_y, ")") for r in eachrow(results_df)],
        "distance" => [round(r.distance, digits=3) for r in eachrow(results_df)],
        "i [to j=(1,1)]" => [string("(", r.i_coord[1], ",", r.i_coord[2], ")") for r in eachrow(results_df)],
        "real part" => results_df.formatted_real,
        "imag part" => results_df.formatted_imag
    )
    
    println(display_df)
    println("----------------------------------------")
end

"""
    print_correlation_results_multi_orbital(results_df, filename, Lx, Ly, orbital_labels, dropmaxmin, startbin, endbin, n_bins)

Print multi-orbital correlation analysis results.

Arguments:
- `results_df`: DataFrame containing analysis results
- `filename`: Name of the correlation file
- `Lx`, `Ly`: Lattice dimensions
- `orbital_labels`: Array of orbital labels
- `dropmaxmin`: Number of maximum and minimum values dropped
- `startbin`, `endbin`: Bin range used
- `n_bins`: Total number of bins
"""
function print_correlation_results_multi_orbital(results_df, filename, Lx, Ly, orbital_labels, dropmaxmin, startbin, endbin, n_bins)
    println("\nMulti-orbital correlation analysis of $filename:")
    println("----------------------------------------")
    println("Lattice size: $(Lx)×$(Ly)")
    println("Orbitals analyzed: $(join(orbital_labels, ", "))")
    
    # Print bin information
    if n_bins > 1
        if startbin > 1 || !isnothing(endbin)
            bin_range = "$startbin-$(isnothing(endbin) ? n_bins : min(endbin, n_bins))"
            println("Bins: using range $bin_range from $n_bins total, dropped $dropmaxmin max/min values")
        else
            println("Bins: dropped $dropmaxmin max/min values from $n_bins total")
        end
    else
        println("Only one bin found - no filtering applied")
    end
    
    # Print detailed results for each orbital
    for orbital in orbital_labels
        println("\n----------------------------------------")
        println("Orbital: $orbital")
        println("----------------------------------------")
        
        # Create display DataFrame
        display_df = DataFrame(
            "imj" => [string("(", r.imj_x, ",", r.imj_y, ")") for r in eachrow(results_df)],
            "distance" => [round(r.distance, digits=3) for r in eachrow(results_df)],
            "i [to j=(1,1)]" => [string("(", r.i_coord[1], ",", r.i_coord[2], ")") for r in eachrow(results_df)],
            "real part" => results_df[:, "$(orbital)_formatted_real"],
            "imag part" => results_df[:, "$(orbital)_formatted_imag"]
        )
        
        println(display_df)
    end
    
    # Print compact summary (real parts only)
    println("\n----------------------------------------")
    println("Compact summary (real parts only):")
    println("----------------------------------------")
    
    # Build compact table
    compact_df = DataFrame("imj" => [string("(", r.imj_x, ",", r.imj_y, ")") for r in eachrow(results_df)])
    
    for orbital in orbital_labels
        compact_df[!, orbital] = results_df[:, "$(orbital)_formatted_real"]
    end
    
    println(compact_df)
    println("----------------------------------------")
end

# ------------------ Main functions for analysis correlation ----------------- #

"""
    CorrelationAnalysis(filename="spsm_r.bin", filedir=pwd();
                       startbin=2, endbin=nothing, dropmaxmin=0,
                       real_column=3, imag_column=4,
                       auto_digits=true,
                       verbose=true)

Analyze correlation function data files like `spsm_r.bin` or `nn_r.bin`, where the first two columns
represent imj coordinates (distance vector between points i and j), and calculate Monte Carlo averages.

Arguments:
- `filename`: Name of the correlation file (default: "spsm_r.bin")
- `filedir`: Directory containing the file (default: current directory)
- `startbin`: First bin to include in analysis (default: 2)
- `endbin`: Last bin to include in analysis (default: all bins)
- `dropmaxmin`: Number of maximum and minimum values to drop (default: 0)
- `real_column`: Column index containing the real part of correlation values (default: 3)
- `imag_column`: Column index containing the imaginary part of correlation values (default: 4)
- `auto_digits`: Whether to automatically determine significant digits based on error of error (default: true)
- `verbose`: Whether to print results to console (default: true)

Returns:
- DataFrame containing the following columns:
  - `coord`: Tuple of (i,j) coordinates
  - `imj_x`, `imj_y`: x and y components of the distance vector
  - `mean_real`, `mean_imag`: Mean values of real and imaginary parts
  - `err_real`, `err_imag`: Errors of real and imaginary parts
  - `formatted_real`, `formatted_imag`: Formatted strings of values with errors
  - `distance`: Euclidean distance from origin
  - `i_coord`: Absolute coordinates relative to (1,1)

Example:
```julia
# Analyze spsm_r.bin in current directory
results = CorrelationAnalysis()

# Analyze nn_r.bin with custom columns
results = CorrelationAnalysis("nn_r.bin", real_column=3, imag_column=4)

# Access results
real_part = results[1, :mean_real]      # mean of real part
imag_part = results[1, :mean_imag]      # mean of imaginary part
real_err = results[1, :err_real]        # error of real part
imag_err = results[1, :err_imag]        # error of imaginary part

# Get coordinates and values
for row in eachrow(results)
    coord = row.coord
    real_val = row.mean_real
    imag_val = row.mean_imag
    real_e = row.err_real
    imag_e = row.err_imag
    # Format as complex number for display
    complex_val = "(\$real_val ± \$real_e) + (\$imag_val ± \$imag_e)im"
    println("At \$coord: \$complex_val")
end
```
"""
function CorrelationAnalysis(filename::String="spsm_r.bin", filedir::String=pwd();
                            startbin::Int=2, endbin::Union{Int,Nothing}=nothing, dropmaxmin::Int=0,
                            real_column::Int=3, imag_column::Int=4,
                            auto_digits::Bool=true,
                            verbose::Bool=true)
    # 文件读取和基本验证
    if !endswith(filename, ".bin")
        filename = filename * ".bin"
    end
    
    filepath = joinpath(filedir, filename)
    if !isfile(filepath)
        error("File not found: $filepath")
    end
    
    # 读取数据
    data_array = readdlm(filepath, Float64)
    
    # 创建DataFrame
    df = DataFrame(
        :imj_x => Int.(data_array[:, 1]),
        :imj_y => Int.(data_array[:, 2]),
        :real_val => data_array[:, real_column],
        :imag_val => data_array[:, imag_column]
    )
    
    # 推断格子大小和bin结构
    Lx = Int(maximum(df.imj_x))
    Ly = Int(maximum(df.imj_y))
    
    # 创建坐标元组列用于分组
    df.coord = [(x, y) for (x, y) in zip(df.imj_x, df.imj_y)]
    
    # 计算bin信息
    unique_coords = unique(df.coord)
    n_coords = length(unique_coords)
    n_bins = nrow(df) ÷ n_coords
    
    # 添加bin列
    df.bin = repeat(1:n_bins, inner=n_coords)
    
    # 添加复数值列
    df.complex_val = Complex.(df.real_val, df.imag_val)
    
    # 输出数据结构信息
    if verbose
        println("File structure analysis:")
        println("  Total rows: $(nrow(df))")
        println("  Unique coordinates: $n_coords")
        println("  Estimated bins: $n_bins")
    end
    
    # 应用bin筛选，对于单轨道数据使用 real_val 和 imag_val 列
    if dropmaxmin > 0
        filtered_df = filter_bins(df, startbin, endbin, dropmaxmin, n_bins, verbose, 
                               value_columns=[:real_val, :imag_val])
    else
        filtered_df = filter_bins(df, startbin, endbin, dropmaxmin, n_bins, verbose)
    end
    
    # 计算统计结果
    results_df = calculate_statistics(filtered_df, auto_digits)
    
    # 计算额外信息（距离和i坐标）
    results_df.distance = [euclidean_distance(x, y) for (x, y) in zip(results_df.imj_x, results_df.imj_y)]
    results_df.i_coord = [calculate_i_coord(x, y, Lx, Ly) for (x, y) in zip(results_df.imj_x, results_df.imj_y)]
    
    # 按距离排序
    sort!(results_df, [:distance, :imj_x, :imj_y])
    
    # 打印结果
    if verbose
        print_correlation_results(results_df, filename, Lx, Ly, dropmaxmin, startbin, endbin, n_bins)
    end
    
    return results_df
end


"""
    MultiOrbitalCorrelationAnalysis(filename="nn_r.bin", filedir=pwd();
                                  startbin=2, endbin=nothing, dropmaxmin=0,
                                  orbital_columns=[(3,4), (5,6), (7,8), (9,10)],
                                  orbital_labels=["AA", "AB", "BA", "BB"],
                                  auto_digits=true,
                                  verbose=true)

Analyze multi-orbital correlation function data files, where each orbital pair has its own
real and imaginary columns.

Arguments:
- `filename`: Name of the correlation file (default: "nn_r.bin")
- `filedir`: Directory containing the file (default: current directory)
- `startbin`: First bin to include in analysis (default: 2)
- `endbin`: Last bin to include in analysis (default: all bins)
- `dropmaxmin`: Number of maximum and minimum values to drop (default: 0)
- `orbital_columns`: Array of (real_column, imag_column) tuples for each orbital (default: [(3,4), (5,6), (7,8), (9,10)])
- `orbital_labels`: Array of labels for each orbital (default: ["AA", "AB", "BA", "BB"])
- `auto_digits`: Whether to automatically determine significant digits based on error of error (default: true)
- `verbose`: Whether to print results to console (default: true)

Returns:
- DataFrame containing statistics for each coordinate and orbital

Example:
```julia
# Analyze nn_r.bin in current directory with default orbital columns
results = MultiOrbitalCorrelationAnalysis()

# Access results
corr_1_1_AA_real = results[1, :AA_mean_real]  # mean real value at first row for AA orbital
corr_1_1_AB_imag = results[1, :AB_mean_imag]  # mean imag value at first row for AB orbital
```
"""
function MultiOrbitalCorrelationAnalysis(filename::String="nn_r.bin", filedir::String=pwd();
                                      startbin::Int=2, endbin::Union{Int,Nothing}=nothing, dropmaxmin::Int=0,
                                      orbital_columns::Vector{Tuple{Int,Int}}=[(3,4), (5,6), (7,8), (9,10)],
                                      orbital_labels::Vector{String}=["AA", "AB", "BA", "BB"],
                                      auto_digits::Bool=true,
                                      verbose::Bool=true)
    # Validate inputs
    if length(orbital_columns) != length(orbital_labels)
        error("Number of orbital columns ($(length(orbital_columns))) must match number of orbital labels ($(length(orbital_labels)))")
    end
    
    # File reading and validation
    if !endswith(filename, ".bin")
        filename = filename * ".bin"
    end
    
    filepath = joinpath(filedir, filename)
    if !isfile(filepath)
        error("File not found: $filepath")
    end
    
    # Read data
    data_array = readdlm(filepath, Float64)
    
    # Create base DataFrame
    df = DataFrame(
        :imj_x => Int.(data_array[:, 1]),
        :imj_y => Int.(data_array[:, 2])
    )
    
    # Add columns for each orbital
    for (i, (real_col, imag_col)) in enumerate(orbital_columns)
        orbital = orbital_labels[i]
        df[!, Symbol("$(orbital)_real")] = data_array[:, real_col]
        df[!, Symbol("$(orbital)_imag")] = data_array[:, imag_col]
    end
    
    # Add coordinate and bin information
    Lx = Int(maximum(df.imj_x))
    Ly = Int(maximum(df.imj_y))
    df.coord = [(x, y) for (x, y) in zip(df.imj_x, df.imj_y)]
    
    # Calculate bin information
    unique_coords = unique(df.coord)
    n_coords = length(unique_coords)
    n_bins = nrow(df) ÷ n_coords
    df.bin = repeat(1:n_bins, inner=n_coords)
    
    # Print data structure info
    if verbose
        println("File structure analysis:")
        println("  Total rows: $(nrow(df))")
        println("  Unique coordinates: $n_coords")
        println("  Estimated bins: $n_bins")
    end
    
    # Apply bin filtering with all orbital columns
    value_columns = Symbol[]
    for orbital in orbital_labels
        push!(value_columns, Symbol("$(orbital)_real"))
        push!(value_columns, Symbol("$(orbital)_imag"))
    end
    
    filtered_df = filter_bins(df, startbin, endbin, dropmaxmin, n_bins, verbose, value_columns=value_columns)
    
    # Calculate statistics
    results_df = calculate_statistics_multi_orbital(filtered_df, orbital_labels, auto_digits)
    
    # Add additional information
    results_df.distance = [euclidean_distance(x, y) for (x, y) in zip(results_df.imj_x, results_df.imj_y)]
    results_df.i_coord = [calculate_i_coord(x, y, Lx, Ly) for (x, y) in zip(results_df.imj_x, results_df.imj_y)]
    
    # Sort by distance
    sort!(results_df, [:distance, :imj_x, :imj_y])
    
    # Print results
    if verbose
        print_correlation_results_multi_orbital(results_df, filename, Lx, Ly, orbital_labels, dropmaxmin, startbin, endbin, n_bins)
    end
    
    return results_df
end
