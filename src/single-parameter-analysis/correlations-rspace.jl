#= 
单参数分析模块 (Single Parameter Analysis) 之实空间关联函数分析

此模块提供了一组函数，用于分析单个参数文件夹中的DQMC模拟数据。
通常在模拟输出目录中直接运行，处理该目录下的 .bin 文件。

主要功能:
- 单轨道关联函数分析 (CorrelationAnalysis)
- 多轨道关联函数分析 (MultiOrbitalCorrelationAnalysis)

=#

export CorrelationAnalysis, MultiOrbitalCorrelationAnalysis

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
  - `i_coord`: coordinates of unit cell i relative to j = (1,1)

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
        @error "File not found: $filepath"
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
    filtered_df = filter_bins(df, startbin, endbin, dropmaxmin, n_bins, verbose, 
                                value_columns=[:real_val, :imag_val])
    
    # 计算统计结果
    results_df = calculate_statistics(
        filtered_df, auto_digits;
        input_coord_columns=[:imj_x, :imj_y],
        output_coord_names=[:imj_x, :imj_y]
    )
    
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
                                  orbital_labels=["AA", "BA", "AB", "BB"],
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
- `orbital_labels`: Array of labels for each orbital (default: ["AA", "BA", "AB", "BB"])
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
                                      orbital_labels::Vector{String}=["AA", "BA", "AB", "BB"],
                                      auto_digits::Bool=true,
                                      verbose::Bool=true)
    # Validate inputs
    if length(orbital_columns) != length(orbital_labels)
        @error "Number of orbital columns ($(length(orbital_columns))) must match number of orbital labels ($(length(orbital_labels)))"
    end
    
    # File reading and validation
    if !endswith(filename, ".bin")
        filename = filename * ".bin"
    end
    
    filepath = joinpath(filedir, filename)
    if !isfile(filepath)
        @error "File not found: $filepath"
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
    filtered_df = filter_bins(df, startbin, endbin, dropmaxmin, n_bins, verbose, 
                            value_columns=value_columns)
    
    # 计算统计结果
    results_df = calculate_statistics(
        filtered_df, auto_digits;
        orbital_labels=orbital_labels,
        input_coord_columns=[:imj_x, :imj_y],
        output_coord_names=[:imj_x, :imj_y]
    )
    
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
