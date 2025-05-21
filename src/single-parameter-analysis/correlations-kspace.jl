#= 
单参数分析模块 (Single Parameter Analysis) 之动量空间关联函数分析

此模块提供了一组函数，用于分析单个参数文件夹中的DQMC模拟数据。
通常在模拟输出目录中直接运行，处理该目录下的 .bin 文件。

主要功能:
- 单轨道关联函数分析 (KSpaceCorrelationAnalysis)
- 多轨道关联函数分析 (MultiOrbitalKSpaceCorrelationAnalysis)

=#

export KSpaceCorrelationAnalysis, MultiOrbitalKSpaceCorrelationAnalysis

# ----------------- Helper functions for correlation analysis ---------------- #

"""
    print_kspace_correlation_results(results_df, filename, dropmaxmin, startbin, endbin, n_bins)

打印k空间关联分析结果。
"""
function print_kspace_correlation_results(results_df, filename, dropmaxmin, startbin, endbin, n_bins)
    println("\nK-space correlation analysis of $filename:")
    println("----------------------------------------")
    
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
        "k-point" => [string("(", round(r.kx, digits=3), ",", round(r.ky, digits=3), ")") for r in eachrow(results_df)],
        "k_magnitude" => [round(r.k_magnitude, digits=3) for r in eachrow(results_df)],
        "real part" => results_df.formatted_real,
        "imag part" => results_df.formatted_imag
    )
    
    println(display_df)
    println("----------------------------------------")
end

"""
    print_kspace_correlation_results_multi_orbital(results_df, filename, orbital_labels, dropmaxmin, startbin, endbin, n_bins)

打印多轨道k空间关联分析结果。
"""
function print_kspace_correlation_results_multi_orbital(results_df, filename, orbital_labels, dropmaxmin, startbin, endbin, n_bins)
    println("\nMulti-orbital k-space correlation analysis of $filename:")
    println("----------------------------------------")
    println("Orbitals analyzed: $(join(orbital_labels, ", "))")
    
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
    
    # 为每个轨道打印详细结果
    for orbital in orbital_labels
        println("\n----------------------------------------")
        println("Orbital: $orbital")
        println("----------------------------------------")
        
        # 创建显示用DataFrame
        display_df = DataFrame(
            "k-point" => [string("(", round(r.kx, digits=3), ",", round(r.ky, digits=3), ")") for r in eachrow(results_df)],
            "k_magnitude" => [round(r.k_magnitude, digits=3) for r in eachrow(results_df)],
            "real part" => results_df[:, "$(orbital)_formatted_real"],
            "imag part" => results_df[:, "$(orbital)_formatted_imag"]
        )
        
        println(display_df)
    end
    
    # 打印简洁摘要（仅实部）
    println("\n----------------------------------------")
    println("Compact summary (real parts only):")
    println("----------------------------------------")
    
    # 构建简洁表格
    compact_df = DataFrame("k-point" => [string("(", round(r.kx, digits=3), ",", round(r.ky, digits=3), ")") for r in eachrow(results_df)])
    compact_df[!, "k_magnitude"] = [round(r.k_magnitude, digits=3) for r in eachrow(results_df)]
    
    for orbital in orbital_labels
        compact_df[!, orbital] = results_df[:, "$(orbital)_formatted_real"]
    end
    
    println(compact_df)
    println("----------------------------------------")
end


# ------------------ Main functions for analysis correlation ----------------- #

"""
    KSpaceCorrelationAnalysis(filename="nn_k.bin", filedir=pwd();
                             startbin=2, endbin=nothing, dropmaxmin=0,
                             real_column=3, imag_column=4,
                             auto_digits=true,
                             verbose=true)

Analyze k-space correlation function data files like `nn_k.bin` or `spsm_k.bin`, where the first two columns
represent k-space coordinates, and calculate Monte Carlo averages.

Arguments:
- `filename`: Name of the correlation file (default: "nn_k.bin")
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
  - `coord`: Tuple of k-point coordinates
  - `kx`, `ky`: x and y components of the k-point
  - `mean_real`, `mean_imag`: Mean values of real and imaginary parts
  - `err_real`, `err_imag`: Errors of real and imaginary parts
  - `formatted_real`, `formatted_imag`: Formatted strings of values with errors
  - `k_magnitude`: Magnitude of the k-point

Example:
```julia
# Analyze nn_k.bin in current directory
results = KSpaceCorrelationAnalysis()

# Analyze spsm_k.bin with custom columns
results = KSpaceCorrelationAnalysis("spsm_k.bin", real_column=3, imag_column=4)

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
function KSpaceCorrelationAnalysis(filename::String="nn_k.bin", filedir::String=pwd();
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
        :kx => data_array[:, 1],
        :ky => data_array[:, 2],
        :real_val => data_array[:, real_column],
        :imag_val => data_array[:, imag_column]
    )
    
    # 创建坐标元组列用于分组
    df.coord = [(kx, ky) for (kx, ky) in zip(df.kx, df.ky)]
    
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
        println("  Unique k-points: $n_coords")
        println("  Estimated bins: $n_bins")
    end
    
    # 应用bin筛选
    filtered_df = filter_bins(df, startbin, endbin, dropmaxmin, n_bins, verbose,
                          value_columns=[:real_val, :imag_val])
    
    # 计算统计结果
    results_df = calculate_statistics(
        filtered_df, auto_digits;
        input_coord_columns=[:kx, :ky],
        output_coord_names=[:kx, :ky]
    )
    
    # 计算k点的模长
    results_df.k_magnitude = [sqrt(kx^2 + ky^2) for (kx, ky) in zip(results_df.kx, results_df.ky)]
    
    # 按k点模长排序
    sort!(results_df, [:k_magnitude, :kx, :ky])
    
    # 打印结果
    if verbose
        print_kspace_correlation_results(results_df, filename, dropmaxmin, startbin, endbin, n_bins)
    end
    
    return results_df
end

"""
    MultiOrbitalKSpaceCorrelationAnalysis(filename="nn_k.bin", filedir=pwd();
                                        startbin=2, endbin=nothing, dropmaxmin=0,
                                        orbital_columns=[(3,4), (5,6), (7,8), (9,10)],
                                        orbital_labels=["AA", "AB", "BA", "BB"],
                                        auto_digits=true,
                                        verbose=true)

Analyze multi-orbital k-space correlation function data files, where each orbital pair has its own
real and imaginary columns.

Arguments:
- `filename`: Name of the correlation file (default: "nn_k.bin")
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
# Analyze nn_k.bin in current directory with default orbital columns
results = MultiOrbitalKSpaceCorrelationAnalysis()

# Access results
corr_1_1_AA_real = results[1, :AA_mean_real]  # mean real value at first row for AA orbital
corr_1_1_AB_imag = results[1, :AB_mean_imag]  # mean imag value at first row for AB orbital
```
"""
function MultiOrbitalKSpaceCorrelationAnalysis(filename::String="nn_k.bin", filedir::String=pwd();
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
        :kx => data_array[:, 1],
        :ky => data_array[:, 2]
    )
    
    # Add columns for each orbital
    for (i, (real_col, imag_col)) in enumerate(orbital_columns)
        orbital = orbital_labels[i]
        df[!, Symbol("$(orbital)_real")] = data_array[:, real_col]
        df[!, Symbol("$(orbital)_imag")] = data_array[:, imag_col]
    end
    
    # Add coordinate and bin information
    df.coord = [(kx, ky) for (kx, ky) in zip(df.kx, df.ky)]
    
    # Calculate bin information
    unique_coords = unique(df.coord)
    n_coords = length(unique_coords)
    n_bins = nrow(df) ÷ n_coords
    df.bin = repeat(1:n_bins, inner=n_coords)
    
    # Print data structure info
    if verbose
        println("File structure analysis:")
        println("  Total rows: $(nrow(df))")
        println("  Unique k-points: $n_coords")
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
        input_coord_columns=[:kx, :ky],
        output_coord_names=[:kx, :ky]
    )
    
    # Add additional information
    results_df.k_magnitude = [sqrt(kx^2 + ky^2) for (kx, ky) in zip(results_df.kx, results_df.ky)]
    
    # Sort by k-point magnitude
    sort!(results_df, [:k_magnitude, :kx, :ky])
    
    # Print results
    if verbose
        print_kspace_correlation_results_multi_orbital(results_df, filename, orbital_labels, dropmaxmin, startbin, endbin, n_bins)
    end
    
    return results_df
end