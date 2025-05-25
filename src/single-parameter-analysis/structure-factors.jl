#= 
单参数分析模块 (Single Parameter Analysis) 之结构因子分析

此模块提供了一组函数，用于从k空间关联函数数据中提取特定动量点的结构因子。
支持单轨道和多轨道系统，并提供计算反铁磁结构因子的专用函数。

主要功能:
- 单轨道结构因子分析 (StructureFactorAnalysis)
- 多轨道结构因子分析 (MultiOrbitalStructureFactorAnalysis)
- 反铁磁结构因子分析 (AFMStructureFactor)
- 反铁磁相关比率分析 (AFMCorrelationRatio)

=#

export StructureFactorAnalysis, MultiOrbitalStructureFactorAnalysis
export AFMStructureFactor, CDWStructureFactor, AFMCorrelationRatio

# ----------------- Helper functions for structure factor analysis ---------------- #

"""
    find_closest_k_point(k_points, target_k, tolerance)

在k点集合中找到最接近目标k点的坐标。

参数:
- `k_points`: 包含k点坐标的元组构成的向量 [(kx1, ky1), (kx2, ky2), ...]
- `target_k`: 目标k点 (kx, ky)
- `tolerance`: 匹配k点时的容差

返回:
- 最接近的k点坐标 (kx, ky)
- 是否为精确匹配
- 距离目标k点的距离
"""
function find_closest_k_point(k_points::Vector{<:Tuple{<:Real,<:Real}}, target_k::Tuple{<:Real,<:Real}, tolerance::Float64)
    # 计算每个k点与目标k点的距离
    distances = [(kx, ky, sqrt((kx - target_k[1])^2 + (ky - target_k[2])^2)) for (kx, ky) in k_points]
    
    # 按距离排序
    sort!(distances, by = x -> x[3])
    
    # 检查最接近的点是否在容差范围内
    closest_k = (distances[1][1], distances[1][2])
    exact_match = distances[1][3] <= tolerance
    
    # 返回最接近的k点、是否精确匹配、距离
    return closest_k, exact_match, distances[1][3]
end

"""
    print_structure_factor_result(k_point, mean_real, mean_imag, err_real, err_imag, 
                                 formatted_real, formatted_imag, filename, orbital="")

打印结构因子分析结果。

参数:
- `k_point`: 分析的k点
- `mean_real`, `mean_imag`: 平均值（实部和虚部）
- `err_real`, `err_imag`: 误差（实部和虚部）
- `formatted_real`, `formatted_imag`: 格式化后的结果字符串
- `filename`: 分析的文件名
- `orbital`: 轨道标签（可选，用于多轨道分析）
"""
function print_structure_factor_result(k_point, mean_real, mean_imag, err_real, err_imag, 
                                      formatted_real, formatted_imag, filename, orbital="")
    println("\nStructure Factor Analysis of $filename:")
    println("----------------------------------------")
    
    k_str = string("(", round(k_point[1], digits=3), ",", round(k_point[2], digits=3), ")")
    
    if !isempty(orbital)
        println("Orbital: $orbital")
    end
    
    println("k-point: $k_str")
    println("Real part: $formatted_real")
    println("Imag part: $formatted_imag")
    println("----------------------------------------")
end

"""
    print_af_structure_factor_result(k_point, mean_value, err_value, formatted_value, filename)

打印反铁磁结构因子分析结果。

参数:
- `k_point`: 分析的k点
- `mean_value`: 平均值
- `err_value`: 误差
- `formatted_value`: 格式化后的结果字符串
- `filename`: 分析的文件名
"""
function print_af_structure_factor_result(k_point, mean_real, mean_imag, err_real, err_imag, 
                                formatted_real, formatted_imag, filename)
    println("\nAntiferromagnetic Structure Factor Analysis of $filename:")
    println("----------------------------------------")
    
    k_str = string("(", round(k_point[1], digits=3), ",", round(k_point[2], digits=3), ")")
    
    println("k-point: $k_str")
    println("S_AF (real): $formatted_real")
    println("S_AF (imag): $formatted_imag")
    println("Formula: S_AF = [spsm_k(0,A,A) + spsm_k(0,B,B) - spsm_k(0,A,B) - spsm_k(0,B,A)]")
    println("----------------------------------------")
end

# ------------------ Main functions for structure factor analysis ----------------- #

"""
    StructureFactorAnalysis(k_point, filename="spsm_k.bin", filedir=pwd();
                           startbin=2, endbin=nothing, dropmaxmin=0,
                           real_column=3, imag_column=4,
                           auto_digits=true, tolerance=1e-6, verbose=true)

从k空间关联函数数据中提取特定动量点的结构因子。

参数:
- `k_point`: 目标动量点，格式为 (kx, ky) 的元组
- `filename`: 关联函数文件名 (默认: "spsm_k.bin")
- `filedir`: 文件目录 (默认: 当前目录)
- `startbin`: 起始bin (默认: 2)
- `endbin`: 结束bin (默认: 所有bin)
- `dropmaxmin`: 丢弃的最大最小值数量 (默认: 0)
- `real_column`: 实部列索引 (默认: 3)
- `imag_column`: 虚部列索引 (默认: 4)
- `auto_digits`: 是否自动确定有效数字 (默认: true)
- `tolerance`: 匹配k点时的容差 (默认: 1e-6)
- `verbose`: 是否输出详细信息 (默认: true)

返回:
- 包含以下字段的命名元组:
  - `k_point`: 实际使用的k点坐标
  - `mean_real`, `mean_imag`: 平均值（实部和虚部）
  - `err_real`, `err_imag`: 误差（实部和虚部）
  - `formatted_real`, `formatted_imag`: 格式化后的结果字符串

示例:
```julia
# 分析 (π,π) 点的结构因子
result = StructureFactorAnalysis((π, π), "spsm_k.bin")
println("Real part: \$(result.mean_real) ± \$(result.err_real)")
println("Imag part: \$(result.mean_imag) ± \$(result.err_imag)")
```
"""
function StructureFactorAnalysis(k_point::Tuple{<:Real,<:Real}, filename::String="spsm_k.bin", filedir::String=pwd();
                                startbin::Int=2, endbin::Union{Int,Nothing}=nothing, dropmaxmin::Int=0,
                                real_column::Int=3, imag_column::Int=4,
                                auto_digits::Bool=true, tolerance::Float64=1e-6, verbose::Bool=true)
    
    # 构建单轨道配置
    orbital_columns = [(real_column, imag_column)]
    orbital_labels = ["default"]
    
    # 调用多轨道函数，但禁用verbose以避免重复输出
    result = MultiOrbitalStructureFactorAnalysis(
        k_point, "default", filename, filedir;
        startbin=startbin, endbin=endbin, dropmaxmin=dropmaxmin,
        orbital_columns=orbital_columns, orbital_labels=orbital_labels,
        auto_digits=auto_digits, tolerance=tolerance, 
        verbose=false  # 禁用verbose，我们将在下面自己处理输出
    )
    
    # 如果需要打印，使用原始的单轨道格式
    if verbose
        print_structure_factor_result(
            result.k_point, result.mean_real, result.mean_imag, 
            result.err_real, result.err_imag,
            result.formatted_real, result.formatted_imag, filename
        )
    end
    
    # 返回与原始函数相同格式的结果（不包含orbital字段）
    return (
        k_point = result.k_point,
        mean_real = result.mean_real,
        mean_imag = result.mean_imag,
        err_real = result.err_real,
        err_imag = result.err_imag,
        formatted_real = result.formatted_real,
        formatted_imag = result.formatted_imag
    )
end

"""
    MultiOrbitalStructureFactorAnalysis(k_point, orbital_pair, filename="spsm_k.bin", filedir=pwd();
                                       startbin=2, endbin=nothing, dropmaxmin=0,
                                       orbital_columns=[(3,4), (5,6), (7,8), (9,10)],
                                       orbital_labels=["AA", "AB", "BA", "BB"],
                                       auto_digits=true, tolerance=1e-6, verbose=true)

从多轨道k空间关联函数数据中提取特定动量点和轨道对的结构因子。

参数:
- `k_point`: 目标动量点，格式为 (kx, ky) 的元组
- `orbital_pair`: 轨道对标签 (如 "AA", "AB" 等)
- `filename`: 关联函数文件名 (默认: "spsm_k.bin")
- `filedir`: 文件目录 (默认: 当前目录)
- `startbin`: 起始bin (默认: 2)
- `endbin`: 结束bin (默认: 所有bin)
- `dropmaxmin`: 丢弃的最大最小值数量 (默认: 0)
- `orbital_columns`: 轨道列索引数组 (默认: [(3,4), (5,6), (7,8), (9,10)])
- `orbital_labels`: 轨道标签数组 (默认: ["AA", "AB", "BA", "BB"])
- `auto_digits`: 是否自动确定有效数字 (默认: true)
- `tolerance`: 匹配k点时的容差 (默认: 1e-6)
- `verbose`: 是否输出详细信息 (默认: true)

返回:
- 包含以下字段的命名元组:
  - `k_point`: 实际使用的k点坐标
  - `orbital`: 轨道对标签
  - `mean_real`, `mean_imag`: 平均值（实部和虚部）
  - `err_real`, `err_imag`: 误差（实部和虚部）
  - `formatted_real`, `formatted_imag`: 格式化后的结果字符串

示例:
```julia
# 分析 (0,0) 点的AA轨道结构因子
result = MultiOrbitalStructureFactorAnalysis((0.0, 0.0), "AA", "spsm_k.bin")
println("AA orbital at Γ point: \$(result.formatted_real)")
```
"""
function MultiOrbitalStructureFactorAnalysis(k_point::Tuple{<:Real,<:Real}, orbital_pair::String, 
                                           filename::String="spsm_k.bin", filedir::String=pwd();
                                           startbin::Int=2, endbin::Union{Int,Nothing}=nothing, dropmaxmin::Int=0,
                                           orbital_columns::Vector{Tuple{Int,Int}}=[(3,4), (5,6), (7,8), (9,10)],
                                           orbital_labels::Vector{String}=["AA", "AB", "BA", "BB"],
                                           auto_digits::Bool=true, tolerance::Float64=1e-6, verbose::Bool=true)
    # Validate inputs
    @assert length(orbital_columns) == length(orbital_labels) "Number of orbital columns ($(length(orbital_columns))) must match number of orbital labels ($(length(orbital_labels)))"
    
    @assert orbital_pair in orbital_labels "Orbital pair '$orbital_pair' not found in orbital_labels: $(join(orbital_labels, ", "))"
    
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
    
    # 找到最接近的k点
    k_points = unique([(data_array[i, 1], data_array[i, 2]) for i in 1:size(data_array, 1)])
    closest_k, exact_match, distance = find_closest_k_point(k_points, k_point, tolerance)
    
    if verbose
        if !exact_match
            println("Requested k-point $(k_point) not found exactly.")
            println("Using closest k-point: $(closest_k) (distance: $(round(distance, digits=6)))")
        else
            println("Exact match found for k-point: $(closest_k)")
        end
    end
    
    # 直接过滤出目标k点的数据
    k_indices = findall(i -> isapprox(data_array[i, 1], closest_k[1], atol=tolerance) && 
                           isapprox(data_array[i, 2], closest_k[2], atol=tolerance), 
                        1:size(data_array, 1))
    
    if isempty(k_indices)
        @error "No data found for k-point $(closest_k)."
    end
    
    # 提取该k点的数据
    k_data = data_array[k_indices, :]
    
    # 创建DataFrame，只包含该k点的数据
    k_df = DataFrame(
        :kx => k_data[:, 1],
        :ky => k_data[:, 2],
        :bin => 1:length(k_indices)
    )
    
    # 添加轨道列
    for (i, (real_col, imag_col)) in enumerate(orbital_columns)
        orbital = orbital_labels[i]
        k_df[!, Symbol("$(orbital)_real")] = k_data[:, real_col]
        k_df[!, Symbol("$(orbital)_imag")] = k_data[:, imag_col]
    end
    
    # 计算bin信息
    n_bins = length(k_indices)
    
    # 准备要过滤的值列
    value_columns = Symbol[]
    for orbital in orbital_labels
        push!(value_columns, Symbol("$(orbital)_real"))
        push!(value_columns, Symbol("$(orbital)_imag"))
    end
    
    # 应用bin筛选（只对该k点的数据）
    filtered_k_df = filter_bins(k_df, startbin, endbin, dropmaxmin, n_bins, verbose,
                              value_columns=value_columns)
    
    if isempty(filtered_k_df)
        @error "No data found for k-point $(closest_k) after filtering."
    end
    
    # 获取指定轨道对的列名
    real_col = Symbol("$(orbital_pair)_real")
    imag_col = Symbol("$(orbital_pair)_imag")
    
    # 使用compute_stats函数计算统计量
    stats = compute_stats(filtered_k_df[!, real_col], filtered_k_df[!, imag_col], auto_digits=auto_digits)
    
    # 打印结果
    if verbose
        print_structure_factor_result(
            closest_k, stats.mean_real, stats.mean_imag, stats.err_real, stats.err_imag,
            stats.formatted_real, stats.formatted_imag, filename, orbital_pair
        )
    end
    
    # 返回结果
    return (
        k_point = closest_k,
        orbital = orbital_pair,
        mean_real = stats.mean_real,
        mean_imag = stats.mean_imag,
        err_real = stats.err_real,
        err_imag = stats.err_imag,
        formatted_real = stats.formatted_real,
        formatted_imag = stats.formatted_imag
    )
end

"""    
    AFMStructureFactor(k_point=(0.0, 0.0), filename="afm_sf_k.bin", filedir=pwd();
                     source_file="spsm_k.bin", startbin=2, endbin=nothing, dropmaxmin=0,
                     auto_digits=true, tolerance=1e-6, verbose=true)

Calculate antiferromagnetic structure factor S_AF(L) = [spsm_k(0,A,A) + spsm_k(0,B,B) - spsm_k(0,A,B) - spsm_k(0,B,A)].

The function can use either:
- A pre-processed file (default: "afm_sf_k.bin") containing the structure factor directly
- A source file (e.g., "spsm_k.bin" or "ss_k.bin") to generate the structure factor file if it doesn't exist

Parameters:
- `k_point`: Target momentum point, default (0.0, 0.0)
- `filename`: Target structure factor file name (default: "afm_sf_k.bin")
- `filedir`: Directory containing the files (default: current directory)
- `source_file`: Source file to generate structure factor if needed (default: "spsm_k.bin")
- `startbin`: Starting bin for statistics (default: 2)
- `endbin`: Ending bin for statistics (default: all bins)
- `dropmaxmin`: Number of max/min values to drop (default: 0)
- `auto_digits`: Whether to automatically determine significant digits (default: true)
- `tolerance`: Tolerance for k-point matching (default: 1e-6)
- `verbose`: Whether to output detailed information (default: true)

Returns:
- A named tuple with the following fields:
  - `k_point`: The actual k-point used
  - `mean_real`: Mean value of the real part
  - `mean_imag`: Mean value of the imaginary part
  - `err_real`: Error estimate for the real part
  - `err_imag`: Error estimate for the imaginary part
  - `formatted_real`: Formatted string for the real part
  - `formatted_imag`: Formatted string for the imaginary part

Examples:
```julia
# Use existing afm_sf_k.bin file
result = AFMStructureFactor()

# Generate from spsm_k.bin if afm_sf_k.bin doesn't exist
result = AFMStructureFactor()

# Generate from ss_k.bin if afm_sf_k.bin doesn't exist
result = AFMStructureFactor(source_file="ss_k.bin")

# Specify custom filenames
result = AFMStructureFactor(filename="custom_afm_sf.bin", source_file="custom_ss.bin")
```
"""
function AFMStructureFactor(k_point::Tuple{<:Real,<:Real}=(0.0, 0.0), filename::String="afm_sf_k.bin", filedir::String=pwd();
                          source_file::String="spsm_k.bin", startbin::Int=2, endbin::Union{Int,Nothing}=nothing, dropmaxmin::Int=0,
                          auto_digits::Bool=true, tolerance::Float64=1e-6, verbose::Bool=true)
    # Ensure filename has .bin extension
    if !endswith(filename, ".bin")
        filename = filename * ".bin"
    end
    target_filepath = joinpath(filedir, filename)
    
    # Check if target file exists, if not, generate it from source_file
    if !isfile(target_filepath)
        if verbose
            println("Target file $filename not found. Generating from $source_file...")
        end
        
        # Ensure source_file has .bin extension
        if !endswith(source_file, ".bin")
            source_file = source_file * ".bin"
        end
        
        # Check if source file exists
        source_filepath = joinpath(filedir, source_file)
        if !isfile(source_filepath)
            @error "Source file not found: $source_filepath"
            return nothing
        end
        
        # Generate target file using merge_staggered_components
        # Use the same column indices for all file types
        merge_staggered_components(
            filename,            # Output filename
            source_file,         # Input filename
            filedir,            # Output directory
            filedir;            # Input directory
            real_columns=[3, 9, 5, 7],  # [AA, BB, AB, BA] real columns
            imag_columns=[4, 10, 6, 8], # [AA, BB, AB, BA] imag columns
            verbose=verbose
        )
        
        # Verify the file was successfully generated
        if !isfile(target_filepath)
            @error "Failed to generate $filename from $source_file"
            return nothing
        end
    end
    
    # Check if target file exists now
    if !isfile(target_filepath)
        @error "Target file not found: $target_filepath. Please provide a valid source_file to generate it."
        return nothing
    end
    
    # Analyze the structure factor file
    if verbose
        println("Analyzing $filename for k-point $k_point...")
    end
    
    # Call StructureFactorAnalysis function
    result = StructureFactorAnalysis(
        k_point,
        filename,
        filedir;
        startbin=startbin,
        endbin=endbin,
        dropmaxmin=dropmaxmin,
        auto_digits=auto_digits,
        tolerance=tolerance,
        verbose=verbose
    )
    
    return result
end

"""
    CDWStructureFactor(k_point=(0.0, 0.0), filename="cdwpair_sf_k.bin", filedir=pwd();
                     source_file="cdwpair_k.bin", startbin=2, endbin=nothing, dropmaxmin=0,
                     auto_digits=true, tolerance=1e-6, verbose=true)

Calculate charge density wave structure factor S_CDW(L) = [cdwpair_k(0,A,A) + cdwpair_k(0,B,B) + cdwpair_k(0,A,B) + cdwpair_k(0,B,A)].

The function can use either:
- A pre-processed file (default: "cdwpair_sf_k.bin") containing the structure factor directly
- A source file (e.g., "cdwpair_k.bin") to generate the structure factor file if it doesn't exist

Parameters:
- `k_point`: Target momentum point, default (0.0, 0.0)
- `filename`: Target structure factor file name (default: "cdwpair_sf_k.bin")
- `filedir`: Directory containing the files (default: current directory)
- `source_file`: Source file to generate structure factor if needed (default: "cdwpair_k.bin")
- `startbin`: Starting bin for statistics (default: 2)
- `endbin`: Ending bin for statistics (default: all bins)
- `dropmaxmin`: Number of max/min values to drop (default: 0)
- `auto_digits`: Whether to automatically determine significant digits (default: true)
- `tolerance`: Tolerance for k-point matching (default: 1e-6)
- `verbose`: Whether to output detailed information (default: true)

Returns:
- A named tuple with the following fields:
  - `k_point`: The actual k-point used
  - `mean_real`: Mean value of the real part
  - `mean_imag`: Mean value of the imaginary part
  - `err_real`: Error estimate for the real part
  - `err_imag`: Error estimate for the imaginary part
  - `formatted_real`: Formatted string for the real part
  - `formatted_imag`: Formatted string for the imaginary part

Examples:
```julia
# Use existing cdwpair_sf_k.bin file
result = CDWStructureFactor()

# Generate from cdwpair_k.bin if cdwpair_sf_k.bin doesn't exist
result = CDWStructureFactor()

# Specify custom filenames
result = CDWStructureFactor(filename="custom_cdw_sf.bin", source_file="custom_cdw.bin")
```
"""
function CDWStructureFactor(k_point::Tuple{<:Real,<:Real}=(0.0, 0.0), filename::String="cdwpair_sf_k.bin", filedir::String=pwd();
                          source_file::String="cdwpair_k.bin", startbin::Int=2, endbin::Union{Int,Nothing}=nothing, dropmaxmin::Int=0,
                          auto_digits::Bool=true, tolerance::Float64=1e-6, verbose::Bool=true)
    # Ensure filename has .bin extension
    if !endswith(filename, ".bin")
        filename = filename * ".bin"
    end
    target_filepath = joinpath(filedir, filename)
    
    # Check if target file exists, if not, generate it from source_file
    if !isfile(target_filepath)
        if verbose
            println("Target file $filename not found. Generating from $source_file...")
        end
        
        # Ensure source_file has .bin extension
        if !endswith(source_file, ".bin")
            source_file = source_file * ".bin"
        end
        
        # Check if source file exists
        source_filepath = joinpath(filedir, source_file)
        if !isfile(source_filepath)
            @error "Source file not found: $source_filepath"
            return nothing
        end
        
        # Generate target file using merge_uniform_components
        merge_uniform_components(
            filename,            # Output filename
            source_file,         # Input filename
            filedir,            # Output directory
            filedir;            # Input directory
            real_columns=[3, 9, 5, 7],  # [AA, BB, AB, BA] real columns
            imag_columns=[4, 10, 6, 8], # [AA, BB, AB, BA] imag columns
            verbose=verbose
        )
        
        # Verify the file was successfully generated
        if !isfile(target_filepath)
            @error "Failed to generate $filename from $source_file"
            return nothing
        end
    end
    
    # Check if target file exists now
    if !isfile(target_filepath)
        @error "Target file not found: $target_filepath. Please provide a valid source_file to generate it."
        return nothing
    end
    
    # Analyze the structure factor file
    if verbose
        println("Analyzing $filename for k-point $k_point...")
    end
    
    # Call StructureFactorAnalysis function
    result = StructureFactorAnalysis(
        k_point,
        filename,
        filedir;
        startbin=startbin,
        endbin=endbin,
        dropmaxmin=dropmaxmin,
        auto_digits=auto_digits,
        tolerance=tolerance,
        verbose=verbose
    )
    
    return result
end

"""    AFMCorrelationRatio(shift_point::Tuple{<:Real,<:Real}, Q_point::Tuple{<:Real,<:Real}=(0.0, 0.0),
                     filename::String="spsm_k.bin", filedir::String=pwd();
                     startbin::Int=2, endbin::Union{Int,Nothing}=nothing, dropmaxmin::Int=0,
                     orbital_columns::Vector{Tuple{Int,Int}}=[(3,4), (5,6), (7,8), (9,10)],
                     orbital_labels::Vector{String}=["AA", "AB", "BA", "BB"],
                     auto_digits::Bool=true, tolerance::Float64=1e-6, verbose::Bool=true)

计算自旋结构因子的相关比（Correlation Ratio）R_{m^2}，用于量化无序-有序相变。

公式: R_{m^2} = 1 - S_AFM(Q+δq) / S_AFM(Q)
其中: δq 是在倒格子空间中的小偏移

参数:
- `shift_point`: 动量空间偏移量 (δq_x, δq_y)
- `Q_point`: 反铁磁向量Q (默认: (0.0, 0.0))
- `filename`: 关联函数文件名 (默认: "spsm_k.bin")
- `filedir`: 文件目录 (默认: 当前目录)
- `startbin`: 起始bin (默认: 2)
- `endbin`: 结束bin (默认: 所有bin)
- `dropmaxmin`: 丢弃的最大最小值数量 (默认: 0)
- `orbital_columns`: 轨道列索引数组 (默认: [(3,4), (5,6), (7,8), (9,10)])
- `orbital_labels`: 轨道标签数组 (默认: ["AA", "AB", "BA", "BB"])
- `auto_digits`: 是否自动确定有效数字 (默认: true)
- `tolerance`: 匹配k点时的容差 (默认: 1e-6)
- `verbose`: 是否输出详细信息 (默认: true)

返回:
- 包含以下字段的命名元组:
  - `Q_point`: 反铁磁向量Q
  - `shift_point`: 动量空间偏移量 (δq_x, δq_y)
  - `Q_shifted`: 偏移后的k点
  - `S_AFM_Q`: Q点处的反铁磁结构因子
  - `err_S_AFM_Q`: Q点处的反铁磁结构因子误差
  - `S_AFM_Q_shifted`: 偏移点处的反铁磁结构因子
  - `err_S_AFM_Q_shifted`: 偏移点处的反铁磁结构因子误差
  - `correlation_ratio`: 计算得到的相关比
  - `err_correlation_ratio`: 相关比的误差
  - `formatted_correlation_ratio`: 格式化后的相关比结果

示例:
```julia
# 计算(0.25,0)偏移的相关比
result = AFMCorrelationRatio((0.25, 0.0))
println("Correlation Ratio: \$(result.formatted_correlation_ratio)")
```

"""
function AFMCorrelationRatio(shift_point::Tuple{<:Real,<:Real}, Q_point::Tuple{<:Real,<:Real}=(0.0, 0.0),
                         filename::String="spsm_k.bin", filedir::String=pwd();
                         startbin::Int=2, endbin::Union{Int,Nothing}=nothing, dropmaxmin::Int=0,
                         orbital_columns::Vector{Tuple{Int,Int}}=[(3,4), (5,6), (7,8), (9,10)],
                         orbital_labels::Vector{String}=["AA", "AB", "BA", "BB"],
                         auto_digits::Bool=true, tolerance::Float64=1e-6, verbose::Bool=true)
    
    # 计算Q点处的反铁磁结构因子
    if verbose
        println("计算Q点 $Q_point 处的反铁磁结构因子")
    end
    
    S_AFM_Q = AFMStructureFactor(Q_point, filename, filedir;
                              startbin=startbin, endbin=endbin, dropmaxmin=dropmaxmin,
                              orbital_columns=orbital_columns, orbital_labels=orbital_labels,
                              auto_digits=auto_digits, tolerance=tolerance, verbose=verbose)
    
    # 计算偏移后的k点坐标
    shift_x, shift_y = shift_point
    Q_shifted = (Q_point[1] + shift_x, Q_point[2] + shift_y)
    
    if verbose
        println("\n计算偏移点 $Q_shifted 处的反铁磁结构因子")
        println("偏移量: ($(round(shift_x, digits=6)), $(round(shift_y, digits=6)))")
    end
    
    # 计算偏移点处的反铁磁结构因子
    S_AFM_Q_shifted = AFMStructureFactor(Q_shifted, filename, filedir;
                                      startbin=startbin, endbin=endbin, dropmaxmin=dropmaxmin,
                                      orbital_columns=orbital_columns, orbital_labels=orbital_labels,
                                      auto_digits=auto_digits, tolerance=tolerance, verbose=verbose)
    
    # 使用实部计算相关比（通常反铁磁结构因子主要是实部）
    ratio = S_AFM_Q_shifted.mean_real / S_AFM_Q.mean_real
    correlation_ratio = 1.0 - ratio
    
    # 误差传播 (使用一阶近似)：
    # 如果 R = 1 - S2/S1，则 δR ≈ sqrt((S2/S1^2 * δS1)^2 + (1/S1 * δS2)^2)
    S1 = S_AFM_Q.mean_real
    S2 = S_AFM_Q_shifted.mean_real
    dS1 = S_AFM_Q.err_real
    dS2 = S_AFM_Q_shifted.err_real
    
    err_ratio = sqrt((S2/(S1^2) * dS1)^2 + (1/S1 * dS2)^2)
    
    # 使用round_error函数保留一位有效数字
    rounded_err, _ = round_error(err_ratio, err_ratio/10)
    
    # 使用format_value_error函数格式化结果
    formatted_val, formatted_err = format_value_error(correlation_ratio, rounded_err, 1)
    formatted_correlation_ratio = "$(formatted_val) ± $(formatted_err)"
    
    # 打印结果
    if verbose
        println("\n自旋结构因子相关比（Correlation Ratio）分析:")
        println("---------------------------------------------------")
        println("Q点: $Q_point")
        println("偏移点: $Q_shifted (偏移量: $shift_point)")
        println("S_AFM(Q) = $(S_AFM_Q.formatted_real)")
        println("S_AFM(Q+δq) = $(S_AFM_Q_shifted.formatted_real)")
        println("相关比 R = 1 - S_AFM(Q+δq)/S_AFM(Q) = $formatted_correlation_ratio")
        println("---------------------------------------------------")
    end
    
    # 返回结果
    return (
        Q_point = Q_point,
        shift_point = shift_point,
        Q_shifted = Q_shifted,
        S_AFM_Q = S_AFM_Q.mean_real,
        err_S_AFM_Q = S_AFM_Q.err_real,
        S_AFM_Q_shifted = S_AFM_Q_shifted.mean_real,
        err_S_AFM_Q_shifted = S_AFM_Q_shifted.err_real,
        correlation_ratio = correlation_ratio,
        err_correlation_ratio = err_ratio,
        formatted_correlation_ratio = formatted_correlation_ratio
    )
end
