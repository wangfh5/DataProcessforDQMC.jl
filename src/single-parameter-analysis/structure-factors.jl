#= 
单参数分析模块 (Single Parameter Analysis) 之结构因子分析

此模块提供了一组函数，用于从k空间关联函数数据中提取特定动量点的结构因子。
支持单轨道和多轨道系统，并提供计算反铁磁结构因子的专用函数。

主要功能:
- 单轨道结构因子分析 (StructureFactorAnalysis)
- 多轨道结构因子分析 (MultiOrbitalStructureFactorAnalysis)
- 反铁磁结构因子分析 (AFStructureFactor)

=#

export StructureFactorAnalysis, MultiOrbitalStructureFactorAnalysis, AFStructureFactor

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
        :real_val => k_data[:, real_column],
        :imag_val => k_data[:, imag_column],
        :bin => 1:length(k_indices)
    )
    
    # 计算bin信息
    n_bins = length(k_indices)
    
    # 应用bin筛选（只对该k点的数据）
    filtered_k_df = filter_bins(k_df, startbin, endbin, dropmaxmin, n_bins, verbose,
                              value_columns=[:real_val, :imag_val])
    
    if isempty(filtered_k_df)
        @error "No data found for k-point $(closest_k) after filtering."
    end
    
    # 使用compute_stats函数计算统计量
    stats = compute_stats(filtered_k_df.real_val, filtered_k_df.imag_val, auto_digits=auto_digits)
    
    # 打印结果
    if verbose
        print_structure_factor_result(
            closest_k, stats.mean_real, stats.mean_imag, stats.err_real, stats.err_imag,
            stats.formatted_real, stats.formatted_imag, filename
        )
    end
    
    # 返回结果
    return (
        k_point = closest_k,
        mean_real = stats.mean_real,
        mean_imag = stats.mean_imag,
        err_real = stats.err_real,
        err_imag = stats.err_imag,
        formatted_real = stats.formatted_real,
        formatted_imag = stats.formatted_imag
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
    AFStructureFactor(k_point=(0.0, 0.0), filename="spsm_k.bin", filedir=pwd();
                     startbin=2, endbin=nothing, dropmaxmin=0,
                     orbital_columns=[(3,4), (5,6), (7,8), (9,10)],
                     orbital_labels=["AA", "AB", "BA", "BB"],
                     auto_digits=true, tolerance=1e-6, verbose=true)

计算反铁磁结构因子 S_AF(L) = [spsm_k(0,A,A) + spsm_k(0,B,B) - spsm_k(0,A,B) - spsm_k(0,B,A)]。

参数:
- `k_point`: 目标动量点，默认为 (0.0, 0.0)
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
  - `mean_value`: 反铁磁结构因子平均值
  - `err_value`: 反铁磁结构因子误差
  - `formatted_value`: 格式化后的结果字符串

示例:
```julia
# 计算Γ点的反铁磁结构因子
result = AFStructureFactor()
println("S_AF at Γ point: \$(result.formatted_value)")
```
"""
function AFStructureFactor(k_point::Tuple{<:Real,<:Real}=(0.0, 0.0), filename::String="spsm_k.bin", filedir::String=pwd();
                          startbin::Int=2, endbin::Union{Int,Nothing}=nothing, dropmaxmin::Int=0,
                          orbital_columns::Vector{Tuple{Int,Int}}=[(3,4), (5,6), (7,8), (9,10)],
                          orbital_labels::Vector{String}=["AA", "AB", "BA", "BB"],
                          auto_digits::Bool=true, tolerance::Float64=1e-6, verbose::Bool=true)
    # 验证轨道标签
    required_labels = ["AA", "AB", "BA", "BB"]
    if !all(label in orbital_labels for label in required_labels)
        @error "AFStructureFactor requires all of these orbital labels: $(join(required_labels, ", "))"
    end
    
    # 文件读取和验证
    if !endswith(filename, ".bin")
        filename = filename * ".bin"
    end
    
    filepath = joinpath(filedir, filename)
    if !isfile(filepath)
        @error "File not found: $filepath"
    end
    
    # 读取数据
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
    
    # 添加轨道列（实部和虚部）
    for (i, (real_col, imag_col)) in enumerate(orbital_columns)
        orbital = orbital_labels[i]
        k_df[!, Symbol("$(orbital)_real")] = k_data[:, real_col]
        k_df[!, Symbol("$(orbital)_imag")] = k_data[:, imag_col]
    end
    
    # 计算bin信息
    n_bins = length(k_indices)
    
    # 准备要过滤的值列（实部和虚部）
    value_columns = Symbol[]
    for orbital in orbital_labels
        push!(value_columns, Symbol("$(orbital)_real"))
        push!(value_columns, Symbol("$(orbital)_imag"))
    end
    
    # 应用bin筛选
    filtered_k_df = filter_bins(k_df, startbin, endbin, dropmaxmin, n_bins, verbose,
                           value_columns=value_columns)
    
    if isempty(filtered_k_df)
        @error "No data found for k-point $(closest_k) after filtering."
    end
    
    # 计算每个bin的反铁磁结构因子 S_AF = AA + BB - AB - BA（分别处理实部和虚部）
    af_real_values = [row.AA_real + row.BB_real - row.AB_real - row.BA_real for row in eachrow(filtered_k_df)]
    af_imag_values = [row.AA_imag + row.BB_imag - row.AB_imag - row.BA_imag for row in eachrow(filtered_k_df)]
    
    # 使用compute_stats函数计算复数的统计量
    stats = compute_stats(af_real_values, af_imag_values, auto_digits=auto_digits)
    
    # 打印结果
    if verbose
        print_af_structure_factor_result(
            closest_k, stats.mean_real, stats.mean_imag, stats.err_real, stats.err_imag, 
            stats.formatted_real, stats.formatted_imag, 
            filename
        )
    end
    
    # 返回结果
    return (
        k_point = closest_k,
        mean_real = stats.mean_real,
        mean_imag = stats.mean_imag,
        err_real = stats.err_real,
        err_imag = stats.err_imag,
        formatted_real = stats.formatted_real,
        formatted_imag = stats.formatted_imag
    )
end
