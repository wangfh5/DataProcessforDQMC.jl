#= 
单参数分析模块 (Single Parameter Analysis) 之结构因子分析

此模块提供了一组函数，用于从k空间关联函数数据中提取特定动量点的结构因子。
支持单/多动量点分析，并提供计算反铁磁结构因子的专用函数。

主要功能:
- 结构因子分析 (StructureFactorAnalysis) - 支持单/多动量点
- 反铁磁结构因子分析 (AFMStructureFactor)
- CDW结构因子分析 (CDWStructureFactor)

=#

export StructureFactorAnalysis
export AFMStructureFactor, CDWStructureFactor
export StructureFactorResult

# ==================================================================================== #
#                              结构体定义
# ==================================================================================== #

"""
    StructureFactorResult

单个动量点处结构因子分析结果的结构体。

# 字段
- `k_point::Tuple{Float64, Float64}`: 实际使用的k点坐标
- `mean_real::Float64`: 实部平均值
- `err_real::Float64`: 实部误差
- `mean_imag::Float64`: 虚部平均值
- `err_imag::Float64`: 虚部误差
- `formatted_real::String`: 格式化后的实部字符串
- `formatted_imag::String`: 格式化后的虚部字符串
"""
struct StructureFactorResult
    k_point::Tuple{Float64, Float64}
    mean_real::Float64
    err_real::Float64
    mean_imag::Float64
    err_imag::Float64
    formatted_real::String
    formatted_imag::String
end

# ==================================================================================== #
#                              辅助函数
# ==================================================================================== #

"""
    print_structure_factor_result(result::StructureFactorResult, filename::String, orbital::String="")

打印结构因子分析结果。
"""
function print_structure_factor_result(result::StructureFactorResult, filename::String, orbital::String="")
    println("\nStructure Factor Analysis of $filename:")
    println("----------------------------------------")
    
    k_str = string("(", round(result.k_point[1], digits=3), ",", round(result.k_point[2], digits=3), ")")
    
    if !isempty(orbital)
        println("Orbital: $orbital")
    end
    
    println("k-point: $k_str")
    println("Real part: $(result.formatted_real)")
    println("Imag part: $(result.formatted_imag)")
    println("----------------------------------------")
end

"""
    print_af_structure_factor_result(result::StructureFactorResult, filename::String)

打印反铁磁结构因子分析结果。
"""
function print_af_structure_factor_result(result::StructureFactorResult, filename::String)
    println("\nAntiferromagnetic Structure Factor Analysis of $filename:")
    println("----------------------------------------")
    
    k_str = string("(", round(result.k_point[1], digits=3), ",", round(result.k_point[2], digits=3), ")")
    
    println("k-point: $k_str")
    println("S_AF (real): $(result.formatted_real)")
    println("S_AF (imag): $(result.formatted_imag)")
    println("Formula: S_AF = [spsm_k(0,A,A) + spsm_k(0,B,B) - spsm_k(0,A,B) - spsm_k(0,B,A)]")
    println("----------------------------------------")
end

# ==================================================================================== #
#                              核心实现（多k点，一次读文件）
# ==================================================================================== #

"""
    _multi_k_structure_factor_analysis_core(k_points, filename, filedir;
                                           startbin, endbin, dropmaxmin,
                                           real_column, imag_column,
                                           auto_digits, k_point_tolerance, verbose) -> Vector{StructureFactorResult}

核心实现：一次读文件，处理多个k点。这是内部函数，不对外导出。
"""
function _multi_k_structure_factor_analysis_core(
    k_points::Vector{Tuple{Float64,Float64}},
    filename::String,
    filedir::String;
    startbin::Int=2,
    endbin::Union{Int,Nothing}=nothing,
    outlier_mode::Symbol=:dropmaxmin,
    outlier_param::Real=0,
    real_column::Int=3,
    imag_column::Int=4,
    auto_digits::Bool=true,
    k_point_tolerance::Float64=1e-6,
    verbose::Bool=true
)
    # 文件验证
    if !endswith(filename, ".bin")
        filename = filename * ".bin"
    end
    
    filepath = joinpath(filedir, filename)
    if !isfile(filepath)
        @error "File not found: $filepath"
        return StructureFactorResult[]
    end
    
    # ===== 一次性读取文件 =====
    data_array = readdlm(filepath, Float64)
    
    # 提取文件中所有可用的k点
    available_k_points = unique([(data_array[i, 1], data_array[i, 2]) for i in 1:size(data_array, 1)])
    
    # 结果容器
    results = StructureFactorResult[]
    
    # 处理每个请求的k点
    for target_k in k_points
        # 找到最接近的k点
        closest_k, exact_match, distance = find_closest_k_point(available_k_points, target_k, k_point_tolerance)
    
        if verbose
            if !exact_match
                println("Requested k-point $target_k not found exactly.")
                println("Using closest k-point: $closest_k (distance: $(round(distance, digits=6)))")
            else
                println("Exact match found for k-point: $closest_k")
            end
        end
        
        # 从已读取的数据中提取该k点的数据
        k_indices = findall(i -> isapprox(data_array[i, 1], closest_k[1], atol=k_point_tolerance) && 
                                isapprox(data_array[i, 2], closest_k[2], atol=k_point_tolerance), 
                                1:size(data_array, 1))
    
        if isempty(k_indices)
                @warn "No data found for k-point $closest_k, skipping..."
                continue
        end
    
        # 提取该k点的数据
        k_data = data_array[k_indices, :]

        # 选择 bin 范围
        n_bins = length(k_indices)
        start_idx = max(startbin, 1)
        end_idx = isnothing(endbin) ? n_bins : min(endbin, n_bins)
        if start_idx > end_idx
            @warn "Invalid bin range startbin=$startbin endbin=$endbin for k-point $closest_k, skipping..."
            continue
        end

        real_vals = k_data[start_idx:end_idx, real_column]
        imag_vals = k_data[start_idx:end_idx, imag_column]

        # 生成离群 mask（metric = abs(real + i*imag)）
        metric = abs.(complex.(real_vals, imag_vals))

        mask = outlier_filter(metric, outlier_mode, outlier_param; min_n=5)
        keep = mask.keep

        if verbose && outlier_mode == :iqrfence
            kept = count(keep)
            total = length(metric)
            println("IQR fence outlier summary: removed_min=$(mask.removed_min), removed_max=$(mask.removed_max), kept=$(kept)/$(total)")
        end

        if !any(keep)
            @warn "No finite data found for k-point $closest_k after outlier filtering, skipping..."
            continue
        end

        # 计算统计量
        stats = compute_stats(real_vals[keep], imag_vals[keep], auto_digits=auto_digits)
        
        # 创建结果结构体
        push!(results, StructureFactorResult(
            closest_k,
            stats.mean_real,
            stats.err_real,
            stats.mean_imag,
            stats.err_imag,
            stats.formatted_real,
            stats.formatted_imag
        ))
    end
    
    return results
end

# ==================================================================================== #
#                              公开接口：StructureFactorAnalysis
# ==================================================================================== #

"""
    StructureFactorAnalysis(k_points::Vector{<:Tuple{<:Real,<:Real}}, 
                           filename::String="spsm_k.bin", filedir::String=pwd();
                           startbin::Int=2, endbin::Union{Int,Nothing}=nothing,
                           outlier_mode::Symbol=:dropmaxmin, outlier_param::Real=0,
                           real_column::Int=3, imag_column::Int=4,
                           auto_digits::Bool=true, k_point_tolerance::Float64=1e-6, 
                           verbose::Bool=true) -> Vector{StructureFactorResult}

从k空间关联函数数据中提取多个动量点的结构因子（一次读文件，高效处理）。

# 参数
- `k_points`: 目标动量点列表，格式为 [(kx1, ky1), (kx2, ky2), ...]
- `filename`: 关联函数文件名 (默认: "spsm_k.bin")
- `filedir`: 文件目录 (默认: 当前目录)
- `startbin`: 起始bin (默认: 2)
- `endbin`: 结束bin (默认: 所有bin)
- `outlier_mode`: `:dropmaxmin` 或 `:iqrfence`（默认 `:dropmaxmin`）
- `outlier_param`: 对应模式的参数（默认 0，表示不做离群过滤）
- `real_column`: 实部列索引 (默认: 3)
- `imag_column`: 虚部列索引 (默认: 4)
- `auto_digits`: 是否自动确定有效数字 (默认: true)
- `k_point_tolerance`: k点匹配容差 (默认: 1e-6)
- `verbose`: 是否输出详细信息 (默认: true)

# 返回值
- `Vector{StructureFactorResult}`: 结构因子结果数组

# 示例
```julia
# 分析多个k点
k_list = [(π, π), (0.0, 0.0), (π, 0.0)]
results = StructureFactorAnalysis(k_list, "spsm_k.bin")
for r in results
    println("k = \$(r.k_point): \$(r.formatted_real)")
end
```
"""
function StructureFactorAnalysis(
    k_points::Vector{<:Tuple{<:Real,<:Real}},
    filename::String="spsm_k.bin",
    filedir::String=pwd();
    startbin::Int=2,
    endbin::Union{Int,Nothing}=nothing,
    outlier_mode::Symbol=:dropmaxmin,
    outlier_param::Real=0,
    real_column::Int=3,
    imag_column::Int=4,
    auto_digits::Bool=true,
    k_point_tolerance::Float64=1e-6,
    verbose::Bool=true
)
    # 类型转换
    k_points_converted = [(Float64(k[1]), Float64(k[2])) for k in k_points]
    
    # 调用核心实现
    results = _multi_k_structure_factor_analysis_core(
        k_points_converted, filename, filedir;
        startbin=startbin, endbin=endbin,
        outlier_mode=outlier_mode, outlier_param=outlier_param,
        real_column=real_column, imag_column=imag_column,
        auto_digits=auto_digits, k_point_tolerance=k_point_tolerance, verbose=verbose
    )
    
    # 打印结果
    if verbose && !isempty(results)
        for result in results
            print_structure_factor_result(result, filename)
        end
    end
    
    return results
end

"""
    StructureFactorAnalysis(k_point::Tuple{<:Real,<:Real}, 
                           filename::String="spsm_k.bin", filedir::String=pwd();
                           startbin::Int=2, endbin::Union{Int,Nothing}=nothing,
                           outlier_mode::Symbol=:dropmaxmin, outlier_param::Real=0,
                           real_column::Int=3, imag_column::Int=4,
                           auto_digits::Bool=true, k_point_tolerance::Float64=1e-6, 
                           verbose::Bool=true) -> StructureFactorResult

从k空间关联函数数据中提取单个动量点的结构因子（调用多k点版本，保持兼容性）。

# 参数
- `k_point`: 目标动量点，格式为 (kx, ky) 的元组
- 其他参数同多k点版本

# 返回值
- `StructureFactorResult`: 结构因子结果结构体

# 示例
```julia
# 分析单个k点
result = StructureFactorAnalysis((π, π), "spsm_k.bin")
println("Real part: \$(result.mean_real) ± \$(result.err_real)")
println("Imag part: \$(result.mean_imag) ± \$(result.err_imag)")
```
"""
function StructureFactorAnalysis(
    k_point::Tuple{<:Real,<:Real},
    filename::String="spsm_k.bin",
    filedir::String=pwd();
    startbin::Int=2,
    endbin::Union{Int,Nothing}=nothing,
    outlier_mode::Symbol=:dropmaxmin,
    outlier_param::Real=0,
    real_column::Int=3,
    imag_column::Int=4,
    auto_digits::Bool=true,
    k_point_tolerance::Float64=1e-6,
    verbose::Bool=true
)
    # 调用多k点版本
    results = StructureFactorAnalysis(
        [k_point], filename, filedir;
        startbin=startbin, endbin=endbin,
        outlier_mode=outlier_mode, outlier_param=outlier_param,
        real_column=real_column, imag_column=imag_column,
        auto_digits=auto_digits, k_point_tolerance=k_point_tolerance, verbose=verbose
    )
    
    # 返回第一个结果（如果存在）
    if isempty(results)
        @error "Failed to analyze k-point $k_point"
    end
    
    return results[1]
end

# ==================================================================================== #
#               AFM/CDW StructureFactor (multi-orbital with predefined columns)
# ==================================================================================== #

"""    
    AFMStructureFactor(k_point::Tuple{<:Real,<:Real}=(0.0, 0.0), 
                      filename::String="afm_sf_k.bin", filedir::String=pwd();
                      force_rebuild::Bool=false, source_file::String="spsm_k.bin", 
                      startbin::Int=2, endbin::Union{Int,Nothing}=nothing,
                      outlier_mode::Symbol=:dropmaxmin, outlier_param::Real=0,
                      auto_digits::Bool=true, k_point_tolerance::Float64=1e-6, verbose::Bool=true) -> StructureFactorResult

Calculate antiferromagnetic structure factor S_AF(L) = [spsm_k(0,A,A) + spsm_k(0,B,B) - spsm_k(0,A,B) - spsm_k(0,B,A)].

The function can use either:
- A pre-processed file (default: "afm_sf_k.bin") containing the structure factor directly
- A source file (e.g., "spsm_k.bin" or "ss_k.bin") to generate the structure factor file if it doesn't exist or `force_rebuild` is true

# Parameters
- `k_point`: Target momentum point, default (0.0, 0.0)
- `filename`: Target structure factor file name (default: "afm_sf_k.bin")
- `filedir`: Directory containing the files (default: current directory)
- `force_rebuild`: Whether to force rebuild the structure factor file even if it exists (default: false)
- `source_file`: Source file to generate structure factor if needed (default: "spsm_k.bin")
- `startbin`: Starting bin for statistics (default: 2)
- `endbin`: Ending bin for statistics (default: all bins)
- `outlier_mode`: `:dropmaxmin` or `:iqrfence` (default: `:dropmaxmin`)
- `outlier_param`: parameter for the selected outlier mode (default: 0, meaning no trimming/fence)
- `auto_digits`: Whether to automatically determine significant digits (default: true)
- `k_point_tolerance`: Tolerance for k-point matching (default: 1e-6)
- `verbose`: Whether to output detailed information (default: true)

# Returns
- `StructureFactorResult`: Structure factor result struct

# Examples
```julia
# Use existing afm_sf_k.bin file
result = AFMStructureFactor()

# Force rebuild from source file
result = AFMStructureFactor(force_rebuild=true)

# Specify custom filenames
result = AFMStructureFactor(filename="custom_afm_sf.bin", source_file="custom_afm_source.bin", force_rebuild=true)
```
"""
function AFMStructureFactor(
    k_point::Tuple{<:Real,<:Real}=(0.0, 0.0),
    filename::String="afm_sf_k.bin",
    filedir::String=pwd();
    force_rebuild::Bool=false,
    source_file::String="spsm_k.bin",
    startbin::Int=2,
    endbin::Union{Int,Nothing}=nothing,
    outlier_mode::Symbol=:dropmaxmin,
    outlier_param::Real=0,
    auto_digits::Bool=true,
    k_point_tolerance::Float64=1e-6,
    verbose::Bool=true
)
    # Ensure filename has .bin extension
    if !endswith(filename, ".bin")
        filename = filename * ".bin"
    end
    # Ensure source_file has .bin extension
    if !endswith(source_file, ".bin")
        source_file = source_file * ".bin"
    end
    target_filepath = joinpath(filedir, filename)
    
    # 如果强制重建或文件不存在，则重新生成文件
    if force_rebuild || !isfile(target_filepath)
        if verbose
            println("$(force_rebuild ? "强制重建" : "文件不存在")，正在生成 $filename...")
            println("  源文件: $source_file")
        end
        
        try
            # 如果用户指定 ss_k.bin 作为源文件但文件不存在，则尝试先生成它
            source_path = joinpath(filedir, source_file)
            if !isfile(source_path) && source_file == "ss_k.bin"
                verbose && println("  未找到 ss_k.bin，尝试从 spsm_k.bin + szsz_k.bin 生成 ss_k.bin...")
                combine_ss_components("spsm_k.bin", "szsz_k.bin", "ss_k.bin", filedir, filedir; verbose=false)
            end

            # 生成（或重建）目标结构因子文件：严格使用用户指定的 source_file -> filename
            merge_afm_sf(source_file, filename, filedir, filedir; verbose=false)

            # 检查文件是否成功生成
            if !isfile(target_filepath)
                @error "无法生成目标文件: $target_filepath"
            end
        catch e
            @error "生成文件时出错: $(sprint(showerror, e))"
        end
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
        outlier_mode=outlier_mode,
        outlier_param=outlier_param,
        auto_digits=auto_digits,
        k_point_tolerance=k_point_tolerance,
        verbose=verbose
    )
    
    return result
end

"""    
    CDWStructureFactor(k_point::Tuple{<:Real,<:Real}=(0.0, 0.0), 
                      filename::String="cdwpair_sf_k.bin", filedir::String=pwd();
                      force_rebuild::Bool=false, source_file::String="cdwpair_k.bin", 
                      startbin::Int=2, endbin::Union{Int,Nothing}=nothing,
                      outlier_mode::Symbol=:dropmaxmin, outlier_param::Real=0,
                      auto_digits::Bool=true, k_point_tolerance::Float64=1e-6, verbose::Bool=true) -> StructureFactorResult

Calculate charge density wave structure factor S_CDW(L) = [cdwpair_k(0,A,A) + cdwpair_k(0,B,B) + cdwpair_k(0,A,B) + cdwpair_k(0,B,A)].

The function can use either:
- A pre-processed file (default: "cdwpair_sf_k.bin") containing the structure factor directly
- A source file (e.g., "cdwpair_k.bin") to generate the structure factor file if it doesn't exist or `force_rebuild` is true

# Parameters
- `k_point`: Target momentum point, default (0.0, 0.0)
- `filename`: Target structure factor file name (default: "cdwpair_sf_k.bin")
- `filedir`: Directory containing the files (default: current directory)
- `force_rebuild`: Whether to force rebuild the structure factor file even if it exists (default: false)
- `source_file`: Source file to generate structure factor if needed (default: "cdwpair_k.bin")
- `startbin`: Starting bin for statistics (default: 2)
- `endbin`: Ending bin for statistics (default: all bins)
- `outlier_mode`: `:dropmaxmin` or `:iqrfence` (default: `:dropmaxmin`)
- `outlier_param`: parameter for the selected outlier mode (default: 0, meaning no trimming/fence)
- `auto_digits`: Whether to automatically determine significant digits (default: true)
- `k_point_tolerance`: Tolerance for k-point matching (default: 1e-6)
- `verbose`: Whether to output detailed information (default: true)

# Returns
- `StructureFactorResult`: Structure factor result struct

# Examples
```julia
# Use existing cdwpair_sf_k.bin file
result = CDWStructureFactor()

# Force rebuild from source file
result = CDWStructureFactor(force_rebuild=true)

# Specify custom filenames
result = CDWStructureFactor(filename="custom_cdw_sf.bin", source_file="custom_cdw.bin", force_rebuild=true)
```
"""
function CDWStructureFactor(
    k_point::Tuple{<:Real,<:Real}=(0.0, 0.0),
    filename::String="cdwpair_sf_k.bin",
    filedir::String=pwd();
    force_rebuild::Bool=false,
    source_file::String="cdwpair_k.bin",
    startbin::Int=2,
    endbin::Union{Int,Nothing}=nothing,
    outlier_mode::Symbol=:dropmaxmin,
    outlier_param::Real=0,
    auto_digits::Bool=true,
    k_point_tolerance::Float64=1e-6,
    verbose::Bool=true
)
    # Ensure filename has .bin extension
    if !endswith(filename, ".bin")
        filename = filename * ".bin"
    end
    # Ensure source_file has .bin extension
    if !endswith(source_file, ".bin")
        source_file = source_file * ".bin"
    end
    target_filepath = joinpath(filedir, filename)
    
    # 如果强制重建或文件不存在，则重新生成文件
    if force_rebuild || !isfile(target_filepath)
        if verbose
            println("$(force_rebuild ? "强制重建" : "文件不存在")，正在生成 $filename...")
            println("  源文件: $source_file")
        end
        
        try
            # 如果源文件不存在，尝试生成常见的派生源文件（cdw_k.bin / cdwpair_k.bin）
            source_path = joinpath(filedir, source_file)
            if !isfile(source_path)
                verbose && println("  未找到源文件 $(source_file)，尝试生成相关派生文件...")
                cdwpair_k_files_generation(filedir; verbose=false)
            end

            # 生成（或重建）目标结构因子文件：严格使用用户指定的 source_file -> filename
            merge_cdw_sf(source_file, filename, filedir, filedir; verbose=false)
            
            # 检查文件是否成功生成
            if !isfile(target_filepath)
                @error "无法生成目标文件: $target_filepath"
            end
        catch e
            @error "生成文件时出错: $(sprint(showerror, e))"
        end
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
        outlier_mode=outlier_mode,
        outlier_param=outlier_param,
        auto_digits=auto_digits,
        k_point_tolerance=k_point_tolerance,
        verbose=verbose
    )
    
    return result
end
