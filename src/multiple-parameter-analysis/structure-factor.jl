#=
多参数分析功能 (Multiple Parameter Analysis) 之结构因子分析
=#

# 导出多参数分析函数
export analyze_structure_factor_multi_parameter,
       analyze_AFM_structure_factor_multi_parameter,
       analyze_CDW_structure_factor_multi_parameter


# ==================================================================================== #
#                              私有辅助函数
# ==================================================================================== #

"""
将参数向量转换为字典（包含 :prefix）
"""
function _params_dict(prefix::String, params_vector)::Dict{Symbol,Any}
    params = Dict{Symbol,Any}()
    params[:prefix] = prefix
    for (key, value, _) in params_vector
        params[Symbol(key)] = value
    end
    return params
end

"""
初始化结果 DataFrame 的列结构
"""
function _init_result_dataframe(
    first_params::Dict{Symbol,Any},
    result_columns::Vector{Symbol},
    result_prefix::String;
    include_k_cols::Bool=false
)::DataFrame
    result_df = DataFrame()
    
    # 参数列
    for param_name in keys(first_params)
        result_df[!, param_name] = typeof(first_params[param_name])[]
    end
    
    # k 点坐标列（仅多 k 版本）
    if include_k_cols
        result_df[!, :kx] = Float64[]
        result_df[!, :ky] = Float64[]
    end
    
    # 结果列
    for col in result_columns
        result_df[!, Symbol("$(result_prefix)_$(col)")] = Float64[]
    end
    
    # 格式化列
    result_df[!, Symbol("$(result_prefix)_formatted")] = String[]
    
    return result_df
end

"""
创建结果行并填充数据
"""
function _create_result_row(
    params::Dict{Symbol,Any},
    result,  # StructureFactorResult or NamedTuple
    result_columns::Vector{Symbol},
    result_prefix::String;
    include_k_point::Bool=false
)::Dict{Symbol,Any}
    row = Dict{Symbol, Any}()
    
    # 添加参数列
    for (param_name, param_value) in params
        row[param_name] = param_value
    end
    
    # 添加 k 点坐标（仅多 k 版本）
    if include_k_point
        row[:kx] = result.k_point[1]
        row[:ky] = result.k_point[2]
    end
    
    # 添加结果列
    for col in result_columns
        result_key = Symbol(result_prefix * "_" * string(col))
        row[result_key] = getproperty(result, col)
    end
    
    # 添加格式化列
    row[Symbol("$(result_prefix)_formatted")] = result.formatted_real
    
    return row
end


# ==================================================================================== #
#                    通用多参数结构因子分析 - 单 k 点版本
# ==================================================================================== #

"""
    analyze_structure_factor_multi_parameter(analyzer_function::Function,
                                          k_point::Tuple{<:Real,<:Real},
                                          base_dir::AbstractString=pwd(); 
                                          filename::String="afm_sf_k.bin",
                                          force_rebuild::Bool=false,
                                          source_file::String="spsm_k.bin",
                                          result_columns::Vector{Symbol}=[:mean_real, :err_real, :mean_imag, :err_imag],
                                          result_prefix::String="S",
                                          startbin::Int=2, 
                                          endbin::Union{Int,Nothing}=nothing, 
                                          outlier_mode::Symbol=:dropmaxmin,
                                          outlier_param::Real=0,
                                          auto_digits::Bool=true, 
                                          k_point_tolerance::Float64=1e-6, 
                                          verbose::Bool=false,
                                          filter_options::Union{Dict, NamedTuple}=Dict()) -> DataFrame

通用的多参数结构因子分析函数（单 k 点版本），可用于分析不同类型的结构因子。

返回宽表，每个参数目录一行，不包含 k 点信息列。

# 参数
- `analyzer_function::Function`: 用于分析结构因子的函数，如 AFMStructureFactor 或 CDWStructureFactor
- `k_point::Tuple{<:Real,<:Real}`: 要分析的k点（位置参数）
- `base_dir::AbstractString=pwd()`: 基础目录路径，默认为当前工作目录
- `filename::String="afm_sf_k.bin"`: 要分析的结构因子文件名
- `force_rebuild::Bool=false`: 是否强制重新构建结构因子文件
- `source_file::String="spsm_k.bin"`: 当 filename 不存在时，用于生成结构因子的源文件名
- `result_columns::Vector{Symbol}`: 结果列的名称，默认为 [:mean_real, :err_real, :mean_imag, :err_imag]
- `result_prefix::String`: 结果列名称的前缀，如 "S_AF" 或 "S_CDW"
- `startbin::Int=2`: 起始bin编号
- `endbin::Union{Int,Nothing}=nothing`: 结束bin编号
- `outlier_mode::Symbol=:dropmaxmin`: 离群过滤模式（`:dropmaxmin` / `:iqrfence`）
- `outlier_param::Real=0`: 对应模式的参数（0 表示不做过滤）
- `auto_digits::Bool=true`: 是否自动确定有效数字
- `k_point_tolerance::Float64=1e-6`: k点匹配容差
- `verbose::Bool=false`: 是否显示详细信息
- `filter_options::Union{Dict, NamedTuple}=Dict()`: 目录筛选选项，可包含:prefix、:b、:U等参数

# 返回值
- `DataFrame`: 包含所有参数和分析结果的DataFrame（宽表，每个参数目录一行）
"""
function analyze_structure_factor_multi_parameter(analyzer_function::Function,
                                           k_point::Tuple{<:Real,<:Real},
                                           base_dir::AbstractString=pwd();
                                           filename::String="afm_sf_k.bin",
                                           force_rebuild::Bool=false,
                                           source_file::String="spsm_k.bin",
                                           result_columns::Vector{Symbol}=[:mean_real, :err_real, :mean_imag, :err_imag],
                                           result_prefix::String="S",
                                           startbin::Int=2, 
                                           endbin::Union{Int,Nothing}=nothing, 
                                           outlier_mode::Symbol=:dropmaxmin,
                                           outlier_param::Real=0,
                                           auto_digits::Bool=true, 
                                           k_point_tolerance::Float64=1e-6, 
                                           verbose::Bool=false,
                                           filter_options::Union{Dict, NamedTuple}=Dict())
    # Scan parameter directories
    param_dirs_with_params = scan_parameter_directories(base_dir; filter_options=filter_options, return_params=true)
    
    if isempty(param_dirs_with_params)
        @warn "No parameter directories found in $base_dir"
        return DataFrame()
    end
    
    # Initialize result DataFrame
    first_dir, first_prefix, first_params_vector = first(param_dirs_with_params)
    first_params = _params_dict(first_prefix, first_params_vector)
    result_df = _init_result_dataframe(first_params, result_columns, result_prefix)

    # Process each directory
    total_dirs = length(param_dirs_with_params)
    println("Found $total_dirs parameter directories to analyze...")
    
    for (i, (dir, prefix, params_vector)) in enumerate(param_dirs_with_params)
        dirname = basename(dir)
        params = _params_dict(prefix, params_vector)
        
        # Check if the structure factor file exists
        filepath = joinpath(dir, filename)
        if !isfile(filepath) && !force_rebuild
            if verbose
                println("($i/$total_dirs) Skipping $dirname: file $filename does not exist")
            end
            continue
        end
        
        try
            if verbose
                println("($i/$total_dirs) Analyzing $dirname...")
            end
            
            analysis_result = analyzer_function(
                k_point, filename, dir;
                force_rebuild=force_rebuild,
                source_file=source_file,
                startbin=startbin, 
                endbin=endbin, 
                outlier_mode=outlier_mode,
                outlier_param=outlier_param,
                auto_digits=auto_digits, 
                k_point_tolerance=k_point_tolerance,
                verbose=false
            )
            
            # Create and add result row
            row = _create_result_row(params, analysis_result, result_columns, result_prefix)
            push!(result_df, row)
            
            if verbose
                println("($i/$total_dirs) Analysis of $dirname completed")
                println("   $(result_prefix)($(k_point)) = $(getproperty(analysis_result, result_columns[1])) + $(getproperty(analysis_result, result_columns[3]))i")
            end
        catch e
            println("($i/$total_dirs) Error analyzing $dirname: $e")
        end
    end
    
    # Check if we have results
    if nrow(result_df) == 0
        @warn "No successful analysis results"
        return DataFrame()
    end
    
    # Sort by parameters
    param_cols = sort(collect(keys(first_params)))
    sort!(result_df, param_cols)
    
    return result_df
end


# ==================================================================================== #
#                    通用多参数结构因子分析 - 多 k 点版本
# ==================================================================================== #

"""
    analyze_structure_factor_multi_parameter(analyzer_function::Function,
                                           k_points::Vector{<:Tuple{<:Real,<:Real}},
                                           base_dir::AbstractString=pwd(); 
                                           filename::String="afm_sf_k.bin",
                                           force_rebuild::Bool=false,
                                           source_file::String="spsm_k.bin",
                                           result_columns::Vector{Symbol}=[:mean_real, :err_real, :mean_imag, :err_imag],
                                           result_prefix::String="S",
                                           startbin::Int=2, 
                                           endbin::Union{Int,Nothing}=nothing, 
                                           outlier_mode::Symbol=:dropmaxmin,
                                           outlier_param::Real=0,
                                           auto_digits::Bool=true, 
                                           k_point_tolerance::Float64=1e-6, 
                                           verbose::Bool=false,
                                           filter_options::Union{Dict, NamedTuple}=Dict()) -> DataFrame

通用的多参数结构因子分析函数（多 k 点版本），可用于分析不同类型的结构因子。

返回长表，每个参数目录 × 每个 k 点一行，包含 :kx 和 :ky 列。

# 参数
- `analyzer_function::Function`: 用于分析结构因子的函数，必须支持多 k 点输入
- `k_points::Vector{<:Tuple{<:Real,<:Real}}`: 要分析的k点列表（位置参数）
- `base_dir::AbstractString=pwd()`: 基础目录路径，默认为当前工作目录
- 其他参数同单 k 点版本

# 返回值
- `DataFrame`: 包含所有参数、k点坐标和分析结果的DataFrame（长表）
  新增列：:kx, :ky
"""
function analyze_structure_factor_multi_parameter(analyzer_function::Function,
                                           k_points::Vector{<:Tuple{<:Real,<:Real}},
                                           base_dir::AbstractString=pwd();
                                           filename::String="afm_sf_k.bin",
                                           force_rebuild::Bool=false,
                                           source_file::String="spsm_k.bin",
                                           result_columns::Vector{Symbol}=[:mean_real, :err_real, :mean_imag, :err_imag],
                                           result_prefix::String="S",
                                           startbin::Int=2, 
                                           endbin::Union{Int,Nothing}=nothing, 
                                           outlier_mode::Symbol=:dropmaxmin,
                                           outlier_param::Real=0,
                                           auto_digits::Bool=true, 
                                           k_point_tolerance::Float64=1e-6, 
                                           verbose::Bool=false,
                                           filter_options::Union{Dict, NamedTuple}=Dict())
    # Scan parameter directories
    param_dirs_with_params = scan_parameter_directories(base_dir; filter_options=filter_options, return_params=true)
    
    if isempty(param_dirs_with_params)
        @warn "No parameter directories found in $base_dir"
        return DataFrame()
    end
    
    # Initialize result DataFrame (with k-point columns)
    first_dir, first_prefix, first_params_vector = first(param_dirs_with_params)
    first_params = _params_dict(first_prefix, first_params_vector)
    result_df = _init_result_dataframe(first_params, result_columns, result_prefix; include_k_cols=true)

    # Process each directory
    total_dirs = length(param_dirs_with_params)
    println("Found $total_dirs parameter directories to analyze (multi-k mode)...")
    
    for (i, (dir, prefix, params_vector)) in enumerate(param_dirs_with_params)
        dirname = basename(dir)
        params = _params_dict(prefix, params_vector)
        
        # Check if the structure factor file exists
        filepath = joinpath(dir, filename)
        if !isfile(filepath) && !force_rebuild
            if verbose
                println("($i/$total_dirs) Skipping $dirname: file $filename does not exist")
            end
            continue
        end
        
        try
            if verbose
                println("($i/$total_dirs) Analyzing $dirname ($(length(k_points)) k-points)...")
            end
            
            # Call analyzer function with multi-k support
            analysis_results = analyzer_function(
                k_points, filename, dir;
                force_rebuild=force_rebuild,
                source_file=source_file,
                startbin=startbin, 
                endbin=endbin, 
                outlier_mode=outlier_mode,
                outlier_param=outlier_param,
                auto_digits=auto_digits, 
                k_point_tolerance=k_point_tolerance,
                verbose=false
            )
            
            # Process each k-point result
            for result in analysis_results
                row = _create_result_row(params, result, result_columns, result_prefix; include_k_point=true)
                push!(result_df, row)
            end
            
            if verbose
                println("($i/$total_dirs) Analysis of $dirname completed: $(length(analysis_results)) k-points")
            end
        catch e
            println("($i/$total_dirs) Error analyzing $dirname: $e")
            if verbose
                println("   Stack trace: ")
                for (exc, bt) in Base.catch_stack()
                    showerror(stdout, exc, bt)
                    println()
                end
            end
        end
    end
    
    # Check if we have results
    if nrow(result_df) == 0
        @warn "No successful analysis results"
        return DataFrame()
    end
    
    # Sort by parameters + kx, ky
    param_cols = sort(collect(keys(first_params)))
    sort_cols = vcat(param_cols, [:kx, :ky])
    sort!(result_df, sort_cols)
    
    return result_df
end


# ==================================================================================== #
#                    专用多参数分析函数（仅单 k 点）
# ==================================================================================== #

"""
    analyze_AFM_structure_factor_multi_parameter(base_dir::AbstractString=pwd(); 
                                               k_point::Tuple{<:Real,<:Real}=(0.0, 0.0),
                                               filename::String="afm_sf_k.bin",
                                               force_rebuild::Bool=false,
                                               source_file::String="spsm_k.bin",
                                               startbin::Int=2, 
                                               endbin::Union{Int,Nothing}=nothing, 
                                               outlier_mode::Symbol=:dropmaxmin,
                                               outlier_param::Real=0,
                                               auto_digits::Bool=true, 
                                               k_point_tolerance::Float64=1e-6, 
                                               verbose::Bool=false,
                                               filter_options::Union{Dict, NamedTuple}=Dict()) -> DataFrame

对多个参数目录执行反铁磁结构因子分析，并将结果整合到一个DataFrame中。

仅支持单 k 点分析。

# 参数
- `base_dir::AbstractString=pwd()`: 基础目录路径，默认为当前工作目录
- `k_point::Tuple{<:Real,<:Real}=(0.0, 0.0)`: 要分析的k点
- `filename::String="afm_sf_k.bin"`: 要分析的结构因子文件名
- `force_rebuild::Bool=false`: 是否强制重新构建结构因子文件
- `source_file::String="spsm_k.bin"`: 当 filename 不存在时，用于生成结构因子的源文件名
- `startbin::Int=2`: 起始bin编号
- `endbin::Union{Int,Nothing}=nothing`: 结束bin编号
- `outlier_mode::Symbol=:dropmaxmin`: 离群过滤模式（`:dropmaxmin` / `:iqrfence`）
- `outlier_param::Real=0`: 对应模式的参数（0 表示不做过滤）
- `auto_digits::Bool=true`: 是否自动确定有效数字
- `k_point_tolerance::Float64=1e-6`: k点匹配容差
- `verbose::Bool=false`: 是否显示详细信息
- `filter_options::Union{Dict, NamedTuple}=Dict()`: 目录筛选选项，可包含:prefix、:b、:U等参数

# 返回值
- `DataFrame`: 包含所有参数和分析结果的DataFrame
"""
function analyze_AFM_structure_factor_multi_parameter(base_dir::AbstractString=pwd(); 
                                                   k_point::Tuple{<:Real,<:Real}=(0.0, 0.0),
                                                   filename::String="afm_sf_k.bin",
                                                   force_rebuild::Bool=false,
                                                   source_file::String="spsm_k.bin",
                                                   startbin::Int=2, 
                                                   endbin::Union{Int,Nothing}=nothing, 
                                                   outlier_mode::Symbol=:dropmaxmin,
                                                   outlier_param::Real=0,
                                                   auto_digits::Bool=true, 
                                                   k_point_tolerance::Float64=1e-6, 
                                                   verbose::Bool=false,
                                                   filter_options::Union{Dict, NamedTuple}=Dict())
    # 调用通用的结构因子分析函数，传递 AFMStructureFactor 函数和相应的参数
    return analyze_structure_factor_multi_parameter(
        AFMStructureFactor,
        k_point,          # k_point is now a positional argument
        base_dir;
        filename=filename,
        force_rebuild=force_rebuild,
        source_file=source_file,
        result_columns=[:mean_real, :err_real, :mean_imag, :err_imag],
        result_prefix="S_AF",
        startbin=startbin,
        endbin=endbin,
        outlier_mode=outlier_mode,
        outlier_param=outlier_param,
        auto_digits=auto_digits,
        k_point_tolerance=k_point_tolerance,
        verbose=verbose,
        filter_options=filter_options
    )
end


"""
    analyze_CDW_structure_factor_multi_parameter(base_dir::AbstractString=pwd(); 
                                               k_point::Tuple{<:Real,<:Real}=(0.0, 0.0),
                                               filename::String="cdwpair_sf_k.bin",
                                               force_rebuild::Bool=false,
                                               source_file::String="cdwpair_k.bin",
                                               startbin::Int=2, 
                                               endbin::Union{Int,Nothing}=nothing, 
                                               outlier_mode::Symbol=:dropmaxmin,
                                               outlier_param::Real=0,
                                               auto_digits::Bool=true, 
                                               k_point_tolerance::Float64=1e-6, 
                                               verbose::Bool=false,
                                               filter_options::Union{Dict, NamedTuple}=Dict()) -> DataFrame

对多个参数目录执行电荷密度波结构因子分析，并将结果整合到一个DataFrame中。

仅支持单 k 点分析。

# 参数
- `base_dir::AbstractString=pwd()`: 基础目录路径，默认为当前工作目录
- `k_point::Tuple{<:Real,<:Real}=(0.0, 0.0)`: 要分析的k点
- `filename::String="cdwpair_sf_k.bin"`: 要分析的结构因子文件名
- `force_rebuild::Bool=false`: 是否强制重新构建结构因子文件
- `source_file::String="cdwpair_k.bin"`: 当 filename 不存在时，用于生成结构因子的源文件名
- `startbin::Int=2`: 起始bin编号
- `endbin::Union{Int,Nothing}=nothing`: 结束bin编号
- `outlier_mode::Symbol=:dropmaxmin`: 离群过滤模式（`:dropmaxmin` / `:iqrfence`）
- `outlier_param::Real=0`: 对应模式的参数（0 表示不做过滤）
- `auto_digits::Bool=true`: 是否自动确定有效数字
- `k_point_tolerance::Float64=1e-6`: k点匹配容差
- `verbose::Bool=false`: 是否显示详细信息
- `filter_options::Union{Dict, NamedTuple}=Dict()`: 目录筛选选项，可包含:prefix、:b、:U等参数

# 返回值
- `DataFrame`: 包含所有参数和分析结果的DataFrame
"""
function analyze_CDW_structure_factor_multi_parameter(base_dir::AbstractString=pwd(); 
                                                   k_point::Tuple{<:Real,<:Real}=(0.0, 0.0),
                                                   filename::String="cdwpair_sf_k.bin",
                                                   force_rebuild::Bool=false,
                                                   source_file::String="cdwpair_k.bin",
                                                   startbin::Int=2, 
                                                   endbin::Union{Int,Nothing}=nothing, 
                                                   outlier_mode::Symbol=:dropmaxmin,
                                                   outlier_param::Real=0,
                                                   auto_digits::Bool=true, 
                                                   k_point_tolerance::Float64=1e-6, 
                                                   verbose::Bool=false,
                                                   filter_options::Union{Dict, NamedTuple}=Dict())
    # 调用通用的结构因子分析函数，传递 CDWStructureFactor 函数和相应的参数
    return analyze_structure_factor_multi_parameter(
        CDWStructureFactor,
        k_point,          # k_point is now a positional argument
        base_dir;
        filename=filename,
        force_rebuild=force_rebuild,
        source_file=source_file,
        result_columns=[:mean_real, :err_real, :mean_imag, :err_imag],
        result_prefix="S_CDW",
        startbin=startbin,
        endbin=endbin,
        outlier_mode=outlier_mode,
        outlier_param=outlier_param,
        auto_digits=auto_digits,
        k_point_tolerance=k_point_tolerance,
        verbose=verbose,
        filter_options=filter_options
    )
end
