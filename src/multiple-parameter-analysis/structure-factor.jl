#=
多参数分析功能 (Multiple Parameter Analysis) 之结构因子分析
=#

# 导出多参数分析函数
export analyze_structure_factor_multi_parameter,
       analyze_AFM_structure_factor_multi_parameter,
       analyze_CDW_structure_factor_multi_parameter


# ---------------------------------------------------------------------------- #
#                         Analysis of Structure Factor                         #
# ---------------------------------------------------------------------------- #

"""
    analyze_structure_factor_multi_parameter(analyzer_function::Function,
                                           base_dir::AbstractString=pwd(); 
                                           k_point::Tuple{<:Real,<:Real}=(0.0, 0.0),
                                           filename::String="afm_sf_k.bin",
                                           source_file::String="spsm_k.bin",
                                           result_columns::Vector{Symbol}=[:S_real, :S_real_err, :S_imag, :S_imag_err],
                                           result_prefix::String="S",
                                           startbin::Int=2, 
                                           endbin::Union{Int,Nothing}=nothing, 
                                           dropmaxmin::Int=0,
                                           auto_digits::Bool=true, 
                                           tolerance::Float64=1e-6, 
                                           verbose::Bool=false,
                                           filter_options::Union{Dict, NamedTuple}=Dict(),
                                           pattern::Regex=r"^proj_fft_honeycomb") -> DataFrame

通用的多参数结构因子分析函数，可用于分析不同类型的结构因子。

# 参数
- `analyzer_function::Function`: 用于分析结构因子的函数，如 AFMStructureFactor 或 CDWStructureFactor
- `base_dir::AbstractString=pwd()`: 基础目录路径，默认为当前工作目录
- `k_point::Tuple{<:Real,<:Real}=(0.0, 0.0)`: 要分析的k点
- `filename::String="afm_sf_k.bin"`: 要分析的结构因子文件名
- `source_file::String="spsm_k.bin"`: 当 filename 不存在时，用于生成结构因子的源文件名
- `result_columns::Vector{Symbol}`: 结果列的名称，默认为 [:S_real, :S_real_err, :S_imag, :S_imag_err]
- `result_prefix::String`: 结果列名称的前缀，如 "S_AF" 或 "S_CDW"
- `startbin::Int=2`: 起始bin编号
- `endbin::Union{Int,Nothing}=nothing`: 结束bin编号
- `dropmaxmin::Int=0`: 丢弃的最大/最小值数量
- `auto_digits::Bool=true`: 是否自动确定有效数字
- `tolerance::Float64=1e-6`: k点匹配容差
- `verbose::Bool=false`: 是否显示详细信息
- `filter_options::Union{Dict, NamedTuple}=Dict()`: 目录筛选选项，可包含:prefix、:b、:U等参数
- `pattern::Regex=r"^proj_fft_honeycomb"`: 用于匹配目录名的正则表达式

# 返回值
- `DataFrame`: 包含所有参数和分析结果的DataFrame
"""
function analyze_structure_factor_multi_parameter(analyzer_function::Function,
                                           base_dir::AbstractString=pwd(); 
                                           k_point::Tuple{<:Real,<:Real}=(0.0, 0.0),
                                           filename::String="afm_sf_k.bin",
                                           source_file::String="spsm_k.bin",
                                           result_columns::Vector{Symbol}=[:S_real, :S_real_err, :S_imag, :S_imag_err],
                                           result_prefix::String="S",
                                           startbin::Int=2, 
                                           endbin::Union{Int,Nothing}=nothing, 
                                           dropmaxmin::Int=0,
                                           auto_digits::Bool=true, 
                                           tolerance::Float64=1e-6, 
                                           verbose::Bool=false,
                                           filter_options::Union{Dict, NamedTuple}=Dict(),
                                           pattern::Regex=r"^proj_fft_honeycomb")
    # 扫描参数目录
    param_dirs = scan_parameter_directories(base_dir; filter_options=filter_options, pattern=pattern)
    
    if isempty(param_dirs)
        @warn "未在 $base_dir 中找到参数目录"
        return DataFrame()
    end
    
    # 创建基本列的DataFrame
    result_df = DataFrame(
        # 参数
        :directory => String[],
        :prefix => String[],
        :b => Float64[],
        :U => Float64[],
        :L => Int[],
        :dtau => Float64[],
        # Gutzwiller参数
        :gw => Union{Float64, Missing}[],
        :lprojgw => Union{Bool, Missing, Nothing}[]  # 允许 nothing 值
    )
    
    # 添加动态结果列
    for col in result_columns
        result_col_name = Symbol(result_prefix * "_" * string(col))
        result_df[!, result_col_name] = Float64[]
    end
    
    # Process each directory
    total_dirs = length(param_dirs)
    println("找到 $total_dirs 个参数目录进行分析...")
    
    for (i, dir) in enumerate(param_dirs)
        dirname = basename(dir)
        
        # Extract parameters
        params = extract_parameters_from_dirname(dirname)
        
        # 检查结构因子文件是否存在
        filepath = joinpath(dir, filename)
        
        # 如果文件不存在，直接跳过此目录
        if !isfile(filepath)
            if verbose
                println("($i/$total_dirs) 跳过 $dirname: 文件 $filename 不存在")
            end
            continue
        end
        
        try
            if verbose
                println("($i/$total_dirs) 分析 $dirname...")
            end
            
            # 执行分析
            analysis_result = analyzer_function(
                k_point, filename, dir;
                source_file=source_file,
                startbin=startbin, 
                endbin=endbin, 
                dropmaxmin=dropmaxmin,
                auto_digits=auto_digits, 
                tolerance=tolerance,
                verbose=false  # 自己处理输出
            )
            
            # 创建基本参数字典
            row_dict = Dict(
                :directory => dirname,
                :prefix => params[:prefix],
                :b => params[:b],
                :U => params[:U],
                :L => params[:L],
                :dtau => params[:dtau],
                # Gutzwiller parameters - use missing for non-existent values
                :gw => get(params, :gw, missing),
                :lprojgw => get(params, :lprojgw, nothing)  # 保持原始的 nothing 值
            )
            
            # 添加分析结果
            for col in result_columns
                result_key = Symbol(result_prefix * "_" * string(col))
                row_dict[result_key] = getproperty(analysis_result, col)
            end
            
            # 添加行到DataFrame
            push!(result_df, row_dict)
            
            if verbose
                println("($i/$total_dirs) 完成分析 $dirname")
                println("   $(result_prefix)($(k_point)) = $(getproperty(analysis_result, result_columns[1])) + $(getproperty(analysis_result, result_columns[3]))i")
            end
        catch e
            println("($i/$total_dirs) 分析 $dirname 时出错: $e")
        end
    end
    
    # Check if we have results
    if nrow(result_df) == 0
        @warn "没有成功的分析结果"
        return DataFrame()
    end
    
    # Sort by (b, U, L, dtau)
    sort!(result_df, [:b, :U, :L, :dtau])
    
    # Print column order (for debugging)
    if verbose
        println("DataFrame列顺序: $(names(result_df))")
    end
    
    return result_df
end

"""
    analyze_AFM_structure_factor_multi_parameter(base_dir::AbstractString=pwd(); 
                                               k_point::Tuple{<:Real,<:Real}=(0.0, 0.0),
                                               filename::String="afm_sf_k.bin",
                                               source_file::String="spsm_k.bin",
                                               startbin::Int=2, 
                                               endbin::Union{Int,Nothing}=nothing, 
                                               dropmaxmin::Int=0,
                                               auto_digits::Bool=true, 
                                               tolerance::Float64=1e-6, 
                                               verbose::Bool=false,
                                               filter_options::Union{Dict, NamedTuple}=Dict(),
                                               pattern::Regex=r"^proj_fft_honeycomb") -> DataFrame

对多个参数目录执行反铁磁结构因子分析，并将结果整合到一个DataFrame中。

# 参数
- `base_dir::AbstractString=pwd()`: 基础目录路径，默认为当前工作目录
- `k_point::Tuple{<:Real,<:Real}=(0.0, 0.0)`: 要分析的k点
- `filename::String="afm_sf_k.bin"`: 要分析的结构因子文件名
- `source_file::String="spsm_k.bin"`: 当 filename 不存在时，用于生成结构因子的源文件名
- `startbin::Int=2`: 起始bin编号
- `endbin::Union{Int,Nothing}=nothing`: 结束bin编号
- `dropmaxmin::Int=0`: 丢弃的最大/最小值数量
- `auto_digits::Bool=true`: 是否自动确定有效数字
- `tolerance::Float64=1e-6`: k点匹配容差
- `verbose::Bool=false`: 是否显示详细信息
- `filter_options::Union{Dict, NamedTuple}=Dict()`: 目录筛选选项，可包含:prefix、:b、:U等参数
- `pattern::Regex=r"^proj_fft_honeycomb"`: 用于匹配目录名的正则表达式

# 返回值
- `DataFrame`: 包含所有参数和分析结果的DataFrame
"""
function analyze_AFM_structure_factor_multi_parameter(base_dir::AbstractString=pwd(); 
                                                   k_point::Tuple{<:Real,<:Real}=(0.0, 0.0),
                                                   filename::String="afm_sf_k.bin",
                                                   source_file::String="spsm_k.bin",
                                                   startbin::Int=2, 
                                                   endbin::Union{Int,Nothing}=nothing, 
                                                   dropmaxmin::Int=0,
                                                   auto_digits::Bool=true, 
                                                   tolerance::Float64=1e-6, 
                                                   verbose::Bool=false,
                                                   filter_options::Union{Dict, NamedTuple}=Dict(),
                                                   pattern::Regex=r"^proj_fft_honeycomb")
    # 调用通用的结构因子分析函数，传递 AFMStructureFactor 函数和相应的参数
    return analyze_structure_factor_multi_parameter(
        AFMStructureFactor,
        base_dir;
        k_point=k_point,
        filename=filename,
        source_file=source_file,
        result_columns=[:mean_real, :err_real, :mean_imag, :err_imag],
        result_prefix="S_AF",
        startbin=startbin,
        endbin=endbin,
        dropmaxmin=dropmaxmin,
        auto_digits=auto_digits,
        tolerance=tolerance,
        verbose=verbose,
        filter_options=filter_options,
        pattern=pattern
    )
end


"""
    analyze_CDW_structure_factor_multi_parameter(base_dir::AbstractString=pwd(); 
                                               k_point::Tuple{<:Real,<:Real}=(0.0, 0.0),
                                               filename::String="cdwpair_sf_k.bin",
                                               source_file::String="cdwpair_k.bin",
                                               startbin::Int=2, 
                                               endbin::Union{Int,Nothing}=nothing, 
                                               dropmaxmin::Int=0,
                                               auto_digits::Bool=true, 
                                               tolerance::Float64=1e-6, 
                                               verbose::Bool=false,
                                               filter_options::Union{Dict, NamedTuple}=Dict(),
                                               pattern::Regex=r"^proj_fft_honeycomb") -> DataFrame

对多个参数目录执行电荷密度波结构因子分析，并将结果整合到一个DataFrame中。

# 参数
- `base_dir::AbstractString=pwd()`: 基础目录路径，默认为当前工作目录
- `k_point::Tuple{<:Real,<:Real}=(0.0, 0.0)`: 要分析的k点
- `filename::String="cdwpair_sf_k.bin"`: 要分析的结构因子文件名
- `source_file::String="cdwpair_k.bin"`: 当 filename 不存在时，用于生成结构因子的源文件名
- `startbin::Int=2`: 起始bin编号
- `endbin::Union{Int,Nothing}=nothing`: 结束bin编号
- `dropmaxmin::Int=0`: 丢弃的最大/最小值数量
- `auto_digits::Bool=true`: 是否自动确定有效数字
- `tolerance::Float64=1e-6`: k点匹配容差
- `verbose::Bool=false`: 是否显示详细信息
- `filter_options::Union{Dict, NamedTuple}=Dict()`: 目录筛选选项，可包含:prefix、:b、:U等参数
- `pattern::Regex=r"^proj_fft_honeycomb"`: 用于匹配目录名的正则表达式

# 返回值
- `DataFrame`: 包含所有参数和分析结果的DataFrame
"""
function analyze_CDW_structure_factor_multi_parameter(base_dir::AbstractString=pwd(); 
                                                   k_point::Tuple{<:Real,<:Real}=(0.0, 0.0),
                                                   filename::String="cdwpair_sf_k.bin",
                                                   source_file::String="cdwpair_k.bin",
                                                   startbin::Int=2, 
                                                   endbin::Union{Int,Nothing}=nothing, 
                                                   dropmaxmin::Int=0,
                                                   auto_digits::Bool=true, 
                                                   tolerance::Float64=1e-6, 
                                                   verbose::Bool=false,
                                                   filter_options::Union{Dict, NamedTuple}=Dict(),
                                                   pattern::Regex=r"^proj_fft_honeycomb")
    # 调用通用的结构因子分析函数，传递 CDWStructureFactor 函数和相应的参数
    return analyze_structure_factor_multi_parameter(
        CDWStructureFactor,
        base_dir;
        k_point=k_point,
        filename=filename,
        source_file=source_file,
        result_columns=[:mean_real, :err_real, :mean_imag, :err_imag],
        result_prefix="S_CDW",
        startbin=startbin,
        endbin=endbin,
        dropmaxmin=dropmaxmin,
        auto_digits=auto_digits,
        tolerance=tolerance,
        verbose=verbose,
        filter_options=filter_options,
        pattern=pattern
    )
end