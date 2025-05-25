#=
多参数分析功能 (Multiple Parameter Analysis)

此文件提供了在包含多个参数文件夹的目录中运行分析的功能。
主要功能包括:
1. 从目录名中提取参数
2. 扫描符合条件的参数目录
3. 对多个参数目录执行分析并整合结果
=#

# 导出多参数分析函数
export extract_parameters_from_dirname,
       scan_parameter_directories,
       save_analysis_results
export analyze_structure_factor_multi_parameter,
       analyze_AFM_structure_factor_multi_parameter,
       analyze_CDW_structure_factor_multi_parameter

# ---------------------------------------------------------------------------- #
#               Helper functions for multiple parameter analysis               #
# ---------------------------------------------------------------------------- #

"""
    extract_parameters_from_dirname(dirname::AbstractString) -> Dict{Symbol, Any}

从目录名中提取参数值，包括前缀类型、温度参数b、相互作用U、晶格尺寸L和时间步长dtau。

# 参数
- `dirname::AbstractString`: 目录名，例如 "proj_fft_honeycomb_exact.b8.000.U4.00.L9.dtau0.05"

# 返回值
- `Dict{Symbol, Any}`: 包含提取参数的字典
"""
function extract_parameters_from_dirname(dirname::AbstractString)
    params = Dict{Symbol, Any}()
    
    # 使用正则表达式提取前缀
    prefix_match = match(r"^(proj_fft_honeycomb(?:_[a-zA-Z]+)?)", dirname)
    params[:prefix] = prefix_match !== nothing ? prefix_match.captures[1] : "unknown"
    
    # 提取数值参数
    # 温度参数b
    b_match = match(r"b(\d+\.\d+)", dirname)
    params[:b] = b_match !== nothing ? parse(Float64, b_match.captures[1]) : NaN
    
    # 相互作用强度U
    u_match = match(r"U(-?\d+\.\d+)", dirname)
    params[:U] = u_match !== nothing ? parse(Float64, u_match.captures[1]) : NaN
    
    # Gutzwiller参数gw
    gw_match = match(r"gw(-?\d+\.\d+)", dirname)
    params[:gw] = gw_match !== nothing ? parse(Float64, gw_match.captures[1]) : NaN
    
    # 晶格尺寸L
    l_match = match(r"L(\d+)", dirname)
    params[:L] = l_match !== nothing ? parse(Int, l_match.captures[1]) : -1
    
    # 时间步长dtau
    dtau_match = match(r"dtau(\d+\.\d+)", dirname)
    params[:dtau] = dtau_match !== nothing ? parse(Float64, dtau_match.captures[1]) : NaN
    
    # 求解投影参数lprojgw
    lprojgw_match = match(r"lprojgw([TF])", dirname)
    params[:lprojgw] = lprojgw_match !== nothing ? (lprojgw_match.captures[1] == "T") : nothing
    
    # CDW参数M_CDW
    mcdw_match = match(r"M_CDW(-?\d+\.\d+)", dirname)
    params[:M_CDW] = mcdw_match !== nothing ? parse(Float64, mcdw_match.captures[1]) : NaN
    
    return params
end

"""
    scan_parameter_directories(base_dir::AbstractString=pwd(); 
                              filter_options::Union{Dict, NamedTuple}=Dict(), 
                              pattern::Regex=r"^proj_fft_honeycomb") -> Vector{String}

扫描指定目录下所有符合筛选条件的子目录。

# 参数
- `base_dir::AbstractString=pwd()`: 基础目录路径，默认为当前工作目录
- `filter_options::Union{Dict, NamedTuple}=Dict()`: 筛选选项，可以包含以下键：
  - `:prefix`: 目录前缀（字符串或字符串数组）
  - `:b`, `:U`, `:L`, `:dtau`, `:gw`: 参数范围（可以是单个值、范围或值数组）
  - `:lprojgw`: 布尔值，是否使用lprojgw
- `pattern::Regex=r"^proj_fft_honeycomb"`: 用于额外匹配目录名的正则表达式

# 返回值
- `Vector{String}`: 符合条件的子目录路径列表

# 示例
```julia
# 筛选前缀为 "proj_fft_honeycomb_exact" 的目录
dirs = scan_parameter_directories(filter_options=(prefix="proj_fft_honeycomb_exact",))

# 筛选多个条件
dirs = scan_parameter_directories(filter_options=Dict(
    :prefix => ["proj_fft_honeycomb_exact", "proj_fft_honeycomb"],
    :U => 4.0,
    :L => [6, 9]
))
```
"""
function scan_parameter_directories(base_dir::AbstractString=pwd(); 
                                   filter_options::Union{Dict, NamedTuple}=Dict(), 
                                   pattern::Regex=r"^proj_fft_honeycomb")
    result_dirs = String[]
    
    # 列出基础目录下的所有条目
    entries = readdir(base_dir, join=true)
    
    # 转换NamedTuple为Dict以统一处理
    filter_dict = filter_options isa NamedTuple ? Dict(pairs(filter_options)) : filter_options
    
    # 过滤出目录
    for entry in entries
        if isdir(entry)
            dirname = basename(entry)
            
            # 应用正则表达式筛选（如果提供）
            if !isempty(pattern.pattern) && match(pattern, dirname) === nothing
                continue
            end
            
            # 如果没有筛选选项，则直接添加目录
            if isempty(filter_dict)
                push!(result_dirs, entry)
                continue
            end
            
            # 提取参数
            params = extract_parameters_from_dirname(dirname)
            
            # 检查是否满足所有筛选条件
            matches_all_filters = true
            
            for (key, filter_value) in filter_dict
                # 获取参数值（如果不存在则为nothing）
                param_value = get(params, key, nothing)
                
                # 如果参数不存在，则不满足条件
                if param_value === nothing
                    matches_all_filters = false
                    break
                end
                
                # 根据筛选值类型进行匹配
                if filter_value isa AbstractArray
                    # 数组：检查参数值是否在数组中
                    if !(param_value in filter_value)
                        matches_all_filters = false
                        break
                    end
                elseif filter_value isa Tuple && length(filter_value) == 2
                    # 范围：检查参数值是否在范围内
                    if !(filter_value[1] <= param_value <= filter_value[2])
                        matches_all_filters = false
                        break
                    end
                else
                    # 单个值：检查参数值是否相等
                    if param_value != filter_value
                        matches_all_filters = false
                        break
                    end
                end
            end
            
            # 如果满足所有条件，则添加目录
            if matches_all_filters
                push!(result_dirs, entry)
            end
        end
    end
    
    return result_dirs
end

"""
    save_analysis_results(df::DataFrame, filename::AbstractString="analysis_results.csv")

将分析结果保存到CSV文件。

# 参数
- `df::DataFrame`: 要保存的DataFrame
- `filename::AbstractString="analysis_results.csv"`: 输出文件名
"""
function save_analysis_results(df::DataFrame, filename::AbstractString="analysis_results.csv")
    if isempty(df)
        @warn "DataFrame为空，无内容可保存"
        return
    end
    
    # 定义期望的列顺序
    desired_columns = [
        :directory, :prefix, :b, :U, :L, :dtau, :gw, :lprojgw,  # 参数
        :S_AF_real, :S_AF_real_err,              # 实部及其误差
        :S_AF_imag, :S_AF_imag_err               # 虚部及其误差
    ]
    
    # 过滤出实际存在的列
    available_columns = intersect(desired_columns, names(df))
    # 添加任何其他列
    other_columns = setdiff(names(df), available_columns)
    all_columns = vcat(available_columns, other_columns)
    
    # 使用columns参数明确指定列顺序
    CSV.write(filename, df, columns=all_columns)
    println("结果已保存到 $filename")
    println("保存的列顺序: $all_columns")
end

# ---------------------------------------------------------------------------- #
#                         Analysis of Structure Factor                         #
# ---------------------------------------------------------------------------- #

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
        
        # 检查结构因子文件是否存在，如果不存在，检查源文件是否存在
        filepath = joinpath(dir, filename)
        source_filepath = joinpath(dir, source_file)
        
        # 如果结构因子文件和源文件都不存在，尝试自动合并文件
        if !isfile(filepath) && !isfile(source_filepath)
            # 尝试自动合并文件
            auto_combined = false
            
            # 根据源文件类型尝试不同的自动合并策略
            auto_combined = try_combine_components(dir, source_file, i, total_dirs, verbose)
            
            # 如果自动合并失败，跳过此目录
            if !auto_combined
                if verbose
                    println("($i/$total_dirs) 跳过 $dirname: 文件 $filename 和源文件 $source_file 都不存在或合并失败")
                end
                continue
            end
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
    try_combine_components(dir::AbstractString, source_file::String, i::Int, total_dirs::Int, verbose::Bool, 
                          combiner_function::Function, args...; kwargs...) -> Bool

通用函数，尝试使用指定的合并函数创建目标文件。

参数：
- `dir`: 包含文件的目录
- `source_file`: 要创建的源文件名
- `i`: 当前目录索引（用于日志记录）
- `total_dirs`: 目录总数（用于日志记录）
- `verbose`: 是否输出详细信息
- `combiner_function`: 用于合并文件的函数
- `args...`: 传递给合并函数的位置参数
- `kwargs...`: 传递给合并函数的关键字参数

返回值：
- `Bool`: 合并是否成功
"""
function try_combine_components(dir::AbstractString, source_file::String, i::Int, total_dirs::Int, verbose::Bool)
    if source_file == "ss_k.bin"
        file1 = "spsm_k.bin"
        file2 = "szsz_k.bin"
        combiner_function = combine_ss_components
    elseif source_file == "cdwpair_k.bin"
        file1 = "cdw_k.bin"
        file2 = "pair_onsite_k.bin"
        combiner_function = combine_cdwpair_components
    else
        if verbose
            println("($i/$total_dirs) 不支持的源文件: $source_file")
        end
        return false
    end
    
    if verbose
        println("($i/$total_dirs) 尝试使用数据 $file1 和 $file2 合并 $source_file...")
    end

    try
        # 调用合并函数
        result = combiner_function(source_file,file1,file2,dir,dir;verbose=false)
        
        # 检查结果
        if result != ""
            if verbose
                println("($i/$total_dirs) 成功创建 $source_file")
            end
            return true
        else
            if verbose
                println("($i/$total_dirs) 创建 $source_file 失败：返回空路径")
            end
            return false
        end
    catch e
        if verbose
            println("($i/$total_dirs) 创建 $source_file 失败: $e")
        end
        return false
    end
end