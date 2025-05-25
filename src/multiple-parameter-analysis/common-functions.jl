#=
多参数分析功能 (Multiple Parameter Analysis) 之公共函数

此文件提供了在包含多个参数文件夹的目录中运行分析的功能。
主要功能包括:
1. 从目录名中提取参数
2. 扫描符合条件的参数目录
3. 保存分析结果（dataframe）到本地文件（csv）
=#

export extract_parameters_from_dirname,
       scan_parameter_directories,
       save_analysis_results

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
    
    # AFM参数 Mz_AFM
    mzafm_match = match(r"Mz_AFM(-?\d+\.\d+)", dirname)
    params[:Mz_AFM] = mzafm_match !== nothing ? parse(Float64, mzafm_match.captures[1]) : NaN

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