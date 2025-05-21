#=
多参数分析功能 (Multiple Parameter Analysis)

此文件提供了在包含多个参数文件夹的目录中运行分析的功能。
主要功能包括:
1. 从目录名中提取参数
2. 扫描符合条件的参数目录
3. 对多个参数目录执行分析并整合结果

使用示例:
1. 分析反铁磁结构因子: 
   ```julia
   using DataProcessforDQMC
   results = analyze_af_structure_factor_multi_parameter(
       k_point=(π, π),
       filename="spsm_k.bin",
       startbin=2,
       dropmaxmin=1
   )
   save_analysis_results(results, "af_structure_factor_results.csv")
   ```

2. 自定义分析:
   ```julia
   using DataProcessforDQMC
   # 扫描目录
   dirs = scan_parameter_directories()
   # 对每个目录执行自定义分析
   results = []
   for dir in dirs
       # 提取参数
       params = extract_parameters_from_dirname(basename(dir))
       # 执行分析...
       push!(results, merge(params, analysis_result))
   end
   # 转换为DataFrame
   df = DataFrame(results)
   ```
=#

# 导出多参数分析函数
export extract_parameters_from_dirname,
       scan_parameter_directories,
       analyze_af_structure_factor_multi_parameter,
       save_analysis_results

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
    
    # 提取前缀类型
    if occursin("proj_fft_honeycomb_exact", dirname)
        params[:prefix] = "proj_fft_honeycomb_exact"
    elseif occursin("proj_fft_honeycomb", dirname)
        params[:prefix] = "proj_fft_honeycomb"
    else
        params[:prefix] = "unknown"
    end
    
    # 提取数值参数
    # 温度参数b
    b_match = match(r"b(\d+\.\d+)", dirname)
    params[:b] = b_match !== nothing ? parse(Float64, b_match.captures[1]) : NaN
    
    # 相互作用强度U
    u_match = match(r"U(\d+\.\d+)", dirname)
    params[:U] = u_match !== nothing ? parse(Float64, u_match.captures[1]) : NaN
    
    # Gutzwiller参数gw
    gw_match = match(r"gw(\d+\.\d+)", dirname)
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
    
    return params
end

"""
    scan_parameter_directories(base_dir::AbstractString=pwd(); pattern::Regex=r"^proj_fft_honeycomb") -> Vector{String}

扫描指定目录下所有符合模式的子目录。

# 参数
- `base_dir::AbstractString=pwd()`: 基础目录路径，默认为当前工作目录
- `pattern::Regex=r"^proj_fft_honeycomb"`: 用于匹配目录名的正则表达式

# 返回值
- `Vector{String}`: 符合条件的子目录路径列表
"""
function scan_parameter_directories(base_dir::AbstractString=pwd(); pattern::Regex=r"^proj_fft_honeycomb")
    result_dirs = String[]
    
    # 列出基础目录下的所有条目
    entries = readdir(base_dir, join=true)
    
    # 过滤出目录并匹配模式
    for entry in entries
        if isdir(entry)
            dirname = basename(entry)
            if match(pattern, dirname) !== nothing
                push!(result_dirs, entry)
            end
        end
    end
    
    return result_dirs
end

"""
    analyze_af_structure_factor_multi_parameter(base_dir::AbstractString=pwd(); 
                                               k_point::Tuple{<:Real,<:Real}=(0.0, 0.0),
                                               filename::String="spsm_k.bin",
                                               startbin::Int=2, 
                                               endbin::Union{Int,Nothing}=nothing, 
                                               dropmaxmin::Int=0,
                                               auto_digits::Bool=true, 
                                               tolerance::Float64=1e-6, 
                                               verbose::Bool=false) -> DataFrame

对多个参数目录执行反铁磁结构因子分析，并将结果整合到一个DataFrame中。

# 参数
- `base_dir::AbstractString=pwd()`: 基础目录路径，默认为当前工作目录
- `k_point::Tuple{<:Real,<:Real}=(0.0, 0.0)`: 要分析的k点
- `filename::String="spsm_k.bin"`: 分析的文件名
- `startbin::Int=2`: 起始bin编号
- `endbin::Union{Int,Nothing}=nothing`: 结束bin编号
- `dropmaxmin::Int=0`: 丢弃的最大/最小值数量
- `auto_digits::Bool=true`: 是否自动确定有效数字
- `tolerance::Float64=1e-6`: k点匹配容差
- `verbose::Bool=false`: 是否显示详细信息

# 返回值
- `DataFrame`: 包含所有参数和分析结果的DataFrame
"""
function analyze_af_structure_factor_multi_parameter(base_dir::AbstractString=pwd(); 
                                                   k_point::Tuple{<:Real,<:Real}=(0.0, 0.0),
                                                   filename::String="spsm_k.bin",
                                                   startbin::Int=2, 
                                                   endbin::Union{Int,Nothing}=nothing, 
                                                   dropmaxmin::Int=0,
                                                   auto_digits::Bool=true, 
                                                   tolerance::Float64=1e-6, 
                                                   verbose::Bool=false)
    # 扫描参数目录
    param_dirs = scan_parameter_directories(base_dir)
    
    if isempty(param_dirs)
        @warn "未在 $base_dir 中找到参数目录"
        return DataFrame()
    end
    
    # 结果存储
    results = []
    
    # 对每个目录执行分析
    total_dirs = length(param_dirs)
    println("找到 $total_dirs 个参数目录进行分析...")
    
    for (i, dir) in enumerate(param_dirs)
        dirname = basename(dir)
        
        # 提取参数
        params = extract_parameters_from_dirname(dirname)
        
        # 检查必要文件是否存在
        filepath = joinpath(dir, filename)
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
            analysis_result = AFStructureFactor(
                k_point, filename, dir;
                startbin=startbin, 
                endbin=endbin, 
                dropmaxmin=dropmaxmin,
                auto_digits=auto_digits, 
                tolerance=tolerance,
                verbose=false  # 我们自己处理输出
            )
            
            # 组合参数和结果
            result_entry = Dict(
                # 参数
                :directory => dirname,
                :prefix => params[:prefix],
                :b => params[:b],
                :U => params[:U],
                :L => params[:L],
                :dtau => params[:dtau],
                # Gutzwiller参数
                :gw => get(params, :gw, NaN),  # 如果不存在则返回NaN
                :lprojgw => get(params, :lprojgw, nothing),  # 如果不存在则返回null
                # k点坐标
                :k_point_x => analysis_result.k_point[1],
                :k_point_y => analysis_result.k_point[2],
                # 结果
                :S_AF_real => analysis_result.mean_real,
                :S_AF_real_err => analysis_result.err_real,
                :S_AF_imag => analysis_result.mean_imag,
                :S_AF_imag_err => analysis_result.err_imag
            )
            
            push!(results, result_entry)
            
            if verbose
                println("($i/$total_dirs) 完成分析 $dirname")
                println("   S_AF($(k_point)) = $(analysis_result.formatted_real) + $(analysis_result.formatted_imag)i")
            end
        catch e
            println("($i/$total_dirs) 分析 $dirname 时出错: $e")
        end
    end
    
    # 转换为DataFrame
    if isempty(results)
        @warn "没有成功的分析结果"
        return DataFrame()
    end
    
    result_df = DataFrame(results)
    
    # 排序（按b, U, L, dtau）
    sort!(result_df, [:b, :U, :L, :dtau])
    
    # 调整列顺序，确保值和误差列相邻
    # 定义期望的列顺序
    desired_columns = [
        :directory, :prefix, :b, :U, :L, :dtau,  # 参数
        :k_point_x, :k_point_y,                  # k点坐标
        :S_AF_real, :S_AF_real_err,              # 实部及其误差
        :S_AF_imag, :S_AF_imag_err               # 虚部及其误差
    ]
    
    # 检查所有列是否存在
    existing_columns = intersect(desired_columns, names(result_df))
    other_columns = setdiff(names(result_df), existing_columns)
    
    # 重新排列列顺序
    result_df = result_df[:, vcat(existing_columns, other_columns)]
    
    # 打印列顺序（用于调试）
    if verbose
        println("DataFrame列顺序: $(names(result_df))")
    end
    
    return result_df
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
        :k_point_x, :k_point_y,                  # k点坐标
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

# 文件结束