#=
多参数分析功能 (Multiple Parameter Analysis) 之公共函数

此文件提供了在包含多个参数文件夹的目录中运行分析的功能。
主要功能包括:
1. 从目录名中提取参数
2. 扫描符合条件的参数目录
3. 保存分析结果（dataframe）到本地文件（csv）
4. 批量迁移legacy格式文件夹名
=#

export scan_parameter_directories,
       save_analysis_results,
       batch_migrate_directories

# ---------------------------------------------------------------------------- #
#               Helper functions for multiple parameter analysis               #
# ---------------------------------------------------------------------------- #


"""
    scan_parameter_directories(base_dir::AbstractString=pwd(); 
                              filter_options::Union{Dict, NamedTuple}=Dict(), 
                              return_params::Bool=false) -> Union{Vector{String}, Vector{Tuple{String,String,Vector{Tuple{String,Any,Int}}}}}

扫描指定目录下所有符合筛选条件的子目录。

# 参数
- `base_dir::AbstractString=pwd()`: 基础目录路径，默认为当前工作目录
- `filter_options::Union{Dict, NamedTuple}=Dict()`: 筛选选项，可以包含以下键：
  - `"prefix"` 或 `:prefix`: 目录前缀（字符串或字符串数组）
  - `"b"`, `"U"`, `"L"`, `"dtau"`, `"gw"` 等: 参数范围（可以是单个值、范围或值数组）
  - `"lprojgw"`: 布尔值，是否使用lprojgw
- `return_params::Bool=false`: 是否同时返回解析的参数信息

# 返回值
- 当 `return_params=false`: `Vector{String}` - 符合条件的子目录路径列表
- 当 `return_params=true`: `Vector{Tuple{String,String,Vector{Tuple{String,Any,Int}}}}` - (目录路径, 前缀, 参数列表)的元组列表

# 示例
```julia
# 筛选前缀为 "proj_fft_honeycomb_exact" 的目录
dirs = scan_parameter_directories(filter_options=Dict("prefix" => "proj_fft_honeycomb_exact"))

# 筛选多个条件并返回参数信息
dirs_with_params = scan_parameter_directories(
    filter_options=Dict("U" => 4.0, "L" => [6, 9]),
    return_params=true
)
```
"""
function scan_parameter_directories(base_dir::AbstractString=pwd(); 
                                   filter_options::Union{Dict, NamedTuple}=Dict(), 
                                   return_params::Bool=false)
    result_dirs = String[]
    result_params = Vector{Tuple{String,String,Vector{Tuple{String,Any,Int}}}}()
    
    # 列出基础目录下的所有条目
    entries = readdir(base_dir, join=true)
    
    # 转换NamedTuple为Dict以统一处理
    filter_dict = filter_options isa NamedTuple ? Dict(pairs(filter_options)) : filter_options
    
    # 过滤出目录
    for entry in entries
        if isdir(entry)
            dirname = basename(entry)
            
            # 直接使用parse_jobname解析参数
            prefix = ""
            params_vector = Tuple{String,Any,Int}[]
            try
                prefix, params_vector = parse_jobname(dirname)
            catch
                # 如果解析失败，跳过此目录
                continue
            end
            
            # 如果没有筛选选项，则直接添加目录
            if isempty(filter_dict)
                push!(result_dirs, entry)
                if return_params
                    push!(result_params, (entry, prefix, params_vector))
                end
                continue
            end
            
            # 检查是否满足所有筛选条件
            matches_all_filters = true
            
            # 创建参数字典以便查找
            params_dict = Dict{String,Any}()
            params_dict["prefix"] = prefix
            for (param_key, param_value, _) in params_vector
                params_dict[param_key] = param_value
            end
            
            for (key, filter_value) in filter_dict
                # 将key转换为字符串（支持Symbol和String）
                key_str = string(key)
                # 获取参数值（如果不存在则为nothing）
                param_value = get(params_dict, key_str, nothing)
                
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
                if return_params
                    push!(result_params, (entry, prefix, params_vector))
                end
            end
        end
    end
    
    return return_params ? result_params : result_dirs
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

"""
    batch_migrate_directories(base_dir::AbstractString=pwd(); 
                             pattern::Regex=r"^proj_fft_honeycomb",
                             dry_run::Bool=true,
                             action::String="mv") -> Vector{Tuple{String,String,Bool}}

批量迁移目录下的legacy格式文件夹名到新格式。

# 参数
- `base_dir::AbstractString=pwd()`: 基础目录路径，默认为当前工作目录
- `pattern::Regex=r"^proj_fft_honeycomb"`: 用于匹配需要迁移的目录名的正则表达式
- `dry_run::Bool=true`: 是否为试运行模式（只显示将要进行的操作，不实际执行）
- `action::String="mv"`: 执行的操作类型
  - `"mv"`: 重命名目录（默认）
  - `"cp"`: 复制目录到新名称，保留原目录
  - `"rm"`: 删除原目录（需要先用"cp"创建新目录）

# 返回值
- `Vector{Tuple{String,String,Bool}}`: 迁移结果列表，每个元组包含 (原名, 新名, 是否成功)

# 示例
```julia
# 试运行，查看将要进行的迁移
results = batch_migrate_directories(dry_run=true)

# 复制到新格式，保留原目录
results = batch_migrate_directories(dry_run=false, action="cp")

# 重命名到新格式
results = batch_migrate_directories(dry_run=false, action="mv")
```
"""
function batch_migrate_directories(base_dir::AbstractString=pwd(); 
                                  pattern::Regex=r"^proj_fft_honeycomb",
                                  dry_run::Bool=true,
                                  action::String="mv")
    
    migration_results = Vector{Tuple{String,String,Bool}}()
    
    # 列出基础目录下的所有条目
    entries = readdir(base_dir, join=false)  # 只获取名称，不包含路径
    
    # 验证action参数
    if !(action in ["mv", "cp", "rm"])
        error("无效的action参数: $(action)。支持的值: \"mv\", \"cp\", \"rm\"")
    end
    
    println("=== 批量迁移Legacy格式目录名 ===")
    println("基础目录: $(base_dir)")
    println("匹配模式: $(pattern)")
    println("模式: $(dry_run ? "试运行" : "实际执行")")
    if !dry_run
        println("操作: $(action)")
    end
    println()
    
    # 过滤出符合条件的目录
    legacy_dirs = String[]
    for entry in entries
        full_path = joinpath(base_dir, entry)
        if isdir(full_path) && match(pattern, entry) !== nothing
            # 检查是否为legacy格式（包含点分隔的参数）
            if occursin(r"\.[a-zA-Z]+\d+", entry)
                push!(legacy_dirs, entry)
            end
        end
    end
    
    if isempty(legacy_dirs)
        println("未找到符合条件的legacy格式目录")
        return migration_results
    end
    
    println("找到 $(length(legacy_dirs)) 个legacy格式目录:")
    
    for (i, old_name) in enumerate(legacy_dirs)
        println("\n[$i/$(length(legacy_dirs))] 处理: $old_name")
        
        try
            # 生成新格式名称
            prefix, params = parse_jobname_legacy(old_name)
            new_name = generate_jobname(prefix, params)
            
            # 验证迁移
            is_valid = verify_migration(old_name, new_name)
            
            if !is_valid
                println("  ❌ 迁移验证失败")
                push!(migration_results, (old_name, new_name, false))
                continue
            end
            
            println("  原名: $old_name")
            println("  新名: $new_name")
            println("  验证: ✅")
            
            # 计算压缩比
            compression = round((1 - length(new_name)/length(old_name))*100, digits=1)
            println("  压缩: $compression%")
            
            if dry_run
                println("  状态: 试运行 - 未实际执行")
                push!(migration_results, (old_name, new_name, true))
            else
                # 实际执行迁移
                old_path = joinpath(base_dir, old_name)
                new_path = joinpath(base_dir, new_name)
                
                # 检查新名称是否已存在（除非是rm操作）
                if action != "rm" && isdir(new_path)
                    println("  ❌ 目标目录已存在: $(new_name)")
                    push!(migration_results, (old_name, new_name, false))
                    continue
                end
                
                # 执行相应操作
                try
                    if action == "cp"
                        cp(old_path, new_path)
                        println("  ✅ 复制成功")
                    elseif action == "mv"
                        mv(old_path, new_path)
                        println("  ✅ 重命名成功")
                    elseif action == "rm"
                        if !isdir(old_path)
                            println("  ❌ 原目录不存在: $(old_name)")
                            push!(migration_results, (old_name, new_name, false))
                            continue
                        end
                        rm(old_path, recursive=true)
                        println("  ✅ 删除成功")
                    end
                    push!(migration_results, (old_name, new_name, true))
                catch e
                    println("  ❌ 操作失败: $(e)")
                    push!(migration_results, (old_name, new_name, false))
                end
            end
            
        catch e
            println("  ❌ 处理失败: $e")
            push!(migration_results, (old_name, "", false))
        end
    end
    
    # 总结
    println("\n=== 迁移总结 ===")
    successful = count(x -> x[3], migration_results)
    total = length(migration_results)
    println("总计: $total 个目录")
    println("成功: $successful 个")
    println("失败: $(total - successful) 个")
    
    if dry_run
        println("\n💡 这是试运行结果。要实际执行迁移，请设置 dry_run=false")
    end
    
    return migration_results
end