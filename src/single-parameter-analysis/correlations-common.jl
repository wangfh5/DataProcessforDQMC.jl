
"""
    filter_bins(df, startbin, endbin, dropmaxmin, n_bins, verbose; value_columns=[:real_val, :imag_val])

通用bin过滤函数，支持实空间和k空间数据。

# 参数
- `df`: 输入DataFrame
- `startbin`: 起始bin
- `endbin`: 结束bin
- `dropmaxmin`: 丢弃的最大最小值数量 （bin大小的metric为所有坐标，所有`value_columns`的平均值）
- `n_bins`: 总bin数
- `verbose`: 是否打印详细信息
- `value_columns`: 用于过滤的数值列，默认为[:real_val, :imag_val]

# 过滤逻辑
该函数按 bin 进行过滤，对每个 bin 计算其所有数据点的平均值，然后排序并丢弃最大/最小的 bin。这确保了在保留数据空间完整性的同时，去除了时间上的离群值。

# 返回
过滤后的DataFrame
"""
function filter_bins(df::DataFrame, startbin::Int, endbin::Union{Int, Nothing}, dropmaxmin::Int, n_bins::Int, verbose::Bool; 
                    value_columns::Vector{Symbol}=[:real_val, :imag_val])
    # 选择bin范围
    if !isnothing(endbin)
        filtered_df = df[(df.bin .>= startbin) .& (df.bin .<= endbin), :]
    else
        filtered_df = df[df.bin .>= startbin, :]
    end
    
    selected_bins = unique(filtered_df.bin)
    selected_count = length(selected_bins)
    
    # 应用max/min过滤（如果需要）
    if dropmaxmin > 0 && selected_count > 2*dropmaxmin
        # 检查所有请求的列是否存在
        valid_columns = [col for col in value_columns if col in propertynames(filtered_df)]
        
        if isempty(valid_columns)
            @warn "None of the specified value_columns exist in the DataFrame, do not apply max/min filtering."
        else
            # 首先，对每个bin计算所有列的平均值总和
            bin_avg_df = combine(groupby(filtered_df, :bin)) do group
                total = 0.0
                for col in valid_columns
                    total += mean(group[!, col])
                end
                (bin=first(group.bin), bin_avg=total)
            end
            
            # 对bin排序
            sorted_bins = sort(bin_avg_df, :bin_avg).bin
            valid_bins = sorted_bins[(dropmaxmin+1):(end-dropmaxmin)]

            # 只保留有效bin的数据
            filtered_df = Base.filter(row -> row.bin in valid_bins, filtered_df)
        end
    end
    
    if verbose
        filtered_count = length(unique(filtered_df.bin))
        
        # 详细的bin选择和过滤信息
        if n_bins > 1
            if startbin > 1 || !isnothing(endbin)
                bin_range = "$startbin-$(isnothing(endbin) ? n_bins : min(endbin, n_bins))"
                println("  Selected bins: $selected_count of $n_bins (range: $bin_range)")
            else
                println("  Selected bins: $selected_count of $n_bins")
            end
            
            if dropmaxmin > 0
                println("  Filtered bins: $filtered_count (removed $dropmaxmin max/min values)")
            else
                println("  No max/min filtering applied")
            end
        end
    end
    
    return filtered_df
end


"""
    calculate_statistics(df, auto_digits; input_coord_columns=[:imj_x, :imj_y], output_coord_names=[:imj_x, :imj_y], orbital_labels=nothing)

通用函数，计算DataFrame中每个坐标的统计信息。
该函数可用于实空间和k空间数据，支持单轨道和多轨道数据。

参数:
- `df`: 包含bin信息的DataFrame
- `auto_digits`: 是否使用自动精度进行误差计算
- `input_coord_columns`: 输入DataFrame中提取坐标值的列 (例如, [:imj_x, :imj_y] 或 [:kx, :ky])
- `output_coord_names`: 输出DataFrame中使用的坐标列名称
- `orbital_labels`: 轨道标签数组。如果为nothing，则处理单轨道数据（使用real_val和imag_val列）

返回:
- 包含每个坐标统计信息的DataFrame
"""
function calculate_statistics(df, auto_digits; 
                             input_coord_columns=[:imj_x, :imj_y],
                             output_coord_names=[:imj_x, :imj_y],
                             orbital_labels=nothing)
    # 定义一个内部函数来计算统计量
    function compute_stats(real_values, imag_values)
        # 计算统计量
        mean_real = mean(real_values)
        mean_imag = mean(imag_values)
        
        # 计算误差（完整精度和格式化精度）
        err_real = error(real_values, sigma=1, bessel=true, auto_digits=false)
        err_imag = error(imag_values, sigma=1, bessel=true, auto_digits=false)
        err_real_fmt = error(real_values, sigma=1, bessel=true, auto_digits=auto_digits)
        err_imag_fmt = error(imag_values, sigma=1, bessel=true, auto_digits=auto_digits)
        
        # 格式化显示
        formatted_real, formatted_real_err = format_value_error(mean_real, err_real_fmt)
        formatted_imag, formatted_imag_err = format_value_error(mean_imag, err_imag_fmt)
        
        return (;
            mean_real = mean_real,
            mean_imag = mean_imag,
            err_real = err_real,
            err_imag = err_imag,
            formatted_real = "$(formatted_real) ± $(formatted_real_err)",
            formatted_imag = "$(formatted_imag) ± $(formatted_imag_err)"
        )
    end
    
    # 对每个坐标分组
    grouped_df = groupby(df, :coord)
    
    # 对每个组计算统计量
    results_df = combine(grouped_df) do group_df
        result = NamedTuple()
        
        # 添加坐标信息
        for (i, col) in enumerate(input_coord_columns)
            result = merge(result, (;Symbol(output_coord_names[i]) => first(group_df[!, col])))
        end
        
        # 单轨道模式
        if isnothing(orbital_labels)
            stats = compute_stats(group_df.real_val, group_df.imag_val)
            result = merge(result, stats)
        else
            # 多轨道模式
            for orbital in orbital_labels
                real_col = Symbol("$(orbital)_real")
                imag_col = Symbol("$(orbital)_imag")
                
                stats = compute_stats(group_df[!, real_col], group_df[!, imag_col])
                
                # 添加轨道前缀到统计量名称
                orbital_stats = NamedTuple{
                    (Symbol("$(orbital)_mean_real"),
                     Symbol("$(orbital)_mean_imag"),
                     Symbol("$(orbital)_err_real"),
                     Symbol("$(orbital)_err_imag"),
                     Symbol("$(orbital)_formatted_real"),
                     Symbol("$(orbital)_formatted_imag"))
                }((
                    stats.mean_real,
                    stats.mean_imag,
                    stats.err_real,
                    stats.err_imag,
                    stats.formatted_real,
                    stats.formatted_imag
                ))
                
                result = merge(result, orbital_stats)
            end
        end
        
        return result
    end
    
    # 按坐标排序
    sort!(results_df, output_coord_names)
    
    return results_df
end