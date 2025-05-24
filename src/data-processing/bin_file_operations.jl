#=
二进制文件操作功能 (Binary File Operations)

此文件提供了处理DQMC输出的二进制数据文件的各种操作函数。
主要功能包括:
1. 合并多个二进制文件
2. 为物理量计算提供便捷函数

使用示例:
1. 合并自旋分量: 
   ```julia
   # 合并XY和Z分量的自旋相关函数
   combined_file = combine_ss_components(
       "ss_k.bin",
       "spsm_k.bin",
       "szsz_k.bin",
       1.0,  # XY分量权重
       1.0   # Z分量权重
   )
   
   # 使用合并后的文件进行结构因子分析
   result = AFMStructureFactor(filename=combined_file)
   ```

2. 自定义文件合并:
   ```julia
   # 合并多个文件，自定义权重和保留列
   combined_file = combine_bin_files(
       "combined.bin",
       [("file1.bin", 0.5), ("file2.bin", 1.0), ("file3.bin", 0.3)],
       preserve_columns=1:3  # 保留前三列不变
   )
   ```
=#

# 导出函数
export combine_bin_files, combine_ss_components

"""
    combine_bin_files(
        output_filename::String,
        input_files::Vector{Tuple{String, Real}},
        output_dir::String=pwd(),
        input_dir::String=pwd();
        preserve_columns::Union{Vector{Int}, UnitRange{Int}}=1:2,
        operation::Function=+,
        verbose::Bool=false
    )

将多个具有类似结构的二进制数据文件按权重组合成一个新文件。

# 参数
- `output_filename::String`: 输出文件名
- `input_files::Vector{Tuple{String, Real}}`: 输入文件名和权重的元组数组，如 [("spsm_k.bin", 0.5), ("szsz_k.bin", 1.0)]
- `output_dir::String=pwd()`: 输出文件目录
- `input_dir::String=pwd()`: 输入文件目录
- `preserve_columns::Union{Vector{Int}, UnitRange{Int}}=1:2`: 保持不变的列索引（默认前两列，通常是坐标）
- `operation::Function=+`: 组合操作函数（默认为加法）
- `verbose::Bool=false`: 是否输出详细信息

# 返回值
- `String`: 生成的文件完整路径

# 示例
```julia
# 合并 spsm_k.bin 和 szsz_k.bin 文件，权重分别为 0.5 和 1.0
combined_file = combine_bin_files(
    "ss_k.bin",
    [("spsm_k.bin", 1.0), ("szsz_k.bin", 1.0)]
)

# 合并实空间关联函数文件，保持前3列不变
combined_file = combine_bin_files(
    "ss_r.bin", 
    [("spsm_r.bin", 1.0), ("szsz_r.bin", 1.0)],
    preserve_columns=1:3
)

# 使用乘法而不是加法合并文件
combined_file = combine_bin_files(
    "product.bin",
    [("file1.bin", 1.0), ("file2.bin", 1.0)],
    operation=*
)
```
"""
function combine_bin_files(
    output_filename::String,
    input_files::Vector{<:Tuple{<:AbstractString, <:Real}},
    output_dir::String=pwd(),
    input_dir::String=pwd();
    preserve_columns::Union{Vector{Int}, UnitRange{Int}}=1:2,
    operation::Function=+,
    verbose::Bool=false
)
    # 确保输出文件名以 .bin 结尾
    if !endswith(output_filename, ".bin")
        output_filename = output_filename * ".bin"
    end
    
    output_path = joinpath(output_dir, output_filename)
    
    # 读取第一个文件以获取基本结构信息
    first_file_path = joinpath(input_dir, first(input_files)[1])
    if !isfile(first_file_path)
        @error "找不到输入文件: $first_file_path"
        return ""
    end
    
    # 读取二进制文件
    first_file_data = readdlm(first_file_path)
    
    # 创建结果数据结构，初始化为第一个文件的数据
    # 这样保证保留列的数据来自第一个文件
    result_data = copy(first_file_data)
    
    # 设置组合列的初始值为第一个文件的值乘以它的权重
    combine_columns = setdiff(1:size(result_data, 2), preserve_columns)
    result_data[:, combine_columns] .*= first(input_files)[2]
    
    if verbose
        println("正在合并文件:")
        println("  保持不变的列: $preserve_columns")
        println("  要组合的列: $combine_columns")
        println("  使用的操作: $operation")
    end
    
    # 处理其余输入文件
    for (i, (filename, weight)) in enumerate(input_files)
        # 跳过第一个文件，因为已经处理过了
        if i == 1
            if verbose
                println("  处理文件 $i: $filename (权重: $weight) (已初始化)")
            end
            continue
        end
        
        if verbose
            println("  处理文件 $i: $filename (权重: $weight)")
        end
        
        file_path = joinpath(input_dir, filename)
        if !isfile(file_path)
            @error "找不到输入文件: $file_path"
            return ""
        end
        
        # 读取文件数据
        data = readdlm(file_path)
        
        # 检查数据结构是否兼容
        if size(data) != size(result_data)
            @error "文件 $filename 的数据结构与其他文件不兼容"
            return ""
        end
        
        # 检查保留列是否匹配
        for col in preserve_columns
            if !all(isapprox.(data[:, col], first_file_data[:, col]))
                @warn "文件 $filename 的保留列 $col 与第一个文件不完全匹配"
            end
        end
        
        # 按权重组合数据（只处理需要组合的列）
        for col in combine_columns
            result_data[:, col] = operation.(result_data[:, col], weight .* data[:, col])
        end
    end
    
    # 写入结果到输出文件
    # 使用与原始文件相同的格式
    open(output_path, "w") do io
        for i in 1:size(result_data, 1)
            for j in 1:size(result_data, 2)
                # 对数值使用科学计数法格式化，确保与原始文件一致
                val = result_data[i, j]
                if isa(val, AbstractFloat)
                    # 使用与 Fortran e16.8 相匹配的格式
                    # 确保数字宽度为 16，小数点后 8 位
                    # 这将保证负数和正数对齐
                    formatted = @sprintf("%16.8e", val)
                    write(io, formatted)
                else
                    # 对于非浮点数，使用固定宽度格式
                    formatted = @sprintf("%16s", string(val))
                    write(io, formatted)
                end
            end
            write(io, "\n")
        end
    end
    
    if verbose
        println("合并后的文件已保存到: $output_path")
    end
    
    return output_path
end

"""
    combine_ss_components(
        output_filename::String="ss_k.bin",
        spsm_filename::String="spsm_k.bin",
        szsz_filename::String="szsz_k.bin",
        xy_weight::Real=1.0,
        zz_weight::Real=1.0,
        output_dir::String=pwd(),
        input_dir::String=pwd();
        preserve_columns::Union{Vector{Int}, UnitRange{Int}}=1:2,
        verbose::Bool=false
    )

合并 XY 和 Z 分量的自旋相关函数数据文件。

# 参数
- `output_filename::String="ss_k.bin"`: 输出文件名
- `spsm_filename::String="spsm_k.bin"`: XY分量文件名
- `szsz_filename::String="szsz_k.bin"`: Z分量文件名
- `xy_weight::Real=1.0`: XY分量的权重
- `zz_weight::Real=1.0`: Z分量的权重
- `output_dir::String=pwd()`: 输出文件目录
- `input_dir::String=pwd()`: 输入文件目录
- `preserve_columns::Union{Vector{Int}, UnitRange{Int}}=1:2`: 保持不变的列索引
- `verbose::Bool=false`: 是否输出详细信息

# 返回值
- `String`: 生成的文件完整路径

# 示例
```julia
# 合并 k 空间的 XY 和 Z 分量
combined_file = combine_ss_components(
    "ss_k.bin",
    "spsm_k.bin",
    "szsz_k.bin"
)

# 合并实空间关联函数，保持前3列不变
combined_file = combine_ss_components(
    "ss_r.bin",
    "spsm_r.bin",
    "szsz_r.bin",
    0.5, 1.0,
    pwd(), pwd(),
    preserve_columns=1:3
)
```
"""
function combine_ss_components(
    output_filename::String="ss_k.bin",
    spsm_filename::String="spsm_k.bin",
    szsz_filename::String="szsz_k.bin",
    xy_weight::Real=1.0,
    zz_weight::Real=1.0,
    output_dir::String=pwd(),
    input_dir::String=pwd();
    preserve_columns::Union{Vector{Int}, UnitRange{Int}}=1:2,
    verbose::Bool=false
)
    return combine_bin_files(
        output_filename,
        [(spsm_filename, xy_weight), (szsz_filename, zz_weight)],
        output_dir,
        input_dir;
        preserve_columns=preserve_columns,
        verbose=verbose
    )
end
