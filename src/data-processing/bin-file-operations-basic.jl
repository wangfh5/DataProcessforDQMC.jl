#=
数据文件操作功能 (Bin File Operations)

此文件提供了处理DQMC输出的数据文件的各种基础操作函数。
主要功能包括:
1. 合并多个数据文件 （combine）
2. 对数据文件的不同列进行标量缩放 （scale）
3. 以线性叠加的方式合并数据文件的不同列为单列 （merge）
    3.1 两轨道交错合并（merge_staggered_components）
    3.2 两轨道均匀合并（merge_uniform_components）
4. 提取特定坐标的所有bin数据 （filter）
=#

# 导出函数
export combine_bin_files, scale_bin_columns, merge_bin_columns, filter_bin_file
export merge_staggered_components, merge_uniform_components

# ---------------------------------------------------------------------------- #
#                         Basic operations on bin files                        #
# ---------------------------------------------------------------------------- #

"""
    combine_bin_files(output_filename::String,
                      input_files::Vector{<:Tuple{<:AbstractString, <:Real}},
                      output_dir::String=pwd(),
                      input_dir::String=pwd();
                      preserve_columns::Union{Vector{Int}, UnitRange{Int}}=1:2,
                      operation::Function=+,
                      verbose::Bool=false)

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
function combine_bin_files(output_filename::String,
                           input_files::Vector{<:Tuple{<:AbstractString, <:Real}},
                           output_dir::String=pwd(),
                           input_dir::String=pwd();
                           preserve_columns::Union{Vector{Int}, UnitRange{Int}}=1:2,
                           operation::Function=+,
                           verbose::Bool=false)
    # 确保输出文件名以 .bin 结尾
    if !endswith(output_filename, ".bin")
        output_filename = output_filename * ".bin"
    end
    
    output_path = joinpath(output_dir, output_filename)
    
    # 读取第一个文件以获取基本结构信息
    first_file_path = joinpath(input_dir, first(input_files)[1])
    if !isfile(first_file_path)
        @error "找不到输入文件: $first_file_path, 无法合并得到 $output_filename"
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
    scale_bin_columns(
        input_filename::String,
        output_filename::String,
        column_factors::Union{Vector{<:Real}, Dict{Int, <:Real}},
        input_dir::String=pwd(),
        output_dir::String=pwd();
        preserve_columns::Union{Vector{Int}, UnitRange{Int}, Nothing}=nothing,
        default_factor::Real=1.0,
        verbose::Bool=false
    )

对二进制数据文件的指定列应用标量因子。

# 参数
- `input_filename::String`: 输入文件名
- `output_filename::String`: 输出文件名
- `column_factors::Union{Vector{<:Real}, Dict{Int, <:Real}}`: 列因子，可以是向量（按顺序应用）或字典（指定列索引）
- `input_dir::String=pwd()`: 输入文件目录
- `output_dir::String=pwd()`: 输出文件目录
- `preserve_columns::Union{Vector{Int}, UnitRange{Int}, Nothing}=nothing`: 保持不变的列索引（这些列将不应用因子）
- `default_factor::Real=1.0`: 未指定列的默认因子
- `verbose::Bool=false`: 是否输出详细信息

# 返回值
- `String`: 生成的文件完整路径

# 示例
```julia
# 对第3-10列应用不同的标量因子
scale_bin_columns(
    "input.bin",          # 输入文件名
    "scaled_output.bin",  # 输出文件名
    [1, 1, 0.5, 0.5, -1, -1, 0.25, 0.25],  # 列因子
    preserve_columns=1:2,  # 前两列保持不变
    verbose=true
)

# 使用字典指定特定列的因子
scale_bin_columns(
    "input.bin",
    "scaled_output.bin",
    Dict(3 => 0.5, 4 => 0.5, 5 => -1, 6 => -1),  # 只修改指定列
    default_factor=1.0,  # 其他列保持原值
    verbose=true
)
```
"""
function scale_bin_columns(
    input_filename::String,
    output_filename::String,
    column_factors::Union{Vector{<:Real}, Dict{Int, <:Real}},
    input_dir::String=pwd(),
    output_dir::String=pwd();
    preserve_columns::Union{Vector{Int}, UnitRange{Int}, Nothing}=nothing,
    default_factor::Real=1.0,
    verbose::Bool=false
)
    # 确保输出文件名以 .bin 结尾
    if !endswith(output_filename, ".bin")
        output_filename = output_filename * ".bin"
    end
    
    output_path = joinpath(output_dir, output_filename)
    input_path = joinpath(input_dir, input_filename)
    
    if !isfile(input_path)
        @error "找不到输入文件: $input_path"
        return ""
    end
    
    # 读取输入文件
    data = readdlm(input_path)
    nrows, ncols = size(data)
    
    # 处理列因子
    factors = fill(default_factor, ncols)
    
    if column_factors isa Dict
        for (col, factor) in column_factors
            if 1 <= col <= ncols
                factors[col] = factor
            else
                @warn "列索引 $col 超出范围 (1-$ncols)"
            end
        end
    else
        # 如果是向量，直接赋值
        n = min(length(column_factors), ncols)
        factors[1:n] = column_factors[1:n]
    end
    
    # 处理保留列（不应用因子）
    if !isnothing(preserve_columns)
        for col in preserve_columns
            if 1 <= col <= ncols
                factors[col] = 1.0
            end
        end
    end
    
    if verbose
        println("正在缩放文件列:")
        println("  输入文件: $input_path")
        println("  输出文件: $output_path")
        println("  列因子: $factors")
        if !isnothing(preserve_columns)
            println("  保留列（不缩放）: $preserve_columns")
        end
    end
    
    # 应用列因子
    result_data = copy(data)
    for col in 1:ncols
        if factors[col] != 1.0
            result_data[:, col] .*= factors[col]
        end
    end
    
    # 写入结果到输出文件
    open(output_path, "w") do io
        for i in 1:nrows
            for j in 1:ncols
                val = result_data[i, j]
                if isa(val, AbstractFloat)
                    formatted = @sprintf("%16.8e", val)
                    write(io, formatted)
                else
                    formatted = @sprintf("%16s", string(val))
                    write(io, formatted)
                end
            end
            write(io, "\n")
        end
    end
    
    if verbose
        println("缩放完成: $output_path")
    end
    
    return output_path
end

"""
    function merge_bin_columns(
        input_filename::String,
        output_filename::String,
        real_columns::Vector{Int},
        imag_columns::Vector{Int},
        weights::Vector{Float64}=ones(length(real_columns)),
        input_dir::String=pwd(),
        output_dir::String=pwd();
        preserve_columns::UnitRange{Int}=1:2,
        verbose::Bool=false
    )

将一个多列数据中的指定列合并为单列，生成新的数据文件。

参数：
- `input_filename`: 输入文件名
- `output_filename`: 输出文件名
- `real_columns`: 要合并的实部列索引
- `imag_columns`: 要合并的虚部列索引
- `weights`: 每对列的权重，默认全为1.0
- `input_dir`: 输入目录，默认为当前目录
- `output_dir`: 输出目录，默认为当前目录
- `preserve_columns`: 要保留的列范围，默认为1:2（通常是k点坐标）
- `verbose`: 是否输出详细信息，默认为false

返回：
- 合并后的文件完整路径或者空字符串（如果失败）

示例：
```julia
# 合并AA和BB轨道到单列文件
merged_file = merge_bin_columns(
    "spsm_k.bin",      # 输入文件名
    "merged.bin",       # 输出文件名
    [3, 9],           # 实部列（AA和BB轨道的实部）
    [4, 10],          # 虚部列（AA和BB轨道的虚部）
    [1.0, 1.0]        # 权重
)
```
"""
function merge_bin_columns(
    input_filename::String,
    output_filename::String,
    real_columns::Vector{Int},
    imag_columns::Vector{Int},
    weights::Vector{Float64}=ones(length(real_columns)),
    input_dir::String=pwd(),
    output_dir::String=pwd();
    preserve_columns::UnitRange{Int}=1:2,
    verbose::Bool=false
)
    # 确保输出文件名以 .bin 结尾
    if !endswith(output_filename, ".bin")
        output_filename = output_filename * ".bin"
    end
    
    # 检查参数一致性
    @assert length(real_columns) == length(imag_columns) "实部列和虚部列数量必须相同"
    @assert length(weights) == length(real_columns) "权重数量必须与列对数量相同"
    
    # 构建输入和输出文件路径
    input_path = joinpath(input_dir, input_filename)
    if !isfile(input_path)
        @error "找不到输入文件: $input_path"
        return ""
    end
    
    output_path = joinpath(output_dir, output_filename)
    
    # 读取输入文件
    if verbose
        println("正在读取文件: $input_path")
    end
    
    data = readdlm(input_path)
    
    # 创建输出数据
    n_rows = size(data, 1)
    n_preserved = length(preserve_columns)
    output_data = zeros(n_rows, n_preserved + 2)  # +2 用于合并后的实部和虚部列
    
    # 复制要保留的列
    output_data[:, 1:n_preserved] = data[:, preserve_columns]
    
    # 合并列
    if verbose
        println("正在合并列:")
        println("  实部列: $real_columns")
        println("  虚部列: $imag_columns")
        println("  权重: $weights")
    end
    
    for i in 1:length(real_columns)
        real_col = real_columns[i]
        imag_col = imag_columns[i]
        weight = weights[i]
        
        # 累加实部和虚部，应用权重
        output_data[:, n_preserved + 1] .+= weight .* data[:, real_col]
        output_data[:, n_preserved + 2] .+= weight .* data[:, imag_col]
        
        if verbose
            println("  应用权重 $weight 到列 $real_col 和 $imag_col")
        end
    end
    
    # 写入结果到输出文件，使用与原始文件相同的格式
    if verbose
        println("正在写入合并后的文件: $output_path")
    end
    
    open(output_path, "w") do io
        for i in 1:size(output_data, 1)
            for j in 1:size(output_data, 2)
                val = output_data[i, j]
                if isa(val, AbstractFloat)
                    # 使用与 Fortran e16.8 相匹配的格式
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
        println("列合并完成: $output_path")
    end
    
    return output_path
end

# ---------------------------------------------------------------------------- #
#                      Merge functions for 2-orbital data                      #
# ---------------------------------------------------------------------------- #

"""    
    merge_staggered_components(
        input_filename::String="ss_k.bin",
        output_filename::String="afm_sf_k.bin",
        input_dir::String=pwd(),
        output_dir::String=pwd();
        real_columns::Vector{Int}=[3, 9, 5, 7],  # [AA, BB, AB, BA] 实部列
        imag_columns::Vector{Int}=[4, 10, 6, 8], # [AA, BB, AB, BA] 虚部列
        preserve_columns::UnitRange{Int}=1:2,
        verbose::Bool=false
    )

将多列数据使用交错相位的方式合并，生成新的数据文件。
公式: S = A + B - C - D，其中 A, B, C, D 分别对应数组中的四个列索引。

这种合并方式常用于计算反铁磁结构因子，其中 S_AF = AA + BB - AB - BA。
默认的文件名设置（ss_k.bin -> afm_sf_k.bin）反映了这一应用场景。

# 参数
- `input_filename::String="ss_k.bin"`: 输入文件名
- `output_filename::String="afm_sf_k.bin"`: 输出文件名
- `input_dir::String=pwd()`: 输入目录，默认为当前目录
- `output_dir::String=pwd()`: 输出目录，默认为当前目录
- `real_columns::Vector{Int}=[3, 9, 5, 7]`: 实部列索引数组，顺序为[A, B, C, D]
- `imag_columns::Vector{Int}=[4, 10, 6, 8]`: 虚部列索引数组，顺序为[A, B, C, D]
- `preserve_columns::UnitRange{Int}=1:2`: 要保留的列范围，默认为1:2（通常是k点坐标）
- `verbose::Bool=false`: 是否输出详细信息，默认为false

# 返回值
- 合并后的文件完整路径或者空字符串（如果失败）

# 示例
```julia
# 计算反铁磁结构因子
afm_sf_file = merge_staggered_components(
    "ss_k.bin",
    "afm_sf_k.bin"
)

# 自定义列索引
result = merge_staggered_components(
    "custom_input.bin",
    "custom_output.bin",
    real_columns=[3, 5, 7, 9],
    imag_columns=[4, 6, 8, 10]
)
```
"""
function merge_staggered_components(
    input_filename::String="ss_k.bin",
    output_filename::String="afm_sf_k.bin",
    input_dir::String=pwd(),
    output_dir::String=pwd();
    real_columns::Vector{Int}=[3, 9, 5, 7],  # [A, B, C, D] 实部列
    imag_columns::Vector{Int}=[4, 10, 6, 8], # [A, B, C, D] 虚部列
    preserve_columns::UnitRange{Int}=1:2,
    verbose::Bool=false
)
    # 确保列数组长度正确
    @assert length(real_columns) == 4 "实部列数组必须包含4个元素 [A, B, C, D]"
    @assert length(imag_columns) == 4 "虚部列数组必须包含4个元素 [A, B, C, D]"
    
    # 定义合并公式: S = A + B - C - D
    weights = [1.0, 1.0, -1.0, -1.0]  # 正权重用于A和B，负权重用于C和D
    
    if verbose
        println("使用交错相位合并列: S = A + B - C - D")
        println("A: 列 $(real_columns[1]) (实部), $(imag_columns[1]) (虚部), 权重: 1.0")
        println("B: 列 $(real_columns[2]) (实部), $(imag_columns[2]) (虚部), 权重: 1.0")
        println("C: 列 $(real_columns[3]) (实部), $(imag_columns[3]) (虚部), 权重: -1.0")
        println("D: 列 $(real_columns[4]) (实部), $(imag_columns[4]) (虚部), 权重: -1.0")
    end
    
    # 调用通用的merge_bin_columns函数
    result = merge_bin_columns(
        input_filename,
        output_filename,
        real_columns,
        imag_columns,
        weights,
        input_dir,
        output_dir;
        preserve_columns=preserve_columns,
        verbose=verbose
    )
    
    return result
end


"""    
    merge_uniform_components(
        input_filename::String="cdwpair_k.bin",
        output_filename::String="cdwpair_sf_k.bin",
        input_dir::String=pwd(),
        output_dir::String=pwd();
        real_columns::Vector{Int}=[3, 9, 5, 7],  # [A, B, C, D] 实部列
        imag_columns::Vector{Int}=[4, 10, 6, 8], # [A, B, C, D] 虚部列
        preserve_columns::UnitRange{Int}=1:2,
        verbose::Bool=false
    )

将多列数据使用统一相位的方式合并，生成新的数据文件。
公式: S = A + B + C + D，其中 A, B, C, D 分别对应数组中的四个列索引。

这种合并方式常用于计算电荷结构因子，其中 S_charge = AA + BB + AB + BA。
默认的文件名设置（cdwpair_k.bin -> cdwpair_sf_k.bin）反映了这一应用场景。

# 参数
- `input_filename::String="cdwpair_k.bin"`: 输入文件名
- `output_filename::String="cdwpair_sf_k.bin"`: 输出文件名
- `input_dir::String=pwd()`: 输入目录，默认为当前目录
- `output_dir::String=pwd()`: 输出目录，默认为当前目录
- `real_columns::Vector{Int}=[3, 9, 5, 7]`: 实部列索引数组，顺序为[A, B, C, D]
- `imag_columns::Vector{Int}=[4, 10, 6, 8]`: 虚部列索引数组，顺序为[A, B, C, D]
- `preserve_columns::UnitRange{Int}=1:2`: 要保留的列范围，默认为1:2（通常是k点坐标）
- `verbose::Bool=false`: 是否输出详细信息，默认为false

# 返回值
- 合并后的文件完整路径或者空字符串（如果失败）

# 示例
```julia
# 计算电荷结构因子
cdwpair_sf_file = merge_uniform_components(
    "cdwpair_k.bin",
    "cdwpair_sf_k.bin"
)

# 自定义列索引
result = merge_uniform_components(
    "custom_input.bin",
    "custom_output.bin",
    real_columns=[3, 5, 7, 9],
    imag_columns=[4, 6, 8, 10]
)
```
"""
function merge_uniform_components(
    input_filename::String="cdwpair_k.bin",
    output_filename::String="cdwpair_sf_k.bin",
    input_dir::String=pwd(),
    output_dir::String=pwd();
    real_columns::Vector{Int}=[3, 9, 5, 7],  # [A, B, C, D] 实部列
    imag_columns::Vector{Int}=[4, 10, 6, 8], # [A, B, C, D] 虚部列
    preserve_columns::UnitRange{Int}=1:2,
    verbose::Bool=false
)
    # 确保列数组长度正确
    @assert length(real_columns) == 4 "实部列数组必须包含4个元素 [A, B, C, D]"
    @assert length(imag_columns) == 4 "虚部列数组必须包含4个元素 [A, B, C, D]"
    
    # 定义合并公式: S = A + B + C + D
    weights = [1.0, 1.0, 1.0, 1.0]  # 所有分量都使用正权重
    
    if verbose
        println("使用统一相位合并列: S = A + B + C + D")
        println("A: 列 $(real_columns[1]) (实部), $(imag_columns[1]) (虚部), 权重: 1.0")
        println("B: 列 $(real_columns[2]) (实部), $(imag_columns[2]) (虚部), 权重: 1.0")
        println("C: 列 $(real_columns[3]) (实部), $(imag_columns[3]) (虚部), 权重: 1.0")
        println("D: 列 $(real_columns[4]) (实部), $(imag_columns[4]) (虚部), 权重: 1.0")
    end
    
    # 调用通用的merge_bin_columns函数
    result = merge_bin_columns(
        input_filename,
        output_filename,
        real_columns,
        imag_columns,
        weights,
        input_dir,
        output_dir;
        preserve_columns=preserve_columns,
        verbose=verbose
    )
    
    return result
end

# ---------------------------------------------------------------------------- #
#                      Filter functions for bin files                          #
# ---------------------------------------------------------------------------- #

"""
    filter_bin_file(
        input_filename::String,
        coordinate::Union{Tuple{Int, Int}, Tuple{Float64, Float64}};
        output_filename::Union{String, Nothing}=nothing,
        dir::String=pwd(),
        k_point_tolerance::Float64=1e-6,
        verbose::Bool=false
    )

Extract all bin data for a specific coordinate from a bin file.

This function filters a bin file to extract only the rows corresponding to a specific 
coordinate across all bins. The coordinate can be either:
- k-space (momentum): Tuple{Float64, Float64}, e.g., (0.5, 0.5)
- r-space (real space): Tuple{Int, Int}, e.g., (1, 2)

# Arguments
- `input_filename::String`: Input .bin file (e.g., "nn_k.bin", "nn_r.bin")
- `coordinate::Union{Tuple{Int, Int}, Tuple{Float64, Float64}}`: Target coordinate
- `output_filename::Union{String, Nothing}=nothing`: Output .bin file name. If nothing, auto-generate from coordinate
- `dir::String=pwd()`: Working directory
- `k_point_tolerance::Float64=1e-6`: Tolerance for k-space coordinate matching (ignored for r-space)
- `verbose::Bool=false`: Whether to print verbose information

# Returns
- `String`: Path to output file, or empty string if failed

# Examples
```julia
# Filter k-space data for k=(π,π)
filter_bin_file("nn_k.bin", (0.5, 0.5), dir="/path/to/data")
# Output: nn_k_0.500_0.500.bin

# Filter r-space data for r=(1, 2)
filter_bin_file("nn_r.bin", (1, 2), dir="/path/to/data")
# Output: nn_r_1_2.bin

# Custom output filename
filter_bin_file("afm_sf_k.bin", (0.5, 0.5), output_filename="afm_pi_pi.bin")
```
"""
function filter_bin_file(
    input_filename::String,
    coordinate::Union{Tuple{Int, Int}, Tuple{Float64, Float64}};
    output_filename::Union{String, Nothing}=nothing,
    dir::String=pwd(),
    k_point_tolerance::Float64=1e-6,
    verbose::Bool=false
)
    # Construct input path
    input_path = joinpath(dir, input_filename)
    if !isfile(input_path)
        @error "Input file not found: $input_path"
        return ""
    end
    
    # Read the bin file
    if verbose
        println("Reading file: $input_path")
    end
    data = readdlm(input_path)
    
    # Extract unique coordinates (first two columns)
    coords = unique(data[:, 1:2], dims=1)
    n_coords = size(coords, 1)
    n_bins = size(data, 1) ÷ n_coords
    
    if verbose
        println("File structure:")
        println("  Total rows: $(size(data, 1))")
        println("  Unique coordinates: $n_coords")
        println("  Number of bins: $n_bins")
        println("  Columns: $(size(data, 2))")
    end
    
    # Find the coordinate index
    coord_idx, is_exact = find_coordinate_index(coords, coordinate, k_point_tolerance)
    
    if isnothing(coord_idx)
        @error "Coordinate $coordinate not found in the file"
        return ""
    end
    
    if verbose
        if coordinate isa Tuple{Int, Int}
            println("Found r-space coordinate: $coordinate at index $coord_idx")
        else
            matched_coord = (coords[coord_idx, 1], coords[coord_idx, 2])
            if is_exact
                println("Found exact k-space match: $coordinate at index $coord_idx")
            else
                dist = sqrt((matched_coord[1] - coordinate[1])^2 + (matched_coord[2] - coordinate[2])^2)
                println("Found closest k-space coordinate: $matched_coord (distance: $dist) at index $coord_idx")
                @warn "No exact match found within tolerance $k_point_tolerance, using closest point"
            end
        end
    end
    
    # Extract rows for this coordinate across all bins
    # Row indices: coord_idx, coord_idx + n_coords, coord_idx + 2*n_coords, ...
    filtered_rows = [coord_idx + i * n_coords for i in 0:(n_bins-1)]
    filtered_data = data[filtered_rows, :]
    
    # Generate output filename if not provided
    if isnothing(output_filename)
        base_name = replace(input_filename, ".bin" => "")
        # Use the actual matched coordinate from the file, not the input coordinate
        matched_coord = (coords[coord_idx, 1], coords[coord_idx, 2])
        if coordinate isa Tuple{Int, Int}
            # R-space: use integer format
            output_filename = "$(base_name)_$(matched_coord[1])_$(matched_coord[2]).bin"
        else
            # K-space: use float format with 3 decimal places
            coord_str1 = @sprintf("%.3f", matched_coord[1])
            coord_str2 = @sprintf("%.3f", matched_coord[2])
            output_filename = "$(base_name)_$(coord_str1)_$(coord_str2).bin"
        end
    end
    
    # Ensure .bin extension
    if !endswith(output_filename, ".bin")
        output_filename = output_filename * ".bin"
    end
    
    output_path = joinpath(dir, output_filename)
    
    # Write filtered data to output file
    if verbose
        println("Writing $(size(filtered_data, 1)) rows to: $output_path")
    end
    
    open(output_path, "w") do io
        for i in 1:size(filtered_data, 1)
            for j in 1:size(filtered_data, 2)
                val = filtered_data[i, j]
                if isa(val, AbstractFloat)
                    # Use Fortran-compatible format
                    formatted = @sprintf("%16.8e", val)
                    write(io, formatted)
                else
                    formatted = @sprintf("%16s", string(val))
                    write(io, formatted)
                end
            end
            write(io, "\n")
        end
    end
    
    if verbose
        println("Filter completed successfully")
        println("Output file: $output_path")
        println("Extracted $n_bins bins for coordinate $coordinate")
    end
    
    return output_path
end