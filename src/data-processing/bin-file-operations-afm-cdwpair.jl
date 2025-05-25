#= 
数据文件操作功能 (Bin File Operations)

应用于反铁磁， CDW 和配对关联函数的计算。
=#

export combine_ss_components, combine_cdwpair_components
export create_cdw_from_nn
export merge_staggered_components, merge_uniform_components

"""
    combine_ss_components(
        output_filename::String="ss_k.bin",
        spsm_filename::String="spsm_k.bin",
        szsz_filename::String="szsz_k.bin",
        output_dir::String=pwd(),
        input_dir::String=pwd(),
        xy_weight::Real=1.0,
        zz_weight::Real=1.0;
        preserve_columns::Union{Vector{Int}, UnitRange{Int}}=1:2,
        verbose::Bool=false
    )

合并 XY 和 Z 分量的自旋相关函数数据文件。

# 参数
- `output_filename::String="ss_k.bin"`: 输出文件名
- `spsm_filename::String="spsm_k.bin"`: XY分量文件名
- `szsz_filename::String="szsz_k.bin"`: Z分量文件名
- `output_dir::String=pwd()`: 输出文件目录
- `input_dir::String=pwd()`: 输入文件目录
- `xy_weight::Real=1.0`: XY分量的权重
- `zz_weight::Real=1.0`: Z分量的权重
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
    output_dir::String=pwd(),
    input_dir::String=pwd(),
    xy_weight::Real=1.0,
    zz_weight::Real=1.0;
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


"""
    create_cdw_from_nn(
        output_filename::String="cdw_k.bin",
        input_filename::String="nn_k.bin",
        output_dir::String=pwd(),
        input_dir::String=pwd();
        preserve_columns::Union{Vector{Int}, UnitRange{Int}}=1:2,
        verbose::Bool=false
    )

将密度密度关联数据(nn_k.bin)转换为CDW关联数据。
应用η_α η_β因子到对应的轨道，其中η_α是+1（α=A轨道）或-1（α=B轨道）。

# 参数
- `output_filename::String="cdw_k.bin"`: 输出文件名
- `input_filename::String="nn_k.bin"`: 输入密度-密度关联文件名
- `output_dir::String=pwd()`: 输出文件目录
- `input_dir::String=pwd()`: 输入文件目录
- `preserve_columns::Union{Vector{Int}, UnitRange{Int}}=1:2`: 保持不变的列索引（默认前两列，通常是k点坐标）
- `verbose::Bool=false`: 是否输出详细信息

# 返回值
- `String`: 生成的文件完整路径

# 示例
```julia
# 创建CDW关联文件
cdw_file = create_cdw_from_nn("cdw_k.bin", "nn_k.bin")
```
"""
function create_cdw_from_nn(
    output_filename::String="cdw_k.bin",
    input_filename::String="nn_k.bin",
    output_dir::String=pwd(),
    input_dir::String=pwd();
    preserve_columns::Union{Vector{Int}, UnitRange{Int}}=1:2,
    verbose::Bool=false
)
    if verbose
        println("正在创建CDW关联文件:")
        println("  输入文件: $input_filename")
        println("  输出文件: $output_filename")
        println("  保持不变的列: $preserve_columns")
    end
    
    # 假设数据列的排列是：AA, AB, BA, BB（实部和虚部交替）
    # 应用η_α η_β因子：AA=+1, AB=-1, BA=-1, BB=+1
    
    # 假设每个轨道有两列（实部和虚部）
    # 创建列因子字典
    column_factors = Dict{Int, Float64}()
    
    # 计算列索引和因子
    # AA轨道：因子 +1
    real_col_aa = preserve_columns[end] + 1
    imag_col_aa = real_col_aa + 1
    column_factors[real_col_aa] = 1.0
    column_factors[imag_col_aa] = 1.0
    
    # AB轨道：因子 -1
    real_col_ab = preserve_columns[end] + 3
    imag_col_ab = real_col_ab + 1
    column_factors[real_col_ab] = -1.0
    column_factors[imag_col_ab] = -1.0
    
    # BA轨道：因子 -1
    real_col_ba = preserve_columns[end] + 5
    imag_col_ba = real_col_ba + 1
    column_factors[real_col_ba] = -1.0
    column_factors[imag_col_ba] = -1.0
    
    # BB轨道：因子 +1
    real_col_bb = preserve_columns[end] + 7
    imag_col_bb = real_col_bb + 1
    column_factors[real_col_bb] = 1.0
    column_factors[imag_col_bb] = 1.0
    
    if verbose
        println("  应用因子：")
        println("    AA轨道 (列 $real_col_aa, $imag_col_aa): +1")
        println("    AB轨道 (列 $real_col_ab, $imag_col_ab): -1")
        println("    BA轨道 (列 $real_col_ba, $imag_col_ba): -1")
        println("    BB轨道 (列 $real_col_bb, $imag_col_bb): +1")
    end
    
    # 调用scale_bin_columns函数应用因子
    return scale_bin_columns(
        output_filename,
        input_filename,
        column_factors,
        output_dir,
        input_dir;
        preserve_columns=preserve_columns,
        default_factor=1.0,
        verbose=verbose
    )
end

"""
    combine_cdwpair_components(
        output_filename::String="cdwpair_k.bin",
        cdw_filename::String="cdw_k.bin",
        pair_filename::String="pair_onsite_k.bin",
        output_dir::String=pwd(),
        input_dir::String=pwd(),
        cdw_weight::Real=1.0,
        pair_weight::Real=2.0;
        preserve_columns::Union{Vector{Int}, UnitRange{Int}}=1:2,
        verbose::Bool=false
    )

将CDW关联和配对关联合并为一个整体结构因子（假设了eta-SU（2）对称性）。
如果CDW关联文件不存在，会自动从密度-密度关联文件创建。

# 参数
- `output_filename::String="cdwpair_k.bin"`: 输出文件名
- `cdw_filename::String="cdw_k.bin"`: CDW关联文件名
- `pair_filename::String="pair_onsite_k.bin"`: 配对关联文件名
- `output_dir::String=pwd()`: 输出文件目录
- `input_dir::String=pwd()`: 输入文件目录
- `cdw_weight::Real=1.0`: CDW关联的权重（默认为1.0）
- `pair_weight::Real=2.0`: 配对关联的权重（默认为2.0）
- `preserve_columns::Union{Vector{Int}, UnitRange{Int}}=1:2`: 保持不变的列索引（默认前两列，通常是k点坐标）
- `verbose::Bool=false`: 是否输出详细信息

# 返回值
- `String`: 生成的文件完整路径

# 示例
```julia
# 合并CDW和配对关联
combined_file = combine_cdwpair_components("cdwpair_k.bin")
```
"""
function combine_cdwpair_components(
    output_filename::String="cdwpair_k.bin",
    cdw_filename::String="cdw_k.bin",
    pair_filename::String="pair_onsite_k.bin",
    output_dir::String=pwd(),
    input_dir::String=pwd(),
    cdw_weight::Real=1.0,
    pair_weight::Real=2.0;
    preserve_columns::Union{Vector{Int}, UnitRange{Int}}=1:2,
    verbose::Bool=false
)
    # 检查CDW文件是否存在，如果不存在，尝试从nn_k.bin创建
    cdw_path = joinpath(input_dir, cdw_filename)
    nn_filename = "nn_k.bin"
    
    if !isfile(cdw_path)
        if verbose
            println("CDW文件不存在: $cdw_path")
            println("尝试从密度-密度关联文件创建...")
        end
        
        # 创建CDW文件
        created_cdw_path = create_cdw_from_nn(
            cdw_filename,
            nn_filename,
            output_dir,
            input_dir;
            preserve_columns=preserve_columns,
            verbose=verbose
        )
        
        if created_cdw_path == ""
            @error "无法找到或者创建CDW文件"
            return ""
        end
    end
    
    # 检查配对文件是否存在
    pair_path = joinpath(input_dir, pair_filename)
    if !isfile(pair_path)
        @error "找不到配对关联文件: $pair_path"
        return ""
    end
    
    if verbose
        println("开始合并CDW和配对关联文件:")
        println("  CDW文件: $cdw_filename (权重: $cdw_weight)")
        println("  配对文件: $pair_filename (权重: $pair_weight)")
    end
    
    # 调用combine_bin_files合并文件
    combined_file = combine_bin_files(
        output_filename,
        [(cdw_filename, cdw_weight), (pair_filename, pair_weight)],
        output_dir,
        input_dir;
        preserve_columns=preserve_columns,
        operation=+,
        verbose=verbose
    )
    
    if verbose && combined_file != ""
        println("CDW和配对关联已合并: $combined_file")
    end
    
    return combined_file
end


"""    
    merge_staggered_components(
        output_filename::String="afm_sf_k.bin",
        input_filename::String="ss_k.bin",
        output_dir::String=pwd(),
        input_dir::String=pwd();
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
- `output_filename::String="afm_sf_k.bin"`: 输出文件名
- `input_filename::String="ss_k.bin"`: 输入文件名
- `output_dir::String=pwd()`: 输出目录，默认为当前目录
- `input_dir::String=pwd()`: 输入目录，默认为当前目录
- `real_columns::Vector{Int}=[3, 9, 5, 7]`: 实部列索引数组，顺序为[A, B, C, D]
- `imag_columns::Vector{Int}=[4, 10, 6, 8]`: 虚部列索引数组，顺序为[A, B, C, D]
- `preserve_columns::UnitRange{Int}=1:2`: 要保留的列范围，默认为1:2（通常是k点坐标）
- `verbose::Bool=false`: 是否输出详细信息，默认为false

# 返回值
- 合并后的文件完整路径

# 示例
```julia
# 计算反铁磁结构因子
afm_sf_file = merge_staggered_components(
    "afm_sf_k.bin",
    "ss_k.bin"
)

# 自定义列索引
result = merge_staggered_components(
    "custom_output.bin",
    "custom_input.bin",
    real_columns=[3, 5, 7, 9],
    imag_columns=[4, 6, 8, 10]
)
```
"""
function merge_staggered_components(
    output_filename::String="afm_sf_k.bin",
    input_filename::String="ss_k.bin",
    output_dir::String=pwd(),
    input_dir::String=pwd();
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
        output_filename,
        input_filename,
        real_columns,
        imag_columns,
        weights,
        output_dir,
        input_dir;
        preserve_columns=preserve_columns,
        verbose=verbose
    )
    
    return result
end


"""    
    merge_uniform_components(
        output_filename::String="cdwpair_sf_k.bin",
        input_filename::String="cdwpair_k.bin",
        output_dir::String=pwd(),
        input_dir::String=pwd();
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
- `output_filename::String="cdwpair_sf_k.bin"`: 输出文件名
- `input_filename::String="cdwpair_k.bin"`: 输入文件名
- `output_dir::String=pwd()`: 输出目录，默认为当前目录
- `input_dir::String=pwd()`: 输入目录，默认为当前目录
- `real_columns::Vector{Int}=[3, 9, 5, 7]`: 实部列索引数组，顺序为[A, B, C, D]
- `imag_columns::Vector{Int}=[4, 10, 6, 8]`: 虚部列索引数组，顺序为[A, B, C, D]
- `preserve_columns::UnitRange{Int}=1:2`: 要保留的列范围，默认为1:2（通常是k点坐标）
- `verbose::Bool=false`: 是否输出详细信息，默认为false

# 返回值
- 合并后的文件完整路径

# 示例
```julia
# 计算电荷结构因子
cdwpair_sf_file = merge_uniform_components(
    "cdwpair_sf_k.bin",
    "cdwpair_k.bin"
)

# 自定义列索引
result = merge_uniform_components(
    "custom_output.bin",
    "custom_input.bin",
    real_columns=[3, 5, 7, 9],
    imag_columns=[4, 6, 8, 10]
)
```
"""
function merge_uniform_components(
    output_filename::String="cdwpair_sf_k.bin",
    input_filename::String="cdwpair_k.bin",
    output_dir::String=pwd(),
    input_dir::String=pwd();
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
        output_filename,
        input_filename,
        real_columns,
        imag_columns,
        weights,
        output_dir,
        input_dir;
        preserve_columns=preserve_columns,
        verbose=verbose
    )
    
    return result
end