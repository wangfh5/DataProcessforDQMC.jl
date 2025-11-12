#=
    derived-bin-generation-cdwpair.jl

该文件包含用于生成电荷密度波和配对相关的衍生数据文件的函数。

数据生成链条：
1. 原始数据：
   - nn_k.bin (密度-密度关联函数)
   - pair_onsite_k.bin (配对关联函数)
   
2. 中间数据：
   - cdw_k.bin (使用create_cdw_from_nn从nn_k.bin生成的CDW关联函数)
   - cdwpair_k.bin (使用combine_cdwpair_components合并的CDW和配对关联函数)
   
3. 最终数据：
   - cdwpair_sf_k.bin (使用merge_cdw_sf生成的CDW结构因子) 具有三种不同的源数据选择
       - cdw_k.bin
       - pair_onsite_k.bin
       - cdwpair_k.bin

每个函数都提供了明确的控制，不会自动判断是否需要更新文件。
用户必须明确指示是否要生成文件。
=#

export cdwpair_k_files_generation
export create_cdw_from_nn, combine_cdwpair_components, merge_cdw_sf


"""
    cdwpair_k_files_generation(dir::AbstractString=pwd(); cdwpair_source::String="cdwpair_k.bin", verbose::Bool=true)

生成目录中与电荷密度波和配对相关的衍生数据文件。

# 参数
- `dir::AbstractString`: 数据目录，默认为当前工作目录
- `cdwpair_source::String`: 用于生成CDW配对结构因子的源文件名，默认为"cdwpair_k.bin"
- `verbose::Bool`: 是否输出详细信息，默认为true

# 返回值
- `Dict{String, String}`: 包含生成的文件路径的字典，键为文件名，值为完整路径

# 示例
```julia
# 使用默认设置生成CDW和配对相关文件
files = cdwpair_k_files_generation()

# 指定目录和源文件
files = cdwpair_k_files_generation("/path/to/data", cdwpair_source="custom_cdwpair_k.bin", verbose=true)
```
"""
function cdwpair_k_files_generation(dir::AbstractString=pwd(); cdwpair_source::String="cdwpair_k.bin", verbose::Bool=true)
    result = Dict{String, String}()
    
    try
        verbose && println("===== 开始生成CDW和配对相关文件 =====")
        
        # 1. 生成cdw_k.bin
        verbose && println("正在生成 cdw_k.bin...")
        cdw_path = create_cdw_from_nn("nn_k.bin", "cdw_k.bin", dir, dir; verbose=false)
        if cdw_path != ""
            result["cdw_k.bin"] = cdw_path
            verbose && println("✓ 成功生成 cdw_k.bin")
        else
            verbose && println("⚠  生成 cdw_k.bin 失败")
        end
        
        # 2. 生成cdwpair_k.bin
        verbose && println("\n正在生成 cdwpair_k.bin...")
        cdwpair_path = combine_cdwpair_components("cdw_k.bin", "pair_onsite_k.bin", "cdwpair_k.bin", dir, dir; verbose=false)
        if cdwpair_path != ""
            result["cdwpair_k.bin"] = cdwpair_path
            verbose && println("✓ 成功生成 cdwpair_k.bin")
        else
            verbose && println("⚠  生成 cdwpair_k.bin 失败")
        end
        
        # 3. 生成cdwpair_sf_k.bin
        verbose && println("\n正在生成 cdwpair_sf_k.bin...")
        cdwpair_sf_path = merge_cdw_sf(cdwpair_source, "cdwpair_sf_k.bin", dir, dir; verbose=false)
        if isfile(joinpath(dir, "cdwpair_sf_k.bin"))
            result["cdwpair_sf_k.bin"] = cdwpair_sf_path
            verbose && println("✓ 成功生成 cdwpair_sf_k.bin")
        else
            verbose && println("⚠  生成 cdwpair_sf_k.bin 失败")
        end
        
        verbose && println("\n===== CDW和配对相关文件生成完成 =====")
        verbose && println("成功生成 $(length(result)) 个文件")
        
        return result
    catch e
        error_msg = "生成CDW和配对相关文件时出错: " * sprint(showerror, e)
        verbose && println("\n✗ 错误: ", error_msg)
        rethrow()
    end
end

"""
    create_cdw_from_nn(
        input_filename::String="nn_k.bin",
        output_filename::String="cdw_k.bin",
        input_dir::String=pwd(),
        output_dir::String=pwd();
        preserve_columns::Union{Vector{Int}, UnitRange{Int}}=1:2,
        verbose::Bool=true
    )

将密度密度关联数据(nn_k.bin)转换为CDW关联数据。
应用η_α η_β因子到对应的轨道，其中η_α是+1（α=A轨道）或-1（α=B轨道）。

# 参数
- `input_filename::String="nn_k.bin"`: 输入密度-密度关联文件名
- `output_filename::String="cdw_k.bin"`: 输出文件名
- `input_dir::String=pwd()`: 输入文件目录
- `output_dir::String=pwd()`: 输出文件目录
- `preserve_columns::Union{Vector{Int}, UnitRange{Int}}=1:2`: 保持不变的列索引（默认前两列，通常是k点坐标）
- `verbose::Bool=true`: 是否输出详细信息

# 返回值
- `String`: 生成的文件完整路径

# 示例
```julia
# 创建CDW关联文件
cdw_file = create_cdw_from_nn("nn_k.bin", "cdw_k.bin")
```
"""
function create_cdw_from_nn(
    input_filename::String="nn_k.bin",
    output_filename::String="cdw_k.bin",
    input_dir::String=pwd(),
    output_dir::String=pwd();
    preserve_columns::Union{Vector{Int}, UnitRange{Int}}=1:2,
    verbose::Bool=true
)
    if verbose
        println("正在创建CDW关联文件:")
        println("  输入文件: $input_filename")
        println("  输出文件: $output_filename")
        println("  保持不变的列: $preserve_columns")
    end
    
    # 假设数据列的排列是：AA, BA, AB, BB（实部和虚部交替）
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
    
    # BA轨道：因子 -1（在新顺序中紧随AA）
    real_col_ba = preserve_columns[end] + 3
    imag_col_ba = real_col_ba + 1
    column_factors[real_col_ba] = -1.0
    column_factors[imag_col_ba] = -1.0
    
    # AB轨道：因子 -1（在新顺序中位于第三对）
    real_col_ab = preserve_columns[end] + 5
    imag_col_ab = real_col_ab + 1
    column_factors[real_col_ab] = -1.0
    column_factors[imag_col_ab] = -1.0
    
    # BB轨道：因子 +1
    real_col_bb = preserve_columns[end] + 7
    imag_col_bb = real_col_bb + 1
    column_factors[real_col_bb] = 1.0
    column_factors[imag_col_bb] = 1.0
    
    if verbose
        println("  应用因子：")
        println("    AA轨道 (列 $real_col_aa, $imag_col_aa): +1")
        println("    BA轨道 (列 $real_col_ba, $imag_col_ba): -1")
        println("    AB轨道 (列 $real_col_ab, $imag_col_ab): -1")
        println("    BB轨道 (列 $real_col_bb, $imag_col_bb): +1")
    end
    
    # 调用scale_bin_columns函数应用因子
    return scale_bin_columns(
        input_filename,
        output_filename,
        column_factors,
        input_dir,
        output_dir;
        preserve_columns=preserve_columns,
        default_factor=1.0,
        verbose=false
    )
end

"""
    combine_cdwpair_components(
        cdw_filename::String="cdw_k.bin",
        pair_filename::String="pair_onsite_k.bin",
        output_filename::String="cdwpair_k.bin",
        input_dir::String=pwd(),
        output_dir::String=pwd(),
        cdw_weight::Real=1.0,
        pair_weight::Real=2.0;
        preserve_columns::Union{Vector{Int}, UnitRange{Int}}=1:2,
        verbose::Bool=true
    )

将CDW关联和配对关联合并为一个整体结构因子（假设了eta-SU（2）对称性）。
如果CDW关联文件不存在，会自动从密度-密度关联文件创建。

# 参数
- `cdw_filename::String="cdw_k.bin"`: CDW关联文件名
- `pair_filename::String="pair_onsite_k.bin"`: 配对关联文件名
- `output_filename::String="cdwpair_k.bin"`: 输出文件名
- `input_dir::String=pwd()`: 输入文件目录
- `output_dir::String=pwd()`: 输出文件目录
- `cdw_weight::Real=1.0`: CDW关联的权重（默认为1.0）
- `pair_weight::Real=2.0`: 配对关联的权重（默认为2.0）
- `preserve_columns::Union{Vector{Int}, UnitRange{Int}}=1:2`: 保持不变的列索引（默认前两列，通常是k点坐标）
- `verbose::Bool=true`: 是否输出详细信息

# 返回值
- `String`: 生成的文件完整路径

# 示例
```julia
# 合并CDW和配对关联
combined_file = combine_cdwpair_components("cdw_k.bin", "pair_onsite_k.bin", "cdwpair_k.bin")
```
"""
function combine_cdwpair_components(
    cdw_filename::String="cdw_k.bin",
    pair_filename::String="pair_onsite_k.bin",
    output_filename::String="cdwpair_k.bin",
    input_dir::String=pwd(),
    output_dir::String=pwd(),
    cdw_weight::Real=1.0,
    pair_weight::Real=2.0;
    preserve_columns::Union{Vector{Int}, UnitRange{Int}}=1:2,
    verbose::Bool=true
)
    return combine_bin_files(
        output_filename,
        [(cdw_filename, cdw_weight), (pair_filename, pair_weight)],
        output_dir,
        input_dir;
        preserve_columns=preserve_columns,
        operation=+,
        verbose=verbose
    )
end

"""
    merge_cdw_sf(
        input_file::AbstractString="cdwpair_k.bin",
        output_file::AbstractString="cdwpair_sf_k.bin",
        input_dir::AbstractString=pwd(),
        output_dir::AbstractString=pwd();
        verbose::Bool=true
    )

计算电荷密度波结构因子并保存到文件。

# 参数
- `input_file::AbstractString="cdwpair_k.bin"`: 输入文件名
- `output_file::AbstractString="cdwpair_sf_k.bin"`: 输出文件名
- `input_dir::AbstractString=pwd()`: 输入文件目录
- `output_dir::AbstractString=pwd()`: 输出文件目录
- `verbose::Bool=true`: 是否显示详细信息

# 返回值
- `String`: 输出文件的完整路径或者空字符串

# 示例
```julia
# 使用默认设置计算电荷密度波结构因子
output_path = merge_cdw_sf()

# 指定输入和输出文件名及目录
output_path = merge_cdw_sf("my_input.bin", "my_output.bin", "/input/dir", "/output/dir", verbose=true)
```
"""
function merge_cdw_sf(
    input_file::AbstractString="cdwpair_k.bin",
    output_file::AbstractString="cdwpair_sf_k.bin",
    input_dir::AbstractString=pwd(),
    output_dir::AbstractString=pwd();
    verbose::Bool=true
)
    verbose && println("从 $input_file 生成电荷密度波结构因子 $output_file...")
    
    try
        # 调用merge_uniform_components函数，使用默认的列索引配置
        result = merge_uniform_components(
            input_file,
            output_file,
            input_dir,
            output_dir;
            real_columns=[3, 9, 7, 5],  # [AA, BB, AB, BA] 实部列 (按AA, BA, AB, BB顺序的列位置)
            imag_columns=[4, 10, 8, 6], # [AA, BB, AB, BA] 虚部列 (按AA, BA, AB, BB顺序的列位置)
            verbose=verbose
        )
        
        verbose && println("电荷密度波结构因子计算完成: ", result)
        return result
    catch e
        error_msg = "计算电荷密度波结构因子时出错: " * sprint(showerror, e)
        verbose && println("错误: ", error_msg)
    end
end
