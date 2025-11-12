#=
    derived-bin-generation-afm.jl

该文件包含用于生成反铁磁相关的衍生数据文件的函数。

数据生成链条：
1. 原始数据：
   - spsm_k.bin (自旋-自旋关联函数，XY分量)
   - szsz_k.bin (自旋-自旋关联函数，Z分量)
   
2. 中间数据：
   - ss_k.bin (使用combine_ss_components进行合并的自旋-自旋关联函数，XY+Z分量)
   
3. 最终数据：
   - afm_sf_k.bin (使用merge_afm_sf生成的反铁磁结构因子) 具有三种不同的源数据选择
       - spsm_k.bin
       - szsz_k.bin
       - ss_k.bin

每个函数都提供了明确的控制，不会自动判断是否需要更新文件。
用户必须明确指示是否要生成文件。
=#

export afm_k_files_generation
export combine_ss_components, merge_afm_sf

"""
    afm_k_files_generation(dir::AbstractString; afm_source::String="ss_k.bin", verbose::Bool=true)

重新生成目录中所有与反铁磁相关的衍生数据文件。

# 参数
- `dir::AbstractString`: 数据目录，默认为当前工作目录
- `afm_source::String`: 用于生成反铁磁结构因子的源文件名，默认为"ss_k.bin"
- `verbose::Bool`: 是否输出详细信息，默认为true

# 返回值
- `Dict{String, String}`: 包含生成的文件路径的字典，键为文件名，值为完整路径

# 示例
```julia
# 使用默认设置重新生成AFM相关文件
files = afm_k_files_generation()

# 指定目录和数据源文件
files = afm_k_files_generation("/path/to/data", afm_source="custom_source.bin", verbose=true)
```
"""
function afm_k_files_generation(dir::AbstractString=pwd(); afm_source::String="ss_k.bin", verbose::Bool=true)
    result = Dict{String, String}()
    
    try
        verbose && println("===== 开始生成AFM相关文件 =====")
        
        # 1. 生成ss_k.bin
        verbose && println("正在生成 ss_k.bin...")
        ss_path = combine_ss_components("spsm_k.bin", "szsz_k.bin", "ss_k.bin", dir, dir; verbose=false)
        if ss_path != ""
            result["ss_k.bin"] = ss_path
            verbose && println("✓ 成功生成 ss_k.bin")
        else
            verbose && println("⚠  生成 ss_k.bin 失败")
        end
        
        # 2. 生成afm_sf_k.bin (利用指定文件afm_source)
        verbose && println("\n正在生成 afm_sf_k.bin...")
        afm_sf_path = merge_afm_sf(afm_source, "afm_sf_k.bin", dir, dir; verbose=false)
        if isfile(afm_sf_path)
            result["afm_sf_k.bin"] = afm_sf_path
            verbose && println("✓ 成功生成 afm_sf_k.bin")
        else
            verbose && println("⚠  生成 afm_sf_k.bin 失败")
        end
        
        verbose && println("\n===== AFM相关文件生成完成 =====")
        verbose && println("成功生成 $(length(result)) 个文件")
        
        return result
    catch e
        error_msg = "生成AFM相关文件时出错: " * sprint(showerror, e)
        verbose && println("\n✗ 错误: ", error_msg)
        rethrow()
    end
end

"""
    combine_ss_components(
        spsm_filename::String="spsm_k.bin",
        szsz_filename::String="szsz_k.bin",
        output_filename::String="ss_k.bin",
        input_dir::String=pwd(),
        output_dir::String=pwd(),
        xy_weight::Real=1.0,
        zz_weight::Real=1.0;
        preserve_columns::Union{Vector{Int}, UnitRange{Int}}=1:2,
        verbose::Bool=true
    )

合并 XY 和 Z 分量的自旋相关函数数据文件。

# 参数
- `spsm_filename::String="spsm_k.bin"`: XY分量文件名
- `szsz_filename::String="szsz_k.bin"`: Z分量文件名
- `output_filename::String="ss_k.bin"`: 输出文件名
- `input_dir::String=pwd()`: 输入文件目录
- `output_dir::String=pwd()`: 输出文件目录
- `xy_weight::Real=1.0`: XY分量的权重
- `zz_weight::Real=1.0`: Z分量的权重
- `preserve_columns::Union{Vector{Int}, UnitRange{Int}}=1:2`: 保持不变的列索引
- `verbose::Bool=true`: 是否输出详细信息

# 返回值
- `String`: 生成的文件完整路径或者空字符串

# 示例
```julia
# 合并 k 空间的 XY 和 Z 分量
combined_file = combine_ss_components(
    "spsm_k.bin",
    "szsz_k.bin",
    "ss_k.bin"
)
```
"""
function combine_ss_components(
    spsm_filename::String="spsm_k.bin",
    szsz_filename::String="szsz_k.bin",
    output_filename::String="ss_k.bin",
    input_dir::String=pwd(),
    output_dir::String=pwd(),
    xy_weight::Real=1.0,
    zz_weight::Real=1.0;
    preserve_columns::Union{Vector{Int}, UnitRange{Int}}=1:2,
    verbose::Bool=true
)
    return combine_bin_files(
        output_filename,
        [(spsm_filename, xy_weight), (szsz_filename, zz_weight)],
        output_dir,
        input_dir;
        preserve_columns=preserve_columns,
        operation=+,
        verbose=verbose
    )
end

"""
    merge_afm_sf(
        input_file::AbstractString="ss_k.bin",
        output_file::AbstractString="afm_sf_k.bin",
        input_dir::AbstractString=pwd(),
        output_dir::AbstractString=pwd();
        verbose::Bool=true
    )

计算反铁磁结构因子并保存到文件。

# 参数
- `input_file::AbstractString="ss_k.bin"`: 输入文件名
- `output_file::AbstractString="afm_sf_k.bin"`: 输出文件名
- `input_dir::AbstractString=pwd()`: 输入文件目录
- `output_dir::AbstractString=pwd()`: 输出文件目录
- `verbose::Bool=true`: 是否显示详细信息

# 返回值
- `String`: 输出文件的完整路径或者空字符串

# 示例
```julia
# 使用默认设置计算反铁磁结构因子
output_path = merge_afm_sf()

# 指定输入和输出文件名及目录
output_path = merge_afm_sf("my_ss_k.bin", "my_afm_sf_k.bin", "/input/dir", "/output/dir", verbose=true)
```
"""
function merge_afm_sf(
    input_file::AbstractString="ss_k.bin",
    output_file::AbstractString="afm_sf_k.bin",
    input_dir::AbstractString=pwd(),
    output_dir::AbstractString=pwd();
    verbose::Bool=true
)
    verbose && println("从 $input_file 生成反铁磁结构因子 $output_file...")
    
    try
        # 调用merge_staggered_components函数，使用默认的列索引配置
        result = merge_staggered_components(
            input_file,
            output_file,
            input_dir,
            output_dir;
            real_columns=[3, 9, 7, 5],  # [AA, BB, AB, BA] 实部列
            imag_columns=[4, 10, 8, 6], # [AA, BB, AB, BA] 虚部列
            verbose=false
        )
        
        verbose && println("反铁磁结构因子计算完成: ", result)
        return result
    catch e
        error_msg = "计算反铁磁结构因子时出错: " * sprint(showerror, e)
        verbose && println("错误: ", error_msg)
    end
end
