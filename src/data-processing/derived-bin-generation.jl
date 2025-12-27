#=
    derived-bin-generation.jl

该文件整合了所有衍生数据生成功能，提供统一的接口用于生成和管理衍生数据文件。

包含的数据生成链条：

1. AFM相关数据链：
   - 原始数据: spsm_k.bin, szsz_k.bin
   - 中间数据: ss_k.bin
   - 最终数据: afm_sf_k.bin （可以源于spsm_k.bin, szsz_k.bin, ss_k.bin）

2. CDW和配对相关数据链：
   - 原始数据: nn_k.bin, pair_onsite_k.bin
   - 中间数据: cdw_k.bin, cdwpair_k.bin
   - 最终数据: cdw_sf_k.bin, pair_onsite_sf_k.bin, cdwpair_sf_k.bin

所有函数都提供明确的控制，不会自动判断是否需要更新文件。
用户必须明确指示是否要生成文件。
=#

export generate_all_derived_files, generate_all_derived_files_recursive

"""    
    generate_all_derived_files(dir::AbstractString=pwd();
                               afm_source::String="ss_k.bin", 
                               verbose::Bool=true)

生成目录中所有衍生数据文件，包括反铁磁和电荷密度波/配对相关的文件。

# 参数
- `dir::AbstractString=pwd()`: 数据目录，默认为当前目录
- `afm_source::String="ss_k.bin"`: 用于生成反铁磁结构因子的源文件名，默认为"ss_k.bin"
- `verbose::Bool=true`: 是否输出详细信息，默认为true

# 返回值
- `Dict{String, String}`: 包含生成的文件路径的字典，键为文件名，值为完整路径

# 示例
```julia
# 在当前目录生成所有衍生数据文件
files = generate_all_derived_files(verbose=true)

# 指定目录和源文件
files = generate_all_derived_files(dir="/path/to/data", afm_source="spsm_k.bin", verbose=true)

# 如果只需要生成特定类型的文件，可以直接调用相应的函数
afm_files = afm_k_files_generation(dir, afm_source="spsm_k.bin")
cdw_files = cdwpair_k_files_generation(dir)
```
"""
function generate_all_derived_files(dir::AbstractString=pwd();
                                    afm_source::String="ss_k.bin", 
                                    verbose::Bool=true)
    result = Dict{String, String}()
    
    verbose && println("===== 开始生成衍生数据文件 =====")
    verbose && println("目录: $dir")
    
    # 生成AFM相关文件
    afm_files = afm_k_files_generation(dir, afm_source=afm_source, verbose=verbose)
    merge!(result, afm_files)
    
    # 生成CDW和配对相关文件
    cdw_files = cdwpair_k_files_generation(dir, verbose=verbose)
    merge!(result, cdw_files)
    
    verbose && println("\n===== 衍生数据文件生成完成 =====")
    verbose && println("总共生成 $(length(result)) 个文件")
    
    # 打印生成的文件列表
    if verbose && !isempty(result)
        println("\n生成的文件列表:")
        for (filename, path) in sort(collect(result))
            println("- $filename: $path")
        end
    end
    
    return result
end

"""    generate_all_derived_files_recursive(base_dir::AbstractString=pwd();
                                         pattern::Regex=r".", 
                                         max_depth::Int=1,
                                         afm_source::String="ss_k.bin", 
                                         verbose::Bool=true)

递归地生成多个目录中的所有衍生数据文件。

# 参数
- `base_dir::AbstractString=pwd()`: 基础目录，默认为当前目录
- `pattern::Regex=r"."`: 用于匹配目录名的正则表达式
- `max_depth::Int=1`: 最大递归深度
- `afm_source::String="ss_k.bin"`: 用于生成反铁磁结构因子的源文件名，默认为"ss_k.bin"
- `verbose::Bool=true`: 是否输出详细信息，默认为true

# 返回值
- `Dict{String, Dict{String, String}}`: 包含每个目录生成文件的嵌套字典

# 示例
```julia
# 在当前目录及其子目录中生成所有衍生数据文件
results = generate_all_derived_files_recursive(verbose=true)

# 只处理符合特定模式的目录
results = generate_all_derived_files_recursive(pattern=r"^proj_fft_", verbose=true)

# 递归处理2层目录
results = generate_all_derived_files_recursive(max_depth=2, verbose=true)

# 指定源文件
results = generate_all_derived_files_recursive(
    base_dir="/path/to/data",
    afm_source="spsm_k.bin",
    verbose=true
)
```
"""
function generate_all_derived_files_recursive(base_dir::AbstractString=pwd(); 
                                              pattern::Regex=r".", 
                                              max_depth::Int=1,
                                              afm_source::String="ss_k.bin", 
                                              verbose::Bool=true)
    result = Dict{String, Dict{String, String}}()
    
    # 处理当前目录
    verbose && println("\n===== 处理目录: $base_dir =====")
    
    # 检查当前目录是否有原始数据文件
    has_data = isfile(joinpath(base_dir, "spsm_k.bin")) || 
               isfile(joinpath(base_dir, "ss_k.bin")) || 
               isfile(joinpath(base_dir, "nn_k.bin")) || 
               isfile(joinpath(base_dir, "cdw_k.bin")) || 
               isfile(joinpath(base_dir, "cdwpair_k.bin"))
    
    if has_data
        dir_result = generate_all_derived_files(
            base_dir; 
            afm_source=afm_source, 
            verbose=verbose
        )
        
        if !isempty(dir_result)
            result[base_dir] = dir_result
        end
    else
        verbose && println("没有需要处理的原始数据文件")
    end
    
    # 如果达到最大深度，不再递归
    if max_depth <= 0
        return result
    end
    
    # 递归处理子目录
    for item in readdir(base_dir, join=true)
        if isdir(item) && match(pattern, basename(item)) !== nothing
            sub_results = generate_all_derived_files_recursive(
                item; 
                pattern=pattern, 
                max_depth=max_depth-1,
                afm_source=afm_source, 
                verbose=verbose
            )
            
            merge!(result, sub_results)
        end
    end
    
    return result
end