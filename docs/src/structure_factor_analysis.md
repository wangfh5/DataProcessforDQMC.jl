# [结构因子分析设计](@id structure_factor_design)

本文档介绍结构因子分析模块的核心设计决策。

## 设计原则

### 1. 多 k 点分析作为核心函数

**核心实现只读文件一次：**

```julia
# 内部实现（私有函数）
function _multi_k_structure_factor_analysis_core(k_points::Vector, ...)
    data_array = readdlm(filepath, Float64)  # 只读一次
    
    results = StructureFactorResult[]
    for k in k_points
        # 从同一个 data_array 提取不同 k 点数据
        push!(results, StructureFactorResult(...))
    end
    return results
end
```

### 2. 将统计结果打包为结构体

```julia
struct StructureFactorResult
    k_point::Tuple{Float64, Float64}
    mean_real::Float64
    err_real::Float64
    mean_imag::Float64
    err_imag::Float64
    formatted_real::String
    formatted_imag::String
end
```

**优势：**
- 类型安全，性能更好
- 支持方法扩展（可添加 `Base.show` 等）
- 字段访问与 NamedTuple 完全相同：`result.mean_real`

**向后兼容：**
```julia
# 多参数接口使用 getproperty，对两者都有效
getproperty(result, :mean_real)  # Struct 和 NamedTuple 都支持
```

### 3. 多重分派实现单/多 k 点接口

```julia
# 单 k 点 → 返回单个结果
StructureFactorAnalysis(k::Tuple, ...) -> StructureFactorResult

# 多 k 点 → 返回结果数组
StructureFactorAnalysis(k_list::Vector{Tuple}, ...) -> Vector{StructureFactorResult}
```

同一个函数名，根据参数类型自动选择实现。

## 使用示例

### 单 k 点

```julia
result = StructureFactorAnalysis((π, π), "spsm_k.bin"; real_column=5, imag_column=6)
println("$(result.formatted_real) + $(result.formatted_imag)i")
```

### 多 k 点（高效）

```julia
k_list = [(1/3+1/L, 1/3), (1/3-1/L, 1/3), (1/3, 1/3+1/L)]
results = StructureFactorAnalysis(k_list, "cpcm_k.bin"; real_column=5, imag_column=6)

for r in results
    magnitude = sqrt(r.mean_real^2 + r.mean_imag^2)
    println("k=$(r.k_point): |G|=$magnitude")
end
```

### AFM/CDW 专用函数

```julia
result = AFMStructureFactor((0.0, 0.0), "afm_sf_k.bin")
result = CDWStructureFactor((0.0, 0.0), "cdwpair_sf_k.bin")
```

这些函数内部调用通用的 `StructureFactorAnalysis`。

## 参考

- 源码：`src/single-parameter-analysis/structure-factor.jl`
- API 文档：[API Reference](@ref api_reference)
