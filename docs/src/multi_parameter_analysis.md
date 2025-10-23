# [多参数分析指南](@id multi_parameter_analysis)

## `filter_options`: 多文件夹扫描时的文件夹筛选

本节介绍多参数分析中 `filter_options` 这个关键字参数：谁在用、底层落点，以及最小可用写法。

### 哪些函数包含 `filter_options`

- `src/multiple-parameter-analysis/structure-factor.jl`
  - `analyze_structure_factor_multi_parameter`
  - `analyze_AFM_structure_factor_multi_parameter`
  - `analyze_CDW_structure_factor_multi_parameter`
- `src/multiple-parameter-analysis/correlation-ratio.jl`
  - `analyze_correlation_ratio_multi_parameter`
  - `analyze_AFM_correlation_ratio_multi_parameter`
  - `analyze_CDW_correlation_ratio_multi_parameter`
- `src/multiple-parameter-analysis/common-functions.jl`
  - `scan_parameter_directories`（最底层）

### `filter_options` 的写法（Dict 或 NamedTuple）

- 值类型分三类（多个条件为 AND 关系）：
  - 单值：相等匹配。支持 Number（Int/Float）、Bool、String（例如前缀字符串）。
  - 数组：集合包含匹配。数组元素可为 Number、Bool、String。
  - 二元元组：(min, max) 闭区间匹配，仅用于数值（不用于 Bool 或 String）。
- 为空表示不过滤。
- 键名需与目录解析得到的参数名一致；常见键包括 `"prefix"`, `"b"`, `"U"`, `"L"`, `"dtau"`, `"gw"`, `"lprojgw"` 等。

### 最小示例

```julia
# 结构因子：同时限定多个参数
analyze_AFM_structure_factor_multi_parameter(
    ".";
    filter_options = Dict(
        "U" => (4.0, 8.0),   # 闭区间
        "L" => [9, 12],      # 多选一
        "dtau" => 0.1        # 单值
    )
)
```

```julia
# 关联比率：使用 NamedTuple 写法
analyze_CDW_correlation_ratio_multi_parameter(
    ".";
    filter_options = (U = [4.0, 6.0], L = 12)
)
```

```julia
# 布尔参数示例（例如 lprojgw 为 true）
scan_parameter_directories(
    ".";
    filter_options = Dict("lprojgw" => true)
)
```

```julia
# 直接扫描目录（底层函数），仅示例用法
scan_parameter_directories("."; filter_options = Dict("prefix" => "your_prefix"))
```

要点：条件按 AND 组合；区间为闭区间；空 `filter_options` 不做筛选。

