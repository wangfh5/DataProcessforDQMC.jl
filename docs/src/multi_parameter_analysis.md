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

## 自定义关联函数分析

本节介绍如何在多参数分析框架中使用匿名函数自定义关联函数分析，适用于 `analyze_structure_factor_multi_parameter` 等批处理函数。核心在于理解参数的**自动继承**、**特殊参数处理**和**固定参数**三个机制。

完整示例可参考 `examples/example_custom_structure_factor.jl`。

### 参数分类与处理原则

#### 1. 自动继承参数（无需显式声明）

以下参数由多参数分析框架自动传递给单目录分析器，**在匿名函数中无需提及**：

- `startbin`、`endbin`、`dropmaxmin`：时间序列的 binning 配置
- `auto_digits`、`k_point_tolerance`、`verbose`：通用配置参数

这些参数通过 `kwargs...` 自动转发，保持代码简洁。

#### 2. 特殊参数（仅需过滤声明）

`force_rebuild` 和 `source_file` 是**唯二**需要在匿名函数签名中声明的参数：

```julia
(coordinate, filename, dir; force_rebuild=false, source_file="", kwargs...) -> ...
```

**重要说明**：
- 这两个参数由框架传入但**不会被底层分析器使用**（即 `StructureFactorAnalysis` 不接受它们）， 这就是为什么需要包装`StructureFactorAnalysis`的其中一个重要原因。
- 不建议在调用 `analyze_structure_factor_multi_parameter` 时设置 `force_rebuild`，应保持默认值 `false`
  - 原因：`force_rebuild=false` 允许框架执行文件存在性检查，跳过缺失源文件的目录
  - 详见源码：`src/multiple-parameter-analysis/structure-factor.jl:L177-L184`

#### 3. 固定参数（需显式指定）

`real_column` 和 `imag_column` 是底层单参数分析器特有的参数，多参数分析函数不接受这两个参数， 需在匿名函数中固定：

```julia
(coordinate, filename, dir; force_rebuild=false, source_file="", kwargs...) -> 
    StructureFactorAnalysis(
        coordinate, filename, dir;
        real_column=7,     # Fix: 指定实部列（如 AB 轨道实部）
        imag_column=8,     # Fix: 指定虚部列（如 AB 轨道虚部）
        kwargs...          # Auto-inherit: startbin, dropmaxmin, etc.
    )
```

通过修改 `real_column` 和 `imag_column`, 可以分析不同的轨道对组合。

### 标准模式（推荐写法）

```julia
# 示例：分析 k 空间 G_AB 关联函数（AB 轨道对应 columns 7,8）
analyze_structure_factor_multi_parameter(
    (k_point, filename, dir; force_rebuild=false, source_file="", kwargs...) -> 
        StructureFactorAnalysis(
            k_point, filename, dir;
            real_column=7,     # Fix: AB orbital real part
            imag_column=8,     # Fix: AB orbital imag part
            kwargs...          # Auto-inherit: startbin, dropmaxmin, etc.
        ),
    k_point_list,
    base_dir;
    filename="cpcm_k.bin",
    # 注意：不设置 force_rebuild，保持默认值 false 以启用文件检查
    startbin=2,
    dropmaxmin=0,
    filter_options=(L=24, U=(3.5, 4.0))
)
```

### 实空间与 k 空间的统一处理

该框架对坐标类型（k-point 或 r-point）无感知，实空间和 k 空间使用**完全相同**的调用模式：

```julia
# k 空间示例（浮点坐标）
analyze_structure_factor_multi_parameter(
    (k, filename, dir; force_rebuild=false, source_file="", kwargs...) -> 
        StructureFactorAnalysis(k, filename, dir; real_column=7, imag_column=8, kwargs...),
    (0.0, 0.0),  # Γ 点
    base_dir; filename="cpcm_k.bin", ...
)

# 实空间示例（整数坐标）
analyze_structure_factor_multi_parameter(
    (r, filename, dir; force_rebuild=false, source_file="", kwargs...) -> 
        StructureFactorAnalysis(r, filename, dir; real_column=7, imag_column=8, kwargs...),
    (12, 12),  # 单元格坐标
    base_dir; filename="cpcm_r.bin", ...
)
```

