# [Bin Analysis 与数据处理指南](@id bin_analysis_guide)

## 功能概览

DataProcessforDQMC.jl 提供了完整的bin数据处理和分析工作流。

### 数据处理流程

```
原始数据 (.bin 文件)
    ↓
数据处理 (四大核心操作)
├─ combine: 多文件加权组合
├─ scale:   列缩放
├─ merge:   列合并
└─ filter:  行提取 (特定坐标)
    ↓
数据导出
├─ CSV:  Excel/Python 友好
└─ JLD2: Julia 原生高效
    ↓
可视化分析
```

### 核心数据处理 (`src/data-processing/`)

**基本操作** (`bin-file-operations-basic.jl`):
- `combine_bin_files` - Combine 操作
- `scale_bin_columns` - Scale 操作  
- `merge_bin_columns` - Merge 操作
- `filter_bin_file` - Filter 操作

**导出功能** (`bin-file-export.jl`):
- `export_bin_to_dataframe` - 转换为 DataFrame
- `export_bin_to_csv` - 导出为 CSV
- `export_bin_to_jld2` - 导出为 JLD2
- `export_directory_bins` - 批量导出

**导出量生成** (`derived-bin-generation-*.jl`):
- `afm_k_files_generation` - AFM 相关文件生成工作流
- `cdw_k_files_generation` - CDW 相关文件生成工作流
- `merge_afm_sf` - AFM 结构因子 (应用 merge 操作)
- `combine_ss_components` - 自旋分量合并 (应用 combine 操作)


## 主要功能

### 1. 数据处理的四大核心操作

数据处理模块位于 `src/data-processing/`，提供四种对 `.bin` 文件的基本操作：

#### 1.1 Combine - 多文件加权组合
将多个 bin 文件按权重组合成新文件。

```julia
combine_bin_files("ss_k.bin", [("spsm_k.bin", 1.0), ("szsz_k.bin", 1.0)])
```

**应用场景**: 合并自旋分量（如 (S⁺S⁻ + S⁻S⁺)/2 + SᶻSᶻ）

#### 1.2 Scale - 列缩放
对 bin 文件的指定列进行缩放。

```julia
scale_bin_columns("scaled.bin", "input.bin", [1, 1, 0.5, 0.5, ...])
```

**应用场景**: 归一化、单位换算

#### 1.3 Merge - 列合并
将多个文件的列合并，或对列进行线性组合。

```julia
merge_afm_sf("afm_sf_k.bin", "ss_k.bin")  # S_AFM = AA + BB - AB - BA
```

**应用场景**: 生成导出量（如 AFM/CDW 结构因子）

#### 1.4 Filter - 行提取
提取特定坐标的所有 bin 数据。

```julia
# K空间：提取特定动量点的所有bin
filter_bin_file("nn_k.bin", (-0.458, -0.458), dir="data/")
# 输出: nn_k_-0.458_-0.458.bin (包含该k点的所有bin)

# R空间：提取特定格点的所有bin
filter_bin_file("nn_r.bin", (1, 2), dir="data/")
# 输出: nn_r_1_2.bin (包含该位置的所有bin)
```

**应用场景**: Bin 收敛性分析、特定点的统计涨落分析

---

### 2. 数据导出功能

#### 导出为DataFrame
```julia
# 默认导出（4个轨道：AA, AB, BA, BB）
df = export_bin_to_dataframe("nn_k.bin", dir="data/")
# 返回DataFrame，包含列: bin, kx, ky, AA_real, AA_imag, AB_real, AB_imag, ...

# 自定义轨道对配置（仅导出对角项 AA, BB）
df = export_bin_to_dataframe("nn_k.bin", dir="data/",
    orbital_columns=[(3,4), (9,10)],
    orbital_labels=["AA", "BB"])
# 返回DataFrame，包含列: bin, kx, ky, AA_real, AA_imag, BB_real, BB_imag
# 注：AA, BB等表示轨道对，如AA表示轨道A到轨道A的关联
```

#### 导出为CSV
```julia
export_bin_to_csv("nn_k.bin", dir="data/")
# 输出: data/exported_csv/nn_k.csv
```

#### 导出为JLD2
```julia
export_bin_to_jld2("nn_k.bin", dir="data/")
# 输出: data/exported_jld2/nn_k.jld2 (dataset: "nn_k")
# 读取: load("nn_k.jld2", "nn_k")

# 如需合并多个文件到一个JLD2
export_bin_to_jld2("nn_k.bin", dir="data/", output_file="all_data.jld2")  # dataset: nn_k
export_bin_to_jld2("nn_r.bin", dir="data/", output_file="all_data.jld2")  # dataset: nn_r
# 输出: data/exported_jld2/all_data.jld2 (datasets: "nn_k", "nn_r")
# 读取: load("all_data.jld2", "nn_k") 或 load("all_data.jld2", "nn_r")
```

#### 批量导出
```julia
# 默认导出为CSV和JLD2
export_directory_bins("data/", file_patterns=["*_k.bin", "*_r.bin"])

# 指定输出格式
export_directory_bins("data/", output_format=:csv)   # 仅CSV
export_directory_bins("data/", output_format=:jld2)  # 仅JLD2
export_directory_bins("data/", output_format=:both)  # CSV + JLD2（默认）
```

## 数据格式说明

### Bin文件结构
```
原始bin文件 (如 nn_k.bin):
- 5760行 = 576个k点 × 10个bin
- 按bin分块存储:
  第1-576行:    bin 0 (所有k点)
  第577-1152行:  bin 1 (所有k点)
  ...
  第5185-5760行: bin 9 (所有k点)
```

### CSV输出格式

**多轨道数据** (nn_k.csv):
```csv
bin,kx,ky,AA_real,AA_imag,AB_real,AB_imag,BA_real,BA_imag,BB_real,BB_imag
1,-0.458,-0.458,0.297,-6.62e-10,0.055,0.003,0.055,-0.003,0.289,-8.10e-10
2,-0.458,-0.458,0.292,-5.75e-10,0.055,-0.002,0.055,0.002,0.299,-3.34e-10
...
```

**单列数据** (afm_sf_k.csv):
```csv
bin,kx,ky,value_real,value_imag
1,-0.458,-0.458,0.503,-9.38e-10
2,-0.458,-0.458,0.497,-6.46e-10
...
```

**过滤后的数据** (afm_sf_k_0.500_0.500.csv):
```csv
bin,kx,ky,value_real,value_imag
1,0.5,0.5,0.496,-9.28e-10
2,0.5,0.5,0.497,-6.46e-10
...
10,0.5,0.5,0.496,-1.15e-9
```


## 典型工作流

### 工作流1: 生成导出量（derive 应用）

使用数据处理操作生成物理导出量：

```julia
using DataProcessforDQMC

# 应用示例：生成 AFM 相关文件
# 内部使用 combine + merge 操作
afm_k_files_generation("data/")

# Step 2: 过滤特定k点 (filter 操作)
filter_bin_file("afm_sf_k.bin", (0.5, 0.5), dir="data/")

# Step 3: 导出为CSV用于可视化
export_bin_to_csv("afm_sf_k_0.500_0.500.bin", dir="data/")

# Step 4: 在Python中绘图
# import pandas as pd
# import matplotlib.pyplot as plt
# df = pd.read_csv('afm_sf_k_0.500_0.500.csv')
# plt.plot(df['bin'], df['value_real'])
# plt.xlabel('Bin index')
# plt.ylabel('S_AFM(π,π)')
```

**涉及操作**: derive 工作流 (combine+merge) → filter → export

---

### 工作流2: 批量数据导出

```julia
using DataProcessforDQMC

# 导出目录下所有关联函数（默认为CSV和JLD2）
export_directory_bins(
    "data/",
    file_patterns=["nn_*.bin", "spsm_*.bin", "afm_sf_*.bin"],
    iscorrelation=true
)

# 在Excel中打开：exported_csv/*.csv
# 在Python中批量分析：
# import pandas as pd
# import glob
# for f in glob.glob('exported_csv/*.csv'):
#     df = pd.read_csv(f)
#     # 进行分析...

# 在Julia中批量加载（使用JLD2）：
# using JLD2, DataFrames
# nn_k_df = load("exported_jld2/nn_k.jld2", "nn_k")
# nn_r_df = load("exported_jld2/nn_r.jld2", "nn_r")
# spsm_k_df = load("exported_jld2/spsm_k.jld2", "spsm_k")
```

**涉及操作**: export

---

### 工作流3: 多k点关联函数的Bin涨落分析

**目的**: 检查DQMC模拟是否已经热化并达到统计平衡，通过观察不同bin的涨落判断收敛性。

**分析指标**:
- 前几个bin: 可能处于热化阶段，数值逐渐稳定
- 后续bin: 应呈现统计涨落（上下振荡），围绕平均值波动
- 收敛标志: bin数据在误差范围内上下振荡，无明显趋势

#### 方案1: 先筛选特定k点再导出（推荐用于关注少数k点）

```julia
using DataProcessforDQMC

# 指定关心的k点（通常是高对称点）
k_points = [
    (0.0, 0.0),        # Γ点
    (0.5, 0.0),        # X点  
    (0.5, 0.5)         # M点
]

for k in k_points
    # 1. 筛选特定k点的所有bin数据
    filter_bin_file("nn_k.bin", k, dir="data/")
    
    # 2. 导出为CSV用于Excel分析
    kx_str = @sprintf("%.3f", k[1])
    ky_str = @sprintf("%.3f", k[2])
    filename = "nn_k_$(kx_str)_$(ky_str).bin"
    export_bin_to_csv(filename, dir="data/")
end

# 3. 在Excel中打开CSV文件，观察每个k点的bin演化
#    - 横轴: bin index (1, 2, 3, ...)
#    - 纵轴: 关联函数值
#    - 检查: 数据是否在后期呈现稳定的上下振荡
```

#### 方案2: 直接导出全部数据（适合系统性分析）

```julia
using DataProcessforDQMC

# 直接导出完整的关联函数数据
export_bin_to_csv("nn_k.bin", dir="data/")

# 在Excel中：
# 1. 打开 exported_csv/nn_k.csv
# 2. 使用筛选功能选择特定k点 (kx, ky)
# 3. 绘制 bin vs value_real 曲线
# 4. 观察数据点的演化趋势
```

**收敛性判断**:
- ✅ **良好收敛**: bin 3-10 在平均值附近上下振荡，标准差较小
- ⚠️ **未完全热化**: 前几个bin数值持续增大或减小
- ❌ **未收敛**: 所有bin呈现明显趋势，无稳定振荡

**涉及操作**: filter → export (或直接 export)