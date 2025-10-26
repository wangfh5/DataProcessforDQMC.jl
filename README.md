# DataProcessforDQMC

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://wangfh5.github.io/DataProcessforDQMC.jl/dev/)
[![Build Status](https://github.com/wangfh5/DataProcessforDQMC.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/wangfh5/DataProcessforDQMC.jl/actions/workflows/CI.yml?query=branch%3Amain)

专门用于处理**行列式量子蒙特卡洛（DQMC）**模拟数据的Julia包，提供方便的数据分析工具链。
仅供个人使用, 对数据规范有特定要求. 

## 🎯 核心功能

- **多参数分析**：批量分析不同参数配置下的模拟结果
- **结构因子计算**：AFM/CDW结构因子的多参数扫描分析
- **相关比率分析**：量化相变和临界现象
- **智能目录管理**：自动识别和解析参数目录结构
- **高性能预编译**：支持系统镜像预编译

## 🚀 快速开始

### 🛠️ 安装

```julia
julia> ]
pkg> add https://github.com/wangfh5/DataProcessforDQMC.jl
```

### 基本使用

假设`base_dir`是一个包含DQMC数据的, 且符合[文件夹命名规则](src/JobManage/JobNaming.jl)的文件夹. 

```julia
using DataProcessforDQMC

# 扫描参数目录
dirs = scan_parameter_directories(
    base_dir,
    filter_options=Dict(:prefix => "honeycomb_model", :L => 12)
)

# 多参数结构因子分析
results = analyze_AFM_structure_factor_multi_parameter(
    base_dir,
    k_point=(0.0, 0.0),
    filter_options=Dict(:U => [3.0, 4.0, 5.0])
)
```

## 📚 文档

- 在线文档：**[开发版](https://wangfh5.github.io/DataProcessforDQMC.jl/dev/)**

- ⚡ **[预编译指南](docs/src/precompilation.md)** - 系统镜像构建和使用
- 📖 **[原始数据操作指南](docs/src/bin_analysis.md)** - 对`.bin`文件进行的常用操作
- 🔧 **[参数分析指南](docs/src/multi_parameter_analysis.md)** - 多参数数据分析

## 🧭 TODO

- 统一 `single-parameter-analysis` 中的轨道对参数接口，复用数据导出模块的自动映射逻辑