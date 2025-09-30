# DataProcessforDQMC

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://wangfh5.github.io/DataProcessforDQMC.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://wangfh5.github.io/DataProcessforDQMC.jl/dev/)
[![Build Status](https://github.com/wangfh5/DataProcessforDQMC.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/wangfh5/DataProcessforDQMC.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/wangfh5/DataProcessforDQMC.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/wangfh5/DataProcessforDQMC.jl)

专门用于处理**行列式量子蒙特卡洛（DQMC）**模拟数据的Julia包，提供高效的数据分析工具链。

## 🎯 核心功能

- **多参数分析**：批量分析不同参数配置下的模拟结果
- **结构因子计算**：AFM/CDW结构因子的多参数扫描分析
- **相关比率分析**：量化相变和临界现象
- **智能目录管理**：自动识别和解析参数目录结构
- **高性能预编译**：支持系统镜像预编译，毫秒级启动

## 🚀 快速开始

### 基本使用
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

- 📖 **[完整文档](https://wangfh5.github.io/DataProcessforDQMC.jl/stable/)**
- ⚡ **[预编译指南](docs/src/precompilation.md)** - 系统镜像构建和使用
- 🔧 **[参数分析指南](docs/src/parameter_analysis.md)** - 多参数数据分析

## 🛠️ 安装

```julia
julia> ]
pkg> add https://github.com/wangfh5/DataProcessforDQMC.jl
```