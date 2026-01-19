# DataProcessforDQMC

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://wangfh5.github.io/DataProcessforDQMC.jl/dev/)
[![Build Status](https://github.com/wangfh5/DataProcessforDQMC.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/wangfh5/DataProcessforDQMC.jl/actions/workflows/CI.yml?query=branch%3Amain)

英文版 README: [README.md](README.md)

专门用于处理 **行列式量子蒙特卡洛（DQMC）** 模拟数据的Julia包，提供方便的数据分析工具链。
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

示例脚本集中在 `examples/`（本地路径：`/Users/ssqc/Projects/DataProcessforDQMC.jl/examples`）, 在该目录中准备好符合规范的数据文件夹 `examples/proj_bt_honeycomb_exact_...` 后, 可以直接运行.
- `example_JobManage.jl`：演示任务目录命名解析/转换，以及批量迁移（dry-run/cp/mv/rm）。
- `example_bin_export.jl`：演示 `.bin` 导出 CSV/JLD2、派生量（如 AFM 结构因子）、自定义轨道列与批量导出。
- `example_bin_filter.jl`：演示按特定 k/r 点过滤 bin 数据，并配合派生量做收敛分析。
- `example_custom_structure_factor.jl`：演示用通用接口分析自定义关联函数，覆盖 k/r 空间与多目录批处理。
- `example_multi_k_structure_factor.jl`：演示一次扫描多个 k 点的结构因子分析及多目录批处理。

## 🧭 TODO

- 统一 `single-parameter-analysis` 中的轨道对参数接口，复用数据导出模块的自动映射逻辑
