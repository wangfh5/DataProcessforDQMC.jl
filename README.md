# DataProcessforDQMC

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://wangfh5.github.io/DataProcessforDQMC.jl/dev/)
[![Build Status](https://github.com/wangfh5/DataProcessforDQMC.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/wangfh5/DataProcessforDQMC.jl/actions/workflows/CI.yml?query=branch%3Amain)

[简体中文](README.zh-CN.md)

A Julia package for processing determinantal quantum Monte Carlo (DQMC) simulation data. It targets a specific data layout and is primarily for personal use.

## Key Features

- Multi-parameter analysis across different parameter configurations
- Structure factor calculations for AFM/CDW scans
- Correlation ratio analysis for phase transitions and critical behavior
- Smart directory parsing for parameterized job folders
- High-performance precompilation with custom sysimages

## Quick Start

### Installation

```julia
julia> ]
pkg> add https://github.com/wangfh5/DataProcessforDQMC.jl
```

### Basic Usage

Example scripts live in `examples/`, and can be run with `julia examples/example_xxx.jl` after preparing the data folders `examples/proj_bt_honeycomb_exact_...`.
- `example_JobManage.jl`: parse legacy job folder names, generate the new naming scheme, and batch migrate (dry-run/copy/move/remove).
- `example_bin_export.jl`: export `.bin` to CSV/JLD2, generate derived quantities (e.g., AFM structure factor), customize orbital columns, and batch export.
- `example_bin_filter.jl`: filter bin data at specific k/r points, then derive and filter for convergence checks.
- `example_custom_structure_factor.jl`: generic structure-factor analysis for custom correlators in k- and r-space, plus multi-directory batch analysis.
- `example_multi_k_structure_factor.jl`: multi-k-point structure-factor analysis in a single scan and multi-directory batches.

## TODO

- Unify the orbital-pair parameter interface in `single-parameter-analysis` and reuse the auto-mapping logic from the data export module.
