# Parameter Analysis Guide

This guide explains how to use the parameter analysis functionality in DataProcessforDQMC.

## Table of Contents
- [Introduction](#introduction)
- [Single Parameter Analysis](#single-parameter-analysis)
- [Multi-Parameter Analysis](#multi-parameter-analysis)
- [Command Line Interface](#command-line-interface)
- [Advanced Usage](#advanced-usage)

## Introduction

The parameter analysis module provides tools for analyzing DQMC simulation data across different parameter values. It supports both single-parameter and multi-parameter analysis, making it easy to study how physical quantities depend on simulation parameters.

## Single Parameter Analysis

### Basic Usage

```julia
using DataProcessforDQMC

# Analyze energy data from a single directory
results = EnergyAnalysis("energy.bin", ".")

# Analyze correlation data
results = CorrelationAnalysis("spsm_r.bin", ".")
```

### Available Functions

- `EnergyAnalysis(filename, dir; kwargs...)`: Analyze energy data
- `CorrelationAnalysis(filename, dir; kwargs...)`: Analyze correlation data
- `MultiOrbitalCorrelationAnalysis(filename, dir; kwargs...)`: Analyze multi-orbital correlation data
- `RenyiNegativity(filename, dir; kwargs...)`: Calculate Renyi negativity

### Common Keyword Arguments

- `startbin`: First bin to include in analysis (default: 1)
- `endbin`: Last bin to include (default: end of data)
- `dropmaxmin`: Number of maximum and minimum values to drop (default: 1)
- `verbose`: Print progress information (default: true)

## Multi-Parameter Analysis

### Basic Usage

```julia
using DataProcessforDQMC

# Analyze energy data across multiple parameter values
results = analyze_multiple_parameters(
    ".",                     # Base directory
    "U",                     # Parameter name
    ["2.0", "3.0", "4.0"],  # Parameter values
    file_pattern="energy.bin",
    analysis_type="energy"
)
```

### Available Analysis Types

- `"energy"`: Analyze energy data (default)
- `"correlation"`: Analyze correlation data
- `"multiorbital"`: Analyze multi-orbital correlations

## Command Line Interface

You can also use the parameter analysis functionality from the command line:

```bash
# Single parameter analysis
julia -e 'using DataProcessforDQMC; EnergyAnalysis("energy.bin", ".")'

# Multi-parameter analysis
julia -e 'using DataProcessforDQMC; include("src/parameter_analysis.jl")' -- \
    --parameter U \
    --values 2.0,3.0,4.0 \
    --analyze energy \
    --output results.csv
```

### Command Line Options

- `--parameter`: Name of the parameter to analyze (default: "U")
- `--values`: Comma-separated list of parameter values (default: "1.0,2.0,3.0")
- `--analyze`: Type of analysis to perform (default: "energy")
- `--output`: Output file for results (default: "results.csv")

## Advanced Usage

### Custom Analysis

You can implement custom analysis functions and use them with the parameter analysis tools:

```julia
# Define a custom analysis function
function my_analysis(filename, dir; kwargs...)
    # Your analysis code here
    return result
end

# Use the custom analysis
results = analyze_multiple_parameters(
    ".", 
    "U", 
    ["2.0", "3.0", "4.0"],
    analysis_func=my_analysis
)
```

### Parallel Processing

For large datasets, you can use Julia's built-in parallel processing:

```julia
using Distributed
addprocs(4)  # Add 4 worker processes

@everywhere using DataProcessforDQMC

# The analysis will now run in parallel
results = analyze_multiple_parameters(
    ".", 
    "U", 
    ["2.0", "3.0", "4.0"],
    parallel=true
)
```

## Troubleshooting

### Common Issues

1. **Missing Data Files**: Ensure the data files exist in the expected locations
2. **Incorrect Parameter Values**: Verify that the parameter values match the directory names
3. **Memory Issues**: For large datasets, consider processing in smaller batches

### Getting Help

If you encounter any issues, please open an issue on the [GitHub repository](https://github.com/wangfh5/DataProcessforDQMC.jl).
