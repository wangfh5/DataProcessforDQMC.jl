module DataProcessforDQMC

using Printf
using Statistics
using DelimitedFiles
using DataFrames
using CSV

# 导入 JobManage 子模块
include("JobManage/JobManage.jl")
using .JobManage
# 重新导出 JobManage 的核心函数，让用户直接使用
export generate_jobname, parse_jobname, parse_jobname_legacy
export migrate_legacy_to_new, verify_migration
export SimulationInfo, save_simulation_info

# Package compilation support
include("packagecompile/compile.jl")
export compile, @Algorithm_str, Algorithm

# Basis statistical functions
include("statistics.jl")

# Legacy functions - based on dataframes.jl
include("format-dataframes.jl")
include("read-dataframes.jl")

# Data processing utilities
include("data-processing/bin-file-operations-basic.jl")
include("data-processing/derived-bin-generation.jl")
include("data-processing/derived-bin-generation-afm.jl")
include("data-processing/derived-bin-generation-cdwpair.jl")

# New functions - simple implementation of analysis scripts
include("single-parameter-analysis/scalar-measurements.jl")
include("single-parameter-analysis/correlations-common.jl")
include("single-parameter-analysis/correlations-rspace.jl")
include("single-parameter-analysis/correlations-kspace.jl")
include("single-parameter-analysis/structure-factor.jl")
include("single-parameter-analysis/correlation-ratio.jl")
include("multiple-parameter-analysis/common-functions.jl")
include("multiple-parameter-analysis/structure-factor.jl")
include("multiple-parameter-analysis/correlation-ratio.jl")

end  # module DataProcessforDQMC
