module DataProcessforDQMC

using Printf
using Statistics
using DelimitedFiles
using DataFrames
using CSV

# Utility functions
include("utilities.jl")

# 作业命名与 SimulationInfo 由独立包 DQMCJobManage 提供
using DQMCJobManage
# 过渡别名：保持 `using DataProcessforDQMC.JobManage` 可用
const JobManage = DQMCJobManage
# 重新导出 DQMCJobManage 的核心函数，让用户直接使用
export generate_jobname, parse_jobname, parse_jobname_legacy
export migrate_legacy_to_new, verify_migration
export SimulationInfo, save_simulation_info

# Package compilation support
include("packagecompile/compile.jl")
export compile, @Algorithm_str, Algorithm

# Basis statistical functions
include("statistics.jl")

# Data processing utilities
include("data-processing/bin-file-operations-basic.jl")
include("data-processing/bin-file-export.jl")
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
