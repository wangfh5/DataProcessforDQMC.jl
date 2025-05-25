module DataProcessforDQMC

using Printf
using Statistics
using DelimitedFiles
using DataFrames
using CSV

# Basis statistical functions
include("statistics.jl")

# Legacy functions - based on dataframes.jl
include("format-dataframes.jl")
include("read-dataframes.jl")

# Data processing utilities
include("data-processing/bin-file-operations-basic.jl")
include("data-processing/bin-file-operations-afm-cdwpair.jl")

# New functions - simple implementation of analysis scripts
include("single-parameter-analysis/scalar-measurements.jl")
include("single-parameter-analysis/correlations-common.jl")
include("single-parameter-analysis/correlations-rspace.jl")
include("single-parameter-analysis/correlations-kspace.jl")
include("single-parameter-analysis/structure-factors.jl")
include("multiple-parameter-analysis.jl")


end  # module DataProcessforDQMC
