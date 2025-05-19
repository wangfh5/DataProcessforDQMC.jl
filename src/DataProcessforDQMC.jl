module DataProcessforDQMC

# Basis statistical functions
include("statistics.jl")

# Legacy functions - based on dataframes.jl
include("format-dataframes.jl")
include("read-dataframes.jl")

# New functions - simple implementation of analysis scripts
include("single-parameter-analysis.jl")
# include("multiple-parameter-analysis.jl")

end  # module DataProcessforDQMC
