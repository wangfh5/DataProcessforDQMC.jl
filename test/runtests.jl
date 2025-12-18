using DataProcessforDQMC
using Test

@testset "DataProcessforDQMC.jl" begin
    include("statistics_format_value_error_test.jl")
    include("statistics_iqr_fence_test.jl")
end
